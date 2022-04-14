depression_na <- read.csv("depression_na.csv")
depression_na <- depression_na$x
y_depna <- array(depression_na,dim = c(193,13,15))
rm(depression_na)
dep_gene <- read.csv("depression_name.csv")
dep_gene <- dep_gene$x
#y <- y_depna
y <- aperm(y_depna,c(2,1,3))
rm(y_depna)
library(RcppArmadillo)
Rcpp::sourceCpp('sparse.cpp')

set.seed(1)
{
  c1track <- list()
  c2track <- list()
  c3track <- list()
  #ztrack <- list()
  #b1track <- list()
  #b2track <- list()
  #//initialization
  al=1;bl=3;av=1;bv=10;am=1;bm=1;sigmal=1;mub=-2
  ;sigmab=1;arho=1;brho=1;psi_1=1/3;psi0=1/3;psi1=1/3;as=1;bs=1;sigmamu=1;
  
  #//dimension
  d=dim(y);
  
  #//class number
  r=1;
  #//observation and latent
  #cube z=randi<cube>(d[0],d[1],d[2],distr_param(-1,1));
  #//class
  #c3o[c3o == 1] = 2 * (runif(sum(c3o == 1)) < 0.5) - 1
  c1 = matrix(as.numeric(runif(d[1] * r) < 0.5), d[1], r)
  c2 =matrix(as.numeric(runif(d[2] * r) < 0.5), d[2], r)
  
  c3 = matrix(as.numeric(runif(d[3] * r) < 0.5), d[3], r)
  #c3[c3 == 1] = 2 * (runif(sum(c3 == 1)) < 0.5) - 1
  #c1 <- c1o
  #c2 <- c2o
  #c3 <- c3o
  #z <- array(0,dim = d)
  #mat lambda1(d[2],r,fill::randn);
  #mat lambda2(d[2],r,fill::randn);
  #lambda1 = matrix(lambda1[,1],ncol = 1)
  #lambda2 = matrix(lambda2[,1],ncol = 1)
  #lambda1 <- matrix(rgamma(d[3]*r, al, 1/bl), d[3],r)
  #lambda2 <- matrix(rgamma(d[3]*r, al, 1/bl), d[3],r)
  lambda1 <- matrix(rgamma(d[3]*r,1,1/3),nrow = d[3],byrow = TRUE)
  lambda2 <- lambda1
  #//paramters
  mu1 <- apply(y, 1, mean,na.rm=TRUE)
  #mu1 <- mu1-mu1
  mu2 <-apply(y, 2, mean,na.rm=TRUE)
  mu2 <- mu2-mu2
  mu3 <- apply(y, 3, mean,na.rm=TRUE)
  sigma2 <- rep(1,d[3])
  v1 <- 5*sigma2
  v2 <- v1
  #b1 <- matrix(rnorm(d[2]*d[3]),nrow = d[2],ncol = d[3])
  b1 <- matrix(-2.3,nrow = d[2],ncol = d[3])
  
  b2 <-b1
  #mat mu(d[1],d[2],fill::zeros);
  #mat v1(d[1],d[2],fill::randu);
  #mat v2(d[1],d[2],fill::randu);
  #mat sigma2(d[1],d[2],fill::ones);
  #mat b1(d[1],d[2],fill::randn);
  #mat b2(d[1],d[2],fill::randn);
  m=1;rho=0.3;g=1;
  
  #vec gamma(3,fill::randu);
  gamma <- c(rho/2,1-rho,rho/2)
  #z <- z0
  z <- array(0,dim = d)
  z[which(is.na(y))] <- NA
}
#//MCMC update
set.seed(2)
library(progress)
pb <- progress_bar$new(format = "  complete [:bar] :percent eta: :eta",
                       total = 10000, clear = FALSE, width= 60)
for(it in 1:10000){
  #ztrack[[it]] <- z
  
  
  z <- update_z(c1, c2, c3, lambda1, lambda2, d, z, b1, b2, v1, v2, av, bv, y, mu1, mu2, mu3, sigma2);
  v1 <- update_v1(z,y,mu1,mu2,mu3,v1,v2,d,av,bv,sigma2)
  v2 <- update_v2(z,y,mu1,mu2,mu3,v1,v2,d,av,bv,sigma2)
  li1 <- update_c1(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma, m,g,rho,al,bl)
  c1 <- li1[[1]]
  c2 <- li1[[2]]
  c3 <- li1[[3]]
  lambda1 <- li1[[4]]
  lambda2 <- li1[[5]]
  r <- ncol(c1)
  
  c2 <- update_c2(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,rho);
  
  li2 <-  update_c3(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma,sigmal,al,bl)
  c3 <-li2[[1]]
  lambda1 <- li2[[2]]
  lambda2 <- li2[[3]]
  
  m <- rgamma(1,am+r,bm+har(d[1]))
  #m = randg( distr_param(am+r,bm+har(d[0])));
  #//update $\rho$
  #rho <-  rbeta(1,arho+sum(c2),brho+d[2]*r-sum(c2))
  #//update $\gamma$
  #probs <- c(psi_1+sum(c3==-1),psi0+sum(c3==0),psi1+sum(c3==1))
  #gamma <- rdirichlet(1,probs)
  #//update $b_1,b_2$
  #li5 <- update_b(d,mub,sigmab,b1,b2,c1,c2,c3,lambda1,lambda2,z);
  
  li5 <- update_bunit(d,mub,sigmab,b1,b2,c1,c2,c3,lambda1,lambda2,z)
  b1 <- li5[[1]]
  b2 <- li5[[2]]
  #up6 = update_norm(z,y,sigmamu,mu,v1,v2);
  sigma2 <- update_sigma2(z,y,mu1,mu2,mu3,v1,v2)
  mu1 <- update_mu1(z,y,sigmamu,sigma2,mu1,mu2,mu3)
  #mu2 <- update_mu2(z,y,sigmamu,sigma2,mu1,mu2,mu3)
  mu3 <- update_mu3(z,y,sigmamu,sigma2,mu1,mu2,mu3)
  #result_multi[[it]] <- rbind(c1,c2,c3)
  c1track[[it]] <- c1
  c2track[[it]] <- c2
  c3track[[it]] <- c3
  #b1track[[it]] <- b1
  #b2track[[it]] <- b2
  pb$tick()
  Sys.sleep(1/100)
}
##########################################
kclass <- vector()
for(i in 1:length(c1track)){
  kclass[i] <- ncol(c1track[[i]])
}
plot(kclass)
table(kclass)

kclass <- vector()
for(i in 5001:10000){
  kclass[i] <- ncol(c1track[[i]])
}
plot(kclass)
table(kclass)
##########################################

kl <- 5
c1left <- c1track[which(kclass==kl)]
c1left <- c1left[10*(1:as.integer(length(c1left)/10))]
min_ham = function(A, B) {
  ## list all permutations
  permutation = permutations(ncol(A))
  
  DHamming = rep(NA, nrow(permutation))
  for(l in 1 : nrow(permutation)) {
    NB = as.vector(B[, permutation[l, ]])
    DHamming[l] = sum(abs(as.vector(A) - NB))
  }
  
  return(min(DHamming))
}
library(e1071)
num <- length(c1left)
c1m <- matrix(nrow = num,ncol = num)

for(i in 1:num){
  for(j in 1:i)
    c1m[i,j] <- min_ham(c1left[[i]],c1left[[j]])
}
c1m[upper.tri(c1m)] <- t(c1m[lower.tri(c1m)])
c1s <- colSums(c1m)
c1p <- c1left[[which.min(c1s)]]
##########################################

r <-kl
c1n <- c1p
set.seed(3)
{
  
  c2ntrack <- list()
  c3ntrack <- list()
  c2n =matrix(as.numeric(runif(d[2] * r) < 0.5), d[2], r)
  
  c3n = matrix(as.numeric(runif(d[3] * r) < 0.5), d[3], r)
  #c3[c3 == 1] = 2 * (runif(sum(c3 == 1)) < 0.5) - 1
  #c1 <- c1o
  #c2 <- c2o
  #c3 <- c3o
  #z <- array(0,dim = d)
  #mat lambda1(d[2],r,fill::randn);
  #mat lambda2(d[2],r,fill::randn);
  #lambda1 = matrix(lambda1[,1],ncol = 1)
  #lambda2 = matrix(lambda2[,1],ncol = 1)
  #lambda1 <- matrix(rgamma(d[3]*r, al, 1/bl), d[3],r)
  #lambda2 <- matrix(rgamma(d[3]*r, al, 1/bl), d[3],r)
  lambda1n <- matrix(rgamma(d[3]*r,1,1/3),nrow = d[3],byrow = TRUE)
  lambda2n <- lambda1n
  zn <- z
  b1n <- b1
  b2n <- b2
  v1n <- v1
  v2n <- v2
  mu1n <- mu1
  mu2n <- mu2
  mu3n <- mu3
  sigma2n <- sigma2
}
pb <- progress_bar$new(format = "  complete [:bar] :percent eta: :eta",
                       total = 1000, clear = FALSE, width= 60)
set.seed(4)
for(it in 1:1000){
  #ztrack[[it]] <- z
  
  
  zn <- update_z(c1n, c2n, c3n, lambda1n, lambda2n, d, zn, b1n, b2n, v1n, v2n, av, bv, y, mu1n, mu2n, mu3n, sigma2n);
  v1n <- update_v1(zn,y,mu1n,mu2n,mu3n,v1n,v2n,d,av,bv,sigma2n)
  v2n <- update_v2(zn,y,mu1n,mu2n,mu3n,v1n,v2n,d,av,bv,sigma2n)
  #li1 <- update_c1(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma, m,g,rho,al,bl)
  #c1 <- li1[[1]]
  #c2 <- li1[[2]]
  #c3 <- li1[[3]]
  #lambda1 <- li1[[4]]
  #lambda2 <- li1[[5]]
  #r <- ncol(c1)
  
  c2n <- update_c2(c1n,c2n,c3n,lambda1n,lambda2n,d,zn,b1n,b2n,rho);
  
  li2 <-  update_c3(c1n,c2n,c3n,lambda1n,lambda2n,d,zn,b1n,b2n,gamma,sigmal,al,bl)
  c3n <-li2[[1]]
  lambda1n <- li2[[2]]
  lambda2n <- li2[[3]]
  
  #m <- rgamma(1,am+r,bm+har(d[1]))
  #m = randg( distr_param(am+r,bm+har(d[0])));
  #//update $\rho$
  #rho <-  rbeta(1,arho+sum(c2),brho+d[2]*r-sum(c2))
  #//update $\gamma$
  #probs <- c(psi_1+sum(c3==-1),psi0+sum(c3==0),psi1+sum(c3==1))
  #gamma <- rdirichlet(1,probs)
  #//update $b_1,b_2$
  #li5 <- update_b(d,mub,sigmab,b1,b2,c1,c2,c3,lambda1,lambda2,z);
  
  li5 <- update_bunit(d,mub,sigmab,b1n,b2n,c1n,c2n,c3n,lambda1n,lambda2n,zn)
  b1n <- li5[[1]]
  b2n <- li5[[2]]
  #up6 = update_norm(z,y,sigmamu,mu,v1,v2);
  sigma2n <- update_sigma2(zn,y,mu1n,mu2n,mu3n,v1n,v2n)
  #mu1 <- update_mu1(z,y,sigmamu,sigma2,mu1,mu2,mu3)
  mu2n <- update_mu2(zn,y,sigmamu,sigma2n,mu1n,mu2n,mu3n)
  mu3n <- update_mu3(zn,y,sigmamu,sigma2n,mu1n,mu2n,mu3n)
  #result_multi[[it]] <- rbind(c1,c2,c3)
  #c1track[[it]] <- c1
  c2ntrack[[it]] <- c2n
  c3ntrack[[it]] <- c3n
  #b1track[[it]] <- b1
  #b2track[[it]] <- b2
  pb$tick()
  Sys.sleep(1/100)
}
save.image("multi.RData")


library(reshape2)
library(ggplot2)
#pc2o = melt(t(c1o), varnames = c("Factor", "secondorder"), value.name = "Indicator")
pc1p = melt(t(c1), varnames = c("Factor", "secondorder"), value.name = "Indicator")
pc2p = melt(t(c2), varnames = c("Factor", "secondorder"), value.name = "Indicator")
pc3p = melt(t(c3), varnames = c("Factor", "thirdorder"), value.name = "Indicator")
#ggplot2::ggplot(pc2o, aes(x = secondorder, y = Factor)) + geom_tile(aes(fill = Indicator), color = "black") +
# scale_fill_gradient(low = "black", high = "green") + theme(legend.position = "none") +
#scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + labs(x = "", y = "")
ggplot2::ggplot(pc1p, aes(x = secondorder, y = Factor)) + geom_tile(aes(fill = Indicator), color = "black") +
  scale_fill_gradient(low = "black", high = "green") + theme(legend.position = "none",axis.text.x = element_text(size = 5, angle = 45)) +
  scale_x_discrete(expand = c(0, 0),guide = guide_axis(n.dodge=10)) + scale_y_discrete(expand = c(0, 0)) + labs(x = "", y = "")
ggplot2::ggplot(pc2p, aes(x = secondorder, y = Factor)) + geom_tile(aes(fill = Indicator), color = "black") +
  scale_fill_gradient(low = "black", high = "green") + theme(legend.position = "none",axis.text.x = element_text(size = 5, angle = 45)) +
  scale_x_discrete(expand = c(0, 0),guide = guide_axis(n.dodge=10)) + scale_y_discrete(expand = c(0, 0)) + labs(x = "", y = "")
ggplot2::ggplot(pc3p, aes(x = thirdorder, y = Factor)) + geom_tile(aes(fill = Indicator), color = "black") +
  scale_fill_gradientn(colors = c("red", "black", "green"), breaks = c(-1, 0, 1)) + theme(legend.position = "none") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + labs(x = "", y = "")
