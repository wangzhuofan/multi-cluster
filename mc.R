count_scale <- read.csv("C://ISBD/bicluster code/0228/count_scale.csv")
count_scale <- as.matrix(count_scale[,-1])

kr <- read.csv("C://ISBD/bicluster code/0228/kr.csv")
kr <- kr$x
kr <- na.omit(kr)
kc <- read.csv("C://ISBD/bicluster code/0228/kc.csv")
kc <- kc$x

brain_mat <- count_scale[kr,kc]

brain_mat <- read.csv("C://ISBD/bicluster code/0228/brain_mat.csv")
brain_mat <- as.matrix(brain_mat[,-1])
#brain_mat <- brain_mat[gcd,]

rm(count_scale)

id <- read.table("C://ISBD/bicluster code/id.csv")
id <- id$x
id <- as.character(id)
brain_id <- id[kc]
br_split <- strsplit(brain_id,"-")
don_br <- unlist(lapply(br_split,function(x){return(x[2])}))
don_set <- names(table(don_br))
att <- read.csv("C://ISBD/bicluster code/0228/att.csv")
att <- att[,-1]
library(dplyr)
brain <- filter(att,SMTS=="Brain")

att_br <- brain[(match(brain_id,brain$SAMPID)),]
tis_br <- as.character(att_br$SMTSD)
tis_set <- names(table(tis_br))
adj <- function(x){
  if(length(x)==0) 
    x <- NA 
  if(length(x) >0)  
    x <- min(x)
  return(x)
}
ind_split <- split(1:1259,list(as.factor(don_br),as.factor(tis_br)))
ind_adj <- lapply(ind_split, adj)
ind_vec <- as.numeric(unlist(ind_adj))
ind_mat <- matrix(ind_vec,nrow = 193)
sub_att <- read.csv("C://ISBD/bicluster code/0228/sub_att.csv")
sub_att <- sub_att[,-1]
subsplit  <- strsplit(as.character(sub_att$SUBJID),"-")
subid <- unlist(lapply(subsplit,function(x){return(x[2])}))

subbr_att <- sub_att[match(don_set,subid),]
subbr_att$AGE <- as.character(subbr_att$AGE)

subbr_att$AGE[subbr_att$AGE=="20-29"] <- 1

subbr_att$AGE[subbr_att$AGE=="30-39"] <- 2
subbr_att$AGE[subbr_att$AGE=="40-49"] <- 3

subbr_att$AGE[subbr_att$AGE=="50-59"] <- 4
subbr_att$AGE[subbr_att$AGE=="60-69"] <- 5
subbr_att$AGE[subbr_att$AGE=="70-79"] <- 6

subbr_att <- subbr_att[,-1]
subbr_att$AGE <- as.integer(subbr_att$AGE)


kg <- nrow(brain_mat)
farray <- function(x){
  if(is.na(x)==1)
    return(rep(NA,kg))
  else
    return(brain_mat[,x])
}
brain_array <- apply(ind_mat, c(1,2), farray)
brain_array <- aperm(brain_array,c(2,3,1))

library(caret)
for(i in 1:13){
  for(j in 1:nrow(brain_mat)){
    teid <- which(is.na(brain_array[,i,j]))
    if(length(teid)>0){
      trainx <- subbr_att[-teid,-3]
      trainy <- brain_array[-teid,i,j]
      testx <- subbr_att[teid,-3]
      testy <- knnregTrain(trainx, testx, trainy, k = 10, use.all = TRUE)
      brain_array[teid,i,j] <- testy
    }
    
  }
}

brain_array <- log(brain_array+1)
brain_array <- brain_array-mean(brain_array)

test <-c(1,8,18,23,26,27,29,31)
y <- brain_array[,,gcd]
del <- c(3,5,6,7,9,11,14,15,19,20,25,26,28,30)
del2 <- c(6,13,14)
y <- y[,,test]
y <- y[,,-del]
y <- y[,,-del2]
del3 <- c(4,8,11,17)
y <- y[,,-del3]
y <- y_sig
y <- y_alz
library(RcppArmadillo)
Rcpp::sourceCpp('test1.cpp')

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
  mu1 <- apply(y, 1, mean)
  mu1 <- mu1-mu1
  mu2 <-apply(y, 2, mean)
  mu2 <- mu2-mu2
  mu3 <- apply(y, 3, mean)
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
  m=3;rho=0.3;g=1;
  
  #vec gamma(3,fill::randu);
  gamma <- c(rho/2,1-rho,rho/2)
  #z <- z0
  z <- array(0,dim = d)
}
  #//MCMC update
set.seed(2)
  for(it in 5001:10000){
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
    
    #m <- rgamma(1,am+r,bm+har(d[1]))
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
    #mu1 <- update_mu1(z,y,sigmamu,sigma2,mu1,mu2,mu3)
    #mu2 <- update_mu2(z,y,sigmamu,sigma2,mu1,mu2,mu3)
    mu3 <- update_mu3(z,y,sigmamu,sigma2,mu1,mu2,mu3)
    #result_multi[[it]] <- rbind(c1,c2,c3)
    c1track[[it]] <- c1
    c2track[[it]] <- c2
    c3track[[it]] <- c3
    #b1track[[it]] <- b1
    #b2track[[it]] <- b2
  }


r <-kl
c1 <-c1left[[162]]
c1 <- c1p
{
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
}
for(it in 1:1000){
  #ztrack[[it]] <- z
  
  
  z <- update_z(c1, c2, c3, lambda1, lambda2, d, z, b1, b2, v1, v2, av, bv, y, mu1, mu2, mu3, sigma2);
  v1 <- update_v1(z,y,mu1,mu2,mu3,v1,v2,d,av,bv,sigma2)
  v2 <- update_v2(z,y,mu1,mu2,mu3,v1,v2,d,av,bv,sigma2)
  #li1 <- update_c1(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma, m,g,rho,al,bl)
  #c1 <- li1[[1]]
  #c2 <- li1[[2]]
  #c3 <- li1[[3]]
  #lambda1 <- li1[[4]]
  #lambda2 <- li1[[5]]
  #r <- ncol(c1)
  
  c2 <- update_c2(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,rho);
  
  li2 <-  update_c3(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma,sigmal,al,bl)
  c3 <-li2[[1]]
  lambda1 <- li2[[2]]
  lambda2 <- li2[[3]]
  
  #m <- rgamma(1,am+r,bm+har(d[1]))
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
  #mu1 <- update_mu1(z,y,sigmamu,sigma2,mu1,mu2,mu3)
  #mu2 <- update_mu2(z,y,sigmamu,sigma2,mu1,mu2,mu3)
  mu3 <- update_mu3(z,y,sigmamu,sigma2,mu1,mu2,mu3)
  #result_multi[[it]] <- rbind(c1,c2,c3)
  #c1track[[it]] <- c1
  c2track[[it]] <- c2
  c3track[[it]] <- c3
  #b1track[[it]] <- b1
  #b2track[[it]] <- b2
}
save.image("multi.RData")