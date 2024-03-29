install.packages(c("proxy","cubature","nloptr","pracma","foreach", "doParallel"))
library(proxy)
library(cubature)
library(nloptr)
library(pracma)
library(foreach)
library(doParallel)

                 
maxEval=1000
mxevl=maxEval

eq_cov_mat = function(x1, x2, l1, sigma1, l2, sigma2, maxEval = mxevl) {
  
  cl = makeCluster(detectCores())
  registerDoParallel(cl)
  
  cov = matrix(0, ncol = 2 * nrow(x2), nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
  
  for(i in 1:(length(as.matrix(x1))/2)){
    if(i%%100==0){print(i)}
    for (j in 1:(length(as.matrix(x2))/2)){
      
      repr1 = function(theta1,theta2) {
        sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2)/l1)
      }
      
      repr2 = function(theta1,theta2) {
        sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2)/l2)
      }
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*cos(theta1)*cos(theta2)+ repr2(theta1,theta2)*sin(theta1)*sin(theta2)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*cos(theta1)*sin(theta2)+ repr2(theta1,theta2)*sin(theta1)*cos(theta2)
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*sin(theta1)*cos(theta2)+ repr2(theta1,theta2)*sin(theta2)*cos(theta1)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta1)*cos(theta2)+ repr1(theta1,theta2)*sin(theta1)*sin(theta2)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                      .export = c("l1", "l2", "sigma1", "sigma2",
                                  "repr1", "repr2", "integrand1", 
                                  "integrand2", "integrand3", "integrand4"))%dopar%{
        adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
      }
      
      
      cov[i, j] = results[1]
      cov[nrow(x1) + i, nrow(x2) + j] = results[4]
      cov[nrow(x1) + i, j] = results[3]
      cov[i, nrow(x2) + j] = results[2]
    }
  }
  
  stopImplicitCluster()
  return(cov)
}

cov_mat=function(x1,x2,l1,sigma1,l2,sigma2){
  
  n1=ifelse(length(as.matrix(x1))==2,1,nrow(x1))
  n2=ifelse(length(as.matrix(x2))==2,1,nrow(x2))
  
  if(n1==1){
    dist_mat=dist(t(x1),x2)
  }else{
    dist_mat=dist(x1,x2)
  }
  cov=matrix(0,nrow=2*n1,ncol=2*n2)
  cov[1:n1,1:n2]=sigma1^2*exp(-.5*dist_mat^2/(l1^2))
  cov[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2^2*exp(-.5*dist_mat^2/(l2^2))
  return(cov)
}

log_likelihood=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs,equivariant=FALSE){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  if(equivariant){
    
    Ktr=eq_cov_mat(xtr,xtr,l1,sigma1,l2,sigma2)+sigma_obs*diag(nrow=2*nrow(xtr))
  }else{
    Ktr=cov_mat(xtr,xtr,l1,sigma1,l2,sigma2)+sigma_obs*diag(nrow=2*nrow(xtr)) 
  }
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  print(c(ll))
  return(ll)
  
  
} 



set.seed(1234)




### Vector field 1

Xtr=cbind(
  -c(1.45,1.1,1.5,1.5,1.28,1.15,1.53,1.27)/1.6,
  1/1.2*c(1.13,1,0.88,0.32,0.02,-0.45,-0.88,-1.07)
)
sigma_obs=0.2
noise=rnorm(16,0,sigma_obs)
Ytr=cbind(-Xtr[,2],Xtr[,1])+cbind(noise[1:8],noise[9:16])



ngrid=17

ny=seq(-1,1,l=ngrid)
nx=seq(-1,1,l=ngrid)

Xte=expand.grid(nx,ny)
Yte=cbind(-Xte[,2],Xte[,1])

initial_par=0.5*c(1,1,1,1,0.4)

#better optimizer for this problem
(opt_par=nloptr(initial_par, function(x){
  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 0)
},lb=rep.int(0.05,5),ub=rep.int(3,5),
opts=list(algorithm = "NLOPT_GN_ISRES","xtol_rel"=1.0e-12))$solution)

#(opt_par=optim(initial_par, function(x){
#  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 0)
#},control=(ifault = 2),method = "L-BFGS-B",lower=0.01)$par)


posterior_mean=cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])%*%solve(cov_mat(Xtr,Xtr,
                                                   opt_par[1],
                                                   opt_par[2],
                                                   opt_par[3],
                                                   opt_par[4])+
                                             diag(opt_par[5],nrow=2*nrow(Xtr)),
                                           c(Ytr[,1],Ytr[,2]))
posterior_mean=cbind(posterior_mean[1:289],posterior_mean[290:578])



(rmse=mean(apply((posterior_mean-Yte)^2,1,sum))^.5)

(opt_par=nloptr(initial_par, function(x){
  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1)
},lb=rep.int(0.05,5),ub=rep.int(3,5),
opts=list(algorithm = "NLOPT_GN_ISRES","xtol_rel"=1.0e-12))$solution)



#(opt_par=optim(initial_par, function(x){
#  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1)
#},control=(ifault = 2),method = "L-BFGS-B",lower=0.05)$par)


###################




##############

posterior_mean_eq=eq_cov_mat(Xte,Xtr,opt_par[1],
                             opt_par[2],
                             opt_par[3],
                             opt_par[4],maxEval = maxEval)%*%solve(eq_cov_mat(Xtr,Xtr,
                                                                              opt_par[1],
                                                                              opt_par[2],
                                                                              opt_par[3],
                                                                              opt_par[4],maxEval = maxEval)+
                                                                     diag(opt_par[5],nrow=2*nrow(Xtr)),
                                                                   c(Ytr[,1],Ytr[,2]))
posterior_mean_eq=cbind(posterior_mean_eq[1:(0.5*length(posterior_mean_eq))],
                        posterior_mean_eq[(0.5*length(posterior_mean_eq)+1):(length(posterior_mean_eq))])


(rmse_eq=mean(apply((posterior_mean_eq-Yte)^2,1,sum))^.5)

par(mfrow=c(2,3))

plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")


points(Xtr,col=2,pch=19)
quiver(Xte[,1], Xte[,2], Yte[,1 ],Yte[,2],scale=.2)

quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2], col = 2,scale=.2)

plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")

data=posterior_mean
points(Xtr,col=2,pch=19)

quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2], col = 2,scale=.2)

quiver(Xte[,1], Xte[,2], data[,1 ],data[,2],col = 4,scale=.2)


data=posterior_mean_eq
plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")
points(Xtr,col=2,pch=19)

quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2], col = 2,scale=.2)

quiver(Xte[,1], Xte[,2], data[,1 ],data[,2],col = 4,scale=.2)


### Vector field 2


bd=0.4
nx=seq(-2,2,l=20)
ny=nx

Xte=expand.grid(nx,ny)

Yte=cbind(Xte[,1]/(bd+(Xte[,1]^2+Xte[,2]^2)^2),Xte[,2]/(bd+(Xte[,1]^2+Xte[,2]^2)^2))


#Xtr=cbind(runif(10,-2,2),runif(10,-2,2))

Xtr=matrix(c(
  -1.94149098, -1.4295386,
  1.13248441,  0.1770790,
  -1.64015467, -1.2153014,
  0.07675992,  1.5943220,
  -0.46293325, -0.4420009,
  -1.71979001, -0.7565169,
  -0.71742231, -1.3598853,
  0.67398159,  1.5847434,
  1.70560191, -1.3344249,
  -0.11236112,  1.6016984
), ncol = 2, byrow = TRUE)


Ytr=cbind(Xtr[,1]/(bd+(Xtr[,1]^2+Xtr[,2]^2)^2),Xtr[,2]/(bd+(Xtr[,1]^2+Xtr[,2]^2)^2)) 


sigma_obs=0.1
noise=rnorm(length(Xtr),0,sigma_obs)
Ytr=Ytr+cbind(noise[1:(length(Xtr)/2)],noise[(length(Xtr)/2+1):length(Xtr)])

(opt_par=nloptr(initial_par, function(x){
  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 0)
},lb=rep.int(0.05,5),ub=rep.int(3,5),
opts=list(algorithm = "NLOPT_GN_ISRES","xtol_rel"=1.0e-12))$solution)


#(opt_par=optim(initial_par, function(x){
#  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 0)
#},control=(ifault = 2),method = "L-BFGS-B",lower=0.05,upper = 3)$par)

posterior_mean=cov_mat(Xte,Xtr,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])%*%solve(cov_mat(Xtr,Xtr,
                                                   opt_par[1],
                                                   opt_par[2],
                                                   opt_par[3],
                                                   opt_par[4])+
                                             diag(opt_par[5],nrow=2*nrow(Xtr)),
                                           c(Ytr[,1],Ytr[,2]))
posterior_mean=cbind(posterior_mean[1:(0.5*length(posterior_mean))],
                     posterior_mean[(0.5*length(posterior_mean)+1):(length(posterior_mean))])




(rmse2=mean(apply((posterior_mean-Yte)^2,1,sum))^.5)

(opt_par=nloptr(initial_par, function(x){
  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1)
},lb=rep.int(0.05,5),ub=rep.int(3,5),
opts=list(algorithm = "NLOPT_GN_ISRES","xtol_rel"=1.0e-12))$solution)




#(opt_par=optim(initial_par, function(x){
#  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1)
#},control=(ifault = 2),method = "L-BFGS-B",lower=0.05,upper=3)$par)

posterior_mean_eq=eq_cov_mat(Xte,Xtr,opt_par[1],
                             opt_par[2],
                             opt_par[3],
                             opt_par[4],maxEval = maxEval)%*%solve(eq_cov_mat(Xtr,Xtr,
                                                                              opt_par[1],
                                                                              opt_par[2],
                                                                              opt_par[3],
                                                                              opt_par[4],maxEval = maxEval)+
                                                                     diag(opt_par[5],nrow=2*nrow(Xtr)),
                                                                   c(Ytr[,1],Ytr[,2]))
posterior_mean_eq=cbind(posterior_mean_eq[1:(0.5*length(posterior_mean_eq))],
                        posterior_mean_eq[(0.5*length(posterior_mean_eq)+1):(length(posterior_mean_eq))])



(rmse_eq2=mean(apply((posterior_mean_eq-Yte)^2,1,sum))^.5)





plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")
points(Xtr,col=2, pch=19)
quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2],scale=.5, col = 2)
quiver(Xte[,1], Xte[,2], Yte[,1 ],Yte[,2],scale=.5)

data=posterior_mean

plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")
points(Xtr,col=2, pch=19)
quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2],scale=.5, col = 2)

quiver(Xte[,1], Xte[,2], data[,1 ],data[,2],scale=.5, col = 4)

data=posterior_mean_eq

plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")
points(Xtr,col=2, pch=19)
quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2],scale=.5, col = 2)
quiver(Xte[,1], Xte[,2], data[,1 ],data[,2],scale=.5, col = 4)


