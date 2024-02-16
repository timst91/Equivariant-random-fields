
grad_standard=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2){
   K=cov_mat(xtr,xtr,l1_2^.5,l2_2^.5,sigma1_2^.5,sigma2_2^.5)+
    diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  K_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    if(i%%100==0){print(i)}
    for (j in 1:(length(as.matrix(xtr))/2)){
      n1=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
      n2=ifelse(length(as.matrix(xtr))==2,1,nrow(xtr))
      
      if(n1==1){
        dist_mat=dist(t(xtr),xtr)
      }else{
        dist_mat=dist(xtr,xtr)
      }
      del_cov1=matrix(0,nrow=2*n1,ncol=2*n2)
      del_cov1[1:n1,1:n2]=sigma1_2*exp(-.5*dist_mat^2/(l1_2))*dist_mat^2/(4*l1_2^4)
      del_cov1[(n1+1):(2*n1),(n2+1):(2*n2)]=0
      
      del_cov2=matrix(0,nrow=2*n1,ncol=2*n2)
      del_cov2[1:n1,1:n2]=0
      del_cov2[(n1+1):(2*n1),(n2+1):(2*n2)]=sigma2_2*exp(-.5*dist_mat^2/(l2_2))*dist_mat^2/(4*l2_2^4)
      
      del_cov3=matrix(0,nrow=2*n1,ncol=2*n2)
      del_cov3[1:n1,1:n2]=exp(-.5*dist_mat^2/(l1_2))
      del_cov3[(n1+1):(2*n1),(n2+1):(2*n2)]=0
      
      del_cov4=matrix(0,nrow=2*n1,ncol=2*n2)
      del_cov4[1:n1,1:n2]=0
      del_cov4[(n1+1):(2*n1),(n2+1):(2*n2)]=exp(-.5*dist_mat^2/(l2_2))
      
      
      
      
    }
  }
  inv_K=solve(K)
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov1%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov1),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov3%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%del_cov3),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov2),
         
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%del_cov4%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%del_cov4),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K))
  
  return(grad)
  
}


grad_fund=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2){
  
  K=eq_cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5,fundamental=1)+
    diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  K_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:(length(as.matrix(xtr))/2)){
      r1=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[i,])^2)^.5)
      r2=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[j,])^2)^.5)
      #dist=abs(xtr[i,1]-xtr[j,1])+abs(xtr[i,2]-xtr[j,2])
      theta2=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                    atan2(xtr[j,2],xtr[j,1]))
      
      
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma1_2*exp(-.5*(r1-r2)^2/(l1_2))*(r1-r2)^2/(2*l1_2^2),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,sigma2_2*exp(-.5*(r1-r2)^2/(l2_2))*(r1-r2)^2/(2*l2_2^2)))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(exp(-.5*(r1-r2)^2/(l1_2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,exp(-.5*(r1-r2)^2/(l2_2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      
      K_l1_2[i, j] = results1[1,1]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K_l1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K_l2_2[i, j] = results2[1,1]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K_l2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K_sigma1_2[i, j] = results3[1,1]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K_sigma1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K_sigma2_2[i, j] = results4[1,1]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K_sigma2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
      
      
      
      
    }
  }
  inv_K=solve(K)
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_l1_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_sigma1_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_l2_2),
         
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_sigma2_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K))
  
  return(grad)
  
}

grad_eq=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2,maxEval=mxevl){
  
  
  K=eq_cov_mat(xtr,xtr,l1_2^.5,sigma1_2^.5,l2_2^.5,sigma2_2^.5,fundamental=0,maxEval = mxevl)+
    diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  K_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    if(i%%100==0){print(i)}
    for (j in 1:(length(as.matrix(xtr))/2)){
      if(length(as.matrix(xtr)) == 2){
        repr1 = function(theta1,theta2) {
          sigma1_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr)-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr))^2)/(l1_2))
        }
        
        repr2 = function(theta1,theta2) {
          sigma2_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr)-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr))^2)/(l2_2))
        }
      }else{
        repr1 = function(theta1,theta2) {
          sigma1_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)/(l1_2))
        }
        
        repr2 = function(theta1,theta2) {
          sigma2_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)/(l2_2))
        }
      } 
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*cos(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(4*l1_2^2)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*cos(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(4*l1_2^2)
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*sin(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(4*l1_2^2)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*sin(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(4*l1_2^2)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results1=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
                                   }
      
      if(length(as.matrix(xtr)) == 2){
        repr1 = function(theta1,theta2) {
          exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr)-
                          cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr))^2)/(l1_2))
        }
        
        
      }else{
        repr1 = function(theta1,theta2) {
          exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                          cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)/(l1_2))
        }
        
      } 
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*cos(theta1)*cos(theta2)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*cos(theta1)*sin(theta2)
        
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        -repr1(theta1,theta2)*sin(theta1)*cos(theta2)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr1(theta1,theta2)*sin(theta1)*sin(theta2)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results2=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
                                   }
      
      
      
      if(length(as.matrix(xtr)) == 2){
        
        repr2 = function(theta1,theta2) {
          sigma2_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr)-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr))^2)/(l2_2))
        }
      }else{
        
        
        repr2 = function(theta1,theta2) {
          sigma2_2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                                   cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)/(l2_2))
        }
      } 
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta1)*sin(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(4*l2_2^2)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta2)*sin(theta1)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(4*l2_2^4)
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        -repr2(theta1,theta2)*sin(theta2)*cos(theta1)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(4*l2_2^2)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta1)*cos(theta2)*
          sum((cbind(c(cos(theta1), sin(theta1)), 
                     c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                 cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%
                 as.numeric(xtr[j,]))^2)/(4*l1_2^2)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results3=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
                                   }
      
      
      if(length(as.matrix(xtr)) == 2){ 
        repr2 = function(theta1,theta2) {
          exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr)-
                          cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr))^2)/(l2_2))
        }
        
        
      }else{
        repr2 = function(theta1,theta2) {
          exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(xtr[i,])-
                          cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(xtr[j,]))^2)/(l1_2))
        }
        
      } 
      
      
      integrand1 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta1)*sin(theta2)
      }
      
      integrand2 = function(theta) {
        
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*sin(theta1)*cos(theta2)
        
      }
      
      integrand3 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta1)*sin(theta2)
      }
      
      integrand4 = function(theta) {
        theta1=theta[1]
        theta2=theta[2]
        repr2(theta1,theta2)*cos(theta1)*cos(theta2)
      }
      fun=list(integrand1,integrand2,integrand3,integrand4)
      
      results4=foreach(l=1:4,.packages = c("cubature"), .combine = cbind,
                       .export = c("l1_2", "l2_2", "sigma1_2", "sigma2_2",
                                   "repr1", "repr2", "integrand1","xtr","ytr",
                                   "integrand2", "integrand3", "integrand4"))%dopar%{
                                     adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
                                   }
      
      
      
      
      
      
      
      K_l1_2[i, j] = results1[1]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[4]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[3]
      K_l1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2]
      
      
      
      K_l2_2[i, j] = results2[1]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[4]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[3]
      K_l2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2]
      
      
      
      K_sigma1_2[i, j] = results3[1]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[4]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[3]
      K_sigma1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2]
      
      
      
      K_sigma2_2[i, j] = results4[1]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[4]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[3]
      K_sigma2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2]
      
    }
    
    
  }
  inv_K=solve(K)
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_l1_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_sigma1_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_l2_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_sigma2_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K))
  
  return(grad)
}


eq_cov_mat = function(x1, x2, l1, sigma1, l2, sigma2, maxEval = mxevl,
                      nu1=0,nu2=nu1,fundamental=FALSE) {
  
  
  
  if (fundamental){
    if(nu1==0){
      cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        if(i%%100==0){print(i)}
        theta1=ifelse(length(as.matrix(x1))==2,atan2(x1[2],x1[1]),
                      atan2(x1[i,2],x1[i,1]))
        for (j in 1:(length(as.matrix(x2))/2)){
          r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                    sum((x1[i,])^2)^.5)
          r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                    sum((x2[j,])^2)^.5)
          #dist=abs(x1[i,1]-x2[j,1])+abs(x1[i,2]-x2[j,2])
          theta2=ifelse(length(as.matrix(x2))==2,atan2(x2[2],x2[1]),
                        atan2(x2[j,2],x2[j,1]))
          
          results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
            diag(c(sigma1^2*exp(-.5*(r1-r2)^2/(l1^2)),sigma2^2*exp(-.5*(r1-r2)^2/(l2^2))))%*%
            cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
          
          cov[i, j] = results[1,1]
          cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i,
              ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[2,2]
          cov[ifelse(length(as.matrix(x1))==2,1,nrow(x1)) + i, j] = results[2,1]
          cov[i,ifelse(length(as.matrix(x2))==2,1,nrow(x2)) + j] = results[1,2]
        }
      }
      return(cov)
    }else{
      cov = matrix(0, ncol =2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        if(i%%100==0){print(i)}
        theta1=atan2(x1[i,2],x1[i,1])
        for (j in 1:(length(as.matrix(x2))/2)){
          r1=ifelse(length(as.matrix(x1))==2,sum((x1)^2)^.5,
                    sum((x1[i,])^2)^.5)
          r2=ifelse(length(as.matrix(x2))==2,sum((x2)^2)^.5,
                    sum((x2[j,])^2)^.5)
          dist=(r1-r2)
          theta2=atan2(x2[j,2],x2[j,1])
          results=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
            diag(c(matern.covariance(dist,l1,nu1,sigma1),matern.covariance(dist,l2,nu2,sigma2)))%*%
            cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
          
          cov[i, j] = results[1,1]
          cov[nrow(x1) + i, nrow(x2) + j] = results[2,2]
          cov[nrow(x1) + i, j] = results[2,1]
          cov[i, nrow(x2) + j] = results[1,2]
        }
      }
      return(cov)
      
    }
    
    
  }else{
    
    if(nu1==0){
      cov = matrix(0, ncol = 2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        if(i%%100==0){print(i)}
        for (j in 1:(length(as.matrix(x2))/2)){
          if(length(as.matrix(x1)) == 2){
            repr1 = function(theta1,theta2) {
              sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2)/(l1^2))
            }
            
            repr2 = function(theta1,theta2) {
              sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2)/(l2^2))
            }
          }else{
            repr1 = function(theta1,theta2) {
              sigma1^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2)/(l1^2))
            }
            
            repr2 = function(theta1,theta2) {
              sigma2^2*exp(-0.5*sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2)/(l2^2))
            }
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
                                      "integrand2", "integrand3", "integrand4","nu1","nu2"))%dopar%{
                                        adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
                                      }
          
          
          cov[i, j] = results[1]
          cov[nrow(x1) + i, nrow(x2) + j] = results[4]
          cov[nrow(x1) + i, j] = results[3]
          cov[i, nrow(x2) + j] = results[2]
        }
      }
      
      
      return(cov)
    }else{
      cov = matrix(0, ncol = 2 * ifelse(length(as.matrix(x2)) == 2, 1, nrow(x2)),
                   nrow = 2 * ifelse(length(as.matrix(x1)) == 2, 1, nrow(x1)))
      
      for(i in 1:(length(as.matrix(x1))/2)){
        if(i%%100==0){print(i)}
        for (j in 1:(length(as.matrix(x2))/2)){
          if(length(as.matrix(x1)) == 2){
            repr1 = function(theta1,theta2) {
              matern.covariance(sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2),
                                l1,nu1,sigma1)
            }
            repr2 = function(theta1,theta2) {
              matern.covariance(sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1)-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2))^2),
                                l2,nu2,sigma2)
            }
          }else{
            repr1 = function(theta1,theta2) {
              matern.covariance(sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2),
                                l1,nu1,sigma1)
            }
            
            repr2 = function(theta1,theta2) {
              matern.covariance(sum((cbind(c(cos(theta1), sin(theta1)), c(-sin(theta1), cos(theta1)))%*%as.numeric(x1[i,])-
                                       cbind(c(cos(theta2), sin(theta2)), c(-sin(theta2), cos(theta2)))%*%as.numeric(x2[j,]))^2),
                                l2,nu2,sigma2)
            }
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
          
          results=foreach(l=1:4,.packages = c("cubature","rSPDE"), .combine = cbind,
                          .export = c("l1", "l2", "sigma1", "sigma2",
                                      "repr1", "repr2", "integrand1",
                                      "integrand2", "integrand3", "integrand4","nu1","nu2"))%dopar%{
                                        adaptIntegrate(fun[[l]], lower=c(0,0), upper=c(2 * pi,2*pi), maxEval = maxEval)$integral
                                      }
          
          
          cov[i, j] = results[1]
          cov[nrow(x1) + i, nrow(x2) + j] = results[4]
          cov[nrow(x1) + i, j] = results[3]
          cov[i, nrow(x2) + j] = results[2]
        }
      }
      
      
      return(cov)
    }
    
  }
  
}


log_likelihood_grad=function(ytr,xtr,l1,sigma1,l2,sigma2,sigma_obs,
                             equivariant=FALSE,fundamental=FALSE,nu1=0,nu2=nu1){
  if(ncol(ytr)==2){
    ytr=c(ytr[,1],ytr[,2])
  }
  
  
  if(equivariant){
    if(fundamental){
      Ktr=eq_cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,fundamental = 1,nu1=nu1,nu2=nu2)+
        sigma_obs*diag(nrow=2*nrow(xtr))
    }else{
      Ktr=eq_cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,nu1=nu1,nu2=nu2)+
        sigma_obs*diag(nrow=2*nrow(xtr))
    }
  }else{
    Ktr=cov_mat(xtr,xtr,l1^.5,sigma1^.5,l2^.5,sigma2^.5,nu1=nu1,nu2=nu2)+
      sigma_obs*diag(nrow=2*nrow(xtr))
  }
  
  
  ll=-0.5*sum((ytr)*solve(Ktr,ytr))-0.5*log(det(Ktr))-nrow(xtr)*log(2*pi)
  ll=ifelse(is.nan(ll),10^6,-ll)
  ll=ifelse(ll>0,ll,1e6)
  print(c(ll))
  return(ll)
  
  
}

initial_par=opt_par^2

initial_par=0.5*c(1,1,1,1,0.4)

(opt_par=optim(initial_par, function(x){
  log_likelihood_grad(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1,fundamental = 0)
},gr=function(x){grad_eq(Xtr,Ytr,x[1],x[2],x[3],x[4],x[5])},control=(ifault = 2),
method = "L-BFGS-B",lower=0.005)$par)


(opt_par=optim(initial_par, function(x){
  log_likelihood_grad(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1,fundamental = 1)
},gr=function(x){grad_fund(Xtr,Ytr,x[1],x[2],x[3],x[4],x[5])},
control=list(pgtol=1e-16, maxit=10000),
method = "L-BFGS-B",lower=0.005)$par^.5)



(opt_par=nloptr(initial_par, function(x){
  log_likelihood_grad(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1,fundamental = 1)
},eval_grad_f  =function(x){grad_eq(Xtr,Ytr,x[1],x[2],x[3],x[4],x[5])}, 
lb=rep.int(0.0001,5),ub=c(30,10,30,10,2),
opts=list(algorithm="NLOPT_LD_LBFGS","xtol_rel"=1e-16,maxeval=10000))$solution^.5)



(opt_par=nloptr(10*initial_par, function(x){
  log_likelihood_grad(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1,fundamental = 1)
},eval_grad_f  =function(x){grad_eq(Xtr,Ytr,x[1],x[2],x[3],x[4],x[5])}, 
lb=rep.int(0.0001,5),ub=rep.int(35,5),
opts=list(algorithm="NLOPT_LD_LBFGS","ftol_abs"=1.0e-16))$solution)

opt_par=opt_par^.5



Adam = function(f, grad_f, initial_params, lr = 0.001, beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8, max_iter = 1000, tol = 1e-6) {
  m = rep(0, length(initial_params))  # First moment estimate
  v = rep(0, length(initial_params))  # Second moment estimate
  t = 0
  
  params = initial_params
  
  for (i in 1:max_iter) {
    t = t + 1
    gradient = grad_f(params)
    
    # Update biased first moment estimate
    m = beta1 * m + (1 - beta1) * gradient
    
    # Update biased second moment estimate
    v = beta2 * v + (1 - beta2) * (gradient^2)
    
    # Compute bias-corrected first moment estimate
    m_hat = m / (1 - beta1^t)
    
    # Compute bias-corrected second moment estimate
    v_hat = v / (1 - beta2^t)
    
    # Update parameters
    params = params - lr * m_hat / (sqrt(v_hat) + epsilon)
    print(f(params))
    print(params)
    print(gradient)
    # Check for convergence
    if (sqrt(sum(gradient^2)) < tol) {
      break
    }
  }
  
  return(params)
}


# Adam optimizer with parameter positivity constraint
Adam <- function(f, grad_f, params_init, lr = 0.001, beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8, max_iter = 1000, tol = 1e-12) {
  # Initialize parameters
  params <- log(params_init)  # Apply logarithmic transformation
  
  # Initialize moments
  m <- numeric(length(params))
  v <- numeric(length(params))
  
  # Initialize iteration counter
  iter <- 0
  
  # Main optimization loop
  while (iter < max_iter) {
    iter <- iter + 1
    
    # Compute gradient
    grad <- grad_f(exp(params))  # Transform back to unconstrained space
    
    # Update biased first moment estimate
    m <- beta1 * m + (1 - beta1) * grad
    
    # Update biased second raw moment estimate
    v <- beta2 * v + (1 - beta2) * (grad^2)
    
    # Correct bias in first moment
    m_hat <- m / (1 - beta1^iter)
    
    # Correct bias in second moment
    v_hat <- v / (1 - beta2^iter)
    
    # Update parameters
    params_new <- params - lr * m_hat / (sqrt(v_hat) + epsilon)
    
    # Apply positivity constraint
    params_new <- pmax(params_new, log(tol))  # Ensure parameters are at least tol
    
    # Check convergence
    if (max(abs(params_new - params)) < tol) {
      break
    }
    
    # Update parameters for next iteration
    params <- params_new
    print(f(exp(params)))
    print(exp(params))
    print(grad)
  }
  
  # Return optimized parameters
  return(exp(params))  # Transform back to positive space
}


initial_par=opt_par^2


opt_par=Adam( function(x){log_likelihood_grad(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],
                                              equivariant = 1,fundamental = 1)},
              function(x){grad_fund(Xtr,Ytr,x[1],x[2],x[3],x[4],x[5])},
              initial_par,max_iter = 100000)^.5

# RMSprop optimizer with parameter positivity constraint
rmsprop<- function(f, grad_f, params_init, lr = 0.001, decay_rate = 0.9, epsilon = 1e-8, max_iter = 1000, tol = 1e-6) {
  # Initialize parameters
  params <- log(params_init)  # Apply logarithmic transformation
  
  # Initialize squared gradient
  cache <- numeric(length(params))
  
  # Initialize iteration counter
  iter <- 0
  
  # Main optimization loop
  while (iter < max_iter) {
    iter <- iter + 1
    
    # Compute gradient
    grad <- grad_f(exp(params))  # Transform back to unconstrained space
    
    # Update squared gradient
    cache <- decay_rate * cache + (1 - decay_rate) * grad^2
    
    # Update parameters
    params_new <- params - lr * grad / (sqrt(cache) + epsilon)
    
    # Apply positivity constraint
    params_new <- pmax(params_new, log(tol))  # Ensure parameters are at least tol
    
    # Check convergence
    if (max(abs(params_new - params)) < tol) {
      break
    }
    
    # Update parameters for next iteration
    params <- params_new
    print(f(exp(params)))
    print(exp(params))
    print(grad)
  }
  
  # Return optimized parameters
  return(exp(params))  # Transform back to positive space
}




opt_par=rmsprop( function(x){log_likelihood_grad(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],
                                              equivariant = 1,fundamental = 1)},
              function(x){grad_fund(Xtr,Ytr,x[1],x[2],x[3],x[4],x[5])},
              opt_par,max_iter = 1000)^.5


grad_fund=function(xtr,ytr,l1_2,sigma1_2,l2_2,sigma2_2,sigma_obs_2){
  
  K=eq_cov_mat(xtr,xtr,l1_2,sigma1_2,l2_2,sigma2_2,fundamental=1)+
    diag(sigma_obs_2,nrow=2*nrow(xtr))
  
  K_l1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_l2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                  nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma1_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  K_sigma2_2 = matrix(0, ncol =2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)),
                      nrow = 2 * ifelse(length(as.matrix(xtr)) == 2, 1, nrow(xtr)))
  
  
  for(i in 1:(length(as.matrix(xtr))/2)){
    if(i%%100==0){print(i)}
    theta1=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                  atan2(xtr[i,2],xtr[i,1]))
    for (j in 1:(length(as.matrix(xtr))/2)){
      r1=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[i,])^2)^.5)
      r2=ifelse(length(as.matrix(xtr))==2,sum((xtr)^2)^.5,
                sum((xtr[j,])^2)^.5)
      #dist=abs(xtr[i,1]-xtr[j,1])+abs(xtr[i,2]-xtr[j,2])
      theta2=ifelse(length(as.matrix(xtr))==2,atan2(xtr[2],xtr[1]),
                    atan2(xtr[j,2],xtr[j,1]))
      
      
      results1=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(sigma1_2^2*exp(-.5*(r1-r2)^2/(l1_2^2))*(r1-r2)^2/(l1_2^3),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results2=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,sigma2_2^2*exp(-.5*(r1-r2)^2/(l2_2^2))*(r1-r2)^2/(l2_2^3)))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results3=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(2*sigma1_2*exp(-.5*(r1-r2)^2/(l1_2^2)),0))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      results4=cbind(c(cos(theta1),sin(theta1)),c(-sin(theta1),cos(theta1)))%*%
        diag(c(0,2*sigma2_2*exp(-.5*(r1-r2)^2/(l2_2^2))))%*%
        cbind(c(cos(theta2),-sin(theta2)),c(sin(theta2),cos(theta2)))
      
      
      
      
      K_l1_2[i, j] = results1[1,1]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[2,2]
      K_l1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results1[2,1]
      K_l1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results1[1,2]
      
      
      
      K_l2_2[i, j] = results2[1,1]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
             ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[2,2]
      K_l2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results2[2,1]
      K_l2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results2[1,2]
      
      
      
      K_sigma1_2[i, j] = results3[1,1]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[2,2]
      K_sigma1_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results3[2,1]
      K_sigma1_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results3[1,2]
      
      
      
      K_sigma2_2[i, j] = results4[1,1]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i,
                 ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[2,2]
      K_sigma2_2[ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + i, j] = results4[2,1]
      K_sigma2_2[i,ifelse(length(as.matrix(xtr))==2,1,nrow(xtr)) + j] = results4[1,2]
      
      
      
      
    }
  }
  inv_K=solve(K)
  
  grad=c(-t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_l1_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma1_2%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(inv_K%*%K_sigma1_2),
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_l2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_l2_2),
         
         -t(c(ytr[,1],ytr[,2]))%*%inv_K%*%K_sigma2_2%*%inv_K%*%(c(ytr[,1],ytr[,2]))+Trace(inv_K%*%K_sigma2_2),
         -t(2*sigma_obs_2*c(ytr[,1],ytr[,2]))%*%inv_K%*%inv_K%*%c(ytr[,1],ytr[,2])+Trace(2*sigma_obs_2*inv_K))
  
  return(grad)
  
}

