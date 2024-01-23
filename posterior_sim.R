bd=0.4
nx=seq(-2,2,l=20)
ny=nx

Xte=expand.grid(nx,ny)

Yte=cbind(Xte[,1]/(bd+(Xte[,1]^2+Xte[,2]^2)^2),Xte[,2]/(bd+(Xte[,1]^2+Xte[,2]^2)^2))


#Xtr=cbind(runif(10,-2,2),runif(10,-2,2))

Xtr = matrix(c(
  -1.25977711, -0.967413407,
  1.03224558, -0.660621061,
  0.26712513, -1.465801326,
  1.72869428, -0.001814459,
  0.55477327,  1.208542531,
  0.80299254, -0.651387027,
  -0.08311013,  0.035682461,
  1.40124767, -0.022245753,
  -0.31067730,  1.188211616,
  -1.87443151,  0.267835582
), ncol = 2, byrow = TRUE)

Ytr=cbind(Xtr[,1]/(bd+(Xtr[,1]^2+Xtr[,2]^2)^2),Xtr[,2]/(bd+(Xtr[,1]^2+Xtr[,2]^2)^2)) 


sigma_obs=0.1
noise=rnorm(length(Xtr),0,sigma_obs)
Ytr=Ytr+cbind(noise[1:(length(Xtr)/2)],noise[(length(Xtr)/2+1):length(Xtr)])


(opt_par=nloptr(initial_par, function(x){
  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 0)
},lb=rep.int(0.05,5),ub=rep.int(3,5),
opts=list(algorithm = "NLOPT_GN_ISRES","xtol_rel"=1.0e-12))$solution)


Ktetr=cov_mat(Xte,Xtr,opt_par[1],
              opt_par[2],
              opt_par[3],
              opt_par[4])

Ktetr_trtr_inv=Ktetr%*%solve(cov_mat(Xtr,Xtr,
              opt_par[1],
              opt_par[2],
              opt_par[3],
              opt_par[4])+
        diag(opt_par[5],nrow=2*nrow(Xtr)))
        
posterior_mean=Ktetr_trtr_inv%*%c(Ytr[,1],Ytr[,2])



posterior_cov= cov_mat(Xte,Xte,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])-Ktetr_trtr_inv%*%t(Ktetr)


(opt_par=nloptr(initial_par, function(x){
  log_likelihood(Ytr,Xtr,x[1],x[2],x[3],x[4],x[5],equivariant = 1)
},lb=rep.int(0.05,5),ub=rep.int(3,5),
opts=list(algorithm = "NLOPT_GN_ISRES","xtol_rel"=1.0e-12))$solution)



Ktetr_eq=eq_cov_mat(Xte,Xtr,opt_par[1],
              opt_par[2],
              opt_par[3],
              opt_par[4])

Ktetr_trtr_inv_eq=Ktetr_eq%*%solve(eq_cov_mat(Xtr,Xtr,
                             opt_par[1],
                             opt_par[2],
                             opt_par[3],
                             opt_par[4])+
                       diag(opt_par[5],nrow=2*nrow(Xtr)))

posterior_mean_eq=Ktetr_trtr_inv_eq%*%c(Ytr[,1],Ytr[,2])

posterior_cov_eq= eq_cov_mat(Xte,Xte,opt_par[1],
                       opt_par[2],
                       opt_par[3],
                       opt_par[4])-Ktetr_trtr_inv_eq%*%t(Ktetr_eq)



svd_eq=svd(posterior_cov_eq)
sqrt_eq=svd_eq$u%*%sqrt(diag(svd_eq$d))%*%t(svd_eq$v)

sim_eq=posterior_mean_eq+sqrt_eq%*%rnorm(length(posterior_mean_eq))
sim_eq=cbind(sim_eq[1:(0.5*length(posterior_mean))],
             sim_eq[(0.5*length(posterior_mean)+1):(length(posterior_mean))])

sim=posterior_mean+chol(posterior_cov+1e-12*diag(length(posterior_mean)))%*%rnorm(length(posterior_mean))
sim=cbind(sim[1:(0.5*length(posterior_mean))],
          sim[(0.5*length(posterior_mean)+1):(length(posterior_mean))])

par(mfrow=c(1,2))
plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")
points(Xtr,col=2,pch=19)
quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2],scale=.5, col = 2)

data=sim
quiver(Xte[,1], Xte[,2], data[,1 ],data[,2],scale=.5, col = 4)

data=sim_eq

plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")
points(Xtr,col=2, pch=19)
quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2],scale=.5, col = 2)
quiver(Xte[,1], Xte[,2], data[,1 ],data[,2],scale=.5, col = 4)
