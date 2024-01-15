
bd=0.4
nx=seq(-2,2,l=20)
ny=nx

Xte2=expand.grid(nx,ny)

Yte2=cbind(Xte[,1]/(bd+(Xte2[,1]^2+Xte2[,2]^2)^2),Xte2[,2]/(bd+(Xte2[,1]^2+Xte2[,2]^2)^2))


#Xtr=cbind(runif(10,-2,2),runif(10,-2,2))

Xtr2 = matrix(c(
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

Ytr2=cbind(Xtr2[,1]/(bd+(Xtr2[,1]^2+Xtr2[,2]^2)^2),Xtr2[,2]/(bd+(Xtr2[,1]^2+Xtr2[,2]^2)^2)) 


sigma_obs=0.1
noise=rnorm(length(Xtr2),0,sigma_obs)
Ytr2=Ytr2+cbind(noise[1:(length(Xtr2)/2)],noise[(length(Xtr2)/2+1):length(Xtr2)])



(opt_par2=nloptr(initial_par, function(x){
  log_likelihood(Ytr2,Xtr2,x[1],x[2],x[3],x[4],x[5],equivariant = 0)
},lb=rep.int(0.05,5),ub=rep.int(3,5),
opts=list(algorithm = "NLOPT_GN_ISRES","xtol_rel"=1.0e-12))$solution)


Ktetr2=cov_mat(Xte2,Xtr2,opt_par2[1],
              opt_par2[2],
              opt_par2[3],
              opt_par2[4])

Ktetr_trtr_inv2=Ktetr2%*%solve(cov_mat(Xtr2,Xtr2,
                                     opt_par2[1],
                                     opt_par2[2],
                                     opt_par2[3],
                                     opt_par2[4])+
                               diag(opt_par2[5],nrow=2*nrow(Xtr2)))

posterior_mean2=Ktetr_trtr_inv2%*%c(Ytr2[,1],Ytr2[,2])



posterior_cov2= cov_mat(Xte2,Xte2,opt_par2[1],
                       opt_par2[2],
                       opt_par2[3],
                       opt_par2[4])-Ktetr_trtr_inv2%*%t(Ktetr2)



sim2=posterior_mean2+chol(posterior_cov2+1e-12*diag(length(posterior_mean2)))%*%rnorm(length(posterior_mean2))
sim2=cbind(sim2[1:(0.5*length(posterior_mean2))],
          sim2[(0.5*length(posterior_mean2)+1):(length(posterior_mean2))])

(opt_par2eq=nloptr(initial_par, function(x){
  log_likelihood(Ytr2,Xtr2,x[1],x[2],x[3],x[4],x[5],equivariant = 1)
},lb=rep.int(0.05,5),ub=rep.int(3,5),
opts=list(algorithm = "NLOPT_GN_ISRES","xtol_rel"=1.0e-12))$solution)



Ktetr_eq2=eq_cov_mat(Xte2,Xtr2,opt_par2eq[1],
                    opt_par2eq[2],
                    opt_par2eq[3],
                    opt_par2eq[4])

Ktetr_trtr_inv_eq2=Ktetr_eq2%*%solve(eq_cov_mat(Xtr2,Xtr2,
                                              opt_par2eq[1],
                                              opt_par2eq[2],
                                              opt_par2eq[3],
                                              opt_par2eq[4])+
                                     diag(opt_par2eq[5],nrow=2*nrow(Xtr2)))

posterior_mean_eq2=Ktetr_trtr_inv_eq2%*%c(Ytr2[,1],Ytr2[,2])

posterior_cov_eq2= eq_cov_mat(Xte2,Xte2,opt_par2eq[1],
                             opt_par2eq[2],
                             opt_par2eq[3],
                             opt_par2eq[4])-Ktetr_trtr_inv_eq2%*%t(Ktetr_eq2)



sim_eq2=posterior_mean_eq2+chol(posterior_cov_eq2)%*%rnorm(length(posterior_mean_eq2))
sim_eq2=cbind(sim_eq2[1:(0.5*length(posterior_mean2))],
             sim[(0.5*length(posterior_mean2)+1):(length(posterior_mean2))])


par(mfrow=c(1,2))
plot(rbind(Xtr2,as.matrix(Xte2)),xaxt="n",yaxt="n")
points(Xtr2,col=2,pch=19)
quiver(Xtr2[,1], Xtr2[,2], Ytr2[,1 ],Ytr2[,2],scale=.5, col = 2)

data=sim2
quiver(Xte2[,1], Xte2[,2], data[,1 ],data[,2],scale=.5, col = 4)

data=sim_eq2

plot(rbind(Xtr2,as.matrix(Xte2)),xaxt="n",yaxt="n")
points(Xtr2,col=2, pch=19)
quiver(Xtr2[,1], Xtr2[,2], Ytr2[,1 ],Ytr2[,2],scale=.5, col = 2)
quiver(Xte2[,1], Xte2[,2], data[,1 ],data[,2],scale=.5, col = 4)
