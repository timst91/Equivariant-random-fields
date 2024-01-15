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



sim=posterior_mean+chol(posterior_cov+1e-12*diag(length(posterior_mean)))%*%rnorm(length(posterior_mean))
sim=cbind(sim[1:(0.5*length(posterior_mean))],
                     sim[(0.5*length(posterior_mean)+1):(length(posterior_mean))])

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



sim_eq=posterior_mean_eq+chol(posterior_cov_eq+1e-12*diag(length(posterior_mean_eq)))%*%rnorm(length(posterior_mean_eq))
sim_eq=cbind(sim_eq[1:(0.5*length(posterior_mean))],
          sim[(0.5*length(posterior_mean)+1):(length(posterior_mean))])

par(mfrow=c(1,2))
plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")
points(Xtr,col=2,pch=19)
quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2],scale=.1, col = 2)

data=sim
quiver(Xte[,1], Xte[,2], data[,1 ],data[,2],scale=.1, col = 4)

data=sim_eq

plot(rbind(Xtr,as.matrix(Xte)),xaxt="n",yaxt="n")
points(Xtr,col=2, pch=19)
quiver(Xtr[,1], Xtr[,2], Ytr[,1 ],Ytr[,2],scale=.5, col = 2)
quiver(Xte[,1], Xte[,2], data[,1 ],data[,2],scale=.5, col = 4)
