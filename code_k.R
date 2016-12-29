################################################
# Variational bayes for mixture gaussians
################################################

# Packages
library(dplyr)
library(ggplot2)
library(ggthemes)
library(animation)

# Aux function, normal mixing coefficients
dnorm2<- function(x, mu, sigma, pipi) {
  pipi*dnorm(x, mu, sigma)
}

################################################
# Create data set from 2 gaussians
################################################
set.seed(1234)
data_points=10000
p1=0.3
p2=1-p1
mu1=3
mu2=8
sd1=1.5
sd2=1
unif=sample(1:2,prob=c(p1,p2),size=data_points,replace=TRUE)
mu_12=c(mu1,mu2)
sd_12=c(sd1,sd2)
x_n=rnorm(n=data_points,mean=mu_12[unif],sd=sd_12[unif])

data.frame(x_n) %>% 
  ggplot()+geom_histogram(aes(x=x_n,..density..),bins = 50,col="black", fill="white" )+
  xlim(-1,12)+theme_tufte()

################################################
# algorithm
################################################

#### initial values
niter=1400
k_gauss=4
# for dirichlet q(pi)
alpha0=0.001
alpha=rep(NA,k_gauss)
# for gaussian-gamma q(u,g)
m0=0
b0=0.1
v0=0.1
w0=0.1
m=rep(NA,k_gauss)
b=rep(NA,k_gauss)
v=rep(NA,k_gauss)
w=rep(NA,k_gauss)
# responsabilities
rho_n=matrix(NA,ncol=k_gauss,nrow=data_points)
rn=matrix(NA,ncol=k_gauss,nrow=data_points)
for(i in 1:data_points){
  set.seed(1234%%i)
  aux=runif(k_gauss)
  for(k in 1:k_gauss){
    rn[i,k]=aux[k]/sum(aux)
  }
}

# plots
niter_plot=c(1:10,10*(2:(niter/10)))
n_plot=length(niter_plot)
mean=matrix(NA,nrow=n_plot,ncol=k_gauss)
sigma=matrix(NA,nrow=n_plot,ncol=k_gauss)
pipi=matrix(NA,nrow=n_plot,ncol=k_gauss)
plot_n=1

for(i in 1:niter){
  for(k in 1:k_gauss){
    # for dirichlet q(pi)
    alpha[k]=alpha0+sum(rn[,k])
    # for gaussian-gamma q(u,g)
    b[k]=b0+sum(rn[,k])
    m[k]=(m0*b0+sum(rn[,k]*x_n))/(b0+sum(rn[,k]))
    v[k]=v0+b0*m0^2/2+sum(rn[,k]*x_n^2)/2-(m0*b0+sum(rn[,k]*x_n))^2/(2*(b0+sum(rn[,k])))
    w[k]=w0+sum(rn[,k])/2
  }

    # for q(z)
  for(k in 1:k_gauss){
    for(j in 1:data_points){
     rho_n[j,k]=exp((digamma(w[k])-log(v[k])-log(2*pi))/2+digamma(alpha[k])-digamma(sum(alpha))
                   -(1/b[k]+m[k]^2*w[k]/v[k]-2*x_n[j]*m[k]*w[k]/v[k]+x_n[j]^2*w[k]/v[k])/2)
     
    }
  }
  for(k in 1:k_gauss){
    for(j in 1:data_points){
      rn[j,k]=rho_n[j,k]/sum(rho_n[j,])
    }

  }
  if(i %in% niter_plot){
   # Expected mean gaussian 1
   mean[plot_n,]=m
   # Expected sd gaussian 1
   sigma[plot_n,]=sqrt(v/w)
   # pi
   pipi[plot_n,]=round((alpha0+colSums(rn))/(k_gauss*alpha0+data_points),2)
   plot_n=plot_n+1
  }
  
}

################################################
# Create final plot
################################################
col=rainbow(k_gauss)

saveGIF({
for(ite in 1:n_plot){
  p=data.frame(x_n) %>% 
    ggplot()+geom_histogram(aes(x=x_n,..density..),bins = 50,col="black", fill="white" )+
    xlim(-1,12)+theme_tufte()+annotate("text", x = 0, y = .28, label = paste("Iteration",niter_plot[ite]))
  
  for(t in 1:k_gauss){
    p=p+stat_function(geom = "line", fun = dnorm2,
                      args = list(mu=mean[ite,t], 
                                  sigma=sigma[ite,t], 
                                  pipi = pipi[ite,t]),
                      colour = col[t], lwd = 1)
  }  
  print(p)
}},
interval = 0.08, movie.name = "ggplot2-variational_bayes.gif", ani.width = 600, ani.height = 600)



