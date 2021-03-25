library(sm)

# population level parameters
#Number of age classes from 1 to max(age)
A=10
age=1:A

#Known parameters of the von Bertalanffy growth model
#Linf.mu*(1-exp(-k.mu(age-t0)))
#Linf = asymptotic size
#k = rate approaching Linf
#t0 = theoretical time at age-0

Linf.mu = 650
k.mu = 0.25

# generate mean length at age
lt<-Linf.mu*(1-exp(-k.mu*(age-0)))

# log-likelihood function
ll <- function(lt,age, linf,k){
  predy<-linf*(1-exp(-k*(age-0)))
  sum(dnorm(lt, predy, 5, log=T))
}

# prior densities
plinf <- function(linf){
  dnorm(linf, 0, 100, log=T)
}

pk<-function(k){
  dunif(k,0,1,log=T)
}


# posterior density function (log scale)
post <- function(lt, age,linf,k){
  ll(lt, age, linf,k) + plinf(linf) + pk(k)
}

geninits <- function(){
  list(linf = runif(1,500,750),
       k = runif(1,0,1))
}

jump <- function(lt, dist = .2){
  lt + rnorm(1, 0, dist)
}

iter = 1000000
chains <- 3
posterior <- array(dim = c(chains, 2, iter))
accepted <- array(dim=c(chains, iter - 1))

set.seed(1234)
for (c in 1:chains){
  theta.post <- array(dim=c(2, iter))
  inits <- geninits()
  theta.post[1, 1] <- inits$linf
  theta.post[2, 1] <- inits$k
  for (t in 2:iter){
    
    theta_star <- c(jump(theta.post[1, t-1]), jump(theta.post[2, t-1]))
    pstar <- post(lt, age,linf=theta_star[1],k=theta_star[2])  
    pprev <- post(lt, age, linf=theta.post[1, t-1],k=theta.post[2,t-1])
    lr <- pstar - pprev
    r <- exp(lr)
    
    accept <- rbinom(1, 1, prob = min(r, 1))
    accepted[c, t - 1] <- accept
    if (accept == 1){
      theta.post[, t] <- theta_star
    } else {
      theta.post[, t] <- theta.post[, t-1]
    }
  }
  posterior[c, , ] <- theta.post
}


sequence <- unique(round(exp(seq(0, log(iter), length.out = 150))))

xlims <- c(400, 900)
ylims <- c(0, 1)
library(animation)
ani.record(reset = TRUE)
for (i in sequence){
  par(mfrow=c(1, 2),mar=c(4,4,2,2))
  plot(posterior[1, 1, 1:i], posterior[1, 2, 1:i],
       type="l", xlim=xlims, ylim=ylims, col="blue",
       xlab=expression(paste("L",infinity)), ylab=expression(k), main="Markov chains")
  lines(posterior[2, 1, 1:i], posterior[2, 2, 1:i],
        col="purple")
  lines(posterior[3, 1, 1:i], posterior[3, 2, 1:i],
        col="red")
  
  sm.density(x=cbind(c(posterior[, 1, 1:i]), c(posterior[, 2, 1:i])),
             xlab="Length infinity", ylab="k",
             zlab="", 
             xlim=xlims, ylim=ylims, col="white", 
             verbose=0)
  title("Posterior density")
  ani.record()
}

oopts = ani.options(interval = 0.5)
ani.replay()
setwd("")

saveGIF(ani.replay(), clean=TRUE,interval=0.07,ani.width=800,img.name = "MCMC_Linf_k_plot")


