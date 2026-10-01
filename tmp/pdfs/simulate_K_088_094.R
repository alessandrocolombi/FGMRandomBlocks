source("PriorK_simulation/PriorK_sim_ab_theta.R")
B <- 20000L
set.seed(20261001L)
load_logC()
K <- integer(B)
for (b in seq_len(B)) {
  sigma <- rbeta(1,1,1)
  theta <- rgamma(1,shape=0.88,rate=0.94)-sigma
  gamma <- rbinom(40,1,eta)
  sizes <- diff(c(0L,which(gamma==1L)))
  ns <- unique(sizes)
  pmfs <- lapply(ns,function(n) prior_Kh_pmf(n,theta,sigma)$probability)
  K[b] <- sum(vapply(sizes,function(n) sample.int(n,1,prob=pmfs[[match(n,ns)]]),integer(1)))
  if (b %% 5000L == 0L) message(b, "/", B)
}
stopifnot(all(K>=9),all(K<=40))
prob <- tabulate(K,nbins=40)/B
write.csv(data.frame(k=1:40,probability=prob),"tmp/pdfs/K_088_094_pmf.csv",row.names=FALSE)
write.csv(data.frame(B=B,mean=mean(K),sd=sd(K),q025=quantile(K,.025,type=1),q975=quantile(K,.975,type=1)),"tmp/pdfs/K_088_094_summary.csv",row.names=FALSE)
print(c(mean=mean(K),sd=sd(K),mass=sum(prob)))
