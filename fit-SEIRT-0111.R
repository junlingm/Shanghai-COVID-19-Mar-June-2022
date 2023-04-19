library(R2jags)
library(ggplot2)

#setwd("C:/Users/lenovo/Desktop/20230111????")
data = read.csv(file="pudong_puxi_data_new.csv", header = TRUE, sep = ",")
rec1=subset(data)
plot(rec1$date,rec1$pd_case,type="l",ylim=c(0, max(rec1$px_case)))
lines(rec1$date,rec1$px_case,col="red")
#plot(rec1$date,rec1$wuzhengzhuang + rec1$wuzhengzhuang, log="y", type="l")


#fit model
model.SEIJRT <- function() {
  E0_pd ~ dcat(rep(20000, 50000))
  E0_px ~ dcat(rep(40000, 60000))
  I0_pd ~ dcat(rep(4000, 6000))
  I0_px ~ dcat(rep(4000, 6000))
  beta1_1 ~ dunif(0, 10)#???ڴ???
  beta1_2 ~ dunif(0, 5)#??֮?䴫??
  beta2_1 ~ dunif(0, beta1_1)
  sigma ~ dunif(0.1, 1)
  gamma ~ dunif(0.01, 1)
  tau1 ~ dunif(0, 1-gamma)
  tau2 ~ dunif(tau1, 1-gamma)
  t1 = 10
  t2 = 14
  t3 = 17
  N_pd = 8283100
  N_px = 16588000
  q ~ dunif(0, 1)
  S0_pd = N_pd*q
  S0_px = N_px*q
  
  incidence_pd[1] ~ dpois(beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q))
  incidence_px[1] ~ dpois(beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q))
  onset_pd[1] ~ dbin(sigma, E0_pd)
  onset_px[1] ~ dbin(sigma, E0_px)
  test_pd[1] ~ dbin(tau1, I0_pd)
  test_px[1] ~ dbin(tau1, I0_px)
  recovery_pd[1] ~ dbin(gamma, I0_pd)
  recovery_px[1] ~ dbin(gamma, I0_px)
  S_pd[1] <- S0_pd - incidence_pd[1]
  S_px[1] <- S0_px - incidence_px[1]
  E_pd[1] <- E0_pd + incidence_pd[1] - onset_pd[1]
  E_px[1] <- E0_px + incidence_px[1] - onset_px[1]
  I_pd[1] <- I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
  I_px[1] <- I0_px + onset_px[1] - test_px[1] - recovery_px[1]
  
  
  for (i in 2:9) {
    incidence_pd[i] ~ dpois(beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q))
    incidence_px[i] ~ dpois(beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q))
    onset_pd[i] ~ dbin(sigma, E_pd[i-1])
    onset_px[i] ~ dbin(sigma, E_px[i-1])
    test_pd[i] ~ dbin(tau1, I_pd[i-1])
    test_px[i] ~ dbin(tau1, I_px[i-1])
    recovery_pd[i] ~ dbin(gamma, I_pd[i-1])
    recovery_px[i] ~ dbin(gamma, I_px[i-1])
    S_pd[i] <- S_pd[i-1] - incidence_pd[i]
    S_px[i] <- S_px[i-1] - incidence_px[i]
    E_pd[i] <- E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] <- E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] <- I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] <- I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
  }
  
  
  for (i in 10:13) {
    incidence_pd[i] ~ dpois(beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q)) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd)
    incidence_px[i] ~ dpois(beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q)) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px)
    onset_pd[i] ~ dbin(sigma, E_pd[i-1])
    onset_px[i] ~ dbin(sigma, E_px[i-1])
    test_pd[i] ~ dbin(tau1, I_pd[i-1])
    test_px[i] ~ dbin(tau1, I_px[i-1])
    recovery_pd[i] ~ dbin(gamma, I_pd[i-1])
    recovery_px[i] ~ dbin(gamma, I_px[i-1])
    S_pd[i] <- S_pd[i-1] - incidence_pd[i]
    S_px[i] <- S_px[i-1] - incidence_px[i]
    E_pd[i] <- E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] <- E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] <- I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] <- I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
  }
  
  for (i in 14:16) {
    incidence_pd[i] ~ dpois(beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q)) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd)
    incidence_px[i] ~ dpois(beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q)) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px)
    onset_pd[i] ~ dbin(sigma, E_pd[i-1])
    onset_px[i] ~ dbin(sigma, E_px[i-1])
    test_pd[i] ~ dbin(tau1, I_pd[i-1])
    test_px[i] ~ dbin(tau1, I_px[i-1])
    recovery_pd[i] ~ dbin(gamma, I_pd[i-1])
    recovery_px[i] ~ dbin(gamma, I_px[i-1])
    S_pd[i] <- S_pd[i-1] - incidence_pd[i]
    S_px[i] <- S_px[i-1] - incidence_px[i]
    E_pd[i] <- E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] <- E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] <- I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] <- I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
  }
  
  for (i in 17:n) {
    incidence_pd[i] ~ dpois(beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q)) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd)
    incidence_px[i] ~ dpois(beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q)) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px)
    onset_pd[i] ~ dbin(sigma, E_pd[i-1])
    onset_px[i] ~ dbin(sigma, E_px[i-1])
    recovery_pd[i] ~ dbin(gamma, I_pd[i-1])
    recovery_px[i] ~ dbin(gamma, I_px[i-1])
    test_pd[i] ~ dbin(tau2, I_pd[i-1])
    test_px[i] ~ dbin(tau2, I_px[i-1])
    S_pd[i] <- S_pd[i-1] - incidence_pd[i]
    S_px[i] <- S_px[i-1] - incidence_px[i]
    E_pd[i] <- E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] <- E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] <- I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] <- I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
  }
  
}

fit2 = jags(data=list(n=34, test_pd = rec1$pd_case, test_px = rec1$px_case),
            parameters.to.save = c("E0_pd", "E0_px", "I0_pd", "I0_px","beta1_1","beta1_2","beta2_1",
                                   "sigma","gamma","tau1", "tau2","q"),
            inits =  function()p0,
            n.chains = 1,
            n.iter =20000,
            n.thin = 100,
            n.burnin =10000,
            model.file = model.SEIJRT)
#fit2 = autojags(fit2)

plot.fit <- function(fit, vars) {
  ns = names(vars)
  if (is.null(ns)) ns = vars
  for (i in 1:length(ns)) if (ns[[i]] == "") ns[[i]] = vars[[i]]
  data = data.frame()
  for (i in 1:length(vars)) {
    h = hist(fit$BUGSoutput$sims.list[[ns[[i]]]], 100, plot=FALSE)
    data = rbind(data, data.frame(value=h$mid, density=h$density, parameter=ns[[i]]))
  }
  ggplot(data, aes(x=value, y=density)) + theme_bw() +
    geom_line() + facet_wrap(~parameter, nrow=2, scales = "free")
}

vars=fit2$BUGSoutput$sims.list[c("beta1_1","beta1_2","beta2_1",
                                 "sigma","gamma","tau1", "q")]
plot.fit(fit2,vars)



#Comparison of the predicted daily number of new cases with the real data 
model.sim= function(p){
  N_pd = 8283100
  N_px = 16588000
  q = p$q
  S0_pd = N_pd*q
  S0_px = N_px*q
  E0_pd = p$E0_pd
  E0_px = p$E0_px
  I0_pd = p$I0_pd
  I0_px = p$I0_px
  beta1_1 = p$beta1_1
  beta1_2 = p$beta1_2
  beta2_1 = p$beta2_1
  sigma = p$sigma
  gamma = p$gamma
  tau1 = p$tau1
  tau2 = p$tau2
  S_pd = rep(1,34)
  S_px = rep(1,34)
  E_pd = rep(1,34)
  E_px = rep(1,34)
  I_pd = rep(1,34)
  I_px = rep(1,34)
  incidence_pd = rep(1, 34)
  incidence_px = rep(1, 34)
  onset_pd = rep(1, 34)
  onset_px = rep(1, 34)
  test_pd = rep(1, 34)
  test_px = rep(1, 34)
  recovery_pd = rep(1, 34)
  recovery_px = rep(1, 34)
  return = rep(1, 34)
  
  
  incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
  incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
  onset_pd[1] = sigma* E0_pd
  onset_px[1] = sigma* E0_px
  test_pd[1] = tau1*I0_pd
  test_px[1] = tau1*I0_px
  recovery_pd[1] = gamma*I0_pd
  recovery_px[1] = gamma*I0_px
  S_pd[1] = S0_pd - incidence_pd[1]
  S_px[1] = S0_px - incidence_px[1]
  E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
  E_px[1] = E0_px + incidence_px[1] - onset_px[1]
  I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
  I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
  return[1] = test_pd[1]+test_px[1]
  
  
  for (i in 2:9) {
    incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
    incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  
  for (i in 10:13) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  for (i in 14:16) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  for (i in 17:34) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau2*I_pd[i-1]
    test_px[i] = tau2*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  return (return = return)
}

data = model.sim(fit2$BUGSoutput$mean)

#p0 = list(E0_pd=1513, E0_px=5127, I0_pd=1897, I0_px=4593,
#          beta1_1=0.778, beta1_2=0.209, beta2_1=0.821, 
#          sigma=0.997, gamma=0.459, tau1=0.076, tau2=0.124,  q=0.221)
#data0 = model.sim(p0)



#plot(rec1$date,rec1$pd_case + rec1$px_case)#,ylim=c(0, max(data0)))
#plot(rec1$date,rec1$pd_case + rec1$px_case)#,ylim=c(0, max(data)))
#lines(rec1$date,data, col="red")
#lines(rec1$date,data0, col="blue")

rec1$yuce = data

ggplot(rec1) +
  geom_line(mapping = aes(x=date, y=yuce),color="red")+
  geom_point(mapping = aes(x=date, y=pd_case+px_case),size=1)+
  theme_bw() +
  labs(x="date", y="reported cases")+
  scale_x_continuous(breaks = c(1,8,15,22,29), labels=c("3.19","3.26","4.2","4.9","4.16"))




#library(mcmc)
#library(ggmcmc)
#F=ggs(as.mcmc(fit2))
#ggs_density(F)

#ggs_traceplot(F)
#ggs_autocorrelation(F)




#Comparison of the predicted daily number of new cases with the real data 
#beyond the time period for fitting our model.
model.sim.2= function(p){
  N_pd = 8283100
  N_px = 16588000
  q = p$q
  S0_pd = N_pd*q
  S0_px = N_px*q
  E0_pd = p$E0_pd
  E0_px = p$E0_px
  I0_pd = p$I0_pd
  I0_px = p$I0_px
  beta1_1 = p$beta1_1
  beta1_2 = p$beta1_2
  beta2_1 = p$beta2_1
  sigma = p$sigma
  gamma = p$gamma
  tau1 = p$tau1
  tau2 = p$tau2
  S_pd = rep(1,100)
  S_px = rep(1,100)
  E_pd = rep(1,100)
  E_px = rep(1,100)
  I_pd = rep(1,100)
  I_px = rep(1,100)
  incidence_pd = rep(1, 100)
  incidence_px = rep(1, 100)
  onset_pd = rep(1, 100)
  onset_px = rep(1, 100)
  test_pd = rep(1, 100)
  test_px = rep(1, 100)
  recovery_pd = rep(1, 100)
  recovery_px = rep(1, 100)
  return = rep(1, 100)
  
  
  incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
  incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
  onset_pd[1] = sigma* E0_pd
  onset_px[1] = sigma* E0_px
  test_pd[1] = tau1*I0_pd
  test_px[1] = tau1*I0_px
  recovery_pd[1] = gamma*I0_pd
  recovery_px[1] = gamma*I0_px
  S_pd[1] = S0_pd - incidence_pd[1]
  S_px[1] = S0_px - incidence_px[1]
  E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
  E_px[1] = E0_px + incidence_px[1] - onset_px[1]
  I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
  I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
  return[1] = test_pd[1]+test_px[1]
  
  
  for (i in 2:9) {
    incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
    incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  
  for (i in 10:13) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  for (i in 14:16) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  for (i in 17:100) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau2*I_pd[i-1]
    test_px[i] = tau2*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  
  return (return = return)
}


data_100 = read.csv(file="pudong_puxi_data_100.csv", header = TRUE, sep = ",")
rec2=subset(data_100)
data_100 = model.sim.2(fit2$BUGSoutput$mean)
rec2$yuce = data_100

ggplot(rec2) +
  geom_line(mapping = aes(x=date, y=yuce),color="blue")+
  geom_point(mapping = aes(x=date, y=pd_case+px_case),size=1)+
  theme_bw() +
  labs(x="date", y="reported cases")+
  scale_x_continuous(breaks = c(1,8,15,22,29,36,43,50,57,64,71,78,85,92,99), 
                     labels=c("3.19","3.26","4.2","4.9","4.16","4.23","4.30","5.7","5.14","5.21","5.28","6.4","6.11","6.18","6.25"))+
  geom_vline(xintercept = 34, color="red", linetype = "dashed")




#Calculation of the actual number of infections through the model
model.sim.3= function(p){
  N_pd = 8283100
  N_px = 16588000
  q = p$q
  S0_pd = N_pd*q
  S0_px = N_px*q
  E0_pd = p$E0_pd
  E0_px = p$E0_px
  I0_pd = p$I0_pd
  I0_px = p$I0_px
  beta1_1 = p$beta1_1
  beta1_2 = p$beta1_2
  beta2_1 = p$beta2_1
  sigma = p$sigma
  gamma = p$gamma
  tau1 = p$tau1
  tau2 = p$tau2
  S_pd = rep(1,100)
  S_px = rep(1,100)
  E_pd = rep(1,100)
  E_px = rep(1,100)
  I_pd = rep(1,100)
  I_px = rep(1,100)
  incidence_pd = rep(1, 100)
  incidence_px = rep(1, 100)
  onset_pd = rep(1, 100)
  onset_px = rep(1, 100)
  test_pd = rep(1, 100)
  test_px = rep(1, 100)
  recovery_pd = rep(1, 100)
  recovery_px = rep(1, 100)
  return = rep(1, 100)
  
  
  incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
  incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
  onset_pd[1] = sigma* E0_pd
  onset_px[1] = sigma* E0_px
  test_pd[1] = tau1*I0_pd
  test_px[1] = tau1*I0_px
  recovery_pd[1] = gamma*I0_pd
  recovery_px[1] = gamma*I0_px
  S_pd[1] = S0_pd - incidence_pd[1]
  S_px[1] = S0_px - incidence_px[1]
  E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
  E_px[1] = E0_px + incidence_px[1] - onset_px[1]
  I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
  I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
  return[1] = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
  
  
  for (i in 2:9) {
    incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
    incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
  }
  
  
  for (i in 10:13) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
  }
  
  for (i in 14:16) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
  }
  
  for (i in 17:100) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
    incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau2*I_pd[i-1]
    test_px[i] = tau2*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
  }
  
 
  return (return = return)
}

data_qz = model.sim.3(fit2$BUGSoutput$mean)
guji_qz = sum(data_qz)




#figure6-3D and Contour map
model.sim.4= function(p){
  N_pd = 8283100
  N_px = 16588000
  q = p$q
  S0_pd = N_pd*q
  S0_px = N_px*q
  E0_pd = p$E0_pd
  E0_px = p$E0_px
  I0_pd = p$I0_pd
  I0_px = p$I0_px
  beta1_1 = p$beta1_1
  beta1_2 = p$beta1_2
  beta2_1 = p$beta2_1
  #beta3_1 = p$beta3_1
  #beta3_2 = p$beta3_2
  sigma = p$sigma
  gamma = p$gamma
  tau1 = p$tau1
  tau2 = p$tau2
  S_pd = rep(1,100)
  S_px = rep(1,100)
  E_pd = rep(1,100)
  E_px = rep(1,100)
  I_pd = rep(1,100)
  I_px = rep(1,100)
  incidence_pd = rep(1, 100)
  incidence_px = rep(1, 100)
  onset_pd = rep(1, 100)
  onset_px = rep(1, 100)
  test_pd = rep(1, 100)
  test_px = rep(1, 100)
  recovery_pd = rep(1, 100)
  recovery_px = rep(1, 100)
  return_finall <- matrix(0,74,79)
  
  for (t1 in 1:1){
    for (t2 in 1:1){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) #+ beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) #+ beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau2*I0_pd
      test_px[1] = tau2*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in (t1+1):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 2:2){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) #+ beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) #+ beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in (t1+1):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+2):(t1+3)){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) #+ beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) #+ beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in (t1+1):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t2):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+4):(t1+4)){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) #+ beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) #+ beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in (t1+1):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+5):74){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) #+ beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) #+ beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in (t1+1):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
  }
  
  for (t1 in 2:2){
    for (t2 in 1:1){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau2*I0_pd
      test_px[1] = tau2*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 2:2){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    
    for (t2 in (t1+1):(t1+3)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in t1:(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t2):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+4):(t1+4)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+5):74){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
  }
  
  for (t1 in 3:3){
    for (t2 in 1:1){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau2*I0_pd
      test_px[1] = tau2*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 2:2){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1):(t1)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+1):(t1+3)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t2):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+4):(t1+4)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+5):74){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
  }
  
  for (t1 in 4:74){
    for (t2 in 1:1){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau2*I0_pd
      test_px[1] = tau2*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 2:2){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 3:(t1-1)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t2-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1):(t1)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+1):(t1+3)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t2):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+4):(t1+4)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:74) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+5):74){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) #+ beta1_2*(I_px[i-1]+J_px[i-1])*S_pd[i-1]/N_pd
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) #+ beta1_2*(I_pd[i-1]+J_pd[i-1])*S_px[i-1]/N_px
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
  }
  return (return = return_finall)
}


huatu_qz <- model.sim.4(fit2$BUGSoutput$mean)

huatu_qz = huatu_qz[,1:74]
huatu_qz1 = huatu_qz/1000000
persp(x=seq(1,74,length.out=nrow(huatu_qz1)), y=seq(1,74,length.out=ncol(huatu_qz1)), huatu_qz1,
      zlim=c(2.9,4.8), theta = 45, phi = 20, expand = 0.75, col = "lightblue",
      ltheta = -45, shade = 1, ticktype = "detailed",
      xlab = "Lockdown", ylab = "Blanket PCR testing", zlab = "Cases(million)")
contour(x=seq(1,74,length.out=nrow(huatu_qz1)), y=seq(1,74,length.out=ncol(huatu_qz1)),
        huatu_qz1, xlab = "Lockdown", ylab = "Blanket PCR testing")


min(huatu_qz)



#若封城后依旧存在区间传播，即beta1_2不为0，这100天的每日感染人数
#比例分别设置为0.01、0.02、0.05、0.1
model.sim.21= function(p){
  x = 0.1
  N_pd = 8283100
  N_px = 16588000
  q = p$q
  S0_pd = N_pd*q
  S0_px = N_px*q
  E0_pd = p$E0_pd
  E0_px = p$E0_px
  I0_pd = p$I0_pd
  I0_px = p$I0_px
  beta1_1 = p$beta1_1
  beta1_2 = p$beta1_2
  beta2_1 = p$beta2_1
  sigma = p$sigma
  gamma = p$gamma
  tau1 = p$tau1
  tau2 = p$tau2
  S_pd = rep(1,100)
  S_px = rep(1,100)
  E_pd = rep(1,100)
  E_px = rep(1,100)
  I_pd = rep(1,100)
  I_px = rep(1,100)
  incidence_pd = rep(1, 100)
  incidence_px = rep(1, 100)
  onset_pd = rep(1, 100)
  onset_px = rep(1, 100)
  test_pd = rep(1, 100)
  test_px = rep(1, 100)
  recovery_pd = rep(1, 100)
  recovery_px = rep(1, 100)
  return = rep(1, 100)
  
  
  incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
  incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
  onset_pd[1] = sigma* E0_pd
  onset_px[1] = sigma* E0_px
  test_pd[1] = tau1*I0_pd
  test_px[1] = tau1*I0_px
  recovery_pd[1] = gamma*I0_pd
  recovery_px[1] = gamma*I0_px
  S_pd[1] = S0_pd - incidence_pd[1]
  S_px[1] = S0_px - incidence_px[1]
  E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
  E_px[1] = E0_px + incidence_px[1] - onset_px[1]
  I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
  I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
  return[1] = test_pd[1]+test_px[1]
  
  
  for (i in 2:9) {
    incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
    incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  
  for (i in 10:13) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
    incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  for (i in 14:16) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
    incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau1*I_pd[i-1]
    test_px[i] = tau1*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  for (i in 17:100) {
    incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
    incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
    onset_pd[i] = sigma* E_pd[i-1]
    onset_px[i] = sigma* E_px[i-1]
    test_pd[i] = tau2*I_pd[i-1]
    test_px[i] = tau2*I_px[i-1]
    recovery_pd[i] = gamma*I_pd[i-1]
    recovery_px[i] = gamma*I_px[i-1]
    S_pd[i] = S_pd[i-1] - incidence_pd[i]
    S_px[i] = S_px[i-1] - incidence_px[i]
    E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
    E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
    I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
    I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
    return[i] = test_pd[i]+test_px[i]
  }
  
  
  return (return = return)
}


data_100 = read.csv(file="pudong_puxi_data_100.csv", header = TRUE, sep = ",")
rec2=subset(data_100)
data_100 = model.sim.21(fit2$BUGSoutput$mean)
rec2$yuce = data_100

ggplot(rec2) +
  geom_line(mapping = aes(x=date, y=yuce),color="blue")+
  geom_point(mapping = aes(x=date, y=pd_case+px_case),size=1)+
  theme_bw() +
  labs(x="date", y="reported cases")+
  scale_x_continuous(breaks = c(1,8,15,22,29,36,43,50,57,64,71,78,85,92,99), 
                     labels=c("3.19","3.26","4.2","4.9","4.16","4.23","4.30","5.7","5.14","5.21","5.28","6.4","6.11","6.18","6.25"))+
  geom_vline(xintercept = 34, color="red", linetype = "dashed")



#若封城后依旧存在区间传播，即beta1_2不为0，这100天的每日感染人数
#比例分别设置为0.01、0.02、0.05、0.1
model.sim.22= function(p){
  x = 0.05
  N_pd = 8283100
  N_px = 16588000
  q = p$q
  S0_pd = N_pd*q
  S0_px = N_px*q
  E0_pd = p$E0_pd
  E0_px = p$E0_px
  I0_pd = p$I0_pd
  I0_px = p$I0_px
  beta1_1 = p$beta1_1
  beta1_2 = p$beta1_2
  beta2_1 = p$beta2_1
  #beta3_1 = p$beta3_1
  #beta3_2 = p$beta3_2
  sigma = p$sigma
  gamma = p$gamma
  tau1 = p$tau1
  tau2 = p$tau2
  S_pd = rep(1,100)
  S_px = rep(1,100)
  E_pd = rep(1,100)
  E_px = rep(1,100)
  I_pd = rep(1,100)
  I_px = rep(1,100)
  incidence_pd = rep(1, 100)
  incidence_px = rep(1, 100)
  onset_pd = rep(1, 100)
  onset_px = rep(1, 100)
  test_pd = rep(1, 100)
  test_px = rep(1, 100)
  recovery_pd = rep(1, 100)
  recovery_px = rep(1, 100)
  return_finall <- matrix(0,74,79)
  
  for (t1 in 1:1){
    for (t2 in 1:1){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*x*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*x*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau2*I0_pd
      test_px[1] = tau2*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in (t1+1):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 2:2){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*x*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*x*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in (t1+1):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+2):(t1+3)){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q)  + beta1_2*x*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q)  + beta1_2*x*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in (t1+1):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t2):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+4):(t1+4)){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*x*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*x*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in (t1+1):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+5):74){
      incidence_pd[1] = beta2_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*x*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*x*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in (t1+1):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
  }
  
  for (t1 in 2:2){
    for (t2 in 1:1){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau2*I0_pd
      test_px[1] = tau2*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 2:2){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    
    for (t2 in (t1+1):(t1+3)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in t1:(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t2):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+4):(t1+4)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+5):74){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
  }
  
  for (t1 in 3:3){
    for (t2 in 1:1){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau2*I0_pd
      test_px[1] = tau2*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 2:2){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1):(t1)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+1):(t1+3)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t2):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+4):(t1+4)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+5):74){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
  }
  
  for (t1 in 4:74){
    for (t2 in 1:1){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau2*I0_pd
      test_px[1] = tau2*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 2:2){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in 3:(t1-1)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t2-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1):(t1)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+1):(t1+3)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t2):(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+4):(t1+4)){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:74) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
    
    for (t2 in (t1+5):74){
      incidence_pd[1] = beta1_1*I0_pd*S0_pd/(N_pd*q) + beta1_2*I0_px*S0_pd/(N_pd*q)
      incidence_px[1] = beta1_1*I0_px*S0_px/(N_px*q) + beta1_2*I0_pd*S0_px/(N_px*q)
      onset_pd[1] = sigma* E0_pd
      onset_px[1] = sigma* E0_px
      test_pd[1] = tau1*I0_pd
      test_px[1] = tau1*I0_px
      recovery_pd[1] = gamma*I0_pd
      recovery_px[1] = gamma*I0_px
      S_pd[1] = S0_pd - incidence_pd[1]
      S_px[1] = S0_px - incidence_px[1]
      E_pd[1] = E0_pd + incidence_pd[1] - onset_pd[1]
      E_px[1] = E0_px + incidence_px[1] - onset_px[1]
      I_pd[1] = I0_pd + onset_pd[1] - test_pd[1] - recovery_pd[1]
      I_px[1] = I0_px + onset_px[1] - test_px[1] - recovery_px[1]
      quezhen = test_pd[1]+test_px[1]+recovery_pd[1]+recovery_px[1]
      
      for (i in 2:(t1-1)) {
        incidence_pd[i] = beta1_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      
      for (i in t1:(t1+4-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta1_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in (t1+4):(t2-1)) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau1*I_pd[i-1]
        test_px[i] = tau1*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      
      for (i in t2:100) {
        incidence_pd[i] = beta2_1*I_pd[i-1]*S_pd[i-1]/(N_pd*q) + beta1_2*x*I_px[i-1]*S_pd[i-1]/(N_pd*q)
        incidence_px[i] = beta2_1*I_px[i-1]*S_px[i-1]/(N_px*q) + beta1_2*x*I_pd[i-1]*S_px[i-1]/(N_px*q)
        onset_pd[i] = sigma* E_pd[i-1]
        onset_px[i] = sigma* E_px[i-1]
        test_pd[i] = tau2*I_pd[i-1]
        test_px[i] = tau2*I_px[i-1]
        recovery_pd[i] = gamma*I_pd[i-1]
        recovery_px[i] = gamma*I_px[i-1]
        S_pd[i] = S_pd[i-1] - incidence_pd[i]
        S_px[i] = S_px[i-1] - incidence_px[i]
        E_pd[i] = E_pd[i-1] + incidence_pd[i] - onset_pd[i]
        E_px[i] = E_px[i-1] + incidence_px[i] - onset_px[i]
        I_pd[i] = I_pd[i-1] + onset_pd[i] - test_pd[i] - recovery_pd[i]
        I_px[i] = I_px[i-1] + onset_px[i] - test_px[i] - recovery_px[i]
        quezhen = quezhen+test_pd[i]+test_px[i]+recovery_pd[i]+recovery_px[i]
      }
      return_finall[t1,t2] = quezhen
    }
  }
  return (return = return_finall)
}


huatu_qz <- model.sim.22(fit2$BUGSoutput$mean)

huatu_qz = huatu_qz[,1:74]
huatu_qz1 = huatu_qz/1000000
persp(x=seq(1,74,length.out=nrow(huatu_qz1)), y=seq(1,74,length.out=ncol(huatu_qz1)), huatu_qz1,
      zlim=c(2.9,4.8), theta = 45, phi = 20, expand = 0.75, col = "lightblue",
      ltheta = -45, shade = 1, ticktype = "detailed",
      xlab = "Lockdown", ylab = "Blanket PCR testing", zlab = "Cases(million)")
contour(x=seq(1,74,length.out=nrow(huatu_qz1)), y=seq(1,74,length.out=ncol(huatu_qz1)),
        huatu_qz1, xlab = "Lockdown", ylab = "Blanket PCR testing")


