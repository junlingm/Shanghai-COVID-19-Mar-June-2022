library(R2jags)
library(ggplot2)
library(readr)
library(dplyr)

data <- read_csv("Shanghai-COVID-19-Mar-Jun-2022.csv", 
                 col_types = cols(date = col_date(format = "%Y-%m-%d"), 
                                  pd_case = col_integer(), 
                                  px_case = col_integer()))
data.fit = data %>% filter(date <= as.Date("2022-04-21"))

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

fit2 = jags(data=list(n=34, test_pd = data.fit$pd_case, test_px = data.fit$px_case),
            parameters.to.save = c("E0_pd", "E0_px", "I0_pd", "I0_px","beta1_1","beta1_2","beta2_1",
                                   "sigma","gamma","tau1", "tau2","q"),
            inits =  function()p0,
            n.chains = 1,
            n.iter =20000,
            n.thin = 100,
            n.burnin =10000,
            model.file = model.SEIJRT)
#fit2 = autojags(fit2)

vars=fit2$BUGSoutput$sims.list[c("beta1_1","beta1_2","beta2_1",
                                 "sigma","gamma","tau1", "q")]

#Comparison of the predicted daily number of new cases with the real data 
model.sim= function(date, p){
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
  n = length(date)
  S_pd = rep(1,n)
  S_px = rep(1,n)
  E_pd = rep(1,n)
  E_px = rep(1,n)
  I_pd = rep(1,n)
  I_px = rep(1,n)
  incidence_pd = rep(1, n)
  incidence_px = rep(1, n)
  onset_pd = rep(1, n)
  onset_px = rep(1, n)
  test_pd = rep(1, n)
  test_px = rep(1, n)
  recovery_pd = rep(1, n)
  recovery_px = rep(1, n)

  
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
  }
  
  for (i in 17:n) {
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
  }
  
  data.frame(
    date = date,
    test_pd = test_pd,
    test_px = test_px,
    total = test_pd + test_px
  )
}

data.sim = model.sim(data$date, fit2$BUGSoutput$mean)

data.plot = rbind(
  data.frame(
    date = data$date,
    cases = data$pd_case + data$px_case,
    curve = "data",
    panel = "total"
  ),
  data.frame(
    date = data$date,
    cases = data$pd_case,
    curve = "data",
    panel = "Pudong"
  ),
  data.frame(
    date = data$date,
    cases = data$px_case,
    curve = "data",
    panel = "Puxi"
  ),
  data.frame(
    date = data.sim$date,
    cases = data.sim$total,
    curve = "prediction",
    panel = "total"
  ),
  data.frame(
    date = data.sim$date,
    cases = data.sim$test_pd,
    curve = "prediction",
    panel = "Pudong"
  ),
  data.frame(
    date = data.sim$date,
    cases = data.sim$test_px,
    curve = "prediction",
    panel = "Puxi"
  )
)
ggplot(data=data.plot, aes(x=date, y=cases, color=curve))+
  theme_bw() + geom_line() + facet_grid(panel~.) +
  geom_line() +
  geom_vline(xintercept = as.Date("2022-04-21")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %d", date_minor_breaks = "1 week")

