

###############################################
# All libraries for this file
#library(dplyr)
library(tidyverse)
library("corrplot")
library(rjags)
library(DescTools)


# Load data
load("rda/galaxy.rda")

##############################################
# Prepepare Data

# gal.select <- galaxy[,column.select]
# gal.select <- na.omit(gal.select)
gal.select <- galaxy %>% filter(N_GC > 0) %>%
  filter(!is.na(D)) %>% filter(!is.na(M_V_T))

# Log-Linear plot of globular cluster population
# Without error bars
# gal.select %>%
#   ggplot(aes(M_V_T, N_GC)) +
#   geom_smooth(method = "loess", span = 0.1) +
#   geom_point(aes(color=Type), alpha=0.6) +
#   scale_x_reverse() +
#   scale_y_continuous(trans='log10') +
#   labs(title="Selected Data") +
#   labs(y="Globular Clusters") +
#   labs(x="Absolute Magnitude")


# Boxplot of globular population vs magnitude
gal.select %>%
  ggplot(aes(x=Type, y=log(N_GC))) + 
  geom_boxplot() + 
  geom_point() +
  labs(y="Globular Clusters (log)") +
  labs(x="Type")
  
#boxplot(log(N_GC) ~ Type, data=gal.select)


####################################################
# Some data will be used for the model
# Others will be used to check the model
gal.rows <- sample(seq(from=1, to= length(gal.select$Type), by =1), 100, replace = FALSE)

gal.select <- gal.select %>%
  mutate(Type.num=as.numeric(Type))
#gal.select$Type.num[gal.select$Type.num==4] <- 3

test.sample <- gal.select[gal.rows,]
train.sample<- gal.select[-gal.rows,]
train.sample<- gal.select

#pairs(train.sample)


# Check types
round(prop.table(table(train.sample$Type)),2)
round(prop.table(table(test.sample$Type)),2)


####################################################
# Jags model
mod_string = " model {
  for (i in 1:n) {
    y[i] ~ dt( mu[i], tau[type[i]], df )
    mu[i] = alpha + beta[type[i]]*mvt[i] + gamma*pow(mvt[i],2) + delta*d[i]
  }

  for (j in 1:max(type)) {
    beta[j] ~ dnorm(-0.8, 1000)
    tau[j] ~ dgamma(5/2.0, 5*10/2.0) # tau is close to, but not equal to the precision
    sigt[j] = sqrt( 1.0 / tau[j] * df / (df - 2.0) ) # standard deviation of errors
  }
  
  
  gamma ~ dnorm(0.0, 10/10^6)
  df = nu + 2.0 # we want degrees of freedom > 2 to guarantee existence of mean and variance
  nu ~ dexp(1.0)

  delta~dnorm(0.0, 1.0/10^6)


  alpha ~ dnorm(-9.5, 10)

  
} "

####################################################
# Fit model


# Prepare data to JAGS
JAGS_data <- list(y = log10(train.sample$N_GC),
                  n = nrow(train.sample),
                  mvt = train.sample$M_V_T,
                  type = train.sample$Type.num,
                  d = train.sample$D)

# Parameters to monitor
params = c("alpha", "beta", "gamma", "delta", "tau", "df")

# Initialisation function
inits1 <- function() {
  inits <- list("alpha"=rnorm(1, -9.5, 1),
                "beta"=rnorm(3, -0.88, 0.5),
                "gamma"=rnorm(1, 0.1, 0.1),
                "delta"=rnorm(1, 0, 0.1))
}

# Model initialization 
mod <- jags.model(textConnection(mod_string),
                  data=JAGS_data,
                  inits = inits1,
                  n.chains=5)
# Burn in
update(mod, 5e3)

# Fit the model
mod_sim <- coda.samples(model=mod,
                        variable.names = params, 
                        n.iter = 7e3)

# Group chains
mod_csim <- as.mcmc(do.call(rbind, mod_sim))

summary(mod_csim)

####################################################
# Convergence diagnostics

plot(mod_sim)
# lmod = lm(log(N_GC) ~ M_V_T, data=train.sample)
# summary(lmod)
# plot(x=train.sample$M_V_T, y=log(train.sample$N_GC))

# Trace plot
#coda::traceplot(mod_sim)

# Density plot
par(mfrow=c(2,2))
densplot(mod_csim[,c("beta[1]", "beta[2]", "beta[3]", "gamma")])
#ggsave("figs/DensityPlotCoeff.pdf", plot = last_plot())

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
#autocorr.plot(mod_sim)
effectiveSize(mod_sim)

# Comptute DIC
(dic1 <- dic.samples(mod, n.iter = 1e3))


# #######################################################
# # Analysis of the posterior distribution
# 
# # Compute posterior distribution for each galaxy
# # Return median and standard deviation
# compute_post_mean <- function(n){
#   post.sim <- as.numeric(mod_csim[,"alpha"] +
#                            mod_csim[,1+test.sample$Type.num[n]]*test.sample$M_V_T[n] +
#                            mod_csim[,"gamma"]*test.sample$M_V_T[n]^2 +
#                            mod_csim[,"delta"]*test.sample$D[n])
#   c(median(10^(post.sim)), sd(10^(post.sim)))
# }
# 
# n <- seq(1:length(test.sample$Type.num))
# # sapply permits to perform element-wise operation on any functions
# estimate <- sapply(n, compute_post_mean)
# 
# test.sample <- test.sample %>% mutate(N_GC.est = round(estimate[1,],2))
# 
# estimate_bound <- ((100 - (100-95)/2))/100
# 
# lower <- estimate[1,]-qnorm(estimate_bound)*estimate[2,]




#######################################################
# Analysis of the posterior distribution

compute_posterior_distrib <- function(type, magnitude, distance){
  len.MC <- length(mod_csim[,"gamma"])
  mu.MC <-  as.numeric(mod_csim[,"alpha"] +
                         mod_csim[,1+type]*magnitude +
                         mod_csim[,"gamma"]*magnitude^2 +
                         mod_csim[,"delta"]*distance)
  mu.MC <- 10^mu.MC
  post.MC <- rt(len.MC, df=mod_csim[,"df"]+3, ncp=mu.MC)
}

res.MC <- data.frame(
  mean.MC = numeric(),
  sd.MC=numeric()
)

n <- seq(1:length(test.sample$Type.num))
for (pos in n) {
  post.distrib <- compute_posterior_distrib(type=test.sample$Type.num[pos], magnitude=test.sample$M_V_T[pos], distance=test.sample$D[pos])
  estimates <- c(median(post.distrib), sd(post.distrib))
  res.MC <- rbind(res.MC, data.frame(mean.MC=estimates[1], sd.MC=estimates[2]))
}

test.sample$N_GC.est <- round(res.MC$mean.MC,2)

target <- 95
estimate_bound <- ((100 - (100-target)/2))/100

lower <- res.MC$mean.MC-qnorm(estimate_bound)*res.MC$sd.MC
upper <- res.MC$mean.MC+qnorm(estimate_bound)*res.MC$sd.MC
test.sample$N_GC.est.low <- ifelse((lower<0.1), 0.1, lower)
test.sample$N_GC.est.up <- upper


test.sample_reshape <- stack(data.frame(Observed = test.sample$M_V_T, Estimated=test.sample$M_V_T))
test.sample_reshape <- test.sample_reshape %>%
  mutate(N_GC=stack(data.frame(Observed = test.sample$N_GC,
                               Estimated=test.sample$N_GC.est))$value)
test.sample_reshape <- test.sample_reshape %>%
  mutate(N_GC.obs.low=stack(data.frame(Observed = test.sample$N_GC.obs.low,
                                       Estimated=test.sample$N_GC.est.low))$value)
test.sample_reshape <- test.sample_reshape %>%
  mutate(N_GC.obs.up=stack(data.frame(Observed = test.sample$N_GC.obs.up,
                                      Estimated=test.sample$N_GC.est.up))$value)

#test.sample_reshape <- gather(test.sample, c(M_V_T, N_GC))


test.sample_reshape %>%
  ggplot(aes(x=values, y=N_GC, color = ind)) +
  geom_point(alpha=0.6) +
  geom_errorbar(aes(ymin=N_GC.obs.low, ymax=N_GC.obs.up, color = ind),
                width=0.2, position=position_dodge(0.05)) +
  scale_x_reverse() +
  scale_y_continuous(trans='log10')+
  #labs(title="Number of Global Clusters", subtitle = "with error bars") +
  labs(color = "Globular cluster population :") +
  labs(y="Globular Clusters (log scale)") +
  labs(x="Absolute Magnitude") +
  theme(legend.position = "bottom")


compute_contained <- function(n){
  c(test.sample$N_GC.est.low[n], test.sample$N_GC.est.up[n]) %overlaps% c(test.sample$N_GC.obs.low[n], test.sample$N_GC.obs.up[n])
}
n <- seq(1:length(test.sample$Type.num))
contained <- sapply(n, compute_contained)

# Overlap(x, y)
# 
# Interval(x, y)

mean(contained)
mean(sapply(n, function(n){
  between(test.sample$N_GC[n], test.sample$N_GC.est.low[n], test.sample$N_GC.est.up[n])
}))

mean(sapply(n, function(n){
  between(test.sample$N_GC.est[n], test.sample$N_GC.obs.low[n], test.sample$N_GC.obs.up[n])
}))


############################################################################•
# Compute the statistics for all the data

res.MC <- data.frame(
  mean.MC = numeric(),
  sd.MC=numeric()
)

n <- seq(1:length(train.sample$Type.num))
for (pos in n) {
  post.distrib <- compute_posterior_distrib(type=train.sample$Type.num[pos],
                                            magnitude=train.sample$M_V_T[pos],
                                            distance=train.sample$D[pos])
  estimates <- c(median(post.distrib), sd(post.distrib))
  res.MC <- rbind(res.MC, data.frame(mean.MC=estimates[1], sd.MC=estimates[2]))
}

train.sample$N_GC.est <- round(res.MC$mean.MC,2)

target <- 95
estimate_bound <- ((100 - (100-target)/2))/100

lower <- res.MC$mean.MC-qnorm(estimate_bound)*res.MC$sd.MC
upper <- res.MC$mean.MC+qnorm(estimate_bound)*res.MC$sd.MC
train.sample$N_GC.est.low <- ifelse((lower<0.1), 0.1, lower)
train.sample$N_GC.est.up <- upper


train.sample_reshape <- stack(data.frame(Observed = train.sample$M_V_T,
                                 Estimated=train.sample$M_V_T))

train.sample_reshape <- train.sample_reshape %>%
  mutate(N_GC=stack(data.frame(Observed = train.sample$N_GC,
                               Estimated=train.sample$N_GC.est))$value)
train.sample_reshape <- train.sample_reshape %>%
  mutate(N_GC.obs.low=stack(data.frame(Observed = train.sample$N_GC.obs.low,
                                       Estimated=train.sample$N_GC.est.low))$value)
train.sample_reshape <- train.sample_reshape %>%
  mutate(N_GC.obs.up=stack(data.frame(Observed = train.sample$N_GC.obs.up,
                                      Estimated=train.sample$N_GC.est.up))$value)


train.sample_reshape %>%
  ggplot(aes(x=values, y=N_GC, color = ind)) +
  geom_point(alpha=0.6) +
  geom_errorbar(aes(ymin=N_GC.obs.low, ymax=N_GC.obs.up, color = ind),
                width=0.2, position=position_dodge(0.05)) +
  scale_x_reverse() +
  scale_y_continuous(trans='log10')+
  #labs(title="Number of Global Clusters", subtitle = "with error bars") +
  labs(color = "Globular cluster population :") +
  labs(y="Globular Clusters (log scale)") +
  labs(x="Absolute Magnitude") +
  theme(legend.position = "bottom")


compute_contained <- function(n){
  c(train.sample$N_GC.est.low[n], train.sample$N_GC.est.up[n]) %overlaps%
    c(train.sample$N_GC.obs.low[n], train.sample$N_GC.obs.up[n])
}

contained <- sapply(n, compute_contained)

mean(contained)

mean(sapply(n, function(n){
  between(train.sample$N_GC[n], train.sample$N_GC.est.low[n], train.sample$N_GC.est.up[n])
}))

  
# test.sample <- test.sample %>%
#   mutate(N_GC.est.low = ifelse((lower<0.1), 0.1, lower))
# test.sample <- test.sample %>%
#   mutate(N_GC.est.up = N_GC.est + qnorm(estimate_bound)*estimate[2,])
# 
# 
#  
#   
# compute_contained <- function(n){
#   c(test.sample$N_GC.est.low[n], test.sample$N_GC.est.up[n]) %overlaps% c(test.sample$N_GC.obs.low[n], test.sample$N_GC.obs.up[n])
# }
# n <- seq(1:length(test.sample$Type.num))
# contained <- sapply(n, compute_contained)
# 
# # Overlap(x, y)
# # 
# # Interval(x, y)
# 
# mean(contained)

#######################################################
# Density plot of the posterior distribution
# given the magnitude

# Compute posterior distribution for one galaxy
# Return the distribution
compute_distrib <- function(n){
  post.sim <- as.numeric(mod_csim[,"alpha"] +
                           mod_csim[,1+test.sample$Type.num[n]]*test.sample$M_V_T[n] +
                           mod_csim[,"gamma"]*test.sample$M_V_T[n]^2 +
                           mod_csim[,"delta"]*test.sample$D[n])
  10^(post.sim)
}


post.exp.sim <- data.frame(res.MC=compute_distrib(13))
post.exp.sim %>%
  ggplot(aes(res.MC, fill="g"), color = "grey", show.legend=FALSE) +
  geom_density(alpha = 0.4, show.legend=FALSE) +
  scale_x_continuous(trans='log10') +
  #labs(title="Density ditribution of clusters") +
  labs(y="Density") +
  labs(x="Globular Clusters")
mean(post.exp.sim$res.MC)
(ci.MC <- c(-1,1)*qnorm(0.975)*sd(post.exp.sim$res.MC)+mean(post.exp.sim$res.MC))


#######################################################
# Density plot of the posterior distribution

train.sample %>%
  ggplot(aes(N_GC, fill=Type)) +
  geom_density(alpha = 0.2) +
  scale_x_continuous(trans='log10') +
  labs(title="Density of trained data") +
  labs(y="Globular Clusters") +
  labs(x="Absolute Magnitude")


test.sample %>%
  ggplot(aes(N_GC.est, fill=Type)) +
  geom_density(alpha = 0.2) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(limits=c(0, 1.5)) +
  #labs(title="Density of estimated clusters") +
  labs(y="Globular Clusters") +
  labs(x="Absolute Magnitude")

test.sample %>%
  ggplot(aes(N_GC, fill=Type)) +
  geom_density(alpha = 0.2) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(limits=c(0, 1.5)) +
  #labs(title="Density of tested clusters", subtitle="Observed") +
  labs(y="Globular Clusters") +
  labs(x="Absolute Magnitude")


####################################################
# Linear regression

train.sample %>%
  ggplot(aes(M_V_T, log10(N_GC))) +
  geom_point(aes(color=Type), alpha=0.6) +
  scale_x_reverse() +
  scale_y_continuous() +
  labs(title="Selected Data") +
  labs(y="Globular Clusters") +
  labs(x="Absolute Magnitude")

# point estimates for the coefficient from the combined simulation
# Take the posteror median instead of the posterior mean
# apply to the column (number 2)
(pmed_coef <- apply(mod_csim, 2, median))



(prec<- 1/ as.numeric(pmed_coef["sig"])^2)


I1 <- test.sample$Type=="E"
I2 <- test.sample$Type=="I"
I3 <- test.sample$Type=="S"
#I4 <- test.sample$Type=="S0"


lX_hat <- pmed_coef["alpha"] + 
  I1*pmed_coef["beta[1]"]*test.sample$M_V_T +
  I2*pmed_coef["beta[2]"]*test.sample$M_V_T +
  I3*pmed_coef["beta[3]"]*test.sample$M_V_T +
#  I4*pmed_coef["beta[3]"]*test.sample$M_V_T +
  I1*pmed_coef["gamma"]*test.sample$M_V_T^2 +
  I2*pmed_coef["gamma"]*test.sample$M_V_T^2 +
  I3*pmed_coef["gamma"]*test.sample$M_V_T^2 +
  #I4*pmed_coef["gamma"]*test.sample$M_V_T^2 +
  test.sample$D*pmed_coef["delta"]

# lX_hat <- pmed_coef["alpha"] + 
#   pmed_coef["beta[1]"]*test.sample$M_V_T +
#   pmed_coef["c[1]"]*test.sample$M_V_T^2 +
#   test.sample$D*pmed_coef["delta"]

# inverse of the link function

X_hat <- data.frame(N_GC=10^(lX_hat))
test.sample <- test.sample %>%
  mutate(N_GCe = round(10^(lX_hat)))


# Calculate the residual
test.sample <- test.sample %>%
  mutate(resid = N_GC-N_GCe)

test.sample %>%
  ggplot(aes(M_V_T, resid)) +
  geom_point(aes(color=Type), alpha=0.8) +
  scale_x_reverse() +
  scale_y_continuous() +
  labs(title="Residuals") +
  labs(y="Globular Clusters") +
  labs(x="Absolute Magnitude")


#test.sample[test.sample=="S0"] <- "S"

test.sample %>%
  ggplot(aes(M_V_T, N_GCe)) +
  geom_smooth(method = "loess", span = 0.3) +
  geom_point(aes(M_V_T, N_GC, color=Type), alpha=0.6, show.legend=FALSE) +
  geom_errorbar(aes(ymin=N_GC-N_GC_u, ymax=N_GC+N_GC_u, 
                    color=Type), width=.2, position=position_dodge(0.05), show.legend=FALSE) +
  geom_point(aes(M_V_T, N_GCe), alpha=0.6) +
  scale_x_reverse() +
  scale_y_continuous(trans='log10') +
  #labs(title="Estimated Global Clusters", subtitle="black: estimated, colored: observed") +
  labs(y="Globular Clusters (log)") +
  labs(x="Absolute Magnitude")


########################################################
# Monte Carlo simulation of the posterior predictive
# Plot of the distribution for one galaxy

# Galaxy to plot
galaxy.num <- 33

# Gather posterior distribution
post.distrib <- compute_posterior_distrib(type=test.sample$Type.num[galaxy.num],
                                  magnitude=test.sample$M_V_T[galaxy.num],
                                  distance=test.sample$D[galaxy.num])

densplot(post.distrib)

post.distrib.df <- data.frame(distrib=post.distrib)

post.distrib.df %>%
  ggplot(aes(var1, fill="g"), color = "grey", show.legend=FALSE) +
  geom_density(alpha = 0.4, show.legend=FALSE) +
  scale_x_continuous() +
  #labs(title="Density ditribution of clusters") +
  labs(y="Density") +
  labs(x="Globular Clusters")

post.distrib.df %>%
  ggplot(aes(var1), color = "grey", show.legend=FALSE) +
  geom_histogram(bins = 30) +
  geom_point(aes(x=test.sample$N_GC[galaxy.num], y=1000), color="blue", alpha=1, size=5)+
  geom_errorbarh(aes(y=1000, xmin=test.sample$N_GC.obs.low[galaxy.num],
                     xmax=test.sample$N_GC.obs.up[galaxy.num]), color="blue", height=20) +
  #scale_x_continuous() +
  xlim(0,60) +
  #labs(title="Density ditribution of clusters") +
  labs(y="Frequency") +
  labs(x="Globular Clusters (linear scale)")




mean(post.distrib.df$var1)
(ci.MC <- c(-1,1)*qnorm(0.975)*sd(post.distrib.df$var1)+mean(post.distrib.df$var1))


# compute_posterior_distrib <- function(type, magnitude, distance){
#   mu.MC <-  as.numeric(pmed_coef["alpha"] +
#                          pmed_coef[1+type]*magnitude +
#                          pmed_coef["gamma"]*magnitude^2 +
#                          pmed_coef["delta"]*distance)
#   mu.MC <- 10^mu.MC
#   post.MC <- rt(15000, df=pmed_coef["df"], ncp=mu.MC)
#   c(mean(post.MC), sd(post.MC))
# }



