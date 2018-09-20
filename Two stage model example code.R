
#################################################################################################
## Load required packages ##
#################################################################################################

require(brms) #Used for bayesian hierarchical modeling
require(rstan) #Used for running the sampler
require(logistf) #Used for Firth corrected penalized likelihood logistic regression

#################################################################################################
## Make example dataset where the outcome is binary and has low prevalence##
#################################################################################################

data<-data.frame(outcome = sample(c(0,1), 100, replace=T, prob=c(0.9,0.1)), covariate1=sample(1:9, 100, replace = TRUE)/10, covariate2= sample(1:100, 100, replace=T), time=rep(0:3,25), randomeff1=rep(1:25, each=4), randomeff2=rep(rep(1:20), 5)) #Make sure model outcome is in numeric format
attach(data)

#################################################################################################
## Stage 1: Running Firth corrected penalized likelihood logistic regression ##
#################################################################################################

logf <- logistf(outcome ~ covariate1 + covariate2 + time + covariate1*time + covariate2*time, data = data, firth=T) #note random effects not included
summary(logf)
coef(logf) #Use for prior distribution mean 
diag(vcov(logf)) #Use for prior distribution variance

#################################################################################################
## Stage 2: Running Bayesian hierarchical model ##
#################################################################################################

#To see what priors are needed to identify in specified Baysiean hierarchical model and identify coefficient labels for priors specified below
get_prior(formula = outcome ~ covariate1 + covariate2 + time + covariate1*time + covariate2*time + (1|randomeff1) + (1|randomeff2), data = data, family = binomial)

#To define prior distributions for Bayesian hierarchical model
#We define the prior disribution as normal with mean obtained from estimate in Stage 1 and sigma obtained from variance-covariance matrix in Stage 1
#We use noninformative, flat priors for the random effects as they just need to account for correlation
#You could use a gelman cauchy distribution in prior - set_prior("cauchy(1,2)", class = "b", coef = "sex1Male"
prior <- c(
set_prior("normal(coef(logf)[2],sqrt(diag(vcov(logf))[2]))", class = "b", coef = "covariate1", lb= NULL,ub=NULL),
set_prior("normal(coef(logf)[3],sqrt(diag(vcov(logf))[3]))", class = "b", coef = "covariate2", lb= NULL,ub=NULL),
set_prior("normal(coef(logf)[4],sqrt(diag(vcov(logf))[4]))", class = "b", coef = "time", lb= NULL,ub=NULL),
set_prior("normal(coef(logf)[5],sqrt(diag(vcov(logf))[5]))", class = "b", coef = "covariate1:time",lb= NULL,ub= NULL),
set_prior("normal(coef(logf)[6],sqrt(diag(vcov(logf))[6]))", class = "b", coef = "covariate2:time", lb= NULL,ub=NULL),
set_prior("inv_gamma(0.001,0.001)", class = "sd", coef= "Intercept", group = "randomeff1", lb= NULL,ub=NULL),
set_prior("inv_gamma(0.001,0.001)", class = "sd", coef= "Intercept", group = "randomeff2", lb= NULL,ub=NULL)
)

#To run the Bayesian hierarchical model using our specified prior
#We use a burn-in of 2,000 iterations, iter for 10,000 with 5 chains and thinning of 5
#Whenever you see the warning "There were x divergent transitions after warmup." you can think about increasing adapt_delta. To do this, write control = list(adapt_delta = ), where should usually be value between 0.8 (current default) and 1. Increasing adapt_delta will slow down the sampler but will decrease the number of divergent transitions threatening the validity of your posterior samples.
fit <- brm(formula = outcome ~ covariate1 + covariate2 + time + covariate1*time + covariate2*time + (1|randomeff1) + (1|randomeff2), data = data, family = binomial, prior = prior, warmup = 2000, iter = 10000, chains = 5, control = list(adapt_delta = 0.95), thin=5)

#################################################################################################
## Checking model fit and obtaining credible intervals ##
################################################################################################# 

summary(fit) #check model fit summary
plot(fit) #plot model fit
launch_shinystan(fit) #can use shinystan app to examine model fit
marginal_effects(fit,surface=T)
WAIC(fit, fit1, compare = T) #if you run multiple model, you can compare different brm models by  “Watanabe-Akaike Information Criterion”

stanplot(fit) #plot posterior intervals
posterior_interval(fit)
stanplot(fit,pars="^b_") #only show population-level effects in the plots
stanplot(fit,type="hist") #show histograms of posterior distributions
stanplot(fit,type="neff") #plot diagnostics of the sampler
stanplot(fit,type="rhat")
posterior_samples(fit, "^b") #extract posterior samples of population-level effects 
posterior_samples(fit, "^sd_") #extract posterior samples of group-level standard deviations

#For obtaining credible intervals
samps<-fit$fit
plot(samps)
round(summary(samps, probs=c(.025, .975))$summary, 3) #this examines 2.5% and 97.5% credible intervals

