# MODELS AND PLOTS

# This script has 4 subsections:

# 1. Models and diagnostics
# 2. Collating results
# 3. Model predictions
# 4. Figures


# Load packages for:

library(data.table) # data wrangling
library(glmmTMB) # fitting models
library(DHARMa) # model diagnostics

# Load cleaned data:

data <- fread("NFDS_guppy_data.csv") # main dataset

unions <- fread("NFDS_guppy_unions.csv") # data on unions/matings

##################################################

# 1. THE MODELS

# Here, we are assessing how male rareness influences components of fitness for:
#  - males
#  - females
#  - offspring / grandoffspring

# Each of these models are run using calculations of male rareness at three levels:

# - local (i.e. the pool or riffle)
# - neighborhood (as defined through network analysis)
# - population (everyone that month)

##################################################
# MALE COMPONENTS OF FITNESS

# subset data
male_data <- data[sampling >12 & sex_stage=="M" & 
                    !is.na(rareness_local) & 
                    !is.na(rareness_neigh) & 
                    !is.na(rareness_pop),]

# NB I am removing all NA rarity scores at ALL spatial scales, because the three models
# at different scales need to be fit with the same amount of data to perform AIC comparison later

# OFFSPRING PER MONTH
# First, we will look at the number of offspring recruited per male per month
# This is our general model of how male reproductive success is influenced by rarity and novelty
# N.B. Because this count data is over-dipsersed (mean < variance),
# we use a negative binomial model with a log-link.
# To account for non-independence of data, we include random intercepts for 
# sampling month (i.e. time since start of experiment), fish ID (to account for repeated measures on individuals),
# and standardized location (i.e. where the fish was observed that month)

# local
n_offspring_local <- glmmTMB(n_offspring_conceived ~ rareness_local + loc_arrived + neigh_arrived +
                               (1 | sampling) + (1 | FishID) +
                               (1 | standardized_location),
                             data = male_data, 
                             family = "nbinom2")

# check model diagnostics
# here we simulating the model residuals and 
# performing outlier tests, dispersion tests, Kolgomorov Smirnoff tests
# and testing for trends in the residuals vs predicted values

# (UNCOMMENT THESE LINES TO RUN DIAGNOSTICS)

#simulationOutput <- simulateResiduals(fittedModel = n_offspring_local, plot = T, n=2000)
#testQuantiles(simulationOutput)


# neighborhood
n_offspring_neigh <- glmmTMB(n_offspring_conceived ~ rareness_neigh + loc_arrived + neigh_arrived +
                               (1 | sampling) + (1 | FishID) +
                               (1 | standardized_location), 
                             data = male_data, 
                             family= "nbinom2")
# check model diagnostics
#simulationOutput <- simulateResiduals(fittedModel = n_offspring_neigh, plot = T, n = 2000)
#testQuantiles(simulationOutput)

# population
n_offspring_pop <- glmmTMB(n_offspring_conceived ~ rareness_pop + loc_arrived + neigh_arrived +
                             (1 | sampling) + (1 | FishID) +
                             (1 | standardized_location),
                           data = male_data, 
                           family = "nbinom2")
# check model diagnostics
#simulationOutput <- simulateResiduals(fittedModel = n_offspring_pop, plot = T, n = 2000)
#testQuantiles(simulationOutput)

##################################################
# SURVIVAL
# Now, we will consider how male rareness and novelty influence monthly survival probability
# Because our response is binary, we use logistic regression with a logit-link
# We additionally inlude "season" as a random intercept to improve the model diagnostics

# local
surv_local <- glmmTMB(survived ~ rareness_local + loc_arrived + neigh_arrived + 
                        (1 | sampling) + (1 | FishID) + (1|season)+
                        (1 | standardized_location),
                      data = male_data,
                      family = "binomial")                                   

# check model diagnostics
#simulationOutput <- simulateResiduals(fittedModel = surv_local, plot = T, n=2000)
#testQuantiles(simulationOutput)


# neighborhood
surv_neigh <- glmmTMB(survived ~ rareness_neigh + loc_arrived + neigh_arrived +
                        (1 | sampling) + (1 | FishID) + (1|season)+
                        (1 | standardized_location), 
                      data = male_data, 
                      family= "binomial")
# check model diagnostics
#simulationOutput <- simulateResiduals(fittedModel = surv_neigh, plot = T, n=2000)
#testQuantiles(simulationOutput)

# population
surv_pop <- glmmTMB(survived ~ rareness_pop + loc_arrived + neigh_arrived +
                      (1 | sampling) + (1 | FishID) + (1|season)+
                      (1 | standardized_location), 
                    data = male_data, 
                    family= "binomial")

# check model diagnostics
#simulationOutput <- simulateResiduals(fittedModel = surv_pop, plot = T, n=2000)
#testQuantiles(simulationOutput)

##################################################

# PROBABILITY OF MATING 
# Next, we decompose male reproductive success
# lets look at the probability of mating given rarity and novelty
# N.B. Because this is a binary outcome (mated/did not mate), we use logistic regression with a logit-link

# local rareness
p_mating_local <- glmmTMB(bred ~ rareness_local + loc_arrived + neigh_arrived +
                            (1 | sampling) + (1 | FishID) +
                            (1 | standardized_location), 
                          data = male_data, 
                          family="binomial")

#simulationOutput <- simulateResiduals(fittedModel = p_mating_local, plot = T, n=2000)
#testQuantiles(simulationOutput)

# neighborhood rareness
p_mating_neigh <- glmmTMB(bred ~ rareness_neigh + loc_arrived + neigh_arrived + 
                            (1 | sampling) + (1 | FishID) + 
                            (1 | standardized_location),  
                          data = male_data, 
                          family="binomial")

#simulationOutput <- simulateResiduals(fittedModel = p_mating_neigh, plot = T, n=2000)
#testQuantiles(simulationOutput)

# population rareness
p_mating_pop <- glmmTMB(bred ~ rareness_pop + loc_arrived + neigh_arrived +
                          (1 | sampling) + (1 | FishID) + 
                          (1 | standardized_location),  
                        data = male_data,
                        family="binomial")

#simulationOutput <- simulateResiduals(fittedModel = p_mating_pop, plot = T, n=2000)
#testQuantiles(simulationOutput)

# NUMBER OF MATING PARTNERS PER MALE PER MONTH
# Now, let's look at the number of mating partners a male has per month, given his rarity and novelty
# N.B. because this is count data, but with no zeroes (we are considering the number of mates a male 
# has, given that he mated at least once), we use a truncated Poisson distribution (log-link)
# This models the *additional* number of partners he had, given that he had at least one.
# Same random effects structure as the survival models

# local
n_partners_local<- glmmTMB(n_mates ~ rareness_local + loc_arrived + neigh_arrived + 
                             (1 | sampling) + (1 | FishID) + (1|season)+
                             (1 | standardized_location), 
                           data=male_data[bred==TRUE,], 
                           family="truncated_poisson")

#simulationOutput <- simulateResiduals(fittedModel = n_partners_local, plot = T, n=3000)
#testQuantiles(simulationOutput)

# neighborhood
n_partners_neigh<- glmmTMB(n_mates ~ rareness_neigh + loc_arrived + neigh_arrived + 
                             (1 | sampling) + (1 | FishID) + (1|season)+
                             (1 | standardized_location),  
                           data=male_data[bred==TRUE,], 
                           family="truncated_poisson")

#simulationOutput <- simulateResiduals(fittedModel = n_partners_neigh, plot = T, n=3000)
#testQuantiles(simulationOutput)

# population
n_partners_pop<- glmmTMB(n_mates ~ rareness_pop + loc_arrived + neigh_arrived + 
                           (1 | sampling) + (1 | FishID) + (1|season)+
                           (1 | standardized_location),
                         data=male_data[bred==TRUE,], 
                         family="truncated_poisson")

#simulationOutput <- simulateResiduals(fittedModel = n_partners_pop, plot = T, n=3000)
#testQuantiles(simulationOutput)

##################################################
# In the next section, we look at how many offspring are produced in each mating (union)
# Note that these are models of the outcome of the UNION - results apply to BOTH male and female partners
##################################################


union_data <- unions[sampling>12 &
                       !is.na(dad_rareness_local) &
                       !is.na(dad_rareness_neigh) &
                       !is.na(dad_rareness_pop),]

# OFFSPRING PER UNION

# N.B. each union results in at least one offspring by default (we can only measure unions
# if they result in at least one offspring that is recruited into the population)

# This means our count data lack any zeros, so we use a truncated Poisson (log-link)
# This models the *additional* number of offspring, given that they had at least one.
# We have a different random effects structure here, to account for repeated measures on both females and males
# We fit random intercepts for momID, dadID, and sampling month

# local
offspring_per_union_local <- glmmTMB(off_per_union ~ dad_rareness_local + loc_arrived + neigh_arrived +
                                       (1|dadID) + (1|momID) + 
                                       (1|sampling),
                                     data = union_data, 
                                     family = "truncated_poisson")

#simulationOutput <- simulateResiduals(fittedModel = offspring_per_union_local, plot = T, n=3000)
#testQuantiles(simulationOutput)

# neighborhood
offspring_per_union_neigh <- glmmTMB(off_per_union ~ dad_rareness_neigh + loc_arrived + neigh_arrived +
                                       (1|dadID) + (1|momID) + 
                                       (1|sampling),
                                     data = union_data, 
                                     family = "truncated_poisson")

#simulationOutput <- simulateResiduals(fittedModel = offspring_per_union_neigh, plot = T, n=3000)
#testQuantiles(simulationOutput)

# population
offspring_per_union_pop <- glmmTMB(off_per_union ~ dad_rareness_pop + loc_arrived + neigh_arrived +
                                     (1|dadID) + (1|momID) + 
                                     (1|sampling),
                                   data = union_data, 
                                   family = "truncated_poisson")

#simulationOutput <- simulateResiduals(fittedModel = offspring_per_union_pop, plot = T, n=3000)
#testQuantiles(simulationOutput)

##################################################

# GRANDOFFSPRING PER UNION

# Next we look at the indirect fitness consequences of male rareness
# To do this, we consider the number of grand-offspring that ultimately arise from unions
# The number of grand-offspring is the sum of the LRS of each offspring produced per union

# N.B. Some (many) individuals have zero grand-offspring, despite producing offspring.
# So again, as for the male reproductive success model, we use a negative binomial model

# local
goffspring_per_union_local <- glmmTMB(goff_per_union ~ dad_rareness_local + loc_arrived + neigh_arrived +
                                        (1|dadID) + (1|momID) + 
                                        (1|sampling),
                                      data = union_data, 
                                      family= "nbinom2",
                                      control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))

#simulationOutput <- simulateResiduals(fittedModel = goffspring_per_union_local, plot = T, n=3000)
#testQuantiles(simulationOutput)

# neighborhood
goffspring_per_union_neigh <- glmmTMB(goff_per_union ~ dad_rareness_neigh + loc_arrived + neigh_arrived +
                                        (1|dadID) + (1|momID) + 
                                        (1|sampling),
                                      data = union_data,
                                      family= "nbinom2",
                                      control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))

#simulationOutput <- simulateResiduals(fittedModel = goffspring_per_union_neigh, plot = T, n=3000)
#testQuantiles(simulationOutput)

# population
goffspring_per_union_pop <- glmmTMB(goff_per_union ~ dad_rareness_pop + loc_arrived + neigh_arrived +
                                      (1|dadID) + (1|momID) + 
                                      (1|sampling),
                                    data = union_data,
                                    family= "nbinom2",
                                    control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))

#simulationOutput <- simulateResiduals(fittedModel = goffspring_per_union_pop, plot = T, n=3000)
#testQuantiles(simulationOutput)
##################################################

# REPRODUCTIVE SUCCESS OF OFFSPRING, GIVEN THE RARITY OF THEIR FATHER

# Here, we model the monthly reproductive success of males and females, as a function of 
# the pattern rarity of their fathers and paternal grandfathers
# Again, this is zero-inflated count data, so we use a negative binomial model with a log link

# We use a slightly different random effects structure (earlier versions led to 
# convergence problems / diagnostic problems). We include "conceived" (which groups individuals
# concieved in the same month into a cohort), and "dadID" and "grampsID" 
# to account for paternal and grandpaternal identity as grouping factors.

# Select data
off_data <- data[sampling>12 &
                   sex_stage!="I" & sex_stage!="X" & sex_stage!="" & 
                   !is.na(dad_rareness_local) & 
                   !is.na(dad_rareness_neigh) & 
                   !is.na(dad_rareness_pop) &
                   !is.na(gramps_rareness_local) &
                   !is.na(gramps_rareness_neigh)   ,]

# local
offspring_reproduction_local <- glmmTMB(n_offspring_conceived ~ dad_rareness_local * sex_stage + gramps_rareness_local*sex_stage +
                                          (1 | sampling) + (1 | conceived) + 
                                          (1 | FishID)  + (1| dadID) + (1|grampsID)+
                                          (1 | standardized_location),   
                                        data = off_data, 
                                        family= "nbinom2")


#simulationOutput <- simulateResiduals(fittedModel = offspring_reproduction_local, plot = T, n=3000)
#testQuantiles(simulationOutput)

# neighborhood
offspring_reproduction_neigh <- glmmTMB(n_offspring_conceived ~ dad_rareness_neigh * sex_stage + gramps_rareness_neigh*sex_stage +
                                          (1 | sampling) + (1 | conceived) + 
                                          (1 | FishID)  + (1| dadID) + (1|grampsID)+
                                          (1 | standardized_location),   
                                        data = off_data, 
                                        family= "nbinom2")

#simulationOutput <- simulateResiduals(fittedModel = offspring_reproduction_neigh, plot = T, n=3000)
#testQuantiles(simulationOutput)

# population
offspring_reproduction_pop <- glmmTMB(n_offspring_conceived ~ dad_rareness_pop*sex_stage + gramps_rareness_pop*sex_stage +
                                        (1 | sampling) + (1 | conceived) + 
                                        (1 | FishID)  + (1| dadID) + (1|grampsID)+
                                        (1 | standardized_location),   
                                      data = off_data, 
                                      family= "nbinom2",
                                      control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))

#simulationOutput <- simulateResiduals(fittedModel = offspring_reproduction_pop, plot = T, n=3000)
#testQuantiles(simulationOutput)

##################################################

# 2. Collating results

# Function for summarizing results from each of the 3 models for a given fitness component
# Extracts sample sizes, fixed effects estimates, sd, p_values, calculates delta AIC scores

results_fixef <- function(m1, m2, m3, model){
  
  # calculate delta AIC scores for the 3 models
  # (best fitting model has dAIC=0)
  dAIC <- round(AIC(m1, m2, m3)[2] - min(AIC(m1,m2, m3)[2]), 2)
  q <-length(fixef(m1)$cond)
  
  # this is the MS friendly version
  r <- data.table(Model = c(paste(model),rep("", (q * 3)-1)),
                  Model_type = c(family(m1)[1], rep("", (q * 3)-1)),
                  N_obs = c(m1$modelInfo$nobs, rep("", (q * 3)-1)),
                  Rareness_scale =  c("Local", rep("",q-1 ),
                                      "Neighborhood",rep("", q-1),
                                      "Population",rep("", q-1)),
                  dAIC = c(dAIC[1,1], rep("", q-1),
                           dAIC[2,1], rep("", q-1),
                           dAIC[3,1], rep("", q-1)),
                  Parameter = c(paste(names(fixef(m1)$cond)), paste(names(fixef(m2)$cond)), paste(names(fixef(m3)$cond))),
                  Estimate = round(c(fixef(m1)$cond, fixef(m2)$cond, fixef(m3)$cond),2),
                  SE = round(c(summary(m1)$coefficients$cond[,2],
                               summary(m2)$coefficients$cond[,2],
                               summary(m3)$coefficients$cond[,2]),2),
                  p_value = signif(c(summary(m1)$coefficients$cond[,4],
                                     summary(m2)$coefficients$cond[,4],
                                     summary(m3)$coefficients$cond[,4]), 3)
  )
  
  r$star <- ifelse(r$p_value<0.001, "***", 
                   ifelse(r$p_value<0.01, "**",
                          ifelse(r$p_value<0.05, "*", "")))
  
  r$p_value[which(r$p_value< 2*10^(-16))] <- paste("< 2e-16")
  
  # This is the R-friendly version (no empty cells, no stars, no < 2e-16s)
  # r <- data.table(Model = rep(paste(model), q),
  #                 Model_type = rep(family(m1)[1], q),
  #                 N_obs = c(m1$modelInfo$nobs, rep("", (q * 3)-1)),
  #                 Rareness_scale =  rep(c("Local", 
  #                                         "Neighborhood",
  #                                         "Population"), each=q),
  #                 dAIC = c(rep(dAIC[1,1],q),
  #                          rep(dAIC[2,1], q),
  #                          rep(dAIC[3,1], q)),
  #                 Parameter = c(paste(names(fixef(m1)$cond)), paste(names(fixef(m2)$cond)), paste(names(fixef(m3)$cond))),
  #                 Estimate = round(c(fixef(m1)$cond, fixef(m2)$cond, fixef(m3)$cond),2),
  #                 SE = round(c(summary(m1)$coefficients$cond[,2], 
  #                              summary(m2)$coefficients$cond[,2],
  #                              summary(m3)$coefficients$cond[,2]),2),
  #                 p_value = signif(c(summary(m1)$coefficients$cond[,4], 
  #                                    summary(m2)$coefficients$cond[,4],
  #                                    summary(m3)$coefficients$cond[,4]), 3)
  # )
  
  
  return(r)
  
}

# run the above function for each of the fitness components and collate in a table
results <- rbind(results_fixef(p_mating_local, p_mating_neigh, p_mating_pop, "Male mating probability"),
                 results_fixef(n_partners_local, n_partners_neigh, n_partners_pop, "Number of partners"),
                 results_fixef(n_offspring_local, n_offspring_neigh, n_offspring_pop, "Male reproduction"),
                 results_fixef(surv_local, surv_neigh, surv_pop, "Male survival"),
                 results_fixef(offspring_per_union_local, offspring_per_union_neigh, offspring_per_union_pop, "Offspring per union"),
                 results_fixef(goffspring_per_union_local, goffspring_per_union_neigh, goffspring_per_union_pop, "Grandoffspring per union"),
                 results_fixef(offspring_reproduction_local, offspring_reproduction_neigh, offspring_reproduction_pop, "Offspring reproduction"))

# convert results into .csv friendly format
df <- apply(results,2,as.character)

# save results in a .csv file
write.csv(file="main_results.csv", df, row.names = FALSE)


# Extract random effects estimates:



##################################################
# 3. MODEL PREDICTIONS

# Here, we generate predicted values using the best-fitting models of fitness components

# data frame of predictors

preds <- data.frame(rareness_neigh = seq(from= -1.5, to = 1.5, length.out = 100),
                    rareness_pop = seq(from= -1.5, to = 1.5, length.out = 100),
                    dad_rareness_neigh = seq(from= -1.5, to = 1.5, length.out = 100),
                    dad_rareness_pop = seq(from= -1.5, to = 1.5, length.out = 100),
                    neigh_arrived = TRUE,
                    loc_arrived = TRUE,
                    FishID = NA,
                    momID = NA,
                    standardized_location = NA,
                    dadID = NA,
                    sampling=NA,
                    season=NA,
                    conceived=NA,
                    dad_conceived=NA)
pred1 <-preds
pred1$neigh_arrived <- FALSE
pred1$loc_arrived <- FALSE
preds <- rbind(preds, pred1)

# negative binomial model of monthly offspring
n_off <- as.data.table(predict(n_offspring_neigh, 
                               newdata = preds, 
                               re.form = NA, 
                               se.fit = TRUE))

names(n_off)<-c("n_off_lin", "n_off_lin_se")

# inverse log to convert to observed scale
n_off$n_off <- exp(n_off$n_off_lin)

# calculate lower and upper 95% confidence intervals on the observed scale
n_off$n_off_l <- exp(n_off$n_off_lin - (1.96 *  n_off$n_off_lin_se))
n_off$n_off_u <- exp(n_off$n_off_lin + (1.96 *  n_off$n_off_lin_se))

# predict probability of mating (logistic regression, inv.link = u/(1+u), where u= exp(linear predictor)
p_mating <- as.data.table(predict(p_mating_neigh, 
                                  newdata = preds, 
                                  re.form = NA, 
                                  se.fit = TRUE))

names(p_mating)<-c("p_mating_lin", "p_mating_lin_se")

# take inverse logit
p_mating$p_mating <- 1 / (1 + exp(-(p_mating$p_mating_lin)))

# confidence intervals on observed scale
p_mating$p_mating_l <- 1 / (1 + exp(- ((p_mating$p_mating_lin) - (1.96 * p_mating$p_mating_lin_se))))
p_mating$p_mating_u <- 1 / (1 + exp(- ((p_mating$p_mating_lin) + (1.96 * p_mating$p_mating_lin_se))))

# predict number of partners, given that they mated

n_partners <- as.data.table(predict(n_partners_neigh, 
                                    newdata = preds, 
                                    re.form = NA, 
                                    se.fit = TRUE))

names(n_partners)<-c("n_partners_lin", "n_partners_se")
# log-link for truncated poisson (all individuals have at least one partner)
n_partners$n_partners <- 1 + exp(n_partners$n_partners_lin)
n_partners$n_partners_l <- 1 + exp(n_partners$n_partners_lin - (1.96 *n_partners$n_partners_se))
n_partners$n_partners_u <- 1 + exp (n_partners$n_partners_lin + (1.96 *n_partners$n_partners_se))

# predict number of offspring per partner per month, given that they mated

n_off_per_partner <- as.data.table(predict(offspring_per_union_pop, 
                                           newdata = preds, 
                                           re.form = NA, 
                                           se.fit = TRUE))

names(n_off_per_partner)<-c("n_off_per_partner_lin", "n_off_per_partner_se")

n_off_per_partner$n_off_per_partner <- 1 + exp(n_off_per_partner$n_off_per_partner_lin)
n_off_per_partner$n_off_per_partner_l <- 1+ exp(n_off_per_partner$n_off_per_partner_lin - (1.96 *  n_off_per_partner$n_off_per_partner_se))
n_off_per_partner$n_off_per_partner_u <- 1 + exp(n_off_per_partner$n_off_per_partner_lin + (1.96 *  n_off_per_partner$n_off_per_partner_se))


# negative binomial model of monthly grand_offspring per union
n_goff <- as.data.table(predict(goffspring_per_union_pop, 
                                newdata = preds, 
                                re.form = NA,
                                se.fit = TRUE))

names(n_goff)<-c("n_goff_lin", "n_goff_se")

# inverse log-link
n_goff$ngoff <- exp(n_goff$n_goff_lin)
n_goff$n_goff_l <- exp(n_goff$n_goff_lin - (1.96 *  n_goff$n_goff_se))
n_goff$n_goff_u <- exp(n_goff$n_goff_lin + (1.96 *  n_goff$n_goff_se))

# survival
surv <- as.data.table(predict(surv_neigh, 
                              newdata = preds, 
                              re.form = NA, 
                              se.fit = TRUE))

names(surv)<-c("surv_lin", "surv_lin_se")

# take inverse logit
surv$surv <- 1 / (1 + exp(-(surv$surv_lin)))

surv$surv_l <- 1 / (1 + exp(- ((surv$surv_lin) - (1.96 * surv$surv_lin_se))))
surv$surv_u <- 1 / (1 + exp(- ((surv$surv_lin) + (1.96 * surv$surv_lin_se))))

# Bind the predictions together into a single data.table
pred <- cbind(preds, n_off, p_mating, n_partners, n_off_per_partner, n_goff, surv )
pred <- as.data.table(pred)

pred$novelty <- ifelse(pred$neigh_arrived==TRUE, "new arrival", "resident")
data$novelty <- ifelse(data$neigh_arrived==TRUE, "new arrival", "resident")

# EXAMPLES OF EFFECT SIZES

# In the MS, we compare the fitness components of rare (r_i = log(0.5)) and
# common (r_i = log(2)) males/fathers. We also contrast novel (i.e. new arrivals) and
# non-novel (i.e. resident) males, at r_i = log(1). 

# As well as using our predictions for plotting (below), we also use them to
# describe the effects sizes.

# We want to compare:
#   1. rare and common residents, with rare r_i = log(0.5), common r_i = log(2)
#   2. residents and new arrivals, with r_i = log(1)

# 1. rare v common residents
comp_rarity <- pred[c(128, 173)]

# 2. novel v resident
comp_novelty <- pred[c(51, 151)]


# Male fitness components (rarity):

# n offspring recruited
comp_rarity[, round(n_off[1]/n_off[2], 2)] # 38% greater in rare
# probability of mating
comp_rarity[, round(p_mating[1]/p_mating[2], 2)] # 37% greater in rare
# number of partners
comp_rarity[, round(n_partners[1]/n_partners[2], 2)] # 12% greater in rare
# survival
comp_rarity[, round(surv[1]/surv[2], 2)] # 0% greater in rare (NOT SIGNIFICANT)

# Male fitness components (novelty):

# n offspring recruited
comp_novelty[, round(n_off[1]/n_off[2], 3)] # 50% greater in novel 
# probability of mating
comp_novelty[, round(p_mating[1]/p_mating[2], 3)] # 79% greater in novel
# number of partners
comp_novelty[, round(n_partners[1]/n_partners[2], 2)] # 11% greater in novel (NOT SIGNIFICANT)
# survival
comp_novelty[, round(surv[2]/surv[1], 2)] # 14% greater in RESIDENTS,
# We can describe survival as mortality risk by considering 1/surv probability
comp_novelty[, round((1/surv[1])/(1|surv[2]), 2)] # mortality risk is 45% greater in novel males

# LITTER SIZE - no significant effects of rarity
# GRAND-OFFSPRING PER LITTER - significant effect of rarity only
comp_rarity[, round(ngoff[1]/ngoff[2], 2)] # 48% greater in rare


##################################################
# 4. FIGURES

# required packages
library(ggplot2)
library(patchwork)

# special home-made guppy colour palette! 
# (using a photo of a male guppy and https://imagecolorpicker.com/)
guppy_palette <- c("#D05300", "#DC9415", "#434B5C", "#94875B",
                   "#D1D9EB", "#F1DBCC", "#79A2CF")

# run this code to see the pretty colours
try <- data.table(test = as.factor(c(1:length(guppy_palette))))
ggplot(try, aes(x=test, y=test, colour=test))+geom_point(size=20)+
  scale_colour_manual(values=guppy_palette)+
  theme_classic()+
  theme(legend.position = "none")

##############
# The actual figures

# distribution of male pattern rarity at the neighborhood level
dist_rare_ALL <- ggplot(data[sex_stage=="M" & sampling>12,], aes(x=rareness_neigh))+
  geom_density(alpha=0.5, adjust=2, fill = guppy_palette[4], colour=NA)+
  labs(x=expression(italic(r)[italic(i)]),
       y="Density")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=0.8, xend=log(0.5)), linetype="dashed")+
  geom_segment(aes(x = log(2), y=-Inf, yend=0.8, xend=log(2)), linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dotted")+
  annotate("text", x = c(log(0.5), log(2)), y = 0.9, label = c("rare", "common"), size = 3)+
  xlim(c(-2, 2))+
  theme_classic()


# plot of monthly reproductive success of males as a function of rarity and novelty
off_per_m_plot <- ggplot(pred, aes(x=rareness_neigh, y=n_off, colour=novelty))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=0.6, xend=log(0.5)), linetype="dashed", colour="grey")+
  geom_segment(aes(x = log(2), y=-Inf, yend=0.6, xend=log(2)), linetype="dashed", colour="grey")+
  annotate("text", x = c(log(0.5), log(2)), y = 0.65, label = c("r", "c"), size = 4)+
  geom_ribbon(data=pred, aes(ymin=n_off_l, ymax=n_off_u,fill=novelty), alpha=0.2, colour=NA)+
  geom_line()+
  labs(x=expression(italic(r)[italic(i)]),
       y="Recruited offspring",
       subtitle="C",
       fill="",
       colour="")+
  scale_colour_manual(values=guppy_palette[c(1,7,1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
  ylim(c(0, 0.65))+
  theme_classic()

# plot of monthly probability of mating for males, as function of rarity and novelty
p_mating_plot <- ggplot(pred, aes(x=rareness_neigh, y=p_mating, colour=novelty))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=0.4, xend=log(0.5)), linetype="dashed", colour="grey")+
  geom_segment(aes(x = log(2), y=-Inf, yend=0.4, xend=log(2)), linetype="dashed", colour="grey")+
  annotate("text", x = c(log(0.5), log(2)), y = 0.45, label = c("r", "c"), size = 4)+
  geom_ribbon(data=pred, aes(ymin=p_mating_l, ymax=p_mating_u, fill=novelty),
              colour=NA, alpha=0.2)+ 
  geom_line()+
  labs(x=expression(italic(r)[italic(i)]), 
       y="Mating probability",
       subtitle="A",       
       fill="",
       colour="")+
  ylim(c(0,0.45))+
  scale_colour_manual(values=guppy_palette[c(1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
  theme_classic()

# plot of number of partners per month for each male (given that they mated at least once)
# as a function of rarity and novelty
n_partners_plot <- ggplot(pred, aes(x=rareness_neigh, y=n_partners, colour=novelty))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=2.65, xend=log(0.5)), linetype="dashed", colour="grey")+
  geom_segment(aes(x = log(2), y=-Inf, yend=2.65, xend=log(2)), linetype="dashed", colour="grey")+
  annotate("text", x = c(log(0.5), log(2)), y = 2.75, label = c("r", "c"), size = 4)+
  geom_ribbon(data=pred, aes(ymin=n_partners_l, ymax=n_partners_u, fill=novelty),
              colour=NA, alpha=0.2)+ 
  geom_line()+                  
  labs(x=expression(italic(r)[italic(i)]), 
       y="Number of partners",
       subtitle="B",
       fill="",
       colour="")+
  ylim(c(1, 3))+
  scale_colour_manual(values=guppy_palette[c(1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
  theme_classic()

# survival as a function of rarity and novelty
surv_plot <- ggplot(pred, aes(x=rareness_neigh, y=surv, colour = novelty))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=0.95, xend=log(0.5)), linetype="dashed", colour="grey")+
  geom_segment(aes(x = log(2), y=-Inf, yend=0.95, xend=log(2)), linetype="dashed", colour="grey")+
  annotate("text", x = c(log(0.5), log(2)), y = 1, label = c("r", "c"), size = 4)+
  geom_ribbon(data=pred, aes(ymin=surv_l, ymax=surv_u, fill=novelty),
              colour=NA, alpha=0.2)+
  geom_line()+
  labs(x=expression(italic(r)[italic(i)]),
       y="Survival probability",
       subtitle = "D",
       fill="",
       colour="")+
  ylim(c(0, 1))+
  annotate("text", x = 1.1, y = 0.9, label = "N. S.", size=4)+
  scale_colour_manual(values=guppy_palette[c(1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
  theme_classic()

# plot of mean litter size, as a function of male partner rarity and novelty
off_per_union_plot <- ggplot(pred, aes(x=dad_rareness_neigh, y=n_off_per_partner, colour=novelty))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=1.27, xend=log(0.5)), linetype="dashed", colour="grey")+
  geom_segment(aes(x = log(2), y=-Inf, yend=1.27, xend=log(2)), linetype="dashed", colour="grey")+
  annotate("text", x = c(log(0.5), log(2)), y = 1.3, label = c("r", "c"), size = 4)+
  geom_ribbon(data=pred, aes(ymin=n_off_per_partner_l, ymax=n_off_per_partner_u,fill=novelty),
              colour=NA, alpha=0.2)+ 
  geom_line()+
  labs(x=expression(italic(r)[italic(i)]), 
       y="Offspring per mating",
       subtitle = "A",
       fill="",
       colour="")+
  ylim(c(1, 1.3))+
  annotate("text", x = 1.2, y = 1.2, label = "N. S.", size=4)+
  scale_colour_manual(values=guppy_palette[c(1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
  theme_classic()

# plot of mean total number of grandoffspring that result from a single union,
# as a function of male partner rarity and novelty
goff_per_union_plot <- ggplot(pred, aes(x=dad_rareness_neigh, y=ngoff, colour=novelty))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=2.75, xend=log(0.5)), linetype="dashed", colour="grey")+
  geom_segment(aes(x = log(2), y=-Inf, yend=2.75, xend=log(2)), linetype="dashed", colour="grey")+
  annotate("text", x = c(log(0.5), log(2)), y = 3, label = c("r", "c"), size = 4)+
  geom_ribbon(data=pred, aes(ymin=n_goff_l, ymax=n_goff_u,fill=novelty),
              colour=NA, alpha=0.2)+ 
  geom_line()+
  labs(x=expression(italic(r)[italic(i)]),
       y="Grand-offspring per mating",
       subtitle = "B",
       fill="",
       colour="")+
  ylim(c(0, 3))+
  scale_colour_manual(values=guppy_palette[c(1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
  theme_classic()

# collate these ointo 2 multipanel figures

Fig2 <- (p_mating_plot | n_partners_plot) /
  (off_per_m_plot | surv_plot) +
  plot_layout(guides="collect") & 
  theme(legend.position = 'bottom')

Fig2

ggsave("Fig2.png")

Fig3 <- (off_per_union_plot / goff_per_union_plot)+
  plot_layout(guides="collect") &
  theme(legend.position = 'bottom')

Fig3

ggsave("Fig3.png", width=4,)


# Figure 4 

library(tidyverse) # for its excellent pivot_longer() function

# In this figure, we want to link the rarity score of males
# to those of their sons and grandsons

dad_son <- data[sex_stage=="M",
                .(son = rareness_pop[1],
                  father = dad_rareness_pop[1],
                  dadID =dadID[1]),
                by=FishID]

# lets rejiggle these so we can get father-son-grandson chains going

gfs <- dad_son[,c("FishID", "father", "dadID")]
names(gfs) <- c("dadID", "grandfather", "grandfather_id" )

gfs <- merge(dad_son, gfs, by="dadID", all.x=TRUE)

gfs$is_rare_grandad <- gfs$grandfather <= log(0.5)

names(gfs)<- c("dad_id", "FishID", 
               "f3", "f2", "f1", 
               "grandfather_id", 
               "is_rare_grandad")

gfs <- pivot_longer(gfs, cols=c("f1", "f2", "f3"), 
                    names_to = "generation", 
                    values_to = "rareness")
gfs<- as.data.table(gfs)

# get percentages that are rare (i.e. negative values) in each generation
gfs[is_rare_grandad==TRUE, as.data.table(table(rareness<log(1))), by=generation]

# f1 = 100% rare
# f2 = 128 / (128+57) = 69% rare
# f3 = 87 / (145+87) = 38% rare


Fig4 <-ggplot(gfs[is_rare_grandad==TRUE,], aes(x=generation, y=rareness, group=FishID))+
  geom_line(alpha=0.1)+
  geom_hline(yintercept = 0, colour=guppy_palette[1], linetype="dashed")+
  annotate("text", x = c(1,2,3), y = 1.4, label = c("100%", "69%", "38%"))+
  labs(x= "Generation", y = expression(italic(r)[italic(i)]))+
  theme_classic()

Fig4

ggsave("Fig4.png", width=7.99*0.75, height=5.67*0.75)

#####################
# END OF SCRIPT
