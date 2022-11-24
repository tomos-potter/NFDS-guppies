# MODELS AND PLOTS

# This script has 4 subsections:

# 1. Models and diagnostics
# 2. Collating results
# 3. Model predictions
# 4. Figures


# Load packages for:

library(data.table) # data wrangling
library(tidyr) # more data wrangling
library(glmmTMB) # fitting models
library(DHARMa) # model diagnostics
library(ggplot2) # plotting figures
library(patchwork) # making multi-panel figures

# Load cleaned data:

data <- fread("NFDS_guppy_data.csv") # main dataset

unions <- fread("NFDS_guppy_unions.csv") # data on unions/matings

##################################################

# 1. THE MODELS

##################################################

# A. FITNESS COMPONENTS

# In this section, we assess how male rarity and novelty
# influences components of fitness for:
#  - males
#  - females
#  - offspring / grandoffspring

# B. INBREEDING AVOIDANCE

# In this section, we assess whether female preference for rare/novel males
# acts as an inbreeding avoidance mechanism

# Each of the models in A and B are run using calculations of male rarity
# at three spatial levels:

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


# NUMBER OF MATING PARTNERS PER MALE PER MONTH
# Now, let's look at the number of mating partners a male has per month, given his rarity and novelty
# N.B. Because this count data is over-dipsersed (mean < variance),
# we use a negative binomial model with a log-link.
# Same random effects structure as the survival models

# local
n_partners_local<- glmmTMB(n_mates ~ rareness_local + loc_arrived + neigh_arrived + 
                             (1 | sampling) + (1 | FishID) + (1|cohort)+
                             (1 | standardized_location), 
                           data=male_data, 
                           family="nbinom2")

#simulationOutput <- simulateResiduals(fittedModel = n_partners_local, plot = T, n=3000)
#testOutliers(simulationOutput, type="bootstrap")
#testQuantiles(simulationOutput)

# neighborhood
n_partners_neigh<- glmmTMB(n_mates ~ rareness_neigh + loc_arrived + neigh_arrived + 
                             (1 | sampling) + (1 | FishID) + (1|cohort) +
                             (1 | standardized_location),  
                           data=male_data, 
                           family="nbinom2")

#simulationOutput <- simulateResiduals(fittedModel = n_partners_neigh, plot = T, n=3000)
#testQuantiles(simulationOutput)

# population
n_partners_pop<- glmmTMB(n_mates ~ rareness_pop + loc_arrived + neigh_arrived + 
                           (1 | sampling) + (1 | FishID) + (1|cohort) +
                           (1 | standardized_location),
                         data=male_data, 
                         family="nbinom2")

simulationOutput <- simulateResiduals(fittedModel = n_partners_pop, plot = T, n=3000)
testQuantiles(simulationOutput)

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
# conceived in the same month into a cohort), and "dadID" and "grampsID" 
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
# IS FEMALE PREFERENCE FOR RARE/NOVEL MALES A MECHANISM TO AVOID INBREEDING?

# Here, we model the kinship coefficient of mating partners as a function of 
# the rarity and novelty of the male partner.

# Rationale: if female preference for rare or novel males acts as a mechanism to avoid
# inbreeding, i.e. to avoid mating with males to whom they are closely related,
# then we would expect the kinship coefficient to be lower among mating partners
# when the male is rare or novel.

# Because kinship coefficients have an unusual distribution (the inclusive
# range of values is 0 - 0.5), we employed a "hurdle" type modelling approach.

# First, we modeled the zero and non-zero kinship terms using logistic regression.
# This model asks: is the probability that mating partners are completely unrelated
# influenced by the rarity or novelty of the male partner?

# We then used a linear model for the log-transformed non-zero terms. This model
# asks: is kinship coefficient among mating partners influenced by the rarity 
# or novelty of the male partner?


# Probability that mating partners are unrelated

# local
unrelated_parents_loc <- glmmTMB(unrelated ~ dad_rareness_local + loc_arrived + neigh_arrived,
                                  data=union_data, family="binomial")

#simulationOutput <- simulateResiduals(fittedModel = unrelated_parents_loc, plot = T, n=3000)

# neigh
unrelated_parents_neigh <- glmmTMB(unrelated ~ dad_rareness_neigh + loc_arrived + neigh_arrived ,
                                 data=union_data, family="binomial")

#simulationOutput <- simulateResiduals(fittedModel = unrelated_parents_neigh, plot = T, n=3000)

# pop
unrelated_parents_pop <- glmmTMB(unrelated ~ dad_rareness_pop + loc_arrived + neigh_arrived,
                                   data=union_data, family="binomial")

#simulationOutput <- simulateResiduals(fittedModel = unrelated_parents_pop, plot = T, n=3000)

# Influence of rarity/novelty on parental kinship

# local
parental_kinship_loc <- glmmTMB(ln_parental_kinship ~ dad_rareness_local + loc_arrived + neigh_arrived,
                                data=union_data[unrelated==FALSE,], family="gaussian")

#simulationOutput <- simulateResiduals(fittedModel = parental_kinship_loc, plot = T, n=3000)

# neigh
parental_kinship_neigh <- glmmTMB(ln_parental_kinship ~ dad_rareness_neigh + loc_arrived + neigh_arrived,
                                data=union_data[unrelated==FALSE,], family="gaussian")

#simulationOutput <- simulateResiduals(fittedModel = parental_kinship_neigh, plot = T, n=3000)

# pop
parental_kinship_pop <- glmmTMB(ln_parental_kinship ~ dad_rareness_pop + loc_arrived + neigh_arrived,
                                  data=union_data[unrelated==FALSE,], family="gaussian")

#simulationOutput <- simulateResiduals(fittedModel = parental_kinship_pop, plot = T, n=3000)


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
results <- rbind(#results_fixef(p_mating_local, p_mating_neigh, p_mating_pop, "Male mating probability"),
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


results_kinship <- rbind(results_fixef(unrelated_parents_loc, unrelated_parents_neigh, unrelated_parents_pop, "Unrelated parents"),
                         results_fixef(parental_kinship_loc, parental_kinship_neigh, parental_kinship_pop, "log-kinship parents"))

# convert results into .csv friendly format
df2 <- apply(results_kinship,2,as.character)

# save results in a .csv file
write.csv(file="kinship_results.csv", df2, row.names = FALSE)


##################################################
# 3. MODEL PREDICTIONS

# Here, we generate predicted values using the best-fitting models of fitness components

# data frame of predictors

preds <- data.frame(rareness_neigh = seq(from= -1.5, to = 1.5, length.out = 100),
                    rareness_pop = seq(from= -1.5, to = 1.5, length.out = 100),
                    dad_rareness_local = seq(from= -1.5, to = 1.5, length.out = 100),
                    dad_rareness_neigh = seq(from= -1.5, to = 1.5, length.out = 100),
                    dad_rareness_pop = seq(from= -1.5, to = 1.5, length.out = 100),
                    neigh_arrived = TRUE,
                    loc_arrived = TRUE,
                    FishID = NA,
                    momID = NA,
                    standardized_location = NA,
                    dadID = NA,
                    sampling=NA,
                    cohort=NA,
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


n_partners <- as.data.table(predict(n_partners_neigh, 
                                    newdata = preds, 
                                    re.form = NA, 
                                    se.fit = TRUE))

names(n_partners)<-c("n_partners_lin", "n_partners_se")

n_partners$n_partners <-  exp(n_partners$n_partners_lin)
n_partners$n_partners_l <-  exp(n_partners$n_partners_lin - (1.96 *n_partners$n_partners_se))
n_partners$n_partners_u <-  exp (n_partners$n_partners_lin + (1.96 *n_partners$n_partners_se))

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

# Bind the predictions together into a single data.table

pred <- cbind(preds, n_off, n_partners, n_off_per_partner, n_goff)
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

# number of partners
comp_rarity[, round(n_partners[1]/n_partners[2], 2)] # 36% greater in rare
# n offspring recruited
comp_rarity[, round(n_off[1]/n_off[2], 2)] # 38% greater in rare

# Male fitness components (novelty):

# n offspring recruited
comp_novelty[, round(n_off[1]/n_off[2], 3)] # 50% greater in novel 
# number of partners
comp_novelty[, round(n_partners[1]/n_partners[2], 2)] # 45% greater in novel 

# LITTER SIZE - no significant effects of rarity
# GRAND-OFFSPRING PER LITTER - significant effect of rarity only
comp_rarity[, round(ngoff[1]/ngoff[2], 2)] # 48% greater in rare


##################################################
# 4. FIGURES

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
# The actual figures - MAIN TEXT

# FIGURE 1 shows examples of guppy color patterns and their inheritance and was not made in R

# FIGURE 2
# Composite figure showing predicted fitness functions

# plot of number of partners per month for each male
# as a function of rarity and novelty
n_partners_plot <- ggplot(pred, aes(x=rareness_neigh, y=n_partners, colour=novelty))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=0.6, xend=log(0.5)), linetype="dashed", colour="grey")+
  geom_segment(aes(x = log(2), y=-Inf, yend=0.6, xend=log(2)), linetype="dashed", colour="grey")+
  annotate("text", x = c(log(0.5), log(2)), y = 0.65, label = c("r", "c"), size = 4)+
  geom_ribbon(data=pred, aes(ymin=n_partners_l, ymax=n_partners_u, fill=novelty),
              colour=NA, alpha=0.2)+ 
  geom_line()+                  
  labs(x=expression(italic(r)[italic(i)]), 
       y="Number of partners",
       subtitle="A",
       fill="",
       colour="")+
  ylim(c(0, 0.65))+
  scale_colour_manual(values=guppy_palette[c(1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
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
       subtitle="B",
       fill="",
       colour="")+
  scale_colour_manual(values=guppy_palette[c(1,7,1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
  ylim(c(0, 0.65))+
  theme_classic()


# plot of mean litter size, as a function of male partner rarity and novelty
off_per_union_plot <- ggplot(pred, aes(x=dad_rareness_neigh, y=n_off_per_partner, colour=novelty))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_segment(aes(x = log(0.5), y=-Inf, yend=1.35, xend=log(0.5)), linetype="dashed", colour="grey")+
  geom_segment(aes(x = log(2), y=-Inf, yend=1.35, xend=log(2)), linetype="dashed", colour="grey")+
  annotate("text", x = c(log(0.5), log(2)), y = 1.5, label = c("r", "c"), size = 4)+
  geom_ribbon(data=pred, aes(ymin=n_off_per_partner_l, ymax=n_off_per_partner_u,fill=novelty),
              colour=NA, alpha=0.2)+ 
  geom_line()+
  labs(x=expression(italic(r)[italic(i)]), 
       y="Offspring per mating",
       subtitle = "C",
       fill="",
       colour="")+
  ylim(c(0, 1.5))+
  annotate("text", x = 1.25, y = 0.9, label = "N. S.", size=4)+
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
       subtitle = "D",
       fill="",
       colour="")+
  ylim(c(0, 3))+
  scale_colour_manual(values=guppy_palette[c(1,7)])+
  scale_fill_manual(values=guppy_palette[c(1,7)])+
  theme_classic()

# collate these ointo 2 multipanel figures

Fig2 <- (n_partners_plot | off_per_m_plot) /
  ( off_per_union_plot | goff_per_union_plot) +
  plot_layout(guides="collect") & 
  theme(legend.position = 'bottom')

Fig2

ggsave("Fig2.pdf", width=4.75, height=5.25, units="in")


# Figure 3
# This figure shows the fate of rarity over three generations

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


Fig3 <-ggplot(gfs[is_rare_grandad==TRUE,], aes(x=generation, y=rareness, group=FishID))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(alpha=0.1, size=0.25)+
  annotate("text", x = c(1,2,3), y = 1.4, label = c("100%", "69%", "38%"))+
  labs(x= "Generation", y = expression(italic(r)[italic(i)]))+
  theme_classic()

Fig3

ggsave("Fig3.pdf", width=2.25, height=2.25, units="in")



#-------

# SUPPLEMENTARY FIGURES 

# Fig S1
# proportion of individuals that move between neighborhoods each month
movement_time <-data[ (sex_stage=="M" | sex_stage=="F"), 
                      .(Neighborhood = table(neigh_arrived)[2]/sum(table(neigh_arrived)),
                        Pool = table(loc_arrived)[2]/sum(table(loc_arrived)) ), 
                      by=c("sampling", "sex_stage")]

movement_time <- pivot_longer(movement_time, cols=c("Neighborhood", "Pool"), 
                              names_to = "loc")

movement_time <- as.data.table(movement_time)

move_plot <- ggplot(movement_time[sampling>12], 
                    aes(x=sampling, y=value, colour=sex_stage))+
  geom_line( )+ 
  labs(x="Months since introduction", 
       y= "Proportion of new arrivals", 
       colour="Sex") + 
  scale_colour_manual(values=guppy_palette[c(1,7)])+
  theme_minimal()+
  facet_wrap(~loc)

ggsave("FigS1.png", width = 8.76/1.5, height=6.19/1.5)


# Fig S2

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

# Figure S3

# change in population/neighborhood size over study duration 
pop_size <- data[!is.na(neighborhood), length(unique(FishID)), by=c("neighborhood", "sampling")]
pop_size$neighborhood<- as.factor(pop_size$neighborhood)

tot_pop <- pop_size[, sum(V1), by=sampling]
tot_pop$neighborhood <- "Total population"

pop_size_fig <- ggplot(pop_size[sampling>12,], 
                       aes(x=sampling, y=V1, colour=neighborhood))+
  geom_line()+
  geom_line(data=tot_pop[sampling>12,], aes(x=sampling, y=V1))+
  labs(y="Number of fish", 
       x="Months since experimental introduction", 
       colour="Neighborhood")+
  scale_colour_manual(values=guppy_palette[c(1,2,3,7,4)])+
  theme_minimal()
ggsave("FigS3.pdf")


# Figure S4
# Distribution of parental kinship values

figs4 <- ggplot(union_data, aes(x=parental_kinship))+
            geom_vline(xintercept = 0.25, linetype="dashed", 
                       colour = guppy_palette[1])+
            geom_vline(xintercept = 0.125, linetype="dashed", 
                       colour = guppy_palette[7])+
            geom_vline(xintercept = 0.125/2, linetype="dashed", 
                       colour = guppy_palette[2])+
            geom_histogram(binwidth=0.0025, fill = guppy_palette[4], alpha=0.7)+
            labs(x="Kinship coefficient of mating partners",
                 y="Frequency")+
            theme_minimal()
ggsave("FigS4.pdf")

#####################
# END OF SCRIPT

