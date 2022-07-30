# DATA PREPARATION

# SUB-SECTIONS:

# 1. Importing, cleaning, and merging data files
# 2. Movement between locations (novelty) and neighborhood classification
# 3. Calculating male pattern rarity
# 4. Quantifying reproductive success
# 5. Parent-offspring relationships

##################################################

# REQUIRED PACKAGES

# (use install.packages(NAME_OF_PACKAGE) to download these if you do not already have them)
library(data.table) # for data wrangling
library(tidyr) # a bit more data wrangling (specifically the pivot_longer() function)
library(igraph) # for network analysis 

##################################################

# 1. Importing, cleaning and merging data --------------------------------

# Three datasets:
#   1. Mark-recapture data with pedigree: tells us where each fish is each month, plus IDs of parents
#   2. Color patterns: tells us the colour pattern code for each male
#   3. Standardized locations: standardized names for pools/riffles each month 


# mark-recapture data
mr_data <- fread("LL_guppy_data.csv")

# colour patterns data:
pat_data <- fread("pat_data.csv")

# standardized locations:
loc_data <- fread("standardized_locations.csv")

# Checking and cleaning pattern data
# remove "99/100/101" in pat_data$pattern - these are codes for no pattern defined
# remove code 31: this indicates a fish has a unique colour pattern (v.rare)
pat_data <- pat_data[pattern<31]


# add in a few demographic bits and bobs to the data

birth_death <- mr_data[, .(first_recorded = min(sampling),
                           last_recorded = max(sampling)), by=FishID]

# merge this birth_death into mr_data
mr_data <- merge(mr_data, birth_death, by="FishID", all.x=TRUE)

# this is the number of months since first observation
mr_data$obs <- mr_data$sampling - (mr_data$first_recorded - 1)


# We need to know how rare a male pattern was at the time of breeding
# We don't actually observe mating, so need to make an assumption
# We assume that breeding occurs 3 months prior to first observation
# (guppies spend ~1 month in utero, then take ~2 months to grow to 14mm - the size at first sampling)
# NB this value also acts as a cohort identifier 
# i.e. all individuals conceived in the same month are within the same cohort
mr_data$conceived <- mr_data$first_recorded - 3


# but lets get rid of negative values for founders
mr_data$conceived[which(mr_data$conceived<0)]<- NA

# Given our high recapture rate, we assume that fish died the month
# following their final observation
mr_data$survived <- mr_data$sampling != mr_data$last_recorded


# We need to merge mr_data and pat_data
# resulting data table has a column for 'pattern' and 'cohort' added, merged by 'FishID'

data <- merge(x = mr_data, y = pat_data, 
              by="FishID",
              all=TRUE)

# Add in standardized locations
# NB location naming convention in the field data depends on the length of the reach
# This can change slighty each month - streams are dynamic systems
# But we have standardized the location names so that names are consistent between months

data <- merge(x=data, y=loc_data,
              by = c("sampling", "location"),
              all.x=TRUE)

# add in season
data$season <- ifelse(data$sampling.month>5 & data$sampling.month<12, "rainy", "dry")

data <- data[sampling<61, ]

##################################################

# 2. Novelty: Movement between localities -------------------------------------------------

# In this section, we first tally up movement of individuals between standardized locations
# Then we use network cluster analysis to classify distinct guppy neighborhoods,
# that are determined by patterns of movement of males among standardised locations

# Movement of individuals between standardized locations

#add columns for location in stream in previous month
data[ , prev_loc :=shift(standardized_location, type="lag"), by= FishID]

# add TRUE/FALSE for whether each fish arrived in / left a location in each month
data[ , loc_arrived := standardized_location !=prev_loc]

# NB - if we did not observe the fish the preceding month, 
# we cannot establish if it is a new arrival or not - will be NA

# NETWORK ANALYSIS TO IDENTIFY GUPPY NEIGHBORHOODS

# Create an adjacency matrix
# this is the number of observations of males movement between standardized locations
links_dat <- data[loc_arrived==TRUE & sex_stage=="M" & 
                    standardized_location!="" & prev_loc!="",]
links <- as.matrix(table(links_dat$standardized_location, links_dat$prev_loc))

# Build the network from the adjacency matrix
net <- graph_from_adjacency_matrix(links)

# Identify clusters within the network
neighbors <- cluster_optimal(net)

#List which neighborhood each pool belongs to
neigh_list<- as.list(membership(neighbors))
neighs_frame <- do.call(cbind.data.frame, neigh_list)
neighs <- pivot_longer(neighs_frame, cols = 1:52, 
                       names_to="standardized_location", 
                       values_to = "neighborhood")

# merge the neighborhood classifications into the data

data <- merge(x=data, y=neighs, 
              by="standardized_location",
              all.x=TRUE)

# we want to treat neighborhood as a factor, not an integer
data$neighborhood <- as.factor(data$neighborhood)

# Now lets add in logical data on whether fish are new arrivals to neighborhoods each month

#add columns for neighborhood in stream in previous/next month
data[ , prev_neigh :=shift(neighborhood, type="lag"), by= FishID]

# add TRUE/FALSE for whether each fish arrived in / left a neighborhood in each month
data[ , neigh_arrived := neighborhood !=prev_neigh]

#data$neigh_arrived[which(is.na(data$neigh_arrived))]<- FALSE


##################################################
# 3. Calculating rareness -------------------------------------------------

# We calculate male pattern rareness at 3 different spatial scales (all within a single month):

# A. local (e.g. within a particular pool / riffle)
# B. neighborhood (i.e. within one of four networks of localities with high degrees of movement within)
# C. population (i.e. at the level of the entire population)

# the general formula is: rareness = log((N_i / sum(N)) * P)
# where N_i is the number of individuals with pattern i,
# sum(N) is the total number of individuals,
# and P is the total number of patterns


##################################################


# A. Local

# How "rare" is each male in it's location in each month?

# Need to calculate relative abundance of each pattern in each location/month

# count data: frequency of each pattern in adult males in each month (sampling) and location (location)
rloc1 <- data[sex_stage=="M", as.data.table(table(pattern, sampling, location))]

# calculate total number of patterns seen in each location/month
# calculate total number of males (with recorded pattern) seen in each location/month

rloc2 <- rloc1[N>0, .(n_local_patterns = length(pattern),
                      n_local_males = sum(N)),
               by=c("sampling", "location")]


# merge location/month counts with location/time pattern frequencies 
rloc3 <- merge(x=rloc1, y=rloc2,
               by = c("sampling", "location"),
               all.y = TRUE)

# relative frequency of each pattern in location/time
rloc3$local_rel_freq <- rloc3$N/rloc3$n_local_males

# null expected relative frequency of each pattern in location/time
#t3$null_local_rel_freq <- 1 / t3$n_local_patterns

rloc3$rareness_local <- log(rloc3$local_rel_freq * rloc3$n_local_patterns)

# remove non-finite scores (i.e. where there were O males of a pattern type)
rloc3 <-rloc3[is.finite(rareness_local),]

rloc3$pattern <- as.factor(rloc3$pattern)
rloc3$sampling <- as.integer(rloc3$sampling)
data$pattern <- as.factor(data$pattern)

rloc3<- rloc3[, c("sampling", "pattern", "location", "rareness_local")]

data <- merge(x=data, y=rloc3,
              by = c("pattern", "sampling", "location"),
              all.x=TRUE)

##################################################
# B. neighborhood

# How "rare" is each adult male in it's neighborhood in each month?

rneigh1 <- data[sex_stage=="M", as.data.table(table(pattern, sampling, neighborhood))]

# calculate total number of patterns seen in each neighborhood/month
# calculate total number of males (with recorded pattern) seen in each neighborhood/month

rneigh2 <- rneigh1[N>0, .(n_neigh_patterns = length(pattern),
                          n_neigh_males = sum(N)),
                   by=c("sampling", "neighborhood")]


# merge location/month counts with location/time pattern frequencies 
rneigh3 <- merge(x=rneigh1, y=rneigh2,
                 by = c("sampling", "neighborhood"),
                 all.y = TRUE)

# relative frequency of each pattern in location/time
rneigh3$neigh_rel_freq <- rneigh3$N/rneigh3$n_neigh_males

rneigh3$rareness_neigh <- log(rneigh3$neigh_rel_freq * rneigh3$n_neigh_patterns)

rneigh3 <- rneigh3[is.finite(rareness_neigh),]

rneigh3$pattern <- as.factor(rneigh3$pattern)
rneigh3$sampling <- as.integer(rneigh3$sampling)

rneigh3<- rneigh3[, c("sampling", "pattern", "neighborhood", "rareness_neigh")]

data <- merge(x=data, y=rneigh3,
              by = c("pattern", "sampling", "neighborhood"),
              all.x=TRUE)

##################################################

# C. total population

# How "rare" is each male in the whole population in each month?

# count data: frequency of each pattern in each month (sampling) only
rpop1 <- data[sex_stage=="M", as.data.table(table(pattern, sampling))]

# calculate total number of patterns seen in each location/time
# calculate total number of males (with recorded pattern) seen in each location/time

rpop2 <- rpop1[N>0, .(n_pop_patterns = length(pattern),
                      n_pop_males = sum(N)),
               by=c("sampling")]


# merge location/time counts with location/time pattern frequencies 
rpop3 <- merge(x=rpop1, y=rpop2,
               by = c("sampling"),
               all.y = TRUE)

# relative frequency of each pattern in location/time
rpop3$pop_rel_freq <- rpop3$N/rpop3$n_pop_males

rpop3$rareness_pop <- log(rpop3$pop_rel_freq * rpop3$n_pop_patterns)

rpop3 <- rpop3[is.finite(rareness_pop),]

rpop3$pattern <- as.factor(rpop3$pattern)
rpop3$sampling <- as.integer(rpop3$sampling)

rpop3<- rpop3[, c("sampling", "pattern", "rareness_pop")]

data <- merge(x=data, y=rpop3,
              by = c("pattern", "sampling"),
              all.x=TRUE)

##################################################
# 4. Quantifying reproductive success -------------------------------------------------

# We quantify 5 measures of reproductive success:

# These are specific to individuals:

# 1. Number of viable offspring (i.e. that ultimately recruit into the population) conceived per month
# 2. Whether an individual mated in a given month 
# 3. The number of partners it mated with per month

# These are per union/mating, and so the same result applies to
# both the female and the male in a given union:

# 4. Number of offspring per union
# 5. Number of grandoffspring ultimately produced per union

##################################################

# NUMBER OF OFFSPRING PER MONTH

# first return a datatable with the conception month and momID for each fish
conception <- data[, .(conceived = conceived[1],
                       momID = momID[1]),
                   by=FishID]

# next return a count of number of offspring for each mum in each month
m_rep <- as.data.table(table(conception$momID, conception$conceived)) 

# switcheroo - change "momID" to "FishID" so we can merge it back in
names(m_rep) <- c("FishID", "sampling", "n_offspring_conceived")
m_rep$sampling<- as.numeric(m_rep$sampling)

# now do the same for dadID for each fish
conception <- data[, .(conceived = conceived[1],
                       dadID = dadID[1]),
                   by=FishID]

d_rep <- as.data.table(table(conception$dadID, conception$conceived)) 

# switcheroo - change "dadID" to "FishID", and "conceived" to "sampling" so we can merge it back in
names(d_rep) <- c("FishID", "sampling", "n_offspring_conceived")
d_rep$sampling<- as.numeric(d_rep$sampling)

# reproductive success of males and females
p_rep <- as.data.table(rbind(m_rep, d_rep))

# merge back in to mark recapture data
data <- merge(x=data, y=p_rep,
              by=c("FishID","sampling"),
              all.x=TRUE)

# individuals that did not have any offspring have NA: replace NA for reproduction with zero
data$n_offspring_conceived[which(is.na(data$n_offspring_conceived))] <- 0 

# pedigree only goes up to month 60, which means max conception sampling month = 57

data$n_offspring_conceived[which(data$sampling>57)] <- NA

##################################################

# DID THEY MATE IN A GIVEN MONTH?
#logical indicator - did an individual breed that month? TRUE=yes, FALSE=no
data$bred <- data$n_offspring_conceived>0

##################################################

# NUMBER OF MATING PARTNERS PER MONTH

# this is just the pedigree with the date of conception of individuals
matings <- data[, .(conceived = conceived[1],
                    momID = momID[1],
                    dadID = dadID[1]),
                by=FishID]

# this counts the number of males each female mated with each month
mum_matings <- matings[, length(unique(dadID)), by=c("momID", "conceived")]
mum_matings <- mum_matings[!is.na(momID),]
mum_matings <- mum_matings[momID!="founder",]

names(mum_matings) <- c("FishID", "sampling", "n_mates")

# now same for dads
dad_matings <- matings[, length(unique(momID)), by=c("dadID", "conceived")]
dad_matings <- dad_matings[!is.na(dadID),]
dad_matings <- dad_matings[dadID!="founder",]
names(dad_matings) <- c("FishID", "sampling", "n_mates")

mating_data <- rbind(mum_matings, dad_matings)

data <- merge(data, mating_data,
              by=c("FishID", "sampling"),
              all.x = TRUE)

# set individuals with no mates from NA to 0
data$n_mates[is.na(data$n_mates)]<-0
data$n_mates[which(data$sampling>57)]<-NA

##################################################

# NUMBER OF OFFSPRING & GRANDOFFSPRING PER UNION

# first, get a list of lifetime measurements for each individual
measured_once <- data[, .( dadID = dadID[1],
                           momID = momID[1],
                           LRS =sum(n_offspring_conceived), # we use this to calculate number of grandoffspring later
                           conceived = conceived[1],
                           last_recorded = last_recorded[1],
                           sampling.year = sampling.year[1],
                           sex_stage = sex_stage[1]),
                      by="FishID"]

# N.B. LRS will not be accurate for those born later (i.e. their offspring not yet counted / included in the pedigree)
# Change LRS values to NA for individuals observed at any point in the last three months

measured_once$LRS[which(measured_once$last_recorded>57)] <- NA

# add a "union" identifier, ie. momID + dadID + date
measured_once$union <- paste(measured_once$dadID, 
                             measured_once$momID, 
                             measured_once$conceived, sep="_")

# we can't measure offspring/grandoffspring where we don't know parental identity:
unions <- measured_once[dadID!="founder" & momID!="founder",]
unions <- unions[!is.na(dadID) & !is.na(momID),]

# Now we return the number of offspring, grandoffspring, identity of mating partners,
# and date of mating for each union:
unions<- unions[, .( off_per_union = length(unique(FishID)), # counts number of offspring per union
                     goff_per_union = sum(LRS), # sum of the LRS of each offspring per union
                     momID = momID[1],
                     dadID = dadID[1],
                     union_year = sampling.year[1],
                     sampling = conceived[1]), 
                by=union]

# add in month (1:12) of union
union_months <- data[, sampling.month[1], by=sampling]
names(union_months)<- c("sampling", "union_month")

unions <- merge(x=unions, y=union_months, by="sampling", all.x=TRUE)

# add in dad rareness and novelty at each union
dad_rareness <- data[, c("FishID", "sampling", 
                         "rareness_local", "rareness_neigh", "rareness_pop", "pattern")]

# switch "FishID" to "dadID" so it connects the monthly rarity score of fathers to union
names(dad_rareness)<- c("dadID","sampling", 
                        "dad_rareness_local", "dad_rareness_neigh", "dad_rareness_pop",
                        "dad_pattern")

unions <- merge(x=unions[!is.na(dadID),], y=dad_rareness,
                by=c("dadID", "sampling"),
                all.x=TRUE)

unions<- unions[!is.na(dad_rareness_local),]

##################################################
# 5. Parent-offspring relationships -------------------------------------------------

# We now need to link paternal traits (pattern, rareness) to each offspring.

# Additionally, we calculate the inbreeding coefficient of each individual

##################################################

# LINKING PATERNAL TRAITS TO OFFSPRING


# First, lets link the father's rareness at conception to each offspring
# individuals, when they were conceived, and the identity of their father
dad_rareness1 <- measured_once[, c("FishID","conceived","dadID", "momID")]
dad_rareness1$union <- paste(dad_rareness1$dadID, 
                             dad_rareness1$momID, 
                             dad_rareness1$conceived, sep="_")

# # re-use dad_rareness object from earlier
names(dad_rareness)

# rename columns then merge these together
names(dad_rareness)[2]<- "conceived"


# this contains the rareness at time of mating of each individual's dad 
dad_rareness <- merge(x=dad_rareness, y=dad_rareness1,
                      by=c("dadID", "conceived"),
                      all.x=TRUE)

dad_rareness <- dad_rareness[, c("FishID", "dadID",
                                 "dad_rareness_local", "dad_rareness_neigh","dad_rareness_pop",
                                 "dad_pattern", "union")]

grandad_rareness <- dad_rareness
names(grandad_rareness)<- c("dadID", "grampsID", "gramps_rareness_local", "gramps_rareness_neigh", "gramps_rareness_pop",       
                            "gramps_pattern", "gramps_union") 

dad_gramps_rareness <- merge(x = dad_rareness, y=grandad_rareness, by="dadID", all.x=TRUE)
dad_gramps_rareness <- dad_gramps_rareness[,-"dadID"]

# merge into main data
data <- merge(x=data, y= dad_gramps_rareness,
              by="FishID", all.x = TRUE)


##################################################

# add movement data to unions

dad_moved <- data[sex_stage=="M", c("FishID", "sampling", "loc_arrived", "neigh_arrived")]

names(dad_moved)[1]<- "dadID"

unions <- merge(x = unions, y=dad_moved,
                by=c("dadID", "sampling"),
                all.x=TRUE)

names(dad_moved)<- c("dadID", "conceived", "loc_arrived_dad", "neigh_arrived_dad")

data <- merge(x = data, y=dad_moved,
              by=c("dadID", "conceived"),
              all.x=TRUE)


# add the date of paternal conception to each row

con_by_id <- data[, .(conceived = conceived[1]), by=FishID]
pat_con <- con_by_id
names(pat_con)<- c("dadID", "dad_conceived")

data <- merge(x=data, y=pat_con,
              by="dadID", all.x=TRUE)

##################################################

write.csv(file = "NFDS_guppy_data.csv", data, row.names = FALSE)
write.csv(file = "NFDS_guppy_unions.csv", unions, row.names = FALSE)


