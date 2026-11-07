library(rhdf5)
library(dplyr)

## User inputs----

use_statbel <- F

use_R0 <- T

# disable infectious from certain year
disable_infections <- F
disabling_year <- 2020

## Population data ----
# Reading population data (2021) statbel
# population_full_data <- cbind(h5read("../data/population.h5", "2021/population/full"))
population_flanders_data <- cbind(h5read("./data/population.h5", "2021/population/flanders"))
population_wal_bxl_data <- cbind(h5read("./data/population.h5", "2021/population/wallonia_brussels"))
colnames(population_flanders_data) <- c("age", "2021")
colnames(population_wal_bxl_data) <- c("age", "2021")

## Mortality data ----

mortality_flanders_data <- cbind(h5read("./data/mortality.h5", "statbel/2021/flanders/mortality_extrapolated"))
mortality_wal_data <- cbind(h5read("./data/mortality.h5", "statbel/2021/wallonia/mortality_extrapolated"))
mortality_bxl_data <- cbind(h5read("./data/mortality.h5", "statbel/2021/brussels/mortality_extrapolated"))
mortality_wal_bxl_data <- mortality_wal_data+mortality_bxl_data

mortality_rate_be <- h5read("./data/mortality.h5", "plan/2021") ## age specific mortality rate

region <- NA
mortality_data <- NA

if(use_flanders_data=="Flanders") {
  population_data <- population_flanders_data[,'2021']
  region <- "flanders"
  if(use_statbel) {
    mortality_data <- mortality_flanders_data
  } else {
    mortality_data <- mortality_rate_be
  }
} else if (use_flanders_data=="Wallonia") {
  population_data <- population_wal_bxl_data[,'2021']
  region <- "walbxl"
  if(use_statbel) {
    mortality_data <- mortality_wal_bxl_data
  } else {
    mortality_data <- mortality_rate_be
  }
} else if (use_flanders_data=="Belgium") {
  region <- "belgium"
  population_data <- population_flanders_data[,'2021'] + population_wal_bxl_data[,'2021']
  mortality_data <- mortality_rate_be

} else {
  stop("Misspecified region, check string in use_flanders_data")
}

sum(population_data)
length(population_data)
nrow(mortality_data)

mortality_data <- c(mortality_data[1], # 0-6 months
                    mortality_data[1], # 6-9 months
                    mortality_data[1], # 9-12 months
                    mortality_data[2:100])
population_data <- c(population_data[1], # 0-6 months
                     0, # 6-9 months
                     0, # 9-12 months
                     population_data[c(2:99)], population_data[100]+population_data[101])

sum(population_data)

options(repr.plot.width = 5, repr.plot.height = 5)

plot(population_data, ylab = "number of individuals")

options(repr.plot.width = 5, repr.plot.height = 5)

plot(mortality_data, xlab = "Age group", ylab="Probability of dying (per year)")

## Contact data ----
# Reading contact matrix from Socrates
contact_rates_small <- h5read("./data/contact-data-socrates.h5", "socrates/original/matrices_per_capita/total")
rownames(contact_rates_small) <- h5read("./data/contact-data-socrates.h5", "socrates/original/matrices_per_capita/agegroup")
colnames(contact_rates_small) <- h5read("./data/contact-data-socrates.h5", "socrates/original/matrices_per_capita/contact_agegroup")

contact_rates_full <- h5read("./data/contact-data-socrates.h5", "socrates/full/matrices_per_capita/total")
rownames(contact_rates_full) <- h5read("./data/contact-data-socrates.h5", "socrates/full/matrices_per_capita/agegroup")
colnames(contact_rates_full) <- h5read("./data/contact-data-socrates.h5", "socrates/full/matrices_per_capita/contact_agegroup")

contact_rates_min15m <- h5read("./data/contact-data-socrates.h5", "socrates/min15m/matrices_per_capita/total")
rownames(contact_rates_min15m) <- h5read("./data/contact-data-socrates.h5", "socrates/min15m/matrices_per_capita/agegroup")
colnames(contact_rates_min15m) <- h5read("./data/contact-data-socrates.h5", "socrates/min15m/matrices_per_capita/contact_agegroup")

contact_rates_plus15m <- h5read("./data/contact-data-socrates.h5", "socrates/plus15m/matrices_per_capita/total")
rownames(contact_rates_plus15m) <- h5read("./data/contact-data-socrates.h5", "socrates/plus15m/matrices_per_capita/agegroup")
colnames(contact_rates_plus15m) <- h5read("./data/contact-data-socrates.h5", "socrates/plus15m/matrices_per_capita/contact_agegroup")

contact_rates_plus1h <- h5read("./data/contact-data-socrates.h5", "socrates/plus1h/matrices_per_capita/total")
rownames(contact_rates_plus1h) <- h5read("./data/contact-data-socrates.h5", "socrates/plus1h/matrices_per_capita/agegroup")
colnames(contact_rates_plus1h) <- h5read("./data/contact-data-socrates.h5", "socrates/plus1h/matrices_per_capita/contact_agegroup")

contact_rates_plus4h <- h5read("./data/contact-data-socrates.h5", "socrates/plus4h/matrices_per_capita/total")
rownames(contact_rates_plus4h) <- h5read("./data/contact-data-socrates.h5", "socrates/plus4h/matrices_per_capita/agegroup")
colnames(contact_rates_plus4h) <- h5read("./data/contact-data-socrates.h5", "socrates/plus4h/matrices_per_capita/contact_agegroup")

date_run <- as.integer(format(Sys.Date(), "%Y%m%d"))

date_dir <- NA
cr_type <- NA
if(use_cr_full) {
  contact_rates <- contact_rates_full
  date_dir <- sprintf("%d-%s-pl-full", date_run, region)
  cr_type <- "full"
} else if(use_cr_min15m) {
  contact_rates <- contact_rates_min15m
  date_dir <- sprintf("%d-%s-pl-min15m", date_run, region)
  cr_type <- "min15m"
} else if(use_cr_plus15m) {
  contact_rates <- contact_rates_plus15m
  date_dir <- sprintf("%d-%s-pl-plus15m", date_run, region)
  cr_type <- "plus15m"
} else if(use_cr_plus1h) {
  contact_rates <- contact_rates_plus1h
  date_dir <- sprintf("%d-%s-pl-plus1h", date_run, region)
  cr_type <- "plus1h"
} else if(use_cr_plus4h) {
  contact_rates <- contact_rates_plus4h
  date_dir <- sprintf("%d-%s-pl-plus4h", date_run, region)
  cr_type <- "plus4h"
}

date_dir

head(contact_rates)

tail(contact_rates)

df <- contact_rates

# original row to split
row_vals <- df["[0,1)", , drop = FALSE]

# create two copies
row1 <- row_vals
row2 <- row_vals
row3 <- row_vals
rownames(row1) <- "[0,6months)"
rownames(row2) <- "[6,9months)"
rownames(row3) <- "[9,12months)"

# add them at the top
df <- rbind(row1, 
            row2, 
            row3, 
            df)

# add new columns
df <- cbind("[0,6months)" = df[,"[0,1)"],
            "[6,9months)" = df[,"[0,1)"],
            "[9,12months)" = df[,"[0,1)"],
            df)

# remove original group
df <- df[!(rownames(df) %in% "[0,1)"), !(colnames(df) %in% "[0,1)")]

contact_rates <- df

head(contact_rates)

num_age_groups <- nrow(contact_rates)
num_age_groups

# based on the values, these are per-capita contact rates 
max(contact_rates)

contact_rates[3,10]
contact_rates[10,3]

# the matrix is symmetric
isSymmetric(contact_rates_full)

m <- contact_rates[1:5, 1:5]

# Define colors
colfunc <- colorRampPalette(c("steelblue", "white"))
age_groups <- rownames(m) 
num <- length(age_groups)

image(1:num, 1:num, m,
      col = rev(colfunc(1e2)),
      xlab = "Participant age group",
      ylab = "Contact age group",
      axes = FALSE, main = "Contact matrix", frame = TRUE)
axis(1, 1:num, age_groups)
axis(2, 1:num, age_groups)
for(i in 1:num) for (j in 1:num) text(i, j, round(m[i, j], 6))

m <- contact_rates

# Define colors
colfunc <- colorRampPalette(c("steelblue", "white"))
age_groups <- rownames(m) 
num <- length(age_groups)

image(1:num, 1:num, m,
      col = rev(colfunc(1e2)),
      xlab = "Participant age group",
      ylab = "Contact age group",
      axes = FALSE, main = "Contact matrix", frame = TRUE)
axis(1, 1:num, age_groups)
axis(2, 1:num, age_groups)
#for(i in 1:num_age_groups) for (j in 1:num_age_groups) text(i, j, round(m[i, j], 6))

#res <- read.csv("./data/coverage_glm_fit_2020.csv")
res <- read.csv("./data/coverage_glm_fit_2020_ci_bootstrap.csv")

head(res)

## Vaccination coverage ########

subtitle <- NA

if(use_flanders_data=="Flanders") {

  region_name <- "Flanders"

  dose_name <- "mmr1"
  mmr1_glm <- res %>%
  filter(region == region_name,
         dose == dose_name)
  
  dose_name <- "mmr2"
  mmr2_glm <- res %>%
  filter(region == region_name,
         dose == dose_name)
  
  subtitle <- "Flanders"

} else if (use_flanders_data=="Wallonia") {

  region_name <- "Wallonia"

  dose_name <- "mmr1"
  mmr1_glm <- res %>%
  filter(region == region_name,
         dose == dose_name)
  
  dose_name <- "mmr2"
  mmr2_glm <- res %>%
  filter(region == region_name,
         dose == dose_name)
  
  subtitle <- "Flanders"

  subtitle <- "Wallonia - Brussels"
} else if (use_flanders_data=="Belgium") {

  region_name <- "FlanWall_weighted_avg" 

  dose_name <- "mmr1"
  mmr1_glm <- res %>%
  filter(region == region_name,
         dose == dose_name)
  
  dose_name <- "mmr2"
  mmr2_glm <- res %>%
  filter(region == region_name,
         dose == dose_name)

  subtitle <- "Belgium"

} else {
  stop("Misspecified region, check string in use_flanders_data")
}

cov1_means <- as.numeric(mmr1_glm$coverage)
cov2_means <- as.numeric(mmr2_glm$coverage)

plot(1985+(0:(length(cov1_means)-1)), cov1_means, xlab = "year", ylab = "vax coverage")
points(1995+(0:(length(cov2_means)-1)), cov2_means)

## Parameterisation #####
gamma <- 1/(365*2.8/12) #1/182.5  # M -> S (6 months = 182.5 days -> gamma = 1/#days = 1/182.5)
sigma <- 1/10  # E -> I, VE -> VI    #1/11        # naar WHO kijken, london school
nu <- 1/8     # I -> R, VI -> VR     #1/9        # from parameterization of Robert and Suffel (papers by Kucharski)
#tau_by_age <- list("2nd dose" = 0.96,
#            "1-2y" = 0.96,
#            "9-12m" = 0.85,
#            "6-9m" = 0.76) # S -> VS, S -> VR, R-> VR       # chance of succesfull vaccination
#waning <-  0.009/365 #   VR -> VS # Everything is expressed in rates and days. This states that one keeps an immunity for 85 years. (Life-long immunity)

N <- sum(population_data)
# mortality rate
mu <- NA
if(use_statbel) {
  mu <- mortality_data/(N*365)    #  time dependent death rate, not age specific
} else {
  mu <- mortality_data/365
}

max(Re(eigen(diag(population_data * sigma/(mu+sigma) * 1./(nu+mu)) %*% contact_rates)$values))

cr <- contact_rates / max(Re(eigen(diag(population_data * sigma/(mu+sigma) * 1./(nu+mu)) %*% contact_rates)$values)) 

isSymmetric(cr)
max(cr)

cr_alternative <- contact_rates * c(((mu+sigma)/sigma)*(nu+mu)) / max(Re(eigen(diag(population_data) %*% contact_rates)$values)) 

isTRUE(all.equal(cr, cr_alternative, tolerance = 1e-7))

isTRUE(all.equal(cr, cr_alternative, tolerance = 1e-8))

## Timesteps for model simulation, RAS and vaccination #####

# Years to simulate
yts_begin <- 1940
yts_end <- 2013

# Vaccination day is on 2nd of September, every year
vacc_month <- 9
vacc_day <- 1
first_sept_per_year <- calculate_day_of_vaccination_for_years(yts_begin, yts_end, vacc_month, vacc_day)

vacc_month <- 3
vacc_day <- 1
first_march_per_year <- calculate_day_of_vaccination_for_years(yts_begin, yts_end, vacc_month, vacc_day)

vacc_month <- 6
vacc_day <- 1
first_june_per_year <- calculate_day_of_vaccination_for_years(yts_begin, yts_end, vacc_month, vacc_day)

# NOTE: 2nd of september!! Need to figure out how to have events on the same day...
days_of_vaccination <- cbind(first_sept_per_year+1, yts_begin:yts_end)

days_of_vaccination_6months <- cbind(first_march_per_year+1, yts_begin:yts_end)

days_of_vaccination_9months <- cbind(first_june_per_year+1, yts_begin:yts_end)

# Day everyone shifts/ages (RAS model)
# days_of_year_change <- calculate_jan1_for_years(yts_begin, yts_end)[-1] # remove first entry
days_of_year_change <- first_sept_per_year

days_of_6months_change <- first_march_per_year

days_of_9months_change <- first_june_per_year

# Time series => output waar we intresse in hebben, elke dag status van alles states 
t <- 1:calculate_days_over_years(yts_begin, yts_end)

days_of_year_change
days_of_6months_change
days_of_9months_change

# start of vaccination
years <- 1985:2013

age_cat_to_vacc <- list()
for(i in 1:length(years)) {
  year <- years[i]
  if(year >= 1995) {
    idx <- i - 10

    if (use_flanders_data == "Flanders") {

      # we currently assume that children receive their dose at the start of the school year (5th grade, 10-11 y.o.) therefore at 10 y.o.
      age_cat_to_vacc[[as.character(year)]] <- rbind(c(1,cov1_means[i], 'MMR1'), 
                                                     c(10,cov2_means[idx],'MMR2'))

      if (year == 1995) {
          if (catch_up_d2_1995 == T) {
          age_cat_to_vacc[[as.character(year)]] <- rbind(c(1,cov1_means[i], 'MMR1'), 
                                                         c(10,cov2_means[idx],'MMR2'),
                                                         c(11,cov2_means[idx],'MMR2'))              
          }
      }

    } else if (use_flanders_data == "Wallonia") {

      # we currently assume that children receive their dose at the start of the school year (6th grade, 11-12 y.o.) therefore at 11 y.o.
      age_cat_to_vacc[[as.character(year)]] <- rbind(c(1,cov1_means[i], 'MMR1'), 
                                                     c(11,cov2_means[idx], 'MMR2'))

    } else if (use_flanders_data == "Belgium") {

      # we assume 11 y.o., in betweem the age groups of 10 and 12 (Flanders-Wallonia) 
      age_cat_to_vacc[[as.character(year)]] <- rbind(c(1,cov1_means[i], 'MMR1'), 
                                                     c(11,cov2_means[idx],'MMR2'))    

    } else {
      stop("Misspecified region, check string in use_flanders_data")
    }
    
  } else {
    age_cat_to_vacc[[as.character(year)]] <- c(1,cov1_means[i], 'MMR1')
  }
}

age_cat_to_vacc["1985"]
age_cat_to_vacc["1995"]
age_cat_to_vacc["2013"]

## initial condition

# 100 age categories (0-99 years) - or more if we stratify 0-1 y.o. with finer resolution
# [,c(3, 4, 7, 8, 9)]) # I, R, VI, VR, M ---> [,c(3, 4, 5, 8, 9, 10, 13, 14, 15)]) # I, R, VP, VI, VR, VVP, VVI, VVR, M

S_r <- 1:num_age_groups 
E_r <- num_age_groups + 1:num_age_groups
I_r <- 2*num_age_groups + 1:num_age_groups
R_r <- 3*num_age_groups + 1:num_age_groups

VP_r <- 4*num_age_groups + 1:num_age_groups
VS_r <- 5*num_age_groups + 1:num_age_groups
VE_r <- 6*num_age_groups + 1:num_age_groups
VI_r <- 7*num_age_groups + 1:num_age_groups
VR_r <- 8*num_age_groups + 1:num_age_groups

VVP_r <- 9*num_age_groups + 1:num_age_groups
VVS_r <- 10*num_age_groups + 1:num_age_groups
VVE_r <- 11*num_age_groups + 1:num_age_groups
VVI_r <- 12*num_age_groups + 1:num_age_groups
VVR_r <- 13*num_age_groups + 1:num_age_groups

M_r <- 14*num_age_groups + 1

NewlyInfectious_r <- 14*num_age_groups + 1 + 1:num_age_groups # must be computed during the simulation; has no influence on the simulation itself, it is only for bookkeeping
NewlyInfectious_dose0 <- 15*num_age_groups + 1 + 1:num_age_groups
NewlyInfectious_dose1 <- 16*num_age_groups + 1 + 1:num_age_groups
NewlyInfectious_dose2 <- 17*num_age_groups + 1 + 1:num_age_groups


age_cat_per_compartment <- list(S_r=S_r, E_r=E_r, I_r=I_r, R_r=R_r,
                                VP_r=VP_r, VS_r=VS_r, VE_r=VE_r, VI_r=VI_r, VR_r=VR_r, 
                                VVP_r=VVP_r, VVS_r=VVS_r, VVE_r=VVE_r, VVI_r=VVI_r, VVR_r=VVR_r,
                                M_r=M_r,
                                NewlyInfectious_r=NewlyInfectious_r,
                                NewlyInfectious_dose0=NewlyInfectious_dose0,
                                NewlyInfectious_dose1=NewlyInfectious_dose1, 
                                NewlyInfectious_dose2=NewlyInfectious_dose2)

tot_comp <- age_cat_per_compartment$NewlyInfectious_dose2[length(age_cat_per_compartment$NewlyInfectious_dose2)]
tot_comp

population_data

y_init <- rep(0, tot_comp)
y_init[S_r[-1]] <- population_data[2:num_age_groups] # 1-99 years; S_r[-1] == 2:100
y_init[M_r] <- population_data[1]     # 0 years

# remove expected deaths from 1 sept 1939 to 1 January 1940
y_init[S_r] <- y_init[S_r]*(1-mu*(30*2+31*2))
y_init[M_r] <- y_init[M_r]*(1-mu[1]*(30*2+31*2))

#infected seed
shift <- num_age_groups-100
y_init[S_r[2+shift]] <- y_init[S_r[2+shift]] - 1
y_init[I_r[2+shift]] <- 1 

calculate_hessian <- T

#survey_years <- unique(df_surveys$data_year)
survey_years <- c(2002, 2006, 2013)


sensitivity <- 0.996
specificity <- 1.0
sero_day <- 1
sero_month <- 2 # we compare the data to the model output on february 1st (as an example)

days_of_6months_change
