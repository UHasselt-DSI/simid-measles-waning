library(deSolve)

tau_by_age <- list("2nd dose" = 0.96,
                   "2nd dose, immune blunting after 1st dose at 9 months" = 0.9,
                   "2nd dose, immune blunting after 1st dose at 6 months" = 0.85,
                   "1-2y" = 0.96,
                   "9-12m" = 0.85, 
                   "6-9m" = 0.76)

# the main function simulateModel: contains two types of functions:
# 1) events_function and 2) model.
# to be used in ode 
#
# events_function:
#   1) shifts everyone via RAS == RAS shift
#      integrated over time; t is the moment of integration at specific event times
#      days_of_year_change => on September 1
#      days_of_vaccination => on September 2
#   2) zombies come to life (newborns and deaths)
#   3) vaccination coverage is now year-dependent at age 1 and 12
#
#  

# as a reminder
# initial_params <- c(num_age_groups, gamma, sigma, nu, tau, waning, mu, N, as.vector(cr), R0)
simulateModel <- function(ts, y_init, initial_params, 
                          days_of_6months_change, 
                          days_of_9months_change, 
                          days_of_year_change, 
                          days_of_vaccination_6months,
                          days_of_vaccination_9months,
                          days_of_vaccination,
                          age_cat_to_vacc,
                          verbose=FALSE, 
                          disable_infections=FALSE, t_disabling_year=NA, 
                          method = "ode45", save_newlyInf=FALSE,
                          introductions = FALSE, start_introductions = NA, n_intro = 10, 
                          VE_transm = 0,
                          immune_blunting = FALSE, 
                          gamma_t = NULL) {  
  
  override_infections <- F
  
  ras_shift_6months <- function(t, states, params) { 
    # the function takes as input the current time, the current condition of compartments (states) and the transition params
    if(verbose)
      print(sprintf("shift day: %d", t))

    # those aged 0-6 months move to the next compartment 6-9 months
    S_idx <- S_r[1]
    E_idx <- E_r[1]
    I_idx <- I_r[1]
    R_idx <- R_r[1]
    
    VP_idx <- VP_r[1]
    VS_idx <- VS_r[1]
    VE_idx <- VE_r[1]
    VI_idx <- VI_r[1]
    VR_idx <- VR_r[1]
    
    VVP_idx <- VVP_r[1]
    VVS_idx <- VVS_r[1]
    VVS_idx <- VVS_r[1]
    VVE_idx <- VVE_r[1]
    VVI_idx <- VVI_r[1]
    VVR_idx <- VVR_r[1]

    states[S_idx+1] <- states[S_idx] ### S[2] = S[1], ..., S[99] = S[98] 
    states[E_idx+1] <- states[E_idx]
    states[I_idx+1] <- states[I_idx]
    states[R_idx+1] <- states[R_idx]
    
    states[VP_idx+1] <- states[VP_idx]
    states[VS_idx+1] <- states[VS_idx]
    states[VE_idx+1] <- states[VE_idx]
    states[VI_idx+1] <- states[VI_idx]
    states[VR_idx+1] <- states[VR_idx]
    
    states[VVP_idx+1] <- states[VVP_idx]
    states[VVS_idx+1] <- states[VVS_idx]
    states[VVE_idx+1] <- states[VVE_idx]
    states[VVI_idx+1] <- states[VVI_idx]
    states[VVR_idx+1] <- states[VVR_idx]
    
    states[S_r[2]] <- states[S_r[2]] + states[M_r] 
    #if some are still in the compartment of maternal immunity, move them to the age2 compartment anyway

    states[S_r[1]] <- 0 # newborns compartments are first set all to zeros
    states[E_r[1]] <- 0
    states[I_r[1]] <- 0
    states[R_r[1]] <- 0
    
    states[VP_r[1]] <- 0
    states[VS_r[1]] <- 0
    states[VE_r[1]] <- 0
    states[VI_r[1]] <- 0
    states[VR_r[1]] <- 0
    
    states[VVP_r[1]] <- 0
    states[VVS_r[1]] <- 0
    states[VVE_r[1]] <- 0
    states[VVI_r[1]] <- 0
    states[VVR_r[1]] <- 0
    
    states[M_r] <- 0
        
    return(states)
  }

  ras_shift_9months <- function(t, states, params) { 
    # the function takes as input the current time, the current condition of compartments (states) and the transition params
    if(verbose)
      print(sprintf("shift day: %d", t))

    # those aged 6-9 months move to the next compartment 9-12 months
    S_idx <- S_r[2]
    E_idx <- E_r[2]
    I_idx <- I_r[2]
    R_idx <- R_r[2]
    
    VP_idx <- VP_r[2]
    VS_idx <- VS_r[2]
    VE_idx <- VE_r[2]
    VI_idx <- VI_r[2]
    VR_idx <- VR_r[2]
    
    VVP_idx <- VVP_r[2]
    VVS_idx <- VVS_r[2]
    VVS_idx <- VVS_r[2]
    VVE_idx <- VVE_r[2]
    VVI_idx <- VVI_r[2]
    VVR_idx <- VVR_r[2]

    states[S_idx+1] <- states[S_idx] ### S[2] = S[1], ..., S[99] = S[98] 
    states[E_idx+1] <- states[E_idx]
    states[I_idx+1] <- states[I_idx]
    states[R_idx+1] <- states[R_idx]
    
    states[VP_idx+1] <- states[VP_idx]
    states[VS_idx+1] <- states[VS_idx]
    states[VE_idx+1] <- states[VE_idx]
    states[VI_idx+1] <- states[VI_idx]
    states[VR_idx+1] <- states[VR_idx]
    
    states[VVP_idx+1] <- states[VVP_idx]
    states[VVS_idx+1] <- states[VVS_idx]
    states[VVE_idx+1] <- states[VVE_idx]
    states[VVI_idx+1] <- states[VVI_idx]
    states[VVR_idx+1] <- states[VVR_idx]
    
    states[S_r[2]] <- 0 # newborns compartments are first set all to zeros
    states[E_r[2]] <- 0
    states[I_r[2]] <- 0
    states[R_r[2]] <- 0

    states[VP_r[2]] <- 0
    states[VS_r[2]] <- 0
    states[VE_r[2]] <- 0
    states[VI_r[2]] <- 0
    states[VR_r[2]] <- 0
    
    states[VVP_r[2]] <- 0
    states[VVS_r[2]] <- 0
    states[VVE_r[2]] <- 0
    states[VVI_r[2]] <- 0
    states[VVR_r[2]] <- 0
        
    return(states)
  }
    
  ras_shift <- function(t, states, params) { 
    # the function takes as input the current time, the current condition of compartments (states) and the transition params
    if(verbose)
      print(sprintf("shift day: %d", t))
    
    N <- params[num_age_groups + 7]

    # Everyone becomes one year older
    # Compute indices for people that need to be shifted, for example
    num_age_cats_to_shift <- length(S_r)-1 # 100 - 1 = 99
    # people aged 0 to 98 are shifted to 1 and 99; old 99 move out (deaths)

    S_idx <- S_r[3:num_age_cats_to_shift]
    E_idx <- E_r[3:num_age_cats_to_shift]
    I_idx <- I_r[3:num_age_cats_to_shift]
    R_idx <- R_r[3:num_age_cats_to_shift]
    
    VP_idx <- VP_r[3:num_age_cats_to_shift]
    VS_idx <- VS_r[3:num_age_cats_to_shift]
    VE_idx <- VE_r[3:num_age_cats_to_shift]
    VI_idx <- VI_r[3:num_age_cats_to_shift]
    VR_idx <- VR_r[3:num_age_cats_to_shift]
    
    VVP_idx <- VVP_r[3:num_age_cats_to_shift]
    VVS_idx <- VVS_r[3:num_age_cats_to_shift]
    VVS_idx <- VVS_r[3:num_age_cats_to_shift]
    VVE_idx <- VVE_r[3:num_age_cats_to_shift]
    VVI_idx <- VVI_r[3:num_age_cats_to_shift]
    VVR_idx <- VVR_r[3:num_age_cats_to_shift]

    states[S_idx+1] <- states[S_idx] ### S[2] = S[1], ..., S[99] = S[98] 
    states[E_idx+1] <- states[E_idx]
    states[I_idx+1] <- states[I_idx]
    states[R_idx+1] <- states[R_idx]
    
    states[VP_idx+1] <- states[VP_idx]
    states[VS_idx+1] <- states[VS_idx]
    states[VE_idx+1] <- states[VE_idx]
    states[VI_idx+1] <- states[VI_idx]
    states[VR_idx+1] <- states[VR_idx]
    
    states[VVP_idx+1] <- states[VVP_idx]
    states[VVS_idx+1] <- states[VVS_idx]
    states[VVE_idx+1] <- states[VVE_idx]
    states[VVI_idx+1] <- states[VVI_idx]
    states[VVR_idx+1] <- states[VVR_idx]
    
    #states[S_r[2]] <- states[S_r[2]] + states[M_r] 
    #if some are still in the compartment of maternal immunity, move them to the age2 compartment anyway
    #underlying assumption: maternal immunity only last at maximum 1 year (average is 2.8 months so it is reasonable)

    states[S_r[3]] <- 0 # newborns compartments are first set all to zeros
    states[E_r[3]] <- 0
    states[I_r[3]] <- 0
    states[R_r[3]] <- 0

    states[VP_r[3]] <- 0
    states[VS_r[3]] <- 0
    states[VE_r[3]] <- 0
    states[VI_r[3]] <- 0
    states[VR_r[3]] <- 0
    
    states[VVP_r[3]] <- 0
    states[VVS_r[3]] <- 0
    states[VVE_r[3]] <- 0
    states[VVI_r[3]] <- 0
    states[VVR_r[3]] <- 0
    
    states[M_r] <- 0
    
    if(introductions && t>=start_introductions) { 
      # put introductions equally spaced the 99 oldest groups (i.e. index from 2 to 100, excluding those aged 0-1)
      for (g in round(seq(num_age_groups-98, num_age_groups, length.out = n_intro))) {
        #print(c(g, I_r[g]))
        states[I_r[g]] <- states[I_r[g]] + 1
        states[S_r[g]] <- states[S_r[g]] - 1 
      }
    }

    # Enforce constant population size, aka zombies
    # we compute current total number in all compartments
    sum_states <- sum(states[S_r], states[E_r], states[I_r], states[R_r],
                      states[VP_r], states[VS_r], states[VE_r], states[VI_r], states[VR_r],
                      states[VVP_r], states[VVS_r], states[VVE_r], states[VVI_r], states[VVR_r],
                      states[M_r])
    # we compute the difference with the expected total population (number of deaths)
    dM <- N - sum_states
    # we set the number of births in M_r to be the same as the number of deaths to keep the population constant
    states[M_r] <- dM
    
    # mortality correction
    # compute expected deaths on the upcoming day
    mod_pop <- rep(0,num_age_groups)
    for (i in 1:num_age_groups) {
      mod_pop[i] <- sum(states[S_r[i]], states[E_r[i]], states[I_r[i]], states[R_r[i]],
                        states[VP_r[i]], states[VS_r[i]], states[VE_r[i]], states[VI_r[i]], states[VR_r[i]],
                        states[VVP_r[i]], states[VVS_r[i]], states[VVE_r[i]], states[VVI_r[i]], states[VVR_r[i]])
    }
    
    exp_deaths <- sum(mod_pop*mu)
    
    states[M_r] <- (states[M_r]+exp_deaths)/(1-mu[1])
    
    return(states)
  }
  
  # MMR1 vaccination starting from birth cohort 1984 (so age 1 in 1985)
  mmr1_vaccination <- function(t, states, params, age_category_to_vacc, kappa) { 
    # this function takes also as argument the age category to be vaccinated, 
    # and the vaccination coverage kappa (defined in the events function)
    # this function applies first doses to unvaccinated individuals,
    # it is only used when age_category_to_vacc <= 1 y.o.
    
    # Index calculations, i.e. where we need to vaccinate (which age category)
    shift <- num_age_groups - 100 # this is equal to 0 if we have 100 age groups, 2 if we have 0-6m, 6-9m, 9-12m
	
    S_vacc_idx <- S_r[1] + age_category_to_vacc + shift # age category of 1 y.o. is element 2 because 0–99
                                                        # we also apply a shift based on how many age categories with monthly resolution
                                                        # we have within the age group of 0-1 y.o. 
    E_vacc_idx <- E_r[1] + age_category_to_vacc + shift
    I_vacc_idx <- I_r[1] + age_category_to_vacc + shift
    R_vacc_idx <- R_r[1] + age_category_to_vacc + shift
    VP_vacc_idx <- VP_r[1] + age_category_to_vacc + shift
    VS_vacc_idx <- VS_r[1] + age_category_to_vacc + shift
    VE_vacc_idx <- VE_r[1] + age_category_to_vacc + shift
    VI_vacc_idx <- VI_r[1] + age_category_to_vacc + shift
    VR_vacc_idx <- VR_r[1] + age_category_to_vacc + shift
	
    # Status needed from the compartments from which individuals leave
    S <- states[S_vacc_idx]
    E <- states[E_vacc_idx] 
    I <- states[I_vacc_idx] 
    R <- states[R_vacc_idx]
    
    #tau <- params[5]    # S -> VS, S -> VR, R-> VR # seroconversion probability

    if (age_category_to_vacc == 1) {
        tau = tau_by_age[['1-2y']]
    }
    if (age_category_to_vacc == 0) {
        tau = tau_by_age[['9-12m']]
    }
    if (age_category_to_vacc == -1) {
        tau = tau_by_age[['6-9m']]
    }
    if (age_category_to_vacc %in% c(1,0,-1) == FALSE) {
        # if the first dose is given later (e.g. to 20y - 40 y in a catch up campaign)
        tau = tau_by_age[['1-2y']]
        #stop("ERROR")
    }
    dS <- -S*kappa*tau - S*kappa*(1-tau) # individuals leave here
    dE <- -E*kappa 
    dI <- -I*kappa 
    dR <- -R*kappa
    
    dVP <- S*kappa*tau # those who DO seroconvert enter here (primary vaccination successfull)
    dVS <- S*kappa*(1-tau) # those who DO NOT seroconvert enter here (primary vaccination failure)
    dVE <- E*kappa
    dVI <- I*kappa 
    dVR <- R*kappa  
	
    # Add the changes to obtain the new state
    states[S_vacc_idx] <- states[S_vacc_idx] + dS
    states[E_vacc_idx] <- states[E_vacc_idx] + dE 
    states[I_vacc_idx] <- states[I_vacc_idx] + dI 
    states[R_vacc_idx] <- states[R_vacc_idx] + dR
    
    states[VP_vacc_idx] <- states[VP_vacc_idx] + dVP
    
    states[VS_vacc_idx] <- states[VS_vacc_idx] + dVS
    states[VE_vacc_idx] <- states[VE_vacc_idx] + dVE 
    states[VI_vacc_idx] <- states[VI_vacc_idx] + dVI 
    states[VR_vacc_idx] <- states[VR_vacc_idx] + dVR
    
    if(verbose) {	
      # check MMR1 coverage
      tot_pop_to_vacc <- sum(states[S_vacc_idx], 
                             states[E_vacc_idx], 
                             states[I_vacc_idx], 
                             states[R_vacc_idx],
                             states[VP_vacc_idx],
                             states[VS_vacc_idx], 
                             states[VE_vacc_idx], 
                             states[VI_vacc_idx], 
                             states[VR_vacc_idx])
      tot_pop_mmr1 <- sum(states[VP_vacc_idx],
                          states[VS_vacc_idx], 
                          states[VE_vacc_idx], 
                          states[VI_vacc_idx], 
                          states[VR_vacc_idx])
      print(sprintf("(%2f)", tot_pop_mmr1/tot_pop_to_vacc))
    }

    if (age_category_to_vacc %in% c(1,0,-1) == FALSE) {
         #print out number of doses administered to adults in the catch up campaign
         num_doses <- - (dS+dE+dI+dR)
         print(c(age_category_to_vacc, num_doses))
        }
        
    return(states)
  }
	
  # MMR2 vaccination starting from 1995
  mmr2_vaccination <- function(t, states, params, age_category_to_vacc, kappa, age_first_dose, immune_blunting) {

    # Index calculations, i.e. where we need to vaccinate (which age category)

    shift <- num_age_groups - 100 # this is equal to 0 if we have 100 age groups, 2 if we have 0-6m, 6-9m, 9-12m
	
    # we apply the second dose only to individuals who received the first dose 
    VP_vacc_idx <- VP_r[1] + age_category_to_vacc + shift # age category of 12 is element 13 because 0–99
                                                      # we also apply a shift based on how many age categories with monthly resolution
                                                      # we have within the age group of 0-1 y.o. 
    VS_vacc_idx <- VS_r[1] + age_category_to_vacc + shift
    VE_vacc_idx <- VE_r[1] + age_category_to_vacc + shift
    VI_vacc_idx <- VI_r[1] + age_category_to_vacc + shift 
    VR_vacc_idx <- VR_r[1] + age_category_to_vacc + shift
      
    VVP_vacc_idx <- VVP_r[1] + age_category_to_vacc + shift
    VVS_vacc_idx <- VVS_r[1] + age_category_to_vacc + shift
    VVE_vacc_idx <- VVE_r[1] + age_category_to_vacc + shift 
    VVI_vacc_idx <- VVI_r[1] + age_category_to_vacc + shift 
    VVR_vacc_idx <- VVR_r[1] + age_category_to_vacc + shift
    
    # Status needed from the compartments from which individuals leave
    VP <- states[VP_vacc_idx]
    VS <- states[VS_vacc_idx]
    VE <- states[VE_vacc_idx] 
    VI <- states[VI_vacc_idx] 
    VR <- states[VR_vacc_idx]
    
    # compute kappa_rel relative to the proportion of individuals with 1 dose
    # to translate kappa (coverage in age group) to kappa rel (coverage in age group with 1 dose)
    tot_pop_to_vacc <- sum(states[S_r[1] + age_category_to_vacc + shift], 
                             states[E_r[1] + age_category_to_vacc + shift], 
                             states[I_r[1] + age_category_to_vacc + shift], 
                             states[R_r[1] + age_category_to_vacc + shift],
                             states[VP_vacc_idx],
                             states[VS_vacc_idx], 
                             states[VE_vacc_idx], 
                             states[VI_vacc_idx], 
                             states[VR_vacc_idx])
    tot_pop_mmr1 <- sum(states[VP_vacc_idx],
                          states[VS_vacc_idx], 
                          states[VE_vacc_idx], 
                          states[VI_vacc_idx], 
                          states[VR_vacc_idx])
      
    # compute relative coverage
    kappa_rel <- min(1, kappa*tot_pop_to_vacc/tot_pop_mmr1) 
    # if kappa*tot_pop_to_vacc > tot_pop_mmr1, cap the relative coverage to 1
      
    if(verbose) {	
        # check if cov1>cov2
        print(sprintf("(%2f) (%2f)", tot_pop_to_vacc, tot_pop_mmr1))
        print(sprintf("(%2f) (%2f)", kappa, kappa_rel))
        print(sprintf("(%2f)", tot_pop_mmr1/tot_pop_to_vacc))
    }
      
    # Compute the change
      
    # we give the second dose to VS (susceptible with first dose)

    VE_tau = tau_by_age[["2nd dose"]]

    if (immune_blunting == TRUE) {
        if (age_first_dose == 1) {
            VE_tau = tau_by_age[["1-2y"]]
            }
        if (age_first_dose == 0.75) {
            #VE_tau = tau_by_age[["9-12m"]]
            VE_tau = tau_by_age[["2nd dose, immune blunting after 1st dose at 9 months"]] # scenario for the SHC: 
                                                                                          # much reduces for D1 at 9 months, some reduction for D2
            } 
        if (age_first_dose == 0.5) {
            VE_tau = tau_by_age[["2nd dose, immune blunting after 1st dose at 6 months"]] 
        }
        if (age_first_dose %in% c(0.5, 0.75, 1) == FALSE) {
            stop("error")
        }
    }
    dVS <- -VS*kappa_rel*VE_tau - VS*kappa_rel*(1-VE_tau) # individuals leave this compartment
    dVVP <- VS*kappa_rel*VE_tau # and enter here if they serocovert
    dVVS <- VS*kappa_rel*(1-VE_tau) # or here if they do not serocovert 
    
    # Add the changes to obtain the new state
    states[VS_vacc_idx] <- states[VS_vacc_idx] + dVS
    states[VVP_vacc_idx] <- states[VVP_vacc_idx] + dVVP
    states[VVS_vacc_idx] <- states[VVS_vacc_idx] + dVVS
    
    # we also give the second dose to the other compartments (VE, VI, VR)
    dVP <- -VP*kappa_rel
    dVE <- -VE*kappa_rel
    dVI <- -VI*kappa_rel
    dVR <- -VR*kappa_rel
    
    # Add the changes to obtain the new state
    states[VP_vacc_idx] <- states[VP_vacc_idx] + dVP
    states[VVP_vacc_idx] <- states[VVP_vacc_idx] - dVP   
    
    states[VE_vacc_idx] <- states[VE_vacc_idx] + dVE
    states[VVE_vacc_idx] <- states[VVE_vacc_idx] - dVE
    states[VI_vacc_idx] <- states[VI_vacc_idx] + dVI
    states[VVI_vacc_idx] <- states[VVI_vacc_idx] - dVI
    states[VR_vacc_idx] <- states[VR_vacc_idx] + dVR
    states[VVR_vacc_idx] <- states[VVR_vacc_idx] - dVR

    return(states)
  }
	
  events_function <- function(t, states, params) {

    # for specific time steps, we define events such as aging (ras_shift), vaccination, and initialization of newly infected
    
    if (save_newlyInf==T) {
        # Reset newly infectious
        # This must be reset to zero every day; you want this for each day
        # afterwards they are summed per year and cumulatively
        if(is.element(t, ts)) { # ts is the full list of time indexes
          states[NewlyInfectious_r] <- rep(0, length(NewlyInfectious_r))
          states[NewlyInfectious_dose0] <- rep(0, length(NewlyInfectious_dose0))
          states[NewlyInfectious_dose1] <- rep(0, length(NewlyInfectious_dose1))
          states[NewlyInfectious_dose2] <- rep(0, length(NewlyInfectious_dose2))
        }
    }

    if(disable_infections && t==t_disabling_year) {
      override_infections <<- T #'superassignment' operator
      print("## Overriding infections ##")
    }
    
    # Day everyone ages (RAS model), RAS shift
    # Applied once per year => days_of_year_change contains the days when shifting occurs
    # as a reminderd: days_of_year_change <- calculate_day_of_vaccination_for_years(yts_begin, yts_end, vacc_month, vacc_day)
    # this function computes September 1 for each year; function name is poorly chosen, because it is called vaccination but we are referring to aging
    if(is.element(t, days_of_6months_change)) { 
      states <- ras_shift_6months(t, states, params)
    }

    if(is.element(t, days_of_9months_change)) { 
      states <- ras_shift_9months(t, states, params)
    }
      
    if(is.element(t, days_of_year_change)) { 
      states <- ras_shift(t, states, params)
    }
    
    # Day of vaccination => check whether it is time to vaccinate
    # Same function as RAS, but +1 so it occurs on September 2
    # (here the function name *is* meaningful)
    # recall that: days_of_vaccination <- a two-column dataframe, first column = index t for Sept 2nd (Sept 1st + 1), second column year of vaccination
    if(is.element(t, days_of_vaccination[,1])) { 
      # if t is the date of vaccination, proceed
        
      idx <- which(days_of_vaccination[,1] == t) 
      # Check which calendar year are we dealing with, because it's a matrix. 
      year <- days_of_vaccination[idx,2] # second column contains the calendar year
      if(verbose)
        print(sprintf("vacc day: %d (%d)", t, year))
      
      #tau <- params[5]    # seroconversion probability # time- and age-independent
      
      # This loop performs the actual vaccination
      # age_cat_to_vacc => list we create ourselves in the main script, region-dependent
      # the keys of the list are as.character(year)
      years_to_vacc <- age_cat_to_vacc[[as.character(year)]] # get the target age group and vax coverage for the year of interest
      if(!is.null(years_to_vacc)) { # if it is not empty, continue
        
        if(is.matrix(years_to_vacc)) { 
          # years_to_vacc is a vector from 1985 to 1995, or a matrix from 1995 onwards (more than 1 target age group) 
          # if it is a matrix (i.e., if we are after 1995)
          for(age_category_to_vacc_idx in 1:length(years_to_vacc[,1])) { 
            
            age_category_to_vacc <- as.numeric(years_to_vacc[age_category_to_vacc_idx,1]) # get the age of the vaccinated individuals
            kappa <- as.numeric(years_to_vacc[age_category_to_vacc_idx,2]) # get the coverage
            dose <- years_to_vacc[age_category_to_vacc_idx,3] # get the dose number (MMR1 or MMR2)

            if (age_category_to_vacc >= 1) { 
                # we only assign first dose to at least 1 y.o. on this date (Sept 1)
                if(verbose)
                  print(sprintf("vacc age category: %d (%f) (%d)", age_category_to_vacc, kappa, year))
                
                if(dose == 'MMR1') { # coverage and age group refer to the first dose, give the first dose
                  states <- mmr1_vaccination(t, states, params, age_category_to_vacc, kappa)
                } else { # otherwise, give the second dose
                  
                  # reconstruct where they received the first dose (if at 6, 9, or 12 months)
                  target <- age_category_to_vacc #e.g. 2 y.o. or 9 y.o.
                  year_first_dose <- year - target + 1 #(e.g. 2035 - 2 + 1 = 2034)
                  #print(dose)
                  #print(year_first_dose)
                  #print(age_cat_to_vacc[[as.character(year_first_dose)]])
                  mm <- age_cat_to_vacc[[as.character(year_first_dose)]]
                  if (is.matrix(mm)) {
                  age_first_dose <- as.numeric(age_cat_to_vacc[[as.character(year_first_dose)]][1,1]) # extract first row (dose 1) and first columns (age group)                      
                  } else {
                  age_first_dose <- as.numeric(age_cat_to_vacc[[as.character(year_first_dose)]][1]) # extract first row (dose 1) and first columns (age group)   
                  }
                  states <- mmr2_vaccination(t, states, params, age_category_to_vacc, kappa, age_first_dose, immune_blunting)
                }
              }
          
          }
        
        } else { 
          # 1985-1994 => only one vaccination received, to 1 y.o.
          # MMR2 was only introduced after 1995, so here we work with a vector
          age_category_to_vacc <- as.numeric(years_to_vacc[1])
          kappa <- as.numeric(years_to_vacc[2]) # corresponding coverage
          
          if(verbose)
            print(sprintf("vacc age category: %d (%f) (%d)", age_category_to_vacc, kappa, year))
          # give the first dose
          states <- mmr1_vaccination(t, states, params, age_category_to_vacc, kappa)
        }
      }
    } 

    if(is.element(t, days_of_vaccination_6months[,1])) { 
      # if t is the date of vaccination, proceed
        
      idx <- which(days_of_vaccination_6months[,1] == t) 
      # Check which calendar year are we dealing with, because it's a matrix. 
      year <- days_of_vaccination_6months[idx,2] # second column contains the calendar year
      if(verbose)
        print(sprintf("vacc day: %d (%d)", t, year))
      
      #tau <- params[5]    # seroconversion probability # time- and age-independent
      
      # This loop performs the actual vaccination
      # age_cat_to_vacc => list we create ourselves in the main script, region-dependent
      # the keys of the list are as.character(year)
      years_to_vacc <- age_cat_to_vacc[[as.character(year)]] # get the target age group and vax coverage for the year of interest
      if(!is.null(years_to_vacc)) { # if it is not empty, continue
        
        if(is.matrix(years_to_vacc)) { 
          # years_to_vacc is a vector from 1985 to 1995, or a matrix from 1995 onwards (more than 1 target age group) 
          # if it is a matrix (i.e., if we are after 1995)
          for(age_category_to_vacc_idx in 1:length(years_to_vacc[,1])) { 
            
            age_category_to_vacc <- as.numeric(years_to_vacc[age_category_to_vacc_idx,1]) # get the age of the vaccinated individuals
            kappa <- as.numeric(years_to_vacc[age_category_to_vacc_idx,2]) # get the coverage
            dose <- years_to_vacc[age_category_to_vacc_idx,3] # get the dose number (MMR1 or MMR2)

            if (age_category_to_vacc == 0.5) { 
                
                # if first dose is assigned to age groups with monthly resolution less than 1 y.o.
                # adjust the index
                # target age 6-9 months --> input age_category_to_vacc = 0.5 
                # --> unshifted index -1 --> corrected index -1+2(shift) = 1 --> S_r[1] + 1

                # target age 9-6 months --> input age_category_to_vacc = 0.75 
                # --> unshifted index 0 --> corrected index 0+2(shift) = 2   --> S_r[1] + 2

                age_category_to_vacc <- -1
            
                if(verbose)
                  print(sprintf("vacc age category: %d (%f) (%d)", age_category_to_vacc, kappa, year))
                
                if(dose == 'MMR1') { # coverage and age group refer to the first dose, give the first dose
                  states <- mmr1_vaccination(t, states, params, age_category_to_vacc, kappa)
                } 
              
              }
            }
        }
    }
  }

    if(is.element(t, days_of_vaccination_9months[,1])) { 
      # if t is the date of vaccination, proceed
        
      idx <- which(days_of_vaccination_9months[,1] == t) 
      # Check which calendar year are we dealing with, because it's a matrix. 
      year <- days_of_vaccination_9months[idx,2] # second column contains the calendar year
      if(verbose)
        print(sprintf("vacc day: %d (%d)", t, year))
      
      #tau <- params[5]    # seroconversion probability # time- and age-independent
      
      # This loop performs the actual vaccination
      # age_cat_to_vacc => list we create ourselves in the main script, region-dependent
      # the keys of the list are as.character(year)
      years_to_vacc <- age_cat_to_vacc[[as.character(year)]] # get the target age group and vax coverage for the year of interest
      if(!is.null(years_to_vacc)) { # if it is not empty, continue
        
        if(is.matrix(years_to_vacc)) { 
          # years_to_vacc is a vector from 1985 to 1995, or a matrix from 1995 onwards (more than 1 target age group) 
          # if it is a matrix (i.e., if we are after 1995)
          for(age_category_to_vacc_idx in 1:length(years_to_vacc[,1])) { 
            
            age_category_to_vacc <- as.numeric(years_to_vacc[age_category_to_vacc_idx,1]) # get the age of the vaccinated individuals
            kappa <- as.numeric(years_to_vacc[age_category_to_vacc_idx,2]) # get the coverage
            dose <- years_to_vacc[age_category_to_vacc_idx,3] # get the dose number (MMR1 or MMR2)

            if (age_category_to_vacc == 0.75) { 
                
                # if first dose is assigned to age groups with monthly resolution less than 1 y.o.
                # adjust the index
                # target age 6-9 months --> input age_category_to_vacc = 0.5 
                # --> unshifted index -1 --> corrected index -1+2(shift) = 1 --> S_r[1] + 1

                # target age 9-6 months --> input age_category_to_vacc = 0.75 
                # --> unshifted index 0 --> corrected index 0+2(shift) = 2   --> S_r[1] + 2

                age_category_to_vacc <- 0
            
                if(verbose)
                  print(sprintf("vacc age category: %d (%f) (%d)", age_category_to_vacc, kappa, year))
                
                if(dose == 'MMR1') { # coverage and age group refer to the first dose, give the first dose
                  states <- mmr1_vaccination(t, states, params, age_category_to_vacc, kappa)
                } 
              
              }
            }
        }
    }
  }

  return(states)
  }

  # Define the transmission model
  model <- function(t, states, params) {
    S <- states[S_r]    
    E <- states[E_r]    
    I <- states[I_r]    
    R <- states[R_r]    
    
    VP <- states[VP_r]  
    VS <- states[VS_r]  
    VE <- states[VE_r]  
    VI <- states[VI_r]  
    VR <- states[VR_r]  
    
    VVP <- states[VVP_r]
    VVS <- states[VVS_r]  
    VVE <- states[VVE_r]  
    VVI <- states[VVI_r]  
    VVR <- states[VVR_r]  
    
    M <- states[M_r]    #1 (0-1 years)

    gamma <- params[2]  # M -> S
    sigma <- params[3]  # E -> I, VE -> VI
    nu <- params[4]     # I -> R, VI -> VR
    #tau <- params[5]   # S -> VP or VS, VS -> VVP or VVS
    waning <- params[6] # VP -> VS, VVP -> VVS

    idx_param <- 6 + 1:num_age_groups
    mu <- params[idx_param]  # age-specific death rate
    
    idx_param <- num_age_groups + 7
    N <- params[idx_param]
	
    # beta is the contact matrix divided by a constant (defined in the main script)
    # beta = c(aa)/maxeig(diagN*c(aa))
    idx_param <- num_age_groups + 7 + 1:(num_age_groups**2)
    beta <- matrix(params[idx_param], ncol=num_age_groups)
    
    meanR0 <- params[length(params)]
    if(override_infections) {
      meanR0 <- 0
    }

    # force of infection
    # VE_transm = 0 vaccine effectiveness against transmission by default
    lambda <- meanR0  * beta %*% (I+(1-VE_transm)*(VI+VVI))
    # matrix multiplication %*%: we multiply c(i,j) with I(j) 

    # compute transitions for each state
    # In general: +incoming -outgoing -death

    if(is.null(gamma_t)) {
        dM <- - M*gamma - M*mu[1] # - M*gamma: loss of maternal immunity; - M*mu[1] mortality for the first age group 0-1
    } else {
        if (t < gamma_t[3]) { # timestep of 30 years (avg age of giving birth) after year when vax coverage > 50% 
            dM <- - M*gamma_t[1] - M*mu[1] # duration of maternal immunity from unvax mothers                 
        } else {
            dM <- - M*gamma_t[2] - M*mu[1] # duration of maternal immunity from vaxxed mothers        
        }
    }
    
    dS <- - S*lambda - S*mu # 0-1, 1-2, ..., 99-100 years
    dS[1] <- dS[1] + M*gamma # 0-1 years

    dE <- S*lambda - E*sigma - E*mu
    dI <- E*sigma - I*nu - I*mu
    dR <- I*nu - R*mu
    
    dVP <- -VP*waning - VP*mu
    dVS <- VP*waning - VS*lambda - VS*mu  
    dVE <- VS*lambda - VE*sigma - VE*mu
    dVI <- VE*sigma - VI*nu - VI*mu
    dVR <- VI*nu - VR*mu 
    
    dVVP <- -VVP*waning - VVP*mu
    dVVS <- VVP*waning - VVS*lambda - VVS*mu  
    dVVE <- VVS*lambda - VVE*sigma - VVE*mu
    dVVI <- VVE*sigma - VVI*nu - VVI*mu
    dVVR <- VVI*nu - VVR*mu 
      
    dNewlyInfected <- E*sigma + VE*sigma + VVE*sigma # newly infected individuals
    dNewlyInfected_dose0 <- E*sigma 
    dNewlyInfected_dose1 <- VE*sigma 
    dNewlyInfected_dose2 <- VVE*sigma 
    
    list(c(dS, dE, dI, dR, 
           dVP, dVS, dVE, dVI, dVR, 
           dVVP, dVVS, dVVE, dVVI, dVVR, 
           dM, 
           dNewlyInfected,
           dNewlyInfected_dose0, 
           dNewlyInfected_dose1, 
           dNewlyInfected_dose2))
  }
  
  # Needed to know when to trigger events
  if (save_newlyInf==T) {
    event_ts <- ts # daily events, needed if you want to correctly save the new daily infections
  } else { 
    # events only on Sept 1st and Sept 2nd
    event_ts <- unique(c(days_of_6months_change,
                         days_of_9months_change,
                         days_of_year_change, 
                         days_of_vaccination_6months[,1],
                         days_of_vaccination_9months[,1],
                         days_of_vaccination[,1]))
  }
  
  y <- ode(y_init, ts, model, initial_params, hmax=1.0, method=method,
           events=list(func=events_function,
                       time=event_ts
                       )) 
  
  return(y)
}