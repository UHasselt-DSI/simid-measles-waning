library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(purrr)
library(viridis)  # for viridis palette
library(binom) # to compute Clopper Pearson coefficients for the seroprevalence data


# Helper functions for timesteps and data ######

calculate_days_per_year <- function(year) {
  jan_1 <- as.Date(sprintf("%d/01/01", year), format="%Y/%m/%d")
  dec_31 <- as.Date(sprintf("%d/12/31", year), format="%Y/%m/%d")
  return(as.integer(dec_31-jan_1)+1)
}

calculate_days_per_years <- function(begin_year, end_year) {
  range_years <- begin_year:end_year
  res <- rep(0, length(range_years))
  for(i in 1:length(range_years)) {
    res[i] <- calculate_days_per_year(range_years[i])
  }
  return(res)
}

calculate_days_over_years <- function(begin_year, end_year) {
  return(sum(calculate_days_per_years(begin_year, end_year)))
}

calculate_jan1_for_years <- function(begin_year, end_year) {
  days_per_year <- calculate_days_per_years(begin_year, end_year)
  range_years <- begin_year:end_year
  res <- rep(1, length(range_years))
  for(i in 2:length(range_years)) {
    res[i] <- sum(days_per_year[1:(i-1)])+1
  }
  return(res)
}

calculate_jan1_for_year <- function(begin_year, end_year, year) {
  jan1_years <- calculate_jan1_for_years(begin_year, end_year)
  idx_year <- which(year == begin_year:end_year)
  return(jan1_years[idx_year])
}

calculate_day_of_year <- function(year, month, day) {
  jan_1 <- as.Date(sprintf("%d/01/01", year), format="%Y/%m/%d")
  vacc_day <- as.Date(sprintf("%d/%02d/%02d", year, month, day), format="%Y/%m/%d")
  return(as.integer(vacc_day-jan_1)+1)
}

calculate_day_of_year_for_years <- function(begin_year, end_year, month, day) {
  jan1_years <- calculate_jan1_for_years(begin_year, end_year)
  range_years <- begin_year:end_year
  res <- rep(0, length(range_years))
  for(i in 1:length(range_years)) {
    res[i] <- calculate_day_of_year(range_years[i], month, day)
  }
  return(jan1_years + res - 1)
}

calculate_day_of_vaccination_for_years <- function(begin_year, end_year, month, day) {
  return(calculate_day_of_year_for_years(begin_year, end_year, month, day))
}

# Functions extracting model outputs ######

calculate_y_of_day_for_year <- function(y, begin_year, end_year, year, month, day, age_cat_per_compartment) {
  jan1_for_year <- calculate_jan1_for_year(begin_year, end_year, year)
  day_of_year <- calculate_day_of_year(year, month, day)
  # extract time step of the desired day 
  day_t <- jan1_for_year + day_of_year - 1

  # extract vector with model output for the chosen time step
  y_day <- y[day_t,((age_cat_per_compartment$S_r[1]:age_cat_per_compartment$M_r)+1)]
  # +1 because we want to omit the first column is omitted because, it contains the timestep t

  # we add fake M compartments with values 0 for ages from 2 to 100 
  # then we reshape the vector to a matrix with number of columns = number of epidemiological compartments 
  if ("NewlyInfectious_r" %in% names(age_cat_per_compartment)) {
    y_year_formatted <- matrix(c(y_day, rep(0, num_age_groups-1)), # we add 99 empty cells 
                               ncol=length(age_cat_per_compartment)-4) # reshape the vector into a matrix
                                                                       # the number of epidemiological compartms is
                                                                       # length(age_cat_per_compartment)-1
                                                                       # excluding the compartment for new infected (total and by vax dose)
  }
  else {
    y_year_formatted <- matrix(c(y_day, rep(0, num_age_groups-1)), 
                               ncol=length(age_cat_per_compartment))
  }
  # y_year_formatted has one age group for each row, one epidemiological compartment for each column
  return(y_year_formatted)
}

calculate_y_for_year <- function(y, begin_year, end_year, year, age_cat_per_compartment) {
  # this function compute the average daily prevalence in a given compartment over 1 year
  jan1_for_year <- calculate_jan1_for_year(begin_year, end_year, year)
  days_per_year <- calculate_days_per_year(year)
  idxs <- jan1_for_year:(jan1_for_year+days_per_year-1)
  y_year_days <- y[idxs,((age_cat_per_compartment$S_r[1]:age_cat_per_compartment$M_r)+1)] # +1; first column is omitted because it contains the timestep
  y_year <- colSums(y_year_days)/days_per_year 
  # colSums adds up all days (rows) for each age category-compartment (column), and then divides by the number of days in the year
  # because the model runs on a daily basis, we need a yearly value, so we take the average number for each compartment
  if ("NewlyInfectious_r" %in% names(age_cat_per_compartment)) {
    y_year_formatted <- matrix(c(y_year, rep(0, num_age_groups-1)), 
                               ncol=length(age_cat_per_compartment)-4) 
  }
  else {
    y_year_formatted <- matrix(c(y_year, rep(0, num_age_groups-1)), 
                               ncol=length(age_cat_per_compartment)) 
  }
  
  return(y_year_formatted)
}


## Immunity level######

calculate_immunity_per_day <- function(y, begin_year, end_year, 
                                       year, month, day,
                                       age_cat_per_compartment) { 
  # this function computes the prevalence in the compartments on a given day
  y_formatted <- calculate_y_of_day_for_year(y, begin_year, end_year, year, month, day, age_cat_per_compartment)
  # rows = age groups, columns = compartments
  y_total_per_age_cat <- rowSums(y_formatted) # we sum all compartments, to get the current population by age group
  # we select the compartment of interest, and we compute the sum for each row (each age group) over the columns
  y_immune_per_age_cat <- rowSums(y_formatted[,c(2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 15)]) # E, I, R, VP, VE, VI, VR, VVP, VVE, VVI, VVR, M
  y_immune <- 100*y_immune_per_age_cat/y_total_per_age_cat # we compute the seroprevalence proportion 
  # we correct by sensitivity and specificity of the test
  # y_seroyear_seroprev <- y_seroyear_prev * sensitivity + (1.-y_seroyear_prev) * (1.-specificity)

  y_immune_per_age_cat <- y_formatted[,15] # M
  y_maternal <- 100*y_immune_per_age_cat/y_total_per_age_cat 

  y_immune_per_age_cat <- rowSums(y_formatted[,c(2, 3, 4)]) # E, I, R
  y_natural <- 100*y_immune_per_age_cat/y_total_per_age_cat 

  y_immune_per_age_cat <- y_formatted[,5] # VP
  y_dose1 <- 100*y_immune_per_age_cat/y_total_per_age_cat 

  y_immune_per_age_cat <- y_formatted[,10] # VVP
  y_dose2 <- 100*y_immune_per_age_cat/y_total_per_age_cat 

  y_immune_per_age_cat <- rowSums(y_formatted[,c(7, 8, 9)]) # VE, VI, VR
  y_hybrid_dose1 <- 100*y_immune_per_age_cat/y_total_per_age_cat 

  y_immune_per_age_cat <- rowSums(y_formatted[,c(12, 13, 14)]) # VVE, VVI, VVR
  y_hybrid_dose2 <- 100*y_immune_per_age_cat/y_total_per_age_cat 

  y_pop <- y_total_per_age_cat

  y_susceptible <- 100 - y_maternal - y_natural - y_dose1 - y_dose2 - y_hybrid_dose1 - y_hybrid_dose2

  # extract those with at least one dose
  y_cov1 <- rowSums(y_formatted[,c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14)]) # VP, VS, VE, VI, VR, VVP, VVS, VVE, VVI, VVR 
  y_cov1 <- 100*y_cov1/y_total_per_age_cat 
  # extract those with exactly two doses
  y_cov2 <- rowSums(y_formatted[,c(10, 11, 12, 13, 14)]) # VVP, VVS, VVE, VVI, VVR 
  y_cov2 <- 100*y_cov2/y_total_per_age_cat 

  out <- rbind(y_immune, y_maternal, y_natural, y_dose1, y_dose2, y_hybrid_dose1, y_hybrid_dose2, 
               y_susceptible, y_pop, 
               y_cov1, y_cov2)
    
  out <- as.data.frame(t(out))
    
  return(out)
}

calculate_overall_immunity_per_day_over_years <- function(y, begin_year, end_year, 
                                                          list_years,
                                                           month, day,
                                                           age_cat_per_compartment, selected_ages = NA) { 

   out_fin <- data.frame()

   for (year in list_years) {
      # this function computes the prevalence in the compartments on a given day
      y_formatted <- calculate_y_of_day_for_year(y, begin_year, end_year, year, month, day, age_cat_per_compartment)

      if (any(is.na(selected_ages))) {
         y_formatted <- calculate_y_of_day_for_year(y, begin_year, end_year, year, month, day, age_cat_per_compartment)          
      } else {
         y_formatted <- calculate_y_of_day_for_year(y, begin_year, end_year, year, month, day, age_cat_per_compartment) 
         y_formatted <- y_formatted[selected_ages,]
      }
      # rows = age groups, columns = compartments
      y_total_per_age_cat <- rowSums(y_formatted) # we sum all compartments, to get the current population by age group
      y_total <- sum(y_total_per_age_cat) # we sum all age groups, to get the total population
      # we select the compartment of interest, and we compute the sum for each row (each age group) over the columns
      y_immune_per_age_cat <- rowSums(y_formatted[,c(2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 15)]) # E, I, R, VP, VE, VI, VR, VVP, VVE, VVI, VVR, M
      y_immune_tot <- sum(y_immune_per_age_cat) # we sum all age groups, to get the total 
      y_immune <- 100*y_immune_tot/y_total # we compute the seroprevalence proportion 
      # we correct by sensitivity and specificity of the test
      # y_seroyear_seroprev <- y_seroyear_prev * sensitivity + (1.-y_seroyear_prev) * (1.-specificity)
    
      y_immune_per_age_cat <- y_formatted[,15] # M
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_maternal <- 100*y_immune_tot/y_total
    
      y_immune_per_age_cat <- rowSums(y_formatted[,c(2, 3, 4)]) # E, I, R
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_natural <- 100*y_immune_tot/y_total
    
      y_immune_per_age_cat <- y_formatted[,5] # VP
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_dose1 <- 100*y_immune_tot/y_total
      y_effcov1 <- y_dose1
    
      y_immune_per_age_cat <- y_formatted[,10] # VVP
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_dose2 <- 100*y_immune_tot/y_total
      y_effcov2 <- y_dose2
    
      y_immune_per_age_cat <- rowSums(y_formatted[,c(7, 8, 9)]) # VE, VI, VR
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_hybrid_dose1 <- 100*y_immune_tot/y_total
    
      y_immune_per_age_cat <- rowSums(y_formatted[,c(12, 13, 14)]) # VVE, VVI, VVR
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_hybrid_dose2 <- 100*y_immune_tot/y_total

      y_immune_per_age_cat <- y_formatted[,1] # S
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_susceptible_novax <- 100*y_immune_tot/y_total

      y_immune_per_age_cat <- y_formatted[,6] # VS
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_susceptible_vax1 <- 100*y_immune_tot/y_total
       
      y_immune_per_age_cat <- y_formatted[,11] # VVS
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_susceptible_vax2 <- 100*y_immune_tot/y_total

      y_check_total <- y_maternal  + y_natural +  y_dose1 +  y_dose2 + y_hybrid_dose1 + y_hybrid_dose2 + 
                  y_susceptible_novax + y_susceptible_vax1 + y_susceptible_vax2
       
      y_immune_per_age_cat <- rowSums(y_formatted[,c(5, 6, 7, 8, 9)]) # VP, VS, VE, VI, VR
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_cov1 <- 100*y_immune_tot/y_total
    
      y_immune_per_age_cat <- rowSums(y_formatted[,c(10, 11, 12, 13, 14)]) # VVP, VVS, VVE, VVI, VVR
      y_immune_tot <- sum(y_immune_per_age_cat)
      y_cov2 <- 100*y_immune_tot/y_total

      y_cov <- y_cov1 + y_cov2
      y_effcov <- y_effcov1 + y_effcov2
    
      out <- rbind(y_immune, y_maternal, y_natural, y_dose1, y_dose2, y_hybrid_dose1, y_hybrid_dose2,
                  y_susceptible_novax, y_susceptible_vax1, y_susceptible_vax2, y_check_total, 
                   y_effcov1, y_effcov2, y_effcov, 
                   y_cov1, y_cov2, y_cov)
      out <- as.data.frame(t(out))
    
      out$year <- year
    
      out_fin <- rbind(out_fin, out)

  }
    
  return(out_fin)
}

plot_immunity_level <- function(y, begin_year, end_year, 
                                       year, month, day,
                                       age_cat_per_compartment) {

    out <- calculate_immunity_per_day(y, yts_begin, yts_end, 
                                           year, month, day,
                                           age_cat_per_compartment)
    
    out <- out[,c('y_immune','y_maternal','y_natural','y_dose1','y_dose2','y_hybrid_dose1','y_hybrid_dose2', 'y_cov1', 'y_cov2')]

    out$y_hybrid <- out$y_hybrid_dose1+out$y_hybrid_dose2
    
    age_groups <- rownames(m) 
    
    out$age <- age_groups
    
    df <- out
    
    #df <- head(df, 40)
    
    # Convert to long format
    df_long <- df %>%
      pivot_longer(
        cols = -c(age, y_immune, y_cov1, y_cov2, y_maternal, y_hybrid_dose1, y_hybrid_dose2),
        names_to = "immunity_type",
        values_to = "value"
      )
    
    # Keep age order
    df$age <- factor(df$age, levels = df$age)
    df_long$age <- factor(df_long$age, levels = df$age)
    
    # Rename for nicer legend
    df_long$immunity_type <- recode(df_long$immunity_type,
      #y_maternal = "Maternal",
      y_natural = "Infection (unvaccinated)",
      y_dose1 = "Vaccine (1 dose)",
      y_dose2 = "Vaccine (2 doses)",
      #y_hybrid_dose1 = "Infection (1 dose)",
      #y_hybrid_dose2 = "Infection (2 doses)"
      y_hybrid = "Infection (vaccinated)"
    )
    
    # CUSTOM STACK / LEGEND ORDER
    custom_order <- c(
      "Infection (unvaccinated)",       #"Maternal",
      "Vaccine (1 dose)",
      #"Infection (1 dose)",
      "Vaccine (2 doses)",
      "Infection (vaccinated)"
      #"Infection (2 doses)",
    )
    
    df_long$immunity_type <- factor(
      df_long$immunity_type,
      levels = custom_order
    )

    color_names <- c(
      #"Vaccine (1 dose)" = "#9ecae1",   # light blue
      #"Infection (1 dose)" = "#3182bd", # darker blue
      
      #"Vaccine (2 doses)" = "#a1d99b",              # light green
      #"Infection (2 doses)" = "#31a354", # darker green
      
      #"Maternal" = "#e6550d",            # strong orange (pops)

      "Vaccine (2 doses)" = "#72B7AD",              # light green
      "Vaccine (1 dose)" = "#28756B",              # darker green
      "Infection (vaccinated)" = "#e6550d",   # strong orange (pops)
      "Infection (unvaccinated)" = "#d9d9d9"            # soft neutral grey
    )

    exclude_ages <- c("[0,6months)", "[6,9months)", "[9,12months)", 
      "[81,82)", "[82,83)", "[83,84)", "[84,85)", "[85,86)",
      "[86,87)", "[87,88)", "[88,89)", "[89,90)", "[90,91)",
      "[91,92)", "[92,93)", "[93,94)", "[94,95)", "[95,96)",
      "[96,97)", "[97,98)", "[98,99)", "[99,100)"
    )
    
    df <- df %>%
      filter(!age %in% exclude_ages)
    
    df_long <- df_long %>%
      filter(!age %in% exclude_ages)
    
    # Plot
    plot <- ggplot(df_long, aes(x = factor(age), y = value, fill = immunity_type)) +
                  geom_bar(stat = "identity", width = 1, colour = "grey80"
                  ) +
                  geom_line(data = df, aes(x = factor(age), y = y_immune, group = 1),
                            colour = "black", linewidth = 1, inherit.aes = FALSE
                  ) +
                  geom_point(data = df, aes(x = factor(age), y = y_cov1, colour = "Vaccination coverage"),
                             size = 1, inherit.aes = FALSE
                  ) + 
                  scale_fill_manual(values = color_names
                  ) + 
                  scale_colour_manual(name = "", values = c("Vaccination coverage" = "black")
                  ) +
                  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)
                  ) + 
                  scale_x_discrete(breaks = levels(factor(df_long$age))[seq(1, length(unique(df_long$age)), by = 3)],
                                   expand = expansion(add = 1.5)
                  ) + 
                  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20), expand = c(0,0)
                  ) +
                  labs(title = paste(year), x = "Age group", y = "Immunity level (%)", fill = "Immunity type"
                  ) +
                  theme(panel.background = element_rect(fill="white"),
                        panel.grid.major.y = element_blank(),  # remove horizontal major grid
                        panel.grid.minor.y = element_blank(),  # remove horizontal minor grid 
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                        plot.title = element_text(hjust = 1),
                        text = element_text(size = 16),
                        panel.border = element_rect(colour = "black",fill = NA,linewidth = 1.),
                   )

    return(list("plot" = plot, "dataset" = df))
}

plot_overall_immunity_level <- function(y, begin_year, end_year, 
                                       list_years, month, day,
                                       age_cat_per_compartment, legend_pos = "right", waning_title = NA) {

    out <- calculate_overall_immunity_per_day_over_years(y, yts_begin, yts_end, 
                                           list_years, month, day,
                                           age_cat_per_compartment)
    
    df <- out

    df <- df[,c("y_immune", "y_maternal", "y_natural", "y_dose1", "y_dose2", "y_hybrid_dose1", "y_hybrid_dose2", "year")]

    df$y_hybrid <- out$y_hybrid_dose1+out$y_hybrid_dose2
    
    # Convert to long format
    df_long <- df %>%
      pivot_longer(
        cols = -c(year, y_immune, y_maternal, y_hybrid_dose1, y_hybrid_dose2),
        names_to = "immunity_type",
        values_to = "value"
      )
    
    # Rename for nicer legend
    df_long$immunity_type <- recode(df_long$immunity_type,
      #y_maternal = "Maternal",
      y_natural = "Infection (unvaccinated)",
      y_dose1 = "Vaccine (1 dose)",
      y_dose2 = "Vaccine (2 doses)",
      #y_hybrid_dose1 = "Infection (1 dose)",
      #y_hybrid_dose2 = "Infection (2 doses)"
      y_hybrid = "Infection (vaccinated)"
    )
    
    # CUSTOM STACK / LEGEND ORDER
    custom_order <- c(
      "Infection (unvaccinated)",       #"Maternal",
      "Vaccine (1 dose)",
      #"Infection (1 dose)",
      "Vaccine (2 doses)",
      "Infection (vaccinated)"
      #"Infection (2 doses)",
    )
    
    df_long$immunity_type <- factor(
      df_long$immunity_type,
      levels = custom_order
    )


    color_names <- c(
      "Vaccine (2 doses)" = "#72B7AD",              # light green
      "Vaccine (1 dose)" = "#28756B",              # darker green
      "Infection (vaccinated)" = "#e6550d",    # strong orange (pops)
      #"Dose 1" = "#9ecae1",              # light blue
      #"Infection after waning dose 1" = "#3182bd",  # darker blue
      #"Dose 2" = "#a1d99b",              # light green
      #"Infection after waning dose 2" = "#31a354",  # darker green
      #"Maternal" = "#e6550d",            # strong orange (pops)
      "Infection (unvaccinated)" = "#d9d9d9"            # soft neutral grey
    )
    
    # Plot
    plot <- ggplot(df_long, aes(x = year, y = value, fill = immunity_type)) +
      geom_bar(stat = "identity", width = 1) +
      # grid on top of bars
      #geom_vline(
      #  xintercept = seq(min(df_long$year), max(df_long$year), by = 10),
      #  colour = "grey50",
      #  linewidth = 0.3
      #) +
      #geom_hline(
      #  yintercept = seq(0, 100, by = 10),
      #  colour = "grey50",
      #  linewidth = 0.3
      #) +
      geom_line(
        data = df,
        aes(x = year, y = y_immune, group = 1),
        colour = "black",
        linewidth = 1.25,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = color_names) + 
        scale_x_continuous(
          breaks = seq(min(df_long$year), max(df_long$year), by = 20),
            #expand = c(0, 0)
        ) + 
      scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
          expand = c(0, 0)
      ) +
      labs(
        title = waning_title,
        x = "Year",
        y = "Immunity level (%)",
        fill = "Immunity type"
      ) +
      theme(panel.background = element_rect(fill="white"),
            panel.grid = element_blank(),
            #panel.ontop = TRUE,
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = legend_pos,
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1.),
            text = element_text(size = 16),
            plot.title = element_text(hjust = 1, size = 16),
           )

    return(list("plot" = plot, "dataset" = df_long))
}


## Seroprevalence ######

calculate_seroprevalence <- function(y, begin_year, end_year, year, age_cat_per_compartment, sensitivity, specificity) { 
  # this function computes the average seroprevalence over 1 year
  y_seroyear_formatted <- calculate_y_for_year(y, begin_year, end_year, year, age_cat_per_compartment)
  # rows = age groups, columns = compartments
  # we select the compartment of interest, and we compute the sum for each row (each age group) over the columns
  y_seroyear_seropos_per_age_cat <- rowSums(y_seroyear_formatted[,c(3, 4, 5, 8, 9, 10, 13, 14, 15)]) # I, R, VP, VI, VR, VVP, VVI, VVR, M
  y_seroyear_total_per_age_cat <- rowSums(y_seroyear_formatted) # we sum all compartments, to get the current population by age group
  y_seroyear_prev <- y_seroyear_seropos_per_age_cat/y_seroyear_total_per_age_cat # we compute the seroprevalence proportion 
  # we correct by sensitivity and specificity of the test
  y_seroyear_seroprev <- y_seroyear_prev * sensitivity + (1.-y_seroyear_prev) * (1.-specificity)
  return(y_seroyear_seroprev)
}

calculate_seroprevalence_for_day <- function(y, begin_year, end_year, year, month, day, age_cat_per_compartment,
                                             sensitivity, specificity) { 
  # this function computes the seroprevalence on a specific day
  # this is also used in the likelihood function
  y_seroyear_formatted <- calculate_y_of_day_for_year(y, begin_year, end_year, year, month, day, age_cat_per_compartment)
  y_seroyear_seropos_per_age_cat <- rowSums(y_seroyear_formatted[,c(3, 4, 5, 8, 9, 10, 13, 14, 15)]) # I, R, VP, VI, VR, VVP, VVI, VVR, M
  y_seroyear_total_per_age_cat <- rowSums(y_seroyear_formatted)
  y_seroyear_prev <- y_seroyear_seropos_per_age_cat/y_seroyear_total_per_age_cat
  y_seroyear_seroprev <- y_seroyear_prev * sensitivity + (1.-y_seroyear_prev) * (1.-specificity)
  return(y_seroyear_seroprev)
}


prepare_susceptibility_percentages_normal_vs_shift <- function(y_normal, baseline_name="no shift", 
                                                               y_shifts, y_shifts_names,
                                                               sim_begin_year, sim_end_year,
                                                               plot_begin_year, plot_end_year, month, day,
                                                               age_cat_per_compartment) {
  
  plot_years <- plot_begin_year:plot_end_year
  t_scope <- rep(NA, length(plot_years))
  # calculates the day of the year at which we want to evaluate the model
  for(i in 1:length(plot_years)) {
    plot_year <- plot_years[i]
    jan1_for_year <- calculate_jan1_for_year(sim_begin_year, sim_end_year, plot_year)
    day_of_year <- calculate_day_of_year(plot_year, month, day)
    t_scope[i] <- jan1_for_year + day_of_year - 1
  }
  # selects the relevant columns we need: S and VS (and also VVS in this case)
  # these are computed separately and then combined into a single vector
  # (general function that is then applied to y_normal and y_shift)
  
  if ("NewlyInfectious_r" %in% names(age_cat_per_compartment)) {
    idx_skip <- c(0, age_cat_per_compartment$NewlyInfectious_r[1]:age_cat_per_compartment$NewlyInfectious_dose2[length(age_cat_per_compartment$NewlyInfectious_dose2)])+1
  }
  else {
    idx_skip <- 1
  }
  idx_susceptibles <- ((age_cat_per_compartment$S_r[1]:age_cat_per_compartment$S_r[length(age_cat_per_compartment$S_r)])+1)
  idx_vacc_susceptibles <- ((age_cat_per_compartment$VS_r[1]:age_cat_per_compartment$VS_r[length(age_cat_per_compartment$VS_r)])+1)
  idx_vacc_vacc_susceptibles <- ((age_cat_per_compartment$VVS_r[1]:age_cat_per_compartment$VVS_r[length(age_cat_per_compartment$VVS_r)])+1)
  idx <- c(idx_susceptibles, idx_vacc_susceptibles, idx_vacc_vacc_susceptibles)
  # Percentages are calculated
  y_normal_perc_susceptibles <- 100*rowSums(y_normal[t_scope,idx])/rowSums(y_normal[t_scope,-idx_skip])   # first row of the matrix corresponds to t = 0, so it must be excluded
  
  df_x <- c(plot_years)
  df_y <- c(y_normal_perc_susceptibles)
  df_strategy <- c(rep(baseline_name, length(plot_years)))
  
  for(i in 1:length(y_shifts)) {
    y_shift <- y_shifts[[i]]
    df_x <- c(df_x,plot_years)
    df_y <- c(df_y,100*rowSums(y_shift[t_scope,idx])/rowSums(y_shift[t_scope,-idx_skip]))
    df_strategy <- c(df_strategy, rep(y_shifts_names[i], length(plot_years)))
  }
  
  return(data.frame(x=df_x, y=df_y, scenario=df_strategy))
}

## New Infections ######

calculate_newly_infectious_for_year_per_day <- function(y, begin_year, end_year, year, age_cat_per_compartment,
                                                        selected_columns, selected_ages) {
  jan1_for_year <- calculate_jan1_for_year(begin_year, end_year, year)
  days_per_year <- calculate_days_per_year(year)
  idxs <- jan1_for_year:(jan1_for_year+days_per_year-1)

  y_year_days <- y[idxs,(age_cat_per_compartment[[selected_columns]][selected_ages]+1)] 

   
  
  return(y_year_days)
}

calculate_newly_infectious_for_year <- function(y, begin_year, end_year, year, age_cat_per_compartment, 
                                                selected_columns, selected_ages) {
  y_year_days <- calculate_newly_infectious_for_year_per_day(y, begin_year, end_year, year, age_cat_per_compartment, 
                                                             selected_columns, selected_ages)
  y_year <- colSums(y_year_days) 
  
  return(y_year)
}

# Plotting functions ######

plot_lexis_diagram_seroprevalence <- function(y, begin_year, end_year, 
                                              list_years, month, day,
                                              ages,
                                              age_cat_per_compartment, sensitivity, specificity) {
    
    df <- data.frame(age = names(ages))

    for (year in list_years) {
        gg <- calculate_seroprevalence_for_day(y, begin_year, end_year, year, month, day,
                                              age_cat_per_compartment, sensitivity, specificity)[(unname(unlist(ages)))] 

        df[[as.character(year)]] <- gg
    }
    
    # Convert wide to long
    df_long <- df %>%
      pivot_longer(
        cols = -age,
        names_to = "year",
        values_to = "seroprev"
      )
    
    df_long <- df_long %>%
      mutate(
        year = as.numeric(gsub("^X", "", year))   # remove leading X
      )
    
    # Convert year to numeric
    df_long$year <- as.numeric(df_long$year)
    
    # Convert age ranges to numeric midpoints
    df_long$age_mid <- sapply(
      strsplit(gsub("\\[|\\)|\\s", "", df_long$age), ","), 
      function(x) mean(as.numeric(x))
    )
    
    # Check seroprev values: must be numeric and non-NA
    df_long$seroprev <- as.numeric(df_long$seroprev)
    
    # Plot
    plot <- ggplot(df_long, aes(x = year, y = age_mid, fill = seroprev)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "A", 
                           limits = c(0.7,1.), 
                           name = "Seroprevalence") +
      scale_y_continuous(name = "Age (years)") +
      scale_x_continuous(name = "Year (cohort)") +
      theme_minimal() +
      labs(title = "Lexis diagram of seroprevalence") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

    return(plot)
}

plot_lexis_diagram_seroprevalence_difference <- function(y1, y2,
                                                         begin_year, end_year, 
                                              list_years, month, day,
                                              ages,
                                              age_cat_per_compartment, sensitivity, specificity) {
    
    df <- data.frame(age = names(ages))

    for (year in list_years) {
        gg1 <- calculate_seroprevalence_for_day(y1, begin_year, end_year, year, month, day,
                                              age_cat_per_compartment, sensitivity, specificity)[(unname(unlist(ages)))] 

        gg2 <- calculate_seroprevalence_for_day(y2, begin_year, end_year, year, month, day,
                                              age_cat_per_compartment, sensitivity, specificity)[(unname(unlist(ages)))] 

        gg <- gg1 - gg2

        df[[as.character(year)]] <- gg
    }
    
    # Convert wide to long
    df_long <- df %>%
      pivot_longer(
        cols = -age,
        names_to = "year",
        values_to = "seroprev"
      )
    
    df_long <- df_long %>%
      mutate(
        year = as.numeric(gsub("^X", "", year))   # remove leading X
      )
    
    # Convert year to numeric
    df_long$year <- as.numeric(df_long$year)
    
    # Convert age ranges to numeric midpoints
    df_long$age_mid <- sapply(
      strsplit(gsub("\\[|\\)|\\s", "", df_long$age), ","), 
      function(x) mean(as.numeric(x))
    )
    
    # Check seroprev values: must be numeric and non-NA
    df_long$seroprev <- as.numeric(df_long$seroprev)
    
    # Plot
    plot <- ggplot(df_long, aes(x = year, y = age_mid, fill = seroprev)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "A", 
                           limits = c(0,0.25), 
                           name = "Absolute difference \n(Seroprevalence)") +
      scale_y_continuous(name = "Age (years)") +
      scale_x_continuous(name = "Year (cohort)") +
      theme_minimal() +
      labs(title = "Lexis diagram of seroprevalence (difference)") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

    return(plot)
}


plot_seroprevalence_on_day_for_years <- function(y, y2,
                                                  begin_year, end_year,
                                                  years, month, day,
                                                  age_cat_per_compartment,
                                                  ages,
                                                  sensitivity, specificity,
                                                  subtitle = NA,
                                                  add_labels = TRUE,
                                                  surveys = NA, 
                                                  tags = c("model", "model 2", "data"), 
                                                  uncertainty = uncertainty, 
                                                  list_col) {

ages_idx <- unname(unlist(ages))

if(is.null(y2)) {
  sim_df <- map_dfr(years, function(yr){
    
    vals1 <- calculate_seroprevalence_for_day(y, begin_year, end_year, yr, month, day,
                                              age_cat_per_compartment, sensitivity, specificity)[ages_idx]
    
    data.frame(x = names(ages),
                x_ticks = ages_idx,
                y = vals1,
                ymin = vals1,
                ymax = vals1, 
                year = yr,
                dataset = tags[1]
              )}
    )
  
  if (is.null(uncertainty) == FALSE) {
   
    sim_df <- map_dfr(years, function(yr) {
      
      vals1 <- calculate_seroprevalence_for_day(y, begin_year, end_year, yr, month, day,
                                                age_cat_per_compartment, sensitivity, specificity)[ages_idx]
      
      gg <- uncertainty

    if (nrow(gg[[as.character(yr)]]) > 2) {
      quantiles <- apply(gg[[as.character(yr)]], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
      
      vals_CI_up <- quantiles["97.5%", ]
      
      vals_CI_dw <- quantiles["2.5%", ]   
    }
    else {
      # fast CI by plotting the curves with the boundaries of the CI
        
      vals_CI_up <- gg[[as.character(yr)]][2,]
      
      vals_CI_dw <- gg[[as.character(yr)]][1,]  
    }
      
      data.frame(
        x = names(ages),
        x_ticks = ages_idx,
        ymin = vals_CI_dw,
        ymax = vals_CI_up,
        y = vals1,
        year = yr,
        dataset = tags[1]
      )
    })
  }
} else {
  sim_df <- map_dfr(years, function(yr){
  
    vals1 <- calculate_seroprevalence_for_day(y, begin_year, end_year, yr, month, day,
                                              age_cat_per_compartment, sensitivity, specificity)[ages_idx]
  
    vals2 <- calculate_seroprevalence_for_day(y2, begin_year, end_year, yr, month, day,
                                              age_cat_per_compartment, sensitivity, specificity)[ages_idx]
  
    bind_rows(
      data.frame(x = names(ages),
                 x_ticks = ages_idx,
                 y = vals1,
                 ymin = vals1,
                 ymax = vals1, 
                 year = yr,
                 dataset = tags[1]
      ),
      data.frame(x = names(ages),
                      x_ticks = ages_idx,
                      y = vals2,
                      year = yr,
                      dataset = tags[2]
      )
    )
  })
}

if(is.data.frame(surveys)){

  idxs <- which(surveys$x %in% names(ages))

  survey_df <- data.frame(
    x = surveys$x[idxs],
    x_ticks = surveys$x_ticks[idxs],
    y = surveys$y[idxs],
    ymin = surveys$lower[idxs],
    ymax = surveys$upper[idxs],
    year = surveys$surveyyear[idxs],
    dataset = tags[length(tags)]
  )
  df_res <- bind_rows(sim_df, survey_df)
} else {
  df_res <- sim_df
}

#title <- paste("Seroprevalence", unique(years))
title <- paste(unique(years))
if(is.na(subtitle)) subtitle <- waiver()

y_lim_low <- 0.4
custom_line_width <- 1.5
if (is.null(uncertainty) == FALSE) {
  y_lim_low <- 0.3
  custom_line_width <- 0.75
}

plot <- ggplot(df_res, aes(x = x_ticks, y = y, colour = dataset)) +
            geom_ribbon(data = df_res %>% filter(dataset == tags[1]),
                        aes(x = x_ticks, ymin = ymin, ymax = ymax,
                            group = interaction(dataset, year)),
                        inherit.aes = FALSE,
                        fill = scales::alpha(list_col[1], 0.22),
                        colour = scales::alpha(list_col[1], 0.22),
                        show.legend = FALSE
            ) +
            geom_errorbar(data = df_res %>% filter(dataset == tags[length(tags)]),
                       aes(x = x_ticks, ymin = ymin, ymax = ymax,
                          group = interaction(dataset, year)),
                       width = 0,  # vertical line only, no horizontal caps
                       colour = scales::alpha("black", 0.3),
                       inherit.aes = FALSE,
                       show.legend = FALSE
            ) +
            geom_line(data = df_res %>% filter(dataset != tags[length(tags)]), linewidth = custom_line_width
            ) +
            geom_point(data = df_res %>% filter(dataset == tags[length(tags)]), size = 1.5
            ) +
            geom_line(data = df_res %>% filter(dataset == tags[length(tags)]), linewidth = 0.5   
            ) +
            scale_x_continuous(breaks = unique(df_res$x_ticks)[seq(1, length(unique(df_res$x_ticks)), by = 5)],
                               labels = unique(df_res$x)[seq(1, length(unique(df_res$x_ticks)), by = 5)]
            ) +
            scale_colour_manual(values = setNames(c(list_col[1], list_col[2], "black"),tags), name = " ")

  plot <- plot + labs(x = "Age",
                      y = "Seroprevalence (%)",
                      title = title,
                      #subtitle = subtitle,
                      colour = "Data",
                ) +
                scale_y_continuous(limits = c(y_lim_low, 1), labels = function(x) x * 100
                ) + 
                theme(panel.background = element_rect(fill="white"),
                      axis.text.x = element_text(size = 12, angle=90, vjust=0.5, hjust=1),
                      axis.text.y = element_text(size = 12),
                      legend.position = "none",
                      panel.border = element_rect(colour = "black",fill = NA,linewidth = 1),
                      axis.line = element_blank(),
                      plot.title = element_text(hjust = 1, size = 14),
                      axis.title.x = element_text(size = 14),
                      axis.title.y = element_text(size = 14)
                    )
    
  if(add_labels){
    plot <- plot +     
            theme(
          legend.position = c(0.95, 0.05),        # bottom-right corner
          legend.justification = c(1, 0),    # anchor legend box to that corner
          legend.text = element_text(size = 14),
        ) 
      }
  
  #if(add_labels){
  #  plot <- plot +     
  #    theme(legend.position = "right")
#}

return(plot)
}

plot_seroprevalence_y_vs_surveys <- function(y, y2, 
                                             yts_begin, yts_end, 
                                             survey_years, sero_month, sero_day, 
                                             age_cat_per_compartment,
                                             sensitivity, specificity, 
                                             ages, 
                                             subtitle, use_flanders_data,
                                             tags = c("model", "model 2", "data"), 
                                             uncertainty = NULL, 
                                             list_col = c("#1f78b4", "#e31a1c")) {
    
  df_survey_2002 <- prepare_survey_data(use_flanders_data, add_2002 = T, add_2006 = F, add_2013 = F)

  plot_seros_2002 <- plot_seroprevalence_on_day_for_years(y, y2, yts_begin, yts_end, 
                                                          2002, sero_month, sero_day, age_cat_per_compartment, 
                                                          ages, sensitivity, specificity, subtitle, 
                                                          surveys=df_survey_2002,
                                                          add_labels=FALSE,
                                                          tags = tags, uncertainty = uncertainty, 
                                                          list_col = list_col)
  
  df_survey_2006 <- prepare_survey_data(use_flanders_data, add_2002 = F, add_2006 = T, add_2013 = F)
    
  plot_seros_2006 <- plot_seroprevalence_on_day_for_years(y, y2, yts_begin, yts_end, 
                                                          2006, sero_month, sero_day, age_cat_per_compartment, 
                                                          ages, sensitivity, specificity, subtitle, 
                                                          surveys=df_survey_2006,
                                                          add_labels=FALSE,
                                                          tags = tags, uncertainty = uncertainty, 
                                                          list_col = list_col)
  
  df_survey_2013 <- prepare_survey_data(use_flanders_data, add_2002 = F, add_2006 = F, add_2013 = T) 
    
  plot_seros_2013 <- plot_seroprevalence_on_day_for_years(y, y2, yts_begin, yts_end, 
                                                          2013, sero_month, sero_day, age_cat_per_compartment, 
                                                          ages, sensitivity, specificity, subtitle, 
                                                          surveys=df_survey_2013,
                                                          tags = tags, 
                                                          uncertainty = uncertainty, 
                                                          list_col = list_col)
  
  return(list(plot_2002=plot_seros_2002,
              plot_2006=plot_seros_2006,
              plot_2013=plot_seros_2013))
}

plot_seroprevalence_years <- function(y, begin_year, end_year, 
                                      years, age_cat_per_compartment, ages,
                                      sensitivity, specificity, 
                                      subtitle=NA, add_labels=TRUE, surveys=NA) {
    x_ticks <- c()
    x <- c()
    vals <- c()  
    var <- c()
    is_sim <- c()
    d_size <- c()
    for(i in 1:length(years)) {
    x_ticks <- c(x_ticks, unname(unlist(ages)))
    x <- c(x, names(ages))
    var <- c(var, as.character(rep(years[i], length(ages))))
    vals <- c(vals, calculate_seroprevalence(y, begin_year, end_year, years[i], age_cat_per_compartment, sensitivity, specificity)[(unname(unlist(ages)))]) 
    }
    
    # read seroprevalence data
    # seroprevalence data only have individuals >= 1 y.o.
    df_res <- NA
    dem_years <- years
    if(is.data.frame(surveys)) {
        idxs <- which(surveys$x %in% names(ages))
        is_sim <- as.character(c(rep(1, length(x)), rep(0, length(surveys$x[idxs]))))
        d_size <- as.numeric(c(rep(0, length(x)), rep(1, length(surveys$x[idxs]))))
        x <- c(x, surveys$x[idxs])
        x_ticks <- c(x_ticks, surveys$x_ticks[idxs])
        vals <- c(vals, surveys$y[idxs])
        var <- c(var, surveys$surveyyear[idxs])
        df_res <- data.frame(x=x, x_ticks = x_ticks, y=vals, var=var, is_sim=is_sim)
        dem_years <- c(dem_years, 1:length(unique(surveys$surveyyear[idxs])))
    } else {
        df_res <- data.frame(x=x, x_ticks = x_ticks, y=vals, var=var)
    }
    
      colors_survey <- c("firebrick1", "cornflowerblue", "chartreuse4",  "#FFAE00", "darkorchid4", "#8ECE98", "darkred", "#F007C9", "sienna4", "gray40", "yellow3")
    colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
              '#cab2d6','#6a3d9a','#b15928', "black", "firebrick1", "cornflowerblue", "chartreuse4", "sienna4")
    
    xlab <- "Age"
    ylab <- "Seroprevalence"
    title <- NA
    
    if(length(unique(years)) == 1) {
    title <- paste("Seroprevalence", unique(years))
    } else {
    title <- paste("Seroprevalence", min(years), "-", max(years))
    }
    
    if(is.na(subtitle)) {
    subtitle <- waiver()
    }
    
    plot <- NA
    if(is.data.frame(surveys)) {
    if(length(unique(years)) == 1) {
      plot <- ggplot(df_res, aes(x = x_ticks, y = vals, colour = is_sim, shape = is_sim)) +
        geom_line(aes(linetype=is_sim)) + geom_point(size=d_size, show.legend=FALSE) +
        scale_linetype_manual(values=c("twodash", "solid")) +
        scale_x_continuous(
        breaks = df_res$x_ticks,
        labels = df_res$x) +
        guides(shape="none", linetype="none")
    } else {
      plot <- ggplot(df_res, aes(x = x_ticks, y = vals, colour = var, shape = is_sim)) +
        geom_line(aes(linetype=is_sim)) + geom_point(size=d_size, show.legend=FALSE) +
        scale_linetype_manual(values=c("twodash", "solid")) +
        scale_x_continuous(
        breaks = df_res$x_ticks,
        labels = df_res$x) +
        guides(shape="none", linetype="none")
    }
    } else {
    plot <- ggplot(df_res, aes(x = x_ticks, y = vals, colour = var)) +
      geom_line() + 
      scale_x_continuous(
      breaks = df_res$x_ticks,
      labels = df_res$x)
    }
    
    if(add_labels) {
    plot <- plot +
      labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
      scale_y_continuous(limits = c(0,1)) +
      scale_colour_manual("Years",values=colors[1:length(dem_years)]) +
      theme(legend.position = "right",
            legend.key=element_blank(),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.line.y.left = element_line(colour = "black"),
            axis.line.x.bottom = element_line(colour = "black"))
    } else {
    plot <- plot +
      scale_y_continuous(limits = c(0.7,1)) +
      scale_colour_manual("Years",values=colors[1:length(dem_years)]) +
      theme(legend.position = "none",
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.line.y.left = element_line(colour = "black"),
            axis.line.x.bottom = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    }
    return(plot)
}


plot_susceptibility_percentages_normal_vs_shift <- function(df_perc_raw, subtitle=NA, add_labels=T, 
                                                            plot_year_begin=NA, plot_year_end=NA,
                                                            uncertainty = NA) {
  colors_survey <- c("firebrick1", "cornflowerblue", "darkorchid4", "chartreuse4",  
                     "#FFAE00", "#8ECE98", "darkred", "#F007C9", "sienna4", "gray40", "yellow3")

  xlab <- "Year"
  ylab <- "Immune population (%)"
  #title <- "Percentage of seronegatives - Normal vs. shift strategy"
  
  if(is.na(plot_year_begin)) {
    plot_year_begin <- min(df_perc_raw$x)
  }
  if(is.na(plot_year_end)) {
    plot_year_end <- max(df_perc_raw$x)
  }
  plotting_years <- plot_year_begin:plot_year_end
  idx <- which(df_perc_raw[,1] %in% plotting_years)
  df_perc <- df_perc_raw[idx,]
  
  if(is.na(subtitle)) {
    subtitle <- waiver()
  }
  
  df_perc$scenario <- factor(
    df_perc$scenario,
    levels = unique(df_perc$scenario)
  )

  df_perc$y <- 100 - df_perc$y

  if (is.data.frame(uncertainty)) {
    
    uncertainty_plot <- uncertainty %>%
      filter(x %in% plotting_years) %>%
      pivot_longer(
        cols = starts_with("y"),
        names_to = "sample",
        values_to = "y"
      ) %>%
      group_by(x, scenario) %>%
      summarise(
        lower = quantile(y, probs = 0.025, na.rm = TRUE),
        upper = quantile(y, probs = 0.975, na.rm = TRUE),
        central = quantile(y, probs = 0.5, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        lower = 100 - lower,
        upper = 100 - upper,
        central = 100 - central
      )
    
    # Set same scenario order as main dataframe
    uncertainty_plot$scenario <- factor(
      uncertainty_plot$scenario,
      levels = levels(df_perc$scenario)
    )

     # Replace central estimate while preserving the original
      # x/scenario order of df_perc
      df_perc <- df_perc %>%
        select(-y) %>%
        left_join(
          uncertainty_plot %>%
            select(x, scenario, central),
          by = c("x", "scenario")
        ) %>%
        rename(y = central)
  }
  
  plot <- ggplot(df_perc, aes(x=x, y=y, colour=scenario))

  # Add uncertainty ribbons first, so lines appear on top
  if (is.data.frame(uncertainty)) {
    
    plot <- plot +
      geom_ribbon(
        data = uncertainty_plot,
        aes(
          x = x,
          ymin = lower,
          ymax = upper,
          fill = scenario
        ),
        alpha = 0.2,
        colour = NA,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

 # Add main lines
  plot <- plot +   
    geom_line() + 
    coord_cartesian(ylim = c(75,100))
  
  num_datasets <- length(unique(df_perc$scenario))
  
  if(add_labels) {
    plot <- plot +
      labs(x=xlab, y=ylab
           #title=title, 
           #subtitle=subtitle
          ) +
      scale_colour_manual("",values=colors_survey[1:num_datasets]) +
      scale_fill_manual(
          values = colors_survey[1:num_datasets]
      ) + 
      theme(legend.position = c(0.05, 0.05), legend.justification = c(0, 0), # bottom left
            legend.key=element_blank(),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1.
              ),
            axis.line = element_blank()
            #axis.line.y.left = element_line(colour = "black"),
            #axis.line.x.bottom = element_line(colour = "black")
           )
  } else {
    plot <- plot +
      scale_colour_manual("",values=colors_survey[1:num_datasets]) +
      scale_fill_manual(
          values = colors_survey[1:num_datasets]
      ) + 
      theme(legend.position = "none",
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1.
              ),
            axis.line = element_blank(),
            #axis.line.y.left = element_line(colour = "black"),
            #axis.line.x.bottom = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
  }
  plot <- plot + theme(text = element_text(size = 16))
    
  return(plot)
}

plot_newly_infectious <- function(y, yts_begin, yts_end, age_cat_per_compartment, subtitle, plot_year_begin=NA, plot_year_end=NA,
                                 selected_columns = "NewlyInfectious_r", selected_ages = c(1:num_age_groups)) {
  if(is.na(plot_year_begin)) {
    plot_year_begin <- yts_begin
  }
  if(is.na(plot_year_end)) {
    plot_year_end <- yts_end
  }
  newly_infectious_years <- plot_year_begin:plot_year_end
  res_newly_infectious_years <- matrix(NA, ncol=length(age_cat_per_compartment$NewlyInfectious_r), nrow=length(newly_infectious_years))
  newly_infectious_years_sum <- rep(NA, length(newly_infectious_years))
  for(i in 1:length(newly_infectious_years)) {
    year <- newly_infectious_years[i]
    res_newly_infectious_years[i,] <- calculate_newly_infectious_for_year(y, yts_begin, yts_end, year, age_cat_per_compartment,
                                                                          selected_columns, selected_ages)
    newly_infectious_years_sum[i] <- sum(res_newly_infectious_years[i,])
  }
  newly_infectious_years_cumsum <- cumsum(newly_infectious_years_sum)
  
  df_newly_infectious_sum <- data.frame(x=newly_infectious_years, y=newly_infectious_years_sum)
  df_newly_infectious_cumsum <- data.frame(x=newly_infectious_years, y=newly_infectious_years_cumsum)
  
  df_newly_infectious_sum_no1940 <- data.frame(x=newly_infectious_years[-1], y=newly_infectious_years_sum[-1])
  df_newly_infectious_cumsum_no1940 <- data.frame(x=newly_infectious_years[-1], y=cumsum(newly_infectious_years_sum[-1]))
  
  df_newly_infectious_sum_start2020 <- data.frame(x=newly_infectious_years[-c(1:80)], y=newly_infectious_years_sum[-c(1:80)])
  df_newly_infectious_cumsum_start2020 <- data.frame(x=newly_infectious_years[-c(1:80)], y=cumsum(newly_infectious_years_sum[-c(1:80)]))
  
  xlab <- "Year"
  ylab <- "Incidence"
  title <- paste("Incidence")
  plot_ni <- ggplot(df_newly_infectious_sum, aes(x, y)) +
    geom_line() +
    labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line.y.left = element_line(colour = "black"),
          axis.line.x.bottom = element_line(colour = "black"))
  
  title <- paste("Cumulative incidence")
  plot_ni_cumsum <- ggplot(df_newly_infectious_cumsum, aes(x, y)) +
    geom_line() +
    labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line.y.left = element_line(colour = "black"),
          axis.line.x.bottom = element_line(colour = "black"))
  
  title <- paste("Incidence")
  plot_ni_no1940 <- ggplot(df_newly_infectious_sum_no1940, aes(x, y)) +
    geom_line() +
    labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line.y.left = element_line(colour = "black"),
          axis.line.x.bottom = element_line(colour = "black"))
  
  title <- paste("Cumulative incidence")
  plot_ni_cumsum_no1940 <- ggplot(df_newly_infectious_cumsum_no1940, aes(x, y)) +
    geom_line() +
    labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line.y.left = element_line(colour = "black"),
          axis.line.x.bottom = element_line(colour = "black"))
  
  title <- paste("Incidence")
  plot_ni_start2020 <- ggplot(df_newly_infectious_sum_start2020, aes(x, y)) +
    geom_line() +
    labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line.y.left = element_line(colour = "black"),
          axis.line.x.bottom = element_line(colour = "black"))
  
  title <- paste("Cumulative incidence")
  plot_ni_cumsum_start2020 <- ggplot(df_newly_infectious_cumsum_start2020, aes(x, y)) +
    geom_line() +
    labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line.y.left = element_line(colour = "black"),
          axis.line.x.bottom = element_line(colour = "black"))
  
  return(list(plot_ni=plot_ni, plot_ni_cumsum=plot_ni_cumsum, 
              plot_ni_no1940=plot_ni_no1940, plot_ni_cumsum_no1940=plot_ni_cumsum_no1940,
              plot_ni_start2020=plot_ni_start2020, plot_ni_cumsum_start2020=plot_ni_cumsum_start2020))
}

plot_newly_infectious_per_category_for_year_groups_list <- function(
    y_vals_list,
    dataset_labels,
    yts_begin,
    yts_end,
    year_groups,
    age_cat_grouped,
    age_cat_per_compartment,
    subtitle,
    selected_columns = "NewlyInfectious_r",
    custom_colors = NA, 
    legend_remove = FALSE, 
    show_data = TRUE
) {

  x <- c()
  y <- c()
  age <- c()
  dataset <- c()
  year_label <- c()

  for (d in seq_along(y_vals_list)) {

    y_vals <- y_vals_list[[d]]

    for (k in seq_along(age_cat_grouped)) {

      for (j in seq_len(nrow(year_groups))) {

        years <- year_groups[j,1]:year_groups[j,2]
        yy <- c()

        for (year in years) {

          ni_raw <- calculate_newly_infectious_for_year(
            y_vals,
            yts_begin,
            yts_end,
            year,
            age_cat_per_compartment,
            selected_columns,
            selected_ages = c(1:num_age_groups)
          )

          yy <- c(yy, sum(ni_raw[age_cat_grouped[[k]]]))
        }

        x <- c(x, names(age_cat_grouped)[k])
        y <- c(y, sum(yy))
        age <- c(age, names(age_cat_grouped)[k])
        dataset <- c(dataset, dataset_labels[d])

        year_label_add <- sprintf("%d-%d", year_groups[j,1], year_groups[j,2])
        if (year_groups[j,1] == year_groups[j,2]) {year_label_add <- sprintf("%d", year_groups[j,1])}
        year_label <- c(year_label,
                        year_label_add)
      }
    }
  }

    colors <- c("firebrick1", "cornflowerblue", "darkorchid4", "gray")

    if (show_data == TRUE) {

        if (year_groups[1,1] >=2024 &  year_groups[1,1] <=2025) {
    
            if (year_groups[1,1] == year_groups[1,2]) {
            
    
            sciensano_data <- read.csv("data/AGE_DISTRIBUTION_2018.csv")
            sciensano_data <- sciensano_data %>% filter(Year==year_groups[1,1]) %>% group_by(AgeGroup_UH) %>%
            summarise(cases = sum(COUNT, na.rm = TRUE))
        
            sciensano_data$AgeGroup_UH <- factor(sciensano_data$AgeGroup_UH, levels = age_levels)
            
            total_by_age <- aggregate(cases ~ AgeGroup_UH, data = sciensano_data, sum)
            
            x <- c(x, names(age_cat_grouped))
            y <- c(y, total_by_age$cases)
            age <- c(age, names(age_cat_grouped))
            dataset <- c(dataset, rep("Data", length(total_by_age$cases)))
    
            year_label_add <- sprintf("%d-%d", year_groups[j,1], year_groups[j,2])
            if (year_groups[j,1] == year_groups[j,2]) {year_label_add <- sprintf("%d", year_groups[j,1])}
                
            year_label <- c(year_label,
                            rep(year_label_add, length(total_by_age$cases)))
        
            colors[length(dataset_labels)+1] <- "gray"
    
            }
    
        }

    }

  df <- data.frame(
    x = x,
    y = y,
    dataset = dataset,
    year = year_label
  )

  df$dataset <- factor(
      df$dataset,
      levels = c(dataset_labels, "Data")
    )

  # normalize within (age group + year)
  df <- df %>%
    group_by(year, dataset) %>%
    mutate(y = y / sum(y)) %>%
    ungroup()

  if (!is.na(custom_colors[1])) colors <- custom_colors

  plot <- ggplot(df, aes(x = x, y = y, fill = dataset)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    facet_wrap(~year) +
    scale_x_discrete(limits = names(age_cat_grouped)) +
    scale_y_continuous(labels = function(x) x * 100) + 
    scale_fill_manual(values = colors) +
    labs(
      x = "Age groups",
      y = "Proportion of cases (%)",
      #subtitle = subtitle,
      fill = NULL
    ) +
    theme(
        text = element_text(size = 15),
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      ),
      axis.line = element_blank(),
      strip.text = element_text(size = 16, hjust = 1),
      strip.background = element_blank(),
      axis.text.y = element_text(size = 15, angle = 90, hjust = 0.5),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
      #axis.line = element_line(colour = "black")
      #axis.line.x = element_line(colour = "black"),
      #axis.line.y = element_line(colour = "black")
    )

    if (legend_remove) {
    plot <- plot + theme(legend.position = "none")
  }

    return(plot)
}

plot_newly_infectious_per_category_for_year_groups <- function(y_vals, yts_begin, yts_end, 
                                                               year_groups, age_cat_grouped, age_cat_per_compartment, 
                                                               subtitle, 
                                                               selected_columns = "NewlyInfectious_r", 
                                                               custom_colors = NA) {
  xlab <- "Age groups"
  ylab <- "Proportion of cases (%)"
  title <- ""#"Incidence per age group"
  x <- c()
  y <- c()
  label <- c()
  for(k in 1:length(age_cat_grouped)) {
    for(j in 1:dim(year_groups)[1]) {
      years <- year_groups[j,1]:year_groups[j,2]
      yy <- c()
      for(i in 1:length(years)) {
        year <- years[i]
        ni_raw <- calculate_newly_infectious_for_year(y_vals, yts_begin, yts_end, year, age_cat_per_compartment, 
                                                     selected_columns, selected_ages = c(1:num_age_groups))
        yy <- c(yy, sum(ni_raw[age_cat_grouped[[k]]]))
      }
      x <- c(x, names(age_cat_grouped)[k])
      y <- c(y, sum(yy))

      year_label_add <- sprintf("%d-%d", year_groups[j,1], year_groups[j,2])
      if (year_groups[j,1] == year_groups[j,2]) {year_label_add <- sprintf("%d", year_groups[j,1])}
      label <- c(label, year_label_add)
    }
  }
  
  df_ni_dp <- data.frame(x=x, y=y, label=label)
  df_ni_dp <- df_ni_dp %>%
    group_by(label) %>%
    mutate(y = y / sum(y))
  #print(df_ni_dp)
  #colors_survey <- c("firebrick1", "cornflowerblue", "chartreuse4", "sienna4", "gray40", "yellow3", "darkorchid4", "darkred")
  #colors <- c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7", "gray40", "yellow3", "darkorchid4", "darkred")

  # gradient using cornflowerblue
  colors <- c( "#18365b", "#477db5","cornflowerblue","#F06292", "#C2185B")
  # gradient brown yeallow teal  
  #colors <- c('#7f3b08','#e08214','#fdb863','#80cdc1','#35978f','#01665e','#003c30')
  
  if (is.na(custom_colors) == FALSE) {
   colors <- custom_colors 
  }
  
  plot <- ggplot(df_ni_dp, aes(fill=label, y=y, x=x)) +
    geom_bar(width=0.8, position="dodge", stat="identity") +
    scale_x_discrete(limits = names(age_cat_grouped)) +
    scale_y_continuous(labels = function(x) x * 100) + 
    #coord_flip() +
    labs(x=xlab, y=ylab, 
         title=title, subtitle=subtitle) +
    scale_fill_manual(name = NULL,, values=colors[1:dim(year_groups)[1]]) +
    theme(text = element_text(size = 15),
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.border = element_rect(
            colour = "black",
            fill = NA,
            linewidth = 1
          ),
          axis.line = element_blank(),
          axis.text.y = element_text(size = 15, angle = 90, hjust = 0.5),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
          #axis.line.y.left = element_line(colour = "black"),
          #axis.line.x.bottom = element_line(colour = "black"), 
          plot.subtitle = element_text(hjust = 1))
  
  return(plot)
}

plot_newly_infectious_per_category_for_year_groups_absolute_numbers <- function(y_vals, yts_begin, yts_end, year_groups, age_cat_grouped, age_cat_per_compartment, subtitle, y_lim_max, selected_columns = "NewlyInfectious_r") {
  xlab <- "Age groups"
  ylab <- "Number of cases"
  title <- "Incidence per age group"
  x <- c()
  y <- c()
  label <- c()
  for(k in 1:length(age_cat_grouped)) {
    for(j in 1:dim(year_groups)[1]) {
      years <- year_groups[j,1]:year_groups[j,2]
      yy <- c()
      for(i in 1:length(years)) {
        year <- years[i]
        ni_raw <- calculate_newly_infectious_for_year(y_vals, yts_begin, yts_end, year, age_cat_per_compartment, 
                                                     selected_columns, selected_ages = c(1:num_age_groups))
        yy <- c(yy, sum(ni_raw[age_cat_grouped[[k]]]))
      }
      x <- c(x, names(age_cat_grouped)[k])
      y <- c(y, sum(yy))
      label <- c(label, sprintf("%d-%d", year_groups[j,1], year_groups[j,2]))
    }
  }
  
  df_ni_dp <- data.frame(x=x, y=y, label=label)
  #df_ni_dp <- df_ni_dp %>%
  #  group_by(label) %>%
  #  mutate(y = y / sum(y))
  #print(df_ni_dp)
  colors_survey <- c("firebrick1", "cornflowerblue", "chartreuse4",  "#FFAE00", "darkorchid4", "#8ECE98", "darkred", "#F007C9", "sienna4", "gray40", "yellow3")
  colors <- c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7", "gray40", "yellow3", "darkorchid4", "darkred")
  plot <- ggplot(df_ni_dp, aes(fill=label, y=y, x=x)) +
    geom_bar(width=0.8, position="dodge", stat="identity") +
    scale_x_discrete(limits = names(age_cat_grouped)) +
    coord_cartesian(ylim = c(0, y_lim_max)) +   # set y-axis limits
    labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
    scale_fill_manual("Years",values=colors[1:dim(year_groups)[1]]) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line.y.left = element_line(colour = "black"),
          axis.line.x.bottom = element_line(colour = "black"))
  
  return(plot)
}

# Likelihood ######

prepare_survey_data <- function(use_flanders_data, add_2002=TRUE, add_2006=TRUE, add_2013=TRUE) {
  data2002 <- NA
  data2006 <- NA
  data2013 <- NA

  
  h5_fn <- "./data/seroprevalence-measles.h5"

  if (use_flanders_data == "Flanders") {

    data2002 <- h5read(h5_fn, "survey/2002/flanders")
    data2006 <- h5read(h5_fn, "survey/2006/flanders")
    data2013 <- h5read(h5_fn, "survey/2013/flanders")

  } else if (use_flanders_data == "Wallonia") {

    data2002 <- h5read(h5_fn, "survey/2002/walbxl")
    data2006 <- h5read(h5_fn, "survey/2006/walbxl")
    data2013 <- h5read(h5_fn, "survey/2013/walbxl")

  } else if (use_flanders_data == "Belgium") {

    data2002 <- h5read(h5_fn, "survey/2002/national")
    data2006 <- h5read(h5_fn, "survey/2006/national")
    data2013 <- h5read(h5_fn, "survey/2013/national")

  } else {
    stop("Misspecified region, check string in use_flanders_data")
  }
  
  surveys_x <- c()
  surveys_y <- c()
  surveys_surveyyear <- c()
  surveys_tot <- c()
  surveys_lower <- c()
  surveys_upper <- c()
  
  if(add_2002) {
    surveys_x <- c(surveys_x, data2002$age)
    surveys_y <- c(surveys_y, data2002$prev)
    surveys_surveyyear <- c(surveys_surveyyear, rep(2002, length(data2002$age)))
    surveys_tot <- c(surveys_tot, data2002$tot)
    # add Clopper Pearson binomial confidence interval
    surveys_lower <- c(surveys_lower, binom.confint(data2002$pos, data2002$tot, conf.level = 0.95, methods = "exact")$lower)
    surveys_upper <- c(surveys_upper, binom.confint(data2002$pos, data2002$tot, conf.level = 0.95, methods = "exact")$upper)
  }
  
  if(add_2006) {
    surveys_x <- c(surveys_x, data2006$age)
    surveys_y <- c(surveys_y, data2006$prev)
    surveys_surveyyear <- c(surveys_surveyyear, rep(2006, length(data2006$age)))
    surveys_tot <- c(surveys_tot, data2006$tot)
    # add Clopper Pearson binomial confidence interval
    surveys_lower <- c(surveys_lower, binom.confint(data2006$pos, data2006$tot, conf.level = 0.95, methods = "exact")$lower)
    surveys_upper <- c(surveys_upper, binom.confint(data2006$pos, data2006$tot, conf.level = 0.95, methods = "exact")$upper)
  }
  
  if(add_2013) {
    surveys_x <- c(surveys_x, data2013$age)
    surveys_y <- c(surveys_y, data2013$prev)
    surveys_surveyyear <- c(surveys_surveyyear, rep(2013, length(data2013$age)))
    surveys_tot <- c(surveys_tot, data2013$tot)
    # add Clopper Pearson binomial confidence interval
    surveys_lower <- c(surveys_lower, binom.confint(data2013$pos, data2013$tot, conf.level = 0.95, methods = "exact")$lower)
    surveys_upper <- c(surveys_upper, binom.confint(data2013$pos, data2013$tot, conf.level = 0.95, methods = "exact")$upper)
  }
  
  surveys <- data.frame(x=surveys_x, y=surveys_y, surveyyear=surveys_surveyyear, 
                        tot=surveys_tot, lower=surveys_lower, upper=surveys_upper)
  shift <- num_age_groups - 100 
  surveys$x_ticks <- surveys$x + 1 + shift
  surveys$x <- paste0("[", surveys$x, ",", surveys$x + 1, ")")
  return(surveys)
}

likelihood_y_binomial <- function(y, yts_begin, yts_end, 
                                  survey_years, sero_month, sero_day,
                                  age_cat_per_compartment, 
                                  sensitivity, specificity, num_age_groups, df_surveys,
                                  beta_binomial = FALSE, k_disp = NA) {
  sumsies <- 0
  for(survey_year in survey_years) { 
    # loop over the survey years (2002, 2006 and 2013)
    
    y_formatted_data <- calculate_seroprevalence_for_day(y, yts_begin, yts_end, survey_year,
                                                         sero_month, sero_day, age_cat_per_compartment,
                                                         sensitivity, specificity)
    
    df_subset <- df_surveys[which(df_surveys$data_year == survey_year), ] 

    # selects the rows of the correct survey year; df_surveys contains all of them
    unique_ages <- unique(df_subset$age) 
    # check how many ages there are (because 2002 and 2006 only go up to age 65)
    # ages go from 1 = [1,2) y.o. to [65-66) y.o.

    for(k in 1:length(unique_ages)) { 
      # this loop performs the fitting, it goes over each of the unique ages
      unique_age <- unique_ages[k]

      # Index calculations
      shift <- num_age_groups - 100 # this is equal to 0 if we have 100 age groups, 2 if we have 0-6m, 6-9m, 9-12m
      
      # age category of 1 y.o. is element 2 because 0–99
      # the seroprevalence data start from age 1
      yh <- y_formatted_data[(unique_age+1+shift)]  # yh = model prediction, index+1+shift because ages go from 0–99.

      # if unique_age = 1, then the column of interest is 1+1+shift = 2 if 100 groups, 4 if 102 groups
      
      df_subset_age <- df_subset[which(df_subset$age == unique_age), ] 

      yy <- df_subset_age$prev # yy = empirical value, seroprevalence

      pos_sample <- df_subset_age$pos
      n_sample <- df_subset_age$tot

      if (beta_binomial == FALSE) {
        sumsies <- sumsies + pos_sample*log(yh) + (n_sample - pos_sample)*log(1.-yh) # log binomial   
      }
      else { #if (beta_binomial == TRUE)

        # log likelihood of BetaBinomial
        # parameterization p = alpha/(alpha + beta) and dispersion parameter k = alpha+beta
        # then alpha = k*p and beta = k*(1-p)

        alpha <- k_disp * yh
        beta  <- k_disp * (1 - yh)

        pos <- pos_sample
        n <- n_sample

        sumsies <- sumsies +
          lgamma(n_sample + 1) - lgamma(pos_sample + 1) - lgamma(n_sample - pos_sample + 1) +   # log 
          lgamma(pos_sample + alpha) +
          lgamma(n_sample - pos_sample + beta) -
          lgamma(n_sample + alpha + beta) -
          lgamma(alpha) -
          lgamma(beta) +
          lgamma(alpha + beta)
          
      }
        
    }
  }
  return(sumsies)
}

