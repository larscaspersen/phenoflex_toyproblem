#helper functions

library(ggplot2)
library(tidyr)
require(reshape2)

#qqplot of phenology data

get_qqplot <- function(cal_data, predicted){
  
  plot_df <- data.frame(cal_data, Predicted = predicted,
                        Error = cal_data$pheno - predicted)
  
  return(ggplot(plot_df, aes(x=pheno,y=Predicted)) +
           geom_point() +
           geom_abline(intercept=0, slope=1) +
           theme_bw(base_size = 15) +
           xlab("Observed bloom date (Day of the year)") +
           ylab("Predicted bloom date (Day of the year)") +
           ggtitle("Predicted vs. observed bloom dates"))
}



get_hist_plot <- function(error){
  return(ggplot(data.frame(error = error), aes(error)) + geom_histogram() +
           ggtitle("Distribution of prediction errors"))
} 


######Temperature response####

# make functions to visualize temperature responses during chilling and forcing
apply_const_temp <- function(temp,
                             A0,
                             A1,
                             E0,
                             E1,
                             Tf,
                             slope,
                             portions = 1200,
                             deg_celsius = TRUE){
  temp_vector <- rep(temp, times = portions)
  res <- chillR::DynModel_driver(
    temp = temp_vector,
    A0 = A0,
    A1 = A1,
    E0 = E0,
    E1 = E1,
    Tf = Tf,
    slope = slope,
    deg_celsius = deg_celsius
  )
  return(res$y[length(res$y)])
}

gen_bell <- function(par, temp_values = seq(-5, 20, 0.1)) {
  E0 <- par[5]
  E1 <- par[6]
  A0 <- par[7]
  A1 <- par[8]
  Tf <- par[9]
  slope <- par[12]
  
  y <- c()
  for (i in seq_along(temp_values)) {
    y[i] <- apply_const_temp(
      temp = temp_values[i],
      A0 = A0,
      A1 = A1,
      E0 = E0,
      E1 = E1,
      Tf = Tf,
      slope = slope
    )
  }
  return(invisible(y))
}

GDH_response <- function(T, par)
{
  Tb <- par[11]
  Tu <- par[4]
  Tc <- par[10]
  GDH_weight <- rep(0, length(T))
  GDH_weight[which(T >= Tb & T <= Tu)] <-
    1 / 2 * (1 + cos(pi + pi * (T[which(T >= Tb &
                                          T <= Tu)] - Tb) / (Tu - Tb)))
  GDH_weight[which(T > Tu & T <= Tc)] <-
    (1 + cos(pi / 2 + pi / 2 * (T[which(T >  Tu &
                                          T <= Tc)] - Tu) / (Tc - Tu)))
  return(GDH_weight)
}


get_temp_response_plot <- function(par, temp_values){
  
  temp_response <- data.frame(
    Temperature = temp_values,
    Chill_response = gen_bell(par, temp_values),
    Heat_response = GDH_response(temp_values, par)
  )
  
  melted_response <- reshape2::melt(temp_response, id.vars = "Temperature")
  
  return(ggplot(melted_response, aes(x = Temperature, y = value)) +
           geom_line(size = 2, aes(col = variable)) +
           ylab("Temperature response (arbitrary units)") +
           xlab("Temperature (°C)") +
           facet_wrap(vars(variable),
                      scales = "free",
                      labeller = labeller(variable = c(
                        Chill_response = c("Chill response"),
                        Heat_response = c("Heat response")
                      ))) +
           scale_color_manual(values = c("Chill_response" = "blue", "Heat_response" = "red")) +
           theme_bw(base_size = 15) +
           theme(legend.position = "none"))
}

#function to see development of the parameters
get_par_plot <- function(fit_list){
  
  #extract the parameters from the list
  pars_in <- fit_list %>%
    purrr::map('par.guess') %>%
    unlist()
  
  #parameter names, should be the same order 
  vars <-  rep(c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1', 'Tf', 'Tc','Tb','slope'),length(fit_list))
  
  #extract upper and lower bounds
  upper <- fit_list %>%
    purrr::map('upper') %>%
    unlist()
  
  lower <- fit_list %>%
    purrr::map('lower') %>%
    unlist()
  
  #create run number object
  run <- rep(seq(1:length(fit_list)), each = 12)
  
  #bind to tibble
  par_df <- tibble::tibble(run, vars, pars_in)
  
  boundary_df <- tibble(run,vars, upper, lower)
  
  min_max <- boundary_df %>%
    dplyr::group_by(vars) %>%
    dplyr::summarise(max_up = max(upper),
              min_low = min (lower))
  
  #add / substract 5% of upper and lower
  min_max$max_up_plus <-  min_max$max_up * 1.05
  min_max$min_low_min <-  min_max$min_low * 0.95
  
  boundary_df <- merge(boundary_df, min_max, by = 'vars')
  
  p <- ggplot(par_df, aes(x = run, y = pars_in)) + geom_line() + 
    geom_ribbon(data = boundary_df, aes(x = run, ymin = upper, ymax = max_up_plus),
                fill = 'grey', alpha = 0.7)+
    geom_ribbon(data = boundary_df, aes(x = run, ymin = min_low_min, ymax = lower),
                fill = 'grey', alpha = 0.7)+
    facet_wrap(~vars, scales = 'free_y')
  
  return(p)
  
}
