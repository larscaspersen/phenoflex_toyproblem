library(chillR)
library(MEIGOR)
library(tidyverse)
library(nleqslv)

source('code/04_functions_meigo.R')

#read the training data of all cultivars
files <- list.files('data/training/', full.names = TRUE)
pheno_train <- lapply(files, read.csv)

#read evaluation data
files <- list.files('data/evaluation/', full.names = TRUE)
pheno_eval <- lapply(files, read.csv)

cultivars <- c('Blanquina', 'Clara', 'Collaos', 'Coloradona',
               'DelaRiega', 'Perezosa', 'Perico', 'Raxao',
               'Teorica', 'Verdialona', 'Xuanina')

names(pheno_train) <- names(pheno_eval) <- cultivars

#read temperature data
hourtemps <- read.csv('data/hourtemps.csv')

temp_list <- temp_list_eval <- list()

#generate list of temperatures for each cultivar
for(cult in cultivars){
  temp_list[[cult]] <- chillR::genSeasonList(hourtemps, years = pheno_train[[cult]]$Year)
  temp_list_eval[[cult]] <- chillR::genSeasonList(hourtemps, years = pheno_eval[[cult]]$Year)
}


#      ('yc',          'zc',          's1',      'Tu', 'theta*', 'theta_c', 'Tau(thetha*)', 'pie_c',   'Tf', 'Tc', 'Tb',  'slope')
x_0 <- c(rep(40,11),   rep(190,11),   rep(0.5,11),  25,   279,    286.1,     47.7,             28,        4,   36,     4,    1.60)
x_U <- c(rep(80,11),   rep(500,11),   rep(1.0,11),  30,   281,      287,       48,             28,       10,   40,    10,    5.00)
x_L <- c(rep(20,11),   rep(100,11),   rep(0.1,11),   0,   279,      286,       16,             24,        0,    0,     0,    0.05)

#limits for the inequality constraints
#         #gdh parameters   #q10 for E0 and E1
c_L <- c(  0,   0,   0,     1.5, 1.5)
c_U <- c(Inf, Inf, Inf,     3.5, 3.5)


evaluation_function_meigo_nonlinear_combined <- function(x, 
                                                modelfn,
                                                bloomJDays,
                                                SeasonList,
                                                na_penalty = 365,
                                                return_bloom_days = FALSE){
  
  #innput:
  #         x is the parameters in meigo
  #         modelfn is the function used to calculate the bloomdays
  #         SeasonList contains the weather data for the different years
  #         na_penalty is the value used if the model failed to predict any day with the current set of parameters for a particular year
  
  #output: inequality constraints g
  #        model performance value F
  
  
  
  #instead of A0, A1, E0 and E1 we have now theta*, theta_c, Tau(thetha*) and pie_c
  #--> we need to solve now a non-linear system to calculate A0, A1, E0 and E1 based on the parameters
  #    ('yc', 'zc', 's1', 'Tu', 'theta*', 'theta_c', 'Tau(thetha*)', 'pie_c', 'Tf', 'Tc', 'Tb',  'slope')
  #x<- c(40,   190,   0.5,  25,   279,      287,       28,             26,       4,   36,    4,    1.60)
  
  params<-numeric(4)
  
  params[1] <- x[5+30]   #theta*
  params[2] <- x[6+30]    #theta_c
  params[3] <- x[7+30]    #Tau(thetha*)
  params[4] <- x[8+30]     #pi_c
  
  
  output<-nleqslv(c(500, 15000), solve_nle, jac=NULL, params, xscalm="auto", method="Newton",
                  control=list(trace=0,allowSingular=TRUE))
  
  
  #This is a numerical method which can produce non-convergence. Check this
  if (output$termcd >= 3){
    #if the nle algorithm has stalled just discard this solution
    E0<-NA; E1<-NA; A0<-NA; A1<-NA
    return(list(F=10^6, g=rep(10^6,5)))
    
    #You would add here a flag to let your optimization procedure know
    #That this solution should be ignored by lack of convergence
    
  } else {
    
    E0 <- output$x[1]
    E1 <- output$x[2]
    
    #A1 and A0 can be calculated through Equations 36 and 37
    
    q=1/params[1]-1/params[2]
    
    A1 <- -exp(E1/params[1])/params[3]*log(1-exp((E0-E1)*q))
    A0 <- A1*exp((E0-E1)/params[2])
  }
  
  
  #change the name of the parameters, dont know if necessary
  par <- x
  par[(5:8)+30] <- c(E0, E1, A0, A1)
  
  #split the parameters and reconstruct them for the different cultivar
  
  
  
  #loop over the cultivars, calculate predicted days for each cultivar
  rss <- purrr::map_dbl(1:length(SeasonList), function(i){
    
    #extract the parameters
    par_cult <- par[c(i, 11+i, 22+i, 34:42)]
    #predict the bloom
    pred_bloom <- unlist(lapply(X = SeasonList[[i]], FUN = modelfn, par = par_cult))
    pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
    #calculate for each cultivar the rss
    rss <- sum((pred_bloom - bloomJDays[[i]])^2)
    return(rss)
    
  })
  
  
  #calculate the model performance value
  F <- sum(rss)
  
  
  
  #####
  #inequality constraints
  #####
  
  #this is the vector containing the values for the inequality constraints
  #at first initialize the vector
  g <- rep(0,5)
  
  
  #equality constraints should be always stated before inequality constraints, according to meigo vignette
  #(but we dont have any in this case)
  
  
  #inequality constraints are reformulated as differences
  
  #Tu >= Tb
  g[1] <- par[4+30] - par[11+30]
  #Tx >= Tb
  g[2] <- par[10+30] - par[11+30]
  #Tc >= Tu
  g[3] <- par[10+30] - par[4+30]
  
  
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[4] <- exp((10 * par[5+30]) / (297 * 279))
  
  g[5] <- exp((10 * par[6+30]) / (297 * 279))
  
  
  if(return_bloom_days == FALSE){
    #output
    return(list(F=F, g=g))
  } else{
    return(pred_bloom)
  }
}

bloomJDays <- purrr::map(pheno_train, 'pheno')

evaluation_function_meigo_nonlinear_combined(x = par, 
                                             modelfn = custom_PhenoFlex_GDHwrapper, 
                                             bloomJDays = bloomJDays, 
                                             SeasonList = temp_list)


set.seed(123456789)
res_list <- list()
local_solver <- c('SOLNP')

for(solver in local_solver){
  problem<-list(f="evaluation_function_meigo_nonlinear_combined",
                x_0 = x_0,
                x_L = x_L,
                x_U = x_U,
                c_L = c_L, 
                c_U = c_U)
  
  #specify options for the solver
  opts<-list(#maxeval = 1000,
    maxtime = 60 * 60, 
    local_solver = solver, 
    local_bestx = 1,
    inter_save = 0,
    ndiverse = 1000)
  
  res_list[[solver]]<-MEIGO(problem = problem,
                                 opts,
                                 algorithm="ESS", 
                                 modelfn = custom_PhenoFlex_GDHwrapper,
                                 bloomJDays = bloomJDays,
                                 SeasonList = temp_list)
  
} 

plot(c(res_list[[1]]$f, res_list[[1]]$fbest) ~ c(res_list[[1]]$neval, res_list[[1]]$numeval))

source('code/99_helper_functions.R')
get_temp_response_plot_v2(convert_parameters(res_list[[1]]$xbest[c(1,12,23,34:42)]), temp_values = seq(-10, 40, by = 0.1), hourtemps = hourtemps, type = 2)

performance_df <- data.frame()
for(i in 1:length(pheno_train)){
  par <- convert_parameters(res_list[[1]]$xbest[c(i,i+11,i+22,34:42)])
  pred <- c(return_predicted_days(par = par, modelfn = custom_PhenoFlex_GDHwrapper, SeasonList = SeasonList[[i]]),
            pred_eval <- return_predicted_days(par = par, modelfn = custom_PhenoFlex_GDHwrapper, SeasonList = temp_list_eval[[i]]))
  source <- c(rep('training', 20), rep('validation', 7))
  

  performance_df <- rbind(performance_df, 
                          data.frame(cultivar = cultivars[i],
                                     source = source,
                                     year = c(pheno_train[[i]]$Year, pheno_eval[[i]]$Year),
                                     observed = c(pheno_train[[i]]$pheno, pheno_eval[[i]]$pheno),
                                     predicted = pred))
}

performance_df$error <- performance_df$predicted - performance_df$observed


error_df <- performance_df %>% 
  filter(source == 'training') %>% 
  group_by(cultivar) %>% 
  summarise(RMSE_train = round(RMSEP(predicted, observed), digits = 1) ,
            RPIQ_train = round(RPIQ(predicted, observed), digits = 1) ,
            mean_bias_train = round(mean(predicted - observed), digits = 1))

error_df <- performance_df %>% 
  filter(source == 'validation') %>% 
  group_by(cultivar, source) %>% 
  summarise(RMSE_val = round(RMSEP(predicted, observed), digits = 1) ,
            RPIQ_val = round(RPIQ(predicted, observed), digits = 1) ,
            mean_bias_val = round(mean(predicted - observed), digits = 1)) %>% 
  merge(error_df, by = 'cultivar')



performance_df %>% 
  ggplot(aes(x = observed, y = predicted)) + 
  geom_point(aes(col = source)) +
  #geom_point(data = performance_df[performance_df$year == 1997,], shape = 1, col = 'red', size = 3, stroke = 2)+ 
  geom_abline(slope = 1, linetype = 'dashed') + 
  geom_text(aes(x = 107, y = 155, label = 'Training (Validation)'))+
  geom_text(data = error_df, aes(x = 105, y = 150, label = paste0('RMSE: ', RMSE_train, ' (', RMSE_val, ')')))+
  geom_text(data = error_df, aes(x = 105, y = 145, label = paste0('RPIQ: ', RPIQ_train, ' (', RPIQ_val, ')')))+
  geom_text(data = error_df, aes(x = 105, y = 140, label = paste0('bias: ', mean_bias_train, ' (', mean_bias_val, ')')))+
  #geom_text(data = error_df[error_df$type == 'validation',], aes(x = 115, y = 140, label = paste0('(',RMSE, ')')))+
  facet_wrap(~cultivar) +
  ylab('Predicted Bloom (Day of the year)')+
  xlab('Observed Bloom (Day of the year)')+
  theme_bw()


error_df_all <- performance_df %>% 
  filter(source == 'training') %>% 
  summarise(RMSE_train = round(RMSEP(predicted, observed), digits = 1) ,
            RPIQ_train = round(RPIQ(predicted, observed), digits = 1) ,
            mean_bias_train = round(mean(predicted - observed), digits = 1))

error_df_all <- performance_df %>% 
  filter(source == 'validation') %>% 
  summarise(RMSE_val = round(RMSEP(predicted, observed), digits = 1) ,
            RPIQ_val = round(RPIQ(predicted, observed), digits = 1) ,
            mean_bias_val = round(mean(predicted - observed), digits = 1)) %>% 
  cbind(error_df_all)

performance_df %>% 
  ggplot(aes(x = observed, y = predicted)) + 
  geom_point(aes(col = source)) +
  #geom_point(data = performance_df[performance_df$year == 1997,], shape = 1, col = 'red', size = 3, stroke = 2)+ 
  geom_abline(slope = 1, linetype = 'dashed') + 
  geom_text(aes(x = 107, y = 155, label = 'Training (Validation)'))+
  geom_text(data = error_df_all, aes(x = 105, y = 150, label = paste0('RMSE: ', RMSE_train, ' (', RMSE_val, ')')))+
  geom_text(data = error_df_all, aes(x = 105, y = 145, label = paste0('RPIQ: ', RPIQ_train, ' (', RPIQ_val, ')')))+
  geom_text(data = error_df_all, aes(x = 105, y = 140, label = paste0('bias: ', mean_bias_train, ' (', mean_bias_val, ')')))+
  #geom_text(data = error_df[error_df$type == 'validation',], aes(x = 115, y = 140, label = paste0('(',RMSE, ')')))+
  ylab('Predicted Bloom (Day of the year)')+
  xlab('Observed Bloom (Day of the year)')+
  theme_bw()

#what are the seasons which lead to the strong differences in model performance?
performance_df %>% 
  filter(abs(error) > 10) %>% 
  group_by(year) %>% 
  summarize(n())


years <- c(1997, 1998, 1991, 1990, 2008)

temp_list_invest <- chillR::genSeasonList(hourtemps, years = years)

accumulation_df <- data.frame()
#run for each cultivar the phenoflex
for(i in 1:length(cultivars)){
  par <- convert_parameters(res_list[[1]]$xbest[c(i,i+11,i+22,34:42)])
  
  for(j in 1:length(years)){
    out <- PhenoFlex(temp_list_invest[[j]]$Temp, times = 1:nrow(temp_list_invest[[j]]), 
                     yc = par[1], zc = par[2], s1 = par[3], 
                     Tu = par[4], E0 = par[5], E1 = par[6], 
                     A0 = par[7], A1 = par[8], Tf = par[9],
                     Tc = par[10], Tb = par[11], slope = par[12], basic_output = FALSE, stopatzc=TRUE)
    
    intermed_df <- data.frame(season = years[j],
                              cultivar = cultivars[i],
                              JDay = temp_list_invest[[j]]$JDay,
                              hour = 0:23,
                              year = temp_list_invest[[j]]$Year,
                              fake_year = c(rep(2020, 3672), rep(2021, nrow(temp_list_invest[[j]]) - 3672)),
                              x = out$x,
                              y = out$y,
                              z = out$z)
    
    #remove chill accumulatio observations after bloom
    intermed_df$y[out$bloomindex:nrow(intermed_df)] <- NA
    intermed_df$z[out$bloomindex:nrow(intermed_df)] <- NA
    
    accumulation_df <- rbind(accumulation_df,intermed_df)
    
    #unitl row 3672 always year 1, rest is year 2
  }

  
  
}

accumulation_df$Date <- lubridate::parse_date_time(x = paste(accumulation_df$year, accumulation_df$JDay, accumulation_df$hour), orders = "yjh")
accumulation_df$fake_Date <- lubridate::parse_date_time(x = paste(accumulation_df$fake_year, accumulation_df$JDay, accumulation_df$hour), orders = "yjh")
accumulation_df$id <- paste0(accumulation_df$season, accumulation_df$cultivar)

accumulation_df_long <- reshape2::melt(accumulation_df, measure.vars = c('x', 'y', 'z'))

performance_df$obs_fake_date <- lubridate::parse_date_time(x = paste(2021, performance_df$observed, 12), orders = "yjh")
performance_df$pred_fake_date <- lubridate::parse_date_time(x = paste(2021, floor(performance_df$predicted), 12), orders = "yjh")
performance_df$season <- performance_df$year

accumulation_df_long %>% 
  filter(variable == 'z',
         fake_year == 2021,
         fake_Date > as.Date('2021-03-15'),
         fake_Date < as.Date('2021-06-07')) %>% 
  ggplot(aes(x = fake_Date, y = value, group = id, col = as.factor(season)))+
  geom_line() +
  geom_vline(data = performance_df[performance_df$year %in% years,], aes(xintercept = obs_fake_date, col = as.factor(season), linetype = 'observed')) +
  geom_vline(data = performance_df[performance_df$year %in% years,], aes(xintercept = pred_fake_date, col = as.factor(season), linetype = 'predicted')) +
  facet_grid(cultivar~variable) +
  theme_bw()

res_list[[1]]$xbest

save.image(file = "combined_fitting.RData")
