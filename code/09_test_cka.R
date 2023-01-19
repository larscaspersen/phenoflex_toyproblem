#try out the fitting for the klein altendorf data
library(chillR)
library(MEIGOR)
library(tidyverse)
library(nleqslv)

# #read phenological data
# boskop_pheno <- read.csv('data/cka_data/Boskoop_first_bloom.csv')
# #split into training and evaluation
# 
# n_train <- floor(nrow(boskop_pheno) * 0.75)
# train_years_boskop <- sample(boskop_pheno$Year, replace = FALSE, size = n_train)
# #order the years
# train_years_boskop <- sort(train_years_boskop)
# boskop_train <- boskop_pheno[boskop_pheno$Year %in% train_years_boskop,]
# boskop_eval <- boskop_pheno[!boskop_pheno$Year %in% train_years_boskop,]
# 
# write.csv(boskop_train, 'data/cka_data/Boskop_train.csv', row.names = FALSE)
# write.csv(boskop_eval, 'data/cka_data/Boskop_eval.csv', row.names = FALSE)
# 
# 
# 
# #read phenological data
# alex_pheno <- read.csv('data/cka_data/AlexanderLukas_first_bloom.csv')
# #split into training and evaluation
# 
# n_train <- floor(nrow(alex_pheno) * 0.75)
# train_years_alex <- sample(alex_pheno$Year, replace = FALSE, size = n_train)
# #order the years
# train_years_alex <- sort(train_years_alex)
# alex_train <- boskop_pheno[boskop_pheno$Year %in% train_years_alex,]
# alex_eval <- boskop_pheno[!boskop_pheno$Year %in% train_years_alex,]
# 
# write.csv(alex_train, 'data/cka_data/Alex_train.csv', row.names = FALSE)
# write.csv(alex_eval, 'data/cka_data/Alex_eval.csv', row.names = FALSE)

#read training and evaluation data
boskop_train <- read.csv('data/cka_data/Boskop_train.csv')
boskop_eval <- read.csv('data/cka_data/Boskop_eval.csv')

alex_train <- read.csv('data/cka_data/Alex_train.csv')
alex_eval <- read.csv('data/cka_data/Alex_eval.csv')


#read temperature data
temp <- read.csv('data/cka_data/TMaxTMin1958-2018_patched.csv')

#make weather hourly
hourtemps <- chillR::stack_hourly_temps(weather = temp, latitude = 50.4)$hourtemps

#split weather data by season and training / evaluation
SeasonList_boskop <- chillR::genSeasonList(hourtemps, years = boskop_train$Year)
SeasonList_boskop_eval <- chillR::genSeasonList(hourtemps, years = boskop_eval$Year)

SeasonList_alex <- chillR::genSeasonList(hourtemps, years = alex_train$Year)
SeasonList_alex_eval <- chillR::genSeasonList(hourtemps, years = alex_eval$Year)

#adjust column names
colnames(alex_train) <- colnames(alex_eval) <- colnames(boskop_train) <- colnames(boskop_eval) <- c('Year', 'pheno')


#read the needed fitting functions
source('code/04_functions_meigo.R')

#define boundaries and starting point
#      ('yc', 'zc', 's1', 'Tu', 'theta*', 'theta_c', 'Tau(thetha*)', 'pie_c',   'Tf', 'Tc', 'Tb',  'slope')
x_0 <- c(40,   190,   0.5,  25,   278.1,    285.1,   47.7,           30.6,        4,   36,     4,    1.60)
x_U <- c(80,   500,   1.0,  30,   281,      287,       48,             32,       10,   40,    10,    5.00)
x_L <- c(20,   100,   0.1,   0,   279,      284,       16,             24,        0,    0,     0,    0.05)

#limits for the inequality constraints
#         #gdh parameters   #q10 for E0 and E1
c_L <- c(  0,   0,   0,     1.5, 1.5)
c_U <- c(Inf, Inf, Inf,     3.5, 3.5)


set.seed(123456789)
res_list <- list()
#local_solver <- c('BFGS','SA', 'SOLNP', 'DHC')
local_solver <- c('SA')

#bind pear and apple data together
pheno_train <- list(boskop_train, alex_train)
temp_list <- list(SeasonList_boskop, SeasonList_alex)

#iterate over species / cultivars
for(i in 1:length(pheno_train)){
  
  #imake new element
  res_list[[i]] <- list()
  
  #iterate over solver
  for(solver in local_solver){
    #define problem
    problem<-list(f="evaluation_function_meigo_nonlinear",
                  x_0 = x_0,
                  x_L = x_L,
                  x_U = x_U,
                  c_L = c_L, 
                  c_U = c_U)
    
    #specify options for the solver
    opts<-list(#maxeval = 1000,
      maxtime = 60 * 5, 
      local_solver = solver, 
      local_bestx = 1,
      inter_save = 0)
    
    #save results
    res_list[[i]][[solver]]<-MEIGO(problem,
                                   opts,
                                   algorithm="ESS", 
                                   modelfn = custom_PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_train[[i]]$pheno,
                                   SeasonList = temp_list[[i]])
    
  } 
}


cultivars <- c('Boskop', 'Alexander Lucas')

#element for to save the fitting performance results
performance_df <- data.frame(NULL)

#iterate over cultivars
for(i in 1:length(res_list)){
  
  #iterate over solver
  for(solver in local_solver){
    
    #find out how many intermediate fitting elements in the results (intermediate steps)
    l <-  length(c(res_list[[i]][[solver]]$f, res_list[[i]][[solver]]$fbest))
    #iterate over intermediate steps
    for(j in 1:l){
      
      #if it is the last item of the intermediate saves, then use the best one
      if(j == l){
        par <- res_list[[i]][[solver]]$xbest
        f <- res_list[[i]][[solver]]$fbest
        eval <- res_list[[i]][[solver]]$numeval
      } else {
        par <- res_list[[i]][[solver]]$x[j,]
        f <- res_list[[i]][[solver]]$f[j]
        eval <- res_list[[i]][[solver]]$neval[j]
      }
      
      #make predictions
      pred_train <- return_predicted_days(convert_parameters(par), 
                                          modelfn = custom_PhenoFlex_GDHwrapper,
                                          SeasonList = temp_list[[i]])
      
      pred_eval <- return_predicted_days(convert_parameters(par), 
                                         modelfn = custom_PhenoFlex_GDHwrapper,
                                         SeasonList = temp_list_eval[[i]])
      
      
      
      performance_df <- rbind(performance_df, 
                              data.frame(cultivar = cultivars[i],
                                         solver = solver,
                                         run = 'run_1',
                                         eval = eval,
                                         f = f,
                                         bias_train = mean(pred_train - pheno_train[[i]]$pheno),
                                         bias_eval = mean(pred_eval - pheno_eval[[i]]$pheno),
                                         rmsep_train = RMSEP(pred_train, pheno_train[[i]]$pheno),
                                         rmsep_eval = RMSEP(pred_eval, pheno_eval[[i]]$pheno),
                                         rpiq_train = RPIQ(pred_train, pheno_train[[i]]$pheno),
                                         rpiq_eval = RPIQ(pred_eval, pheno_eval[[i]]$pheno),
                                         yc = par[1],
                                         zc = par[2],
                                         s1 = par[3],
                                         Tu = par[4],
                                         theta_star = par[5],
                                         theta_c = par[6],
                                         tau = par[7],
                                         pie_c = par[8],
                                         Tf = par[9],
                                         Tc = par[10],
                                         Tb = par[11],
                                         slope = par[12]
                              ))
      
      
    }
  }
  
}

performance_df %>% 
  ggplot(aes(x = eval/1000)) +
  geom_point(aes(y = f, col = 'training')) +
  geom_line(aes(y = f, col = 'training')) +
  #geom_point(aes(y = rmsep_eval, col = 'eval')) +
  #geom_line(aes(y = rmsep_eval, col = 'eval')) +
  facet_grid(solver~cultivar)+
  ylab('Residual sum of squares (RSS)') +
  xlab(expression ("Number of Evaluations * "~10^-2))+
  scale_color_manual(breaks = 'training', values = 'Firebrick')+
  #xlab('Number of Evaluations * 10^-2')
  theme_bw()
ggsave('fitting_CKA_rss.jpeg', height = 10, width = 20, units = 'cm',
       device = 'jpeg')


performance_df %>% 
  ggplot(aes(x = eval/1000)) +
  geom_point(aes(y = rmsep_train, col = 'training')) +
  geom_line(aes(y = rmsep_train, col = 'training')) +
  geom_point(aes(y = rmsep_eval, col = 'evaluation')) +
  geom_line(aes(y = rmsep_eval, col = 'evaluation')) +
  facet_grid(solver~cultivar)+
  ylab('RMSE (days)') +
  xlab(expression ("Number of Evaluations * "~10^-2))+
  scale_color_manual(breaks = c('training', 'evaluation'), 
                     values = c('Firebrick', 'Steelblue'))+
  theme_bw()
ggsave('fitting_CKA.jpeg', height = 10, width = 20, units = 'cm',
       device = 'jpeg')

source('code/99_helper_functions.R')
temp_val <- seq(-30,40, by = 0.1)
get_temp_response_plot_v2(par = res_list[[1]][[3]]$xbest, 
                          temp_values = temp_val, 
                          par_type = 'new', hourtemps = hourtemps)

get_temp_response_plot_v2(par = res_list[[1]][[3]]$x[5,], 
                          temp_values = temp_val, 
                          par_type = 'new', hourtemps = hourtemps)

get_temp_response_plot_v2(par = res_list[[2]][[4]]$xbest, 
                          temp_values = temp_val, 
                          par_type = 'new', hourtemps = hourtemps)
