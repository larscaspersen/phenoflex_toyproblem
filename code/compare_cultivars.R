#read the different training and validation splits of the 
#data set by alvaro

library(chillR)
library(MEIGOR)
library(tidyverse)

#read the training data of all cultivars
files <- list.files('data/training/')
files <- paste0('data/training/', files)
pheno_train <- lapply(files, read.csv)

cultivars <- c('Blanquina', 'Clara', 'Collaos', 'Coloradona',
               'DelaRiega', 'Perezosa', 'Perico', 'Raxao',
               'Teorica', 'Verdialona', 'Xuanina')

names(pheno_train) <- cultivars

#read temperature data
hourtemps <- read.csv('data/hourtemps.csv')

temp_list <- list()

#generate list of temperatures for each cultivar
for(i in 1:length(pheno_train)){
  temp_list[[i]] <- chillR::genSeasonList(hourtemps, years = pheno_train[[i]]$Year)
}

#load meigo functions
source('code/04_functions_meigo.R')

#limits and start point of parameters
x_0 <- c(40, 190, 0.5, 25, 3372.8,  9900.3, log(6319.5), log(5.939917e13),  4, 36,  4,  1.60)
x_U <- c(80, 500, 1.0, 30, 6000.0, 12000.0, log(1e15),    log(1e19),       10, 40, 10, 5.00)
x_L <- c(20, 100, 0.1, 0 , 2000.0, 9000.0, log(1e3),      log(1e13),        0,  0,  0,  0.05)

#limits for the inequality constraints
#         #gdh parameters   #q10 for E0 and E1
c_L <- c(  0,   0,   0,     1.5, 1.5)
c_U <- c(Inf, Inf, Inf,     3.5, 3.5)


set.seed(123456789)
res_list <- list()
local_solver <- c('BFGS','SA', 'SOLNP', 'DHC')

for(i in 1:length(pheno_train)){
  
  res_list[[i]] <- list()
  
 for(solver in local_solver){
   problem<-list(f="evaluation_function_meigo",
                 x_0 = x_0,
                 x_L = x_L,
                 x_U = x_U,
                 c_L = c_L, 
                 c_U = c_U)
   
   #specify options for the solver
   opts<-list(#maxeval = 1000,
     maxtime = 60 * 5, 
     local_solver = solver, 
     log_var = 7:8,
     local_bestx = 1,
     inter_save = 0)
   
   res_list[[i]][[solver]]<-MEIGO(problem,
                                  opts,
                                  algorithm="ESS", 
                                  modelfn = custom_PhenoFlex_GDHwrapper,
                                  bloomJDays = pheno_train[[i]]$pheno,
                                  SeasonList = temp_list[[i]])
   
 } 
}

names(res_list) <- cultivars
f_df <- data.frame(NULL)

for(cult in cultivars){
  for(solver in local_solver){
    f_df <- rbind(f_df,
                  data.frame(f =  c(res_list[[cult]][[solver]]$f, res_list[[cult]][[solver]]$fbest),
                             eval = c(res_list[[cult]][[solver]]$neval, res_list[[cult]][[solver]]$numeval),
                             cultivar = cult,
                             solver = solver))
    
  }
}

f_df %>% 
  ggplot(aes(x = eval, y = f, col = solver)) +
  geom_point() +
  geom_line() +
  facet_grid(~cultivar) +
  coord_cartesian(ylim = c(0,1700))


temp_values <- seq(-10, 40, by = 0.1)

source('code/99_helper_functions.R')

res_list[['Perico']]

par <- res_list[['Xuanina']][['SOLNP']]$xbest
#collaos bfgs looks not good
#rteorica sa looks not good, DHC looks really slim, solnp looks okay
#all models kinda bad for Perico, but looks reasonable temperature response plot
#verdialona bfgs looks bad, solnp also looks shifted, chose DHC
#xuanina looks a bit shifted
get_temp_response_plot(par = par, temp_values = temp_values, log_A = TRUE, hourtemps = hourtemps)
#save.image(file = 'cultivar_comparison.RData')

load('cultivar_comparison.RData')

#choose the starting points for all of the cultivars
start_list <- list()
start_list[['Blanquina']] <- res_list[['Blanquina']][['DHC']]$xbest
start_list[['Clara']] <- res_list[['Clara']][['DHC']]$xbest
start_list[['Collaos']] <- res_list[['Collaos']][['SOLNP']]$xbest
start_list[['Coloradona']] <- res_list[['Coloradona']][['SA']]$xbest
start_list[['DelaRiega']] <- res_list[['DelaRiega']][['DHC']]$xbest
start_list[['Perezosa']] <- res_list[['Perezosa']][['DHC']]$xbest
start_list[['Perico']] <- res_list[['Perico']][['SOLNP']]$xbest
start_list[['Raxao']] <- res_list[['Raxao']][['SA']]$xbest
start_list[['Teorica']] <- res_list[['Teorica']][['SOLNP']]$xbest
start_list[['Verdialona']] <- res_list[['Verdialona']][['DHC']]$xbest
start_list[['Xuanina']] <- res_list[['Xuanina']][['SOLNP']]$xbest

names(temp_list) <- cultivars

#start again 
res_list_r2 <- list()
local_solver <- c('BFGS','SA', 'SOLNP', 'DHC')

for(cult in cultivars){
  
  res_list_r2[[cult]] <- list()
  
  x_start <- start_list[[cult]]
  
  for(solver in local_solver){
    problem<-list(f="evaluation_function_meigo",
                  x_0 = x_start,
                  x_L = x_L,
                  x_U = x_U,
                  c_L = c_L, 
                  c_U = c_U)
    
    #specify options for the solver
    opts<-list(#maxeval = 1000,
      maxtime = 60 * 5, 
      local_solver = solver, 
      log_var = 7:8,
      local_bestx = 1,
      inter_save = 0)
    
    res_list_r2[[cult]][[solver]]<-MEIGO(problem,
                                   opts,
                                   algorithm="ESS", 
                                   modelfn = custom_PhenoFlex_GDHwrapper,
                                   bloomJDays = pheno_train[[cult]]$pheno,
                                   SeasonList = temp_list[[cult]])
    
  } 
}


f_df_r2 <- data.frame(NULL)

for(cult in cultivars){
  for(solver in local_solver){
    f_df_r2 <- rbind(f_df_r2,
                  data.frame(f =  c(res_list_r2[[cult]][[solver]]$f, res_list_r2[[cult]][[solver]]$fbest),
                             eval = c(res_list_r2[[cult]][[solver]]$neval, res_list_r2[[cult]][[solver]]$numeval),
                             cultivar = cult,
                             solver = solver))
    
  }
}

f_df_r2 %>% 
  ggplot(aes(x = eval, y = f, col = solver)) +
  geom_point() +
  geom_line() +
  facet_grid(~cultivar) +
  coord_cartesian(ylim = c(0,1700))
