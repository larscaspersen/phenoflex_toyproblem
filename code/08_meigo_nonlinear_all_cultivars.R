library(chillR)
library(MEIGOR)
library(tidyverse)
library(nleqslv)

source('code/04_functions_meigo.R')

#read the training data of all cultivars
files <- list.files('data/training/')
files <- paste0('data/training/', files)
pheno_train <- lapply(files, read.csv)

#read evaluation data
files <- list.files('data/evaluation/')
files <- paste0('data/evaluation/', files)
pheno_eval <- lapply(files, read.csv)

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

temp_list_eval <- list()

#generate list of temperatures for each cultivar
for(i in 1:length(pheno_eval)){
  temp_list_eval[[i]] <- chillR::genSeasonList(hourtemps, years = pheno_eval[[i]]$Year)
}


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
local_solver <- c('DHC')

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

purrr::map(res_list, function(x) x[['DHC']]$fbest)

source('code/99_helper_functions.R')

purrr::map(res_list, function(x){
  
  get_temp_response_plot_v2(par = x[['DHC']]$xbest, temp_values = seq(-10,40, by = 0.1), par_type = 'new', log_A = FALSE, hourtemps = hourtemps)
  
} )

names(res_list) <- cultivars

f_df <- data.frame(NULL)

for(cult in cultivars){
  for(solver in c('SOLNP', 'DHC')){
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



#saveRDS(res_list, 'result_r1_nonlin.RData')


for(i in 1:length(pheno_train)){
  
  for(solver in c('SOLNP')){
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


for(cult in cultivars){
  for(solver in c('SOLNP', 'DHC')){
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



#saveRDS(res_list, 'result_r1_nonlin.RData')
performance_df <- data.frame(NULL)

for(i in 1:length(res_list)){
  
  for(solver in c('SOLNP', 'DHC')){
    #solver <- 'DHC'
    l <-  length(c(res_list[[i]][[solver]]$f, res_list[[i]][[solver]]$fbest))
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
      pred_train <- evaluation_function_meigo(x = par, 
                                              modelfn = custom_PhenoFlex_GDHwrapper,
                                              bloomJDays = pheno_train[[i]]$pheno,
                                              SeasonList = temp_list[[i]], 
                                              return_bloom_days = TRUE)
      
      pred_eval <- evaluation_function_meigo(x = par, 
                                             modelfn = custom_PhenoFlex_GDHwrapper,
                                             bloomJDays = pheno_eval[[i]]$pheno,
                                             SeasonList = temp_list_eval[[i]], 
                                             return_bloom_days = TRUE)
      
      
      
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
  geom_point(aes(y = rmsep_train, col = 'training')) +
  geom_line(aes(y = rmsep_train, col = 'training')) +
  #geom_point(aes(y = rmsep_eval, col = 'eval')) +
  #geom_line(aes(y = rmsep_eval, col = 'eval')) +
  facet_grid(solver~cultivar)+
  ylab('RMSE (days)') +
  xlab(expression ("Number of Evaluations * "~10^-2))+
  scale_color_manual(breaks = 'training', values = 'Firebrick')+
  #xlab('Number of Evaluations * 10^-2')
  theme_bw()
ggsave('fitting_nonlinear_train.jpeg', height = 10, width = 20, units = 'cm',
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
  #xlab('Number of Evaluations * 10^-2')
  theme_bw()
ggsave('fitting_nonlinear_both.jpeg', height = 10, width = 20, units = 'cm',
       device = 'jpeg')

performance_df %>% 
  ggplot(aes(x = eval)) +
  geom_point(aes(y = rpiq_train, col = 'train')) +
  geom_line(aes(y = rpiq_train, col = 'train')) +
  geom_point(aes(y = rpiq_eval, col = 'eval')) +
  geom_line(aes(y = rpiq_eval, col = 'eval')) +
  facet_grid(solver~cultivar)

performance_df %>% 
  ggplot(aes(x = eval)) +
  geom_point(aes(y = bias_train, col = 'train')) +
  geom_line(aes(y = bias_train, col = 'train')) +
  geom_point(aes(y = bias_eval, col = 'eval')) +
  geom_line(aes(y = bias_eval, col = 'eval')) +
  facet_grid(solver~cultivar)

performance_df %>% 
  reshape2::melt(measure.vars = c('yc', 'zc', 's1', 'Tu', 'theta_star', 'theta_c', 'tau', 'pie_c',   'Tf', 'Tc', 'Tb',  'slope')) %>% 
  ggplot(aes(x = eval, y = value, col = cultivar)) +
  geom_point() +
  geom_line() +
  facet_wrap(solver~variable, scales = 'free_y')




#read hajars parameters
files <- list.files('data/hajar_parameters/', full.names = TRUE)
par_list <- lapply(files,read.csv) 



par_df <- data.frame()

for(i in 1:length(res_list)){
  
  par_df <- rbind(par_df ,
                 data.frame(cultivar = cultivars[i],
                            hajar = par_list[[i]]$Last_run,
                            meigo = convert_parameters(res_list[[i]]$DHC$xbest),
                            parameter = c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1',   'Tf', 'Tc', 'Tb',  'slope'))
  )
  
}

par_df$cult_abb <- c('B', 'Cl', 'Coa', 'Col', 'D', 'Pez', 'Per', 'R', 'T', 'V', 'X')
par_df$parameter <- factor(par_df$parameter,  
                          levels = c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1',   'Tf', 'Tc', 'Tb',  'slope'))

par_df %>% 
  ggplot(aes(x = cult_abb)) + 
  geom_point(aes(y = hajar, col = 'hajar')) + 
  geom_point(aes(y = meigo, col = 'meigo')) + 
  facet_wrap(~parameter, scales = 'free_y')+
  theme_bw()
ggsave('compare_hajar_meigo.jpeg', height = 10, width = 20, units = 'cm',
       device = 'jpeg')


#see how the evaluation performance for hajar is

for(i in 1:length(pheno_eval)){
  pred_eval <- return_predicted_days(par = par_list[[i]]$Last_run, 
                                     modelfn = custom_PhenoFlex_GDHwrapper_no_control, 
                                     SeasonList = temp_list_eval[[i]])
  print(cultivars[i])
  print(paste0('Number runs: ', ))
  print(RMSEP(pred_eval, pheno_eval[[i]]$pheno))
  print('----')
}

pred_eval <- return_predicted_days(par = par_list[[3]]$Last_run, 
                                   modelfn = custom_PhenoFlex_GDHwrapper_no_control, 
                                   SeasonList = temp_list_eval[[3]])

RMSEP(pred_eval, pheno_eval[[3]]$pheno)
