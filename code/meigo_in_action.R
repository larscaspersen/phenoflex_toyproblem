# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("MEIGOR")
library(MEIGOR)
library(chillR)
library(tidyverse)


#read the phenology data
bloom <- readxl::read_excel('data/Blanquina.xlsx')
#read temperature data
hourtemps <- read.csv('data/hourtemps.csv')

#remove years without observations
bloom <- na.omit(bloom)
bloom$Year <- as.numeric(bloom$Year)

#split the data into training and evaluation
n_train <- floor(nrow(bloom)*0.75)
years <- sample(bloom$Year, n_train, replace = FALSE)

#extract the days of bloom
bloomJDays <- c()
for(year in years){
  bloomJDays <- c(bloomJDays,bloom[bloom$Year == year,]$pheno)
}

#evalulation
eval_years <-  bloom$Year[!(bloom$Year %in% years)]
eval_bloomJDays <- c()
for(year in eval_years){
  eval_bloomJDays <- c(eval_bloomJDays,bloom[bloom$Year == year,]$pheno)
}

#create a list of the temperature data split by each chill season
SeasonList <- chillR::genSeasonList(hourtemps, years = years)
eval_SeasonList <- chillR::genSeasonList(hourtemps, years = eval_years)

#this function takes the weather data (x) and the model parameters (par) and returns the expected day of bloom 
custom_PhenoFlex_GDHwrapper <- function (x, par){
  #x is one of the elements in season List
  #par are the parameters of the model
  
  #transfrom log to normla value of parameters
  par[7]<-exp(par[7])
  par[8]<-exp(par[8])
  
  #in case the parameters do not make sense, return NA
  if (par[4] <= par[11]){
    return(NA)
  } else if(par[10] <= par[4]){
    return(NA)
  } else if(exp((10 * par[5]) / (297 * 279)) < 1.5 |  exp((10 * par[5]) / (297 * 279)) > 3.5 ){
    #also when the q10 criterion is outside the limits of 1.5 and 3.5
    return(NA)
  } else if(exp((10 * par[6]) / (297 * 279)) < 1.5 |  exp((10 * par[6]) / (297 * 279)) > 3.5 ){
    return(NA)
  } else{
    #calculat the bloom day
    bloomindex <- PhenoFlex(temp = x$Temp, times = seq_along(x$Temp), 
                            yc = par[1], zc = par[2], s1 = par[3], Tu = par[4], E0 = par[5], 
                            E1 = par[6], A0 = par[7], A1 = par[8], Tf = par[9], Tc = par[10], 
                            Tb = par[11], slope = par[12], Imodel = 0L, basic_output = TRUE)$bloomindex
    
    #return values
    if (bloomindex == 0){
      return(NA)
    } 
      
    JDay <- x$JDay[bloomindex]
    JDaylist <- which(x$JDay == JDay)
    n <- length(JDaylist)
    if (n == 1){
      return(JDay)
    } 
    return(JDay + which(JDaylist == bloomindex)/n - 1/(n/ceiling(n/2)))
  }

}


#in meigo the start guess is x_0, the upper range is x_U and the lower range is x_L
#       yc,  zc,  s1, Tu,      E0,      E1,          A0,                A1, Tf, Tc, Tb,  slope
x_0 <- c(40, 190, 0.5, 25, 3372.8,  9900.3, log(6319.5), log(5.939917e13),  4, 36,  4,  1.60)
x_U <- c(80, 500, 1.0, 30, 6000.0, 12000.0, log(1e15),    log(1e19),       10, 40, 10, 5.00)
x_L <- c(20, 100, 0.1, 0 , 2000.0, 9000.0, log(1e3),      log(1e13),        0,  0,  0,  0.05)

#limits for the inequality constraints
#         #gdh parameters   #q10 for E0 and E1
c_L <- c(  0,   0,   0,     1.5, 1.5)
c_U <- c(Inf, Inf, Inf,     3.5, 3.5)


#this function creates the model performance the optimizer "sees"
#it takes the parameters (x), the function to create the predictions (modelfn), the actual observed bloom days (bloomJdays) and the temperatures orginized
#in a list by the seasons (years)
#usually it returns the model performance expressed in RSS but it can also simply return the predicted bloom days if "return_bloom_days" is set TRUE

evaluation_function_meigo <- function(x, 
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
  
  #change the name of the parameters, dont know if necessary
  par <- x
  
  #calculate the predicted flower dates
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, par = par))
  
  #if the model returns no bloom day, then give penalty term instead
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
  
  #calculate the model performance value
  F <- sum((pred_bloom - bloomJDays)^2)
  
  
  
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
  g[1] <- x[4] - x[11]
  #Tx >= Tb
  g[2] <- x[10] - x[11]
  #Tc >= Tu
  g[3] <- x[10] - x[4]
  
  
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[4] <- exp((10 * x[5]) / (297 * 279))
  
  g[5] <- exp((10 * x[6]) / (297 * 279))
  
  
  if(return_bloom_days == FALSE){
    #output
    return(list(F=F, g=g))
  } else{
    return(pred_bloom)
  }

}

#here I mainly follow the way MEIGO needs the optimization problem to be specified
#it needs:
#   - the function that takes the parameters and returns the model performance
#   - first guess of the parameters (but we could leave it blank)
#   - upper (x_U) and and lower (x_L) bounds of the parameters
#   - upper (c_U) and lower (c_L) bounds of the constraints (they are specified in the evaluation function)


#=========================
#PROBLEM SPECIFICATIONS
#=========================
problem<-list(f="evaluation_function_meigo",
              x_0 = x_0,
              x_L = x_L,
              x_U = x_U,
              c_L = c_L, 
              c_U = c_U)

#specify options for the solver
opts<-list(#maxeval = 1000,
  maxtime = 60 * 2, 
  local_solver = 'DHC', 
  log_var = 7:8,
  local_bestx = 1,
  inter_save = 0)

Results_meigo<-MEIGO(problem,
                     opts,
                     algorithm="ESS", 
                     modelfn = custom_PhenoFlex_GDHwrapper,
                     bloomJDays = bloomJDays,
                     SeasonList = SeasonList)


#use evaluation function to return predicted days instead of f
pred_days <- evaluation_function_meigo(x =  Results_meigo$xbest,
                                       modelfn = custom_PhenoFlex_GDHwrapper,
                                       bloomJDays = bloomJDays,
                                       SeasonList = SeasonList, 
                                       return_bloom_days = TRUE)

chillR::RPIQ(predicted = pred_days, observed = bloomJDays)
#baseline 0.32
#everything below 1 is extremely bad, better would be values of around 3
chillR::RMSEP(predicted = pred_days, observed = bloomJDays)
#we aim for RMSE of 3 to 4 days....

#get function to create temperature response plot
source('code/99_helper_functions.R')
get_temp_response_plot(Results_meigo$xbest, temp_values = seq(-10,40, by = 0.1), log_A = TRUE)









#down below I was doing a more intensive screening on different start points, different solvers and option settings




#object in which I want to store the results
results_list <- list()

#set of local optimizer we can choose from in MEIGO
local_solver <- c('NM', 'BFGS', 'CG', 'LBFGSB', 'SA', 'SOLNP', 'DHC')

#JAE perform n_init initial simulations to find nice initial guesses for the parameters 
n_init <- 10000

#Initialize a matrix of initial solutions with zeros and add the first solution
#with our own initial guess
sol_x <- matrix(0, n_init, 12)
sol_x[1,]<- x_0
sol_f <- rep(Inf, n_init)
sol_f[1] <- evaluation_function_meigo(sol_x[1,], modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, SeasonList = SeasonList)$F

#create random sets of parameters combinations
for (i in 2:n_init){
  sol_x[i,] <- runif(12)*(x_U-x_L)+x_L
  sol_f[i] <- evaluation_function_meigo(x = sol_x[i,], modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, SeasonList = SeasonList)$F
}

#take the ten best start points
n_best <- 10
i_best <- order(unlist(sol_f))[1:n_best]
sol_x[i_best,]
sol_f[i_best]

#I want to experiment on different starting points, that is why I put them in a list
start_points <- lapply(i_best, function(i) sol_x[i,])
#add the traditional start point aswell

write.csv(start_points, 'start_points.csv', row.names = FALSE)
write.csv(bloomJDays, 'bloomJDays.csv', row.names = FALSE)
write.csv(years, 'years.csv', row.names = FALSE)


#iterate other the start points
for(i in 5:7){
#for(i in 1:length(start_points)){
  start_parameters <- start_points[[i]]
  
  #=========================
  #PROBLEM SPECIFICATIONS
  #=========================
  problem<-list(f="evaluation_function_meigo",
                x_0 = start_parameters,
                x_L = x_L,
                x_U = x_U,
                c_L = c_L, 
                c_U = c_U)
  
  results_list[[i]] <- list()
  
  #iterate over the solver
  for(solver in local_solver){
    
    evaluation_function_meigo(start_parameters, modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, SeasonList = SeasonList)$g
    
    #specify options for the solver
    opts<-list(#maxeval = 1000,
               maxtime = 60 * 2, 
               local_solver = solver, 
               log_var = 7:8,
               local_bestx = 1,
               inter_save = 0)
    
    Results_meigo<-MEIGO(problem,
                              opts,
                              algorithm="ESS", 
                              modelfn = custom_PhenoFlex_GDHwrapper,
                              bloomJDays = bloomJDays,
                              SeasonList = SeasonList)
    
    results_list[[i]][[solver]] <- Results_meigo
    
    # pred_days <- evaluation_function_meigo(x =  Results_meigo$xbest,
    #                           modelfn = custom_PhenoFlex_GDHwrapper,
    #                           bloomJDays = bloomJDays,
    #                           SeasonList = SeasonList, 
    #                           return_bloom_days = TRUE)
    # 
    # chillR::RPIQ(predicted = pred_days, observed = bloomJDays)
    # chillR::RMSEP(predicted = pred_days, observed = bloomJDays)

    
  }
  

  
}

results_list[[1]]

#evaluation of the fitted parameters
result_list_8to10 <- readRDS("result_run1_8to10.rds")
result_list_combined <- c(results_list, result_list_8to10[8:10])


fbest <- purrr::map(result_list_combined, function(x){
  purrr::map(x, function(y){
    data.frame(
      f = c(y$f, y$fbest),
      neval = c(y$neval, y$numeval)
    )
  }) %>% 
    bind_rows(.id = 'solver')
  
}) %>% 
  bind_rows(.id = 'start_point')


# ggplot(fbest_long, aes(x = solver, y = f)) +
#   geom_point() +
#   facet_grid(~start_point)
# 
# #get rank per run
# fbest_long %>% 
#   group_by(start_point) %>% 
#   mutate(rank = rank(f)) %>% 
#   ungroup() %>% 
#   group_by(solver) %>% 
#   summarise(mean_rank = mean(rank),
#             median_rank = median(rank),
#             sd_rank = sd(rank))


#maybe take the 10 best intermediate results and use that as a new starting point?
top_results <- fbest %>% 
  group_by(solver, start_point) %>% 
  summarise(f = min(f),
            neval = max(neval)) %>% 
  ungroup() %>% 
  slice_min(order_by = f, n = 10)


ggplot(fbest, aes(x = neval, y = f, col = as.factor(start_point))) + 
  geom_point() +
  geom_line() +
  facet_grid(~solver) +
  coord_cartesian(ylim = c(0, 5000))


ggplot(fbest, aes(x = neval/1000, y = f, col = solver)) + 
  geom_point() +
  geom_line() +
#  geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~as.factor(start_point)) +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  theme_bw()

ggsave('moderl_perofrmance_r1_total.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')


ggplot(fbest, aes(x = neval/1000, y = f, col = solver)) + 
  geom_point() +
  geom_line() +
#  geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~as.factor(start_point)) +
  coord_cartesian(ylim = c(0, 2000)) +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  theme_bw()

ggsave('moderl_perofrmance_r1_zoom.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')



ggplot(fbest, aes(x = neval/1000, y = f, col = solver)) + 
  geom_point() +
  geom_line() +
  geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~as.factor(start_point)) +
  coord_cartesian(ylim = c(0, 2000)) +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  theme_bw()

ggsave('moderl_perofrmance_r1_zoom_highlight.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')





saveRDS(result_list_combined, 'results_fit_run1.rds')

#start points for second run of fitting

start_point_r2 <- list()

for(i in seq_along(top_results$f)){
  start <- as.numeric(top_results$start_point[i])
  start_point_r2[[i]] <- result_list_combined[[start]][[top_results$solver[i]]]$xbest
  p <- get_temp_response_plot(result_list_combined[[start]][[top_results$solver[i]]]$xbest, temp_values = seq(-10,40, by = 0.1), log_A = TRUE)
  fname <- paste0('figures/run1_start_point_',i, 'temp_response.jpeg')
  ggsave(plot = p, filename = fname, height = 10, width = 15, units = 'cm', device = 'jpeg')
}


#write the model performance
for(i in seq_along(top_results$f)){
  print(i)
  start <- as.numeric(top_results$start_point[i])
  par <- result_list_combined[[start]][[top_results$solver[i]]]$xbest
  #training
  pred_days <- evaluation_function_meigo(par, modelfn = custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, SeasonList = SeasonList, return_bloom_days = TRUE)
  print(RMSEP(pred_days, bloomJDays))
  #evaluation
  pred_eval_days <- evaluation_function_meigo(par, modelfn = custom_PhenoFlex_GDHwrapper, bloomJDays = eval_bloomJDays, SeasonList = eval_SeasonList, return_bloom_days = TRUE)
  print(RMSEP(pred_eval_days, eval_bloomJDays))
  print('------')
}

saveRDS(start_point_r2, 'start_points_r2.rds')

result_list_r2 <- list()

#iterate other the start points
for(i in 1:5){
  #for(i in 1:length(start_points)){
  start_parameters <- start_point_r2[[i]]
  
  #=========================
  #PROBLEM SPECIFICATIONS
  #=========================
  problem<-list(f="evaluation_function_meigo",
                x_0 = start_parameters,
                x_L = x_L,
                x_U = x_U,
                c_L = c_L, 
                c_U = c_U)
  
  result_list_r2[[i]] <- list()
  
  #iterate over the solver
  for(solver in local_solver){
    
    #evaluation_function_meigo(start_parameters, modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, SeasonList = SeasonList)$g
    
    #specify options for the solver
    opts<-list(#maxeval = 1000,
      maxtime = 60 * 2, 
      local_solver = solver, 
      log_var = 7:8,
      local_bestx = 1,
      inter_save = 0)
    
    Results_meigo<-MEIGO(problem,
                         opts,
                         algorithm="ESS", 
                         modelfn = custom_PhenoFlex_GDHwrapper,
                         bloomJDays = bloomJDays,
                         SeasonList = SeasonList)
    
    result_list_r2[[i]][[solver]] <- Results_meigo
    
    # pred_days <- evaluation_function_meigo(x =  Results_meigo$xbest,
    #                           modelfn = custom_PhenoFlex_GDHwrapper,
    #                           bloomJDays = bloomJDays,
    #                           SeasonList = SeasonList, 
    #                           return_bloom_days = TRUE)
    # 
    # chillR::RPIQ(predicted = pred_days, observed = bloomJDays)
    # chillR::RMSEP(predicted = pred_days, observed = bloomJDays)
    
    
  }
  
  
  
}

result_list_r2_6to10 <- readRDS("result_run2_6to10.rds")
result_list_r2_combined <- c(result_list_r2, result_list_r2_6to10[6:10])

fbest_r2 <- purrr::map(result_list_r2_combined, function(x){
  purrr::map(x, function(y){
    data.frame(
      f = c(y$f, y$fbest),
      neval = c(y$neval, y$numeval)
    )
  }) %>% 
    bind_rows(.id = 'solver')
  
}) %>% 
  bind_rows(.id = 'start_point')


ggplot(fbest_r2, aes(x = neval/1000, y = f, col = solver)) + 
  geom_point() +
  geom_line() +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  #geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~as.factor(start_point)) +
  coord_cartesian(ylim = c(0, 2000)) +
  theme_bw()

ggsave('moderl_perofrmance_r2.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')


#combine fbest2 and fbest, also

neval_r2 <- fbest %>% 
  group_by(solver, start_point) %>% 
  summarise(eval_r1 = max(neval))


start_r2 <- fbest %>% 
  group_by(solver, start_point) %>% 
  summarise(f = min(f),
            neval = max(neval)) %>% 
  mutate(solver = paste0(solver, '_r2'),
         run = 'run_2')

#add another column to top_results
top_results$start_point_r2 <- 1:10
top_results$neval_r1 <- top_results$neval

#safe new start point in seperate column
fbest_r2$start_point_r2 <- fbest_r2$start_point

#change start_point to old start point of r1


fbest$run <- 'run_1'
fbest$start_point_r2 <- NA

#start_point of r2 does not correspond to start point of r1
#need to have two labels?
#start_point_r1 and start_point_r2


test <- fbest_r2 %>% 
  mutate(run = 'run_2',
         start_point = recode(start_point,
                              `1` = 3, `2` = 8, `3` = 1, `4` = 8, `5` = 4,
                              `6` = 6, `7` = 1, `8` = 2, `9` = 1, `10` = 10)) %>% 
  merge(top_results[,c('start_point_r2', 'neval_r1')], by = 'start_point_r2') %>% 
#  merge(neval_r2, by = c('start_point', 'solver')) %>% 
   # mutate(neval = neval + eval_r1,
   #      solver = paste0(solver, '_r2')) %>%
  mutate(neval = neval + neval_r1) %>% 
  select(start_point, solver, f, neval, run, start_point_r2) %>% 
  rbind(fbest) %>% 
  mutate(id = paste0(solver, start_point, start_point_r2))


#  rbind(start_r2) %>% 
#  ggplot(aes(x = neval/1000, y = f, col = solver)) + 
ggplot(test, aes(x = neval/1000, y = f)) + 
  geom_point(aes(col = run)) +
  geom_line(aes(col = run, group = id)) +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  #geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~as.factor(start_point), scales = 'free_x') +
  #geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  coord_cartesian(ylim = c(0, 2000)) +
  theme_bw()

ggsave('moderl_perofrmance_combined.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')

#combine it with the other data.frame



top_r2 <- test %>% 
  filter(run == 'run_2') %>% 
  group_by(solver, start_point_r2, start_point) %>% 
  summarise(f = min(f),
            neval = max(neval)) %>%
  ungroup() %>% 
  distinct(f, .keep_all = TRUE) %>% 
  slice_min(order_by = f, n = 10)





ggplot(test, aes(x = neval/1000, y = f)) + 
  geom_point(aes(col = run)) +
  geom_line(aes(col = run, group = id)) +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  #geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~as.factor(start_point), scales = 'free_x') +
  geom_point(data = top_r2, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  coord_cartesian(ylim = c(0, 2000)) +
  theme_bw()

ggsave('moderl_perofrmance_combined_highlight.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')



for(i in seq_along(top_r2$f)){
  start <- as.numeric(top_r2$start_point_r2[i])
  solver <- top_r2$solver[i]
  par <- result_list_r2_combined[[start]][[solver]]$xbest
  #p <- get_temp_response_plot(result_list_combined[[start]][[top_results$solver[i]]]$xbest, temp_values = seq(-10,40, by = 0.1), log_A = TRUE)
  p <- get_temp_response_plot(result_list_combined[[start]][[top_results$solver[i]]]$xbest, 
                         temp_values = seq(-10,40, by = 0.1), log_A = TRUE, hourtemps = hourtemps[hourtemps$Year %in% years,])
  fname <- paste0('figures/run3_start_point_',i, 'temp_response.jpeg')
  fname <- paste0('figures/run2_start_point_',i, 'temp_response.jpeg')
  ggsave(plot = p, filename = fname, height = 10, width = 15, units = 'cm', device = 'jpeg')

  
  print(i)
  #training
  pred_days <- evaluation_function_meigo(par, modelfn = custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, SeasonList = SeasonList, return_bloom_days = TRUE)
  print(RMSEP(pred_days, bloomJDays))
  #evaluation
  pred_eval_days <- evaluation_function_meigo(par, modelfn = custom_PhenoFlex_GDHwrapper, bloomJDays = eval_bloomJDays, SeasonList = eval_SeasonList, return_bloom_days = TRUE)
  print(RMSEP(pred_eval_days, eval_bloomJDays))
  print('------')
}

#--> looks like DHC has best rank and simular sd
#SA and SOLNP also look decent



par_hajar <- c(4.030065e+01, 1.761458e+02, 1.009330e-01, 2.753515e+01, 3.372800e+03, 9.900298e+03, 6.155735e+03,
               5.939889e+13, 2.104434e+00, 3.520949e+01, 6.656520e+00, 1.167121e+01)
get_temp_response_plot(par_hajar, temp_values = seq(-10,40, by = 0.1), log_A = FALSE)

fname <- paste0('figures/hajar_fit_blanquina.jpeg')
ggsave(plot = p, filename = fname, height = 10, width = 15, units = 'cm', device = 'jpeg')







#third run of fitting

#write the model performance
start_point_r3 <- list()

for(i in seq_along(top_r2$f)){
  start <- as.numeric(top_r2$start_point_r2[i])
  solver <- top_r2$solver[i]
  start_point_r3[[i]] <- result_list_r2_combined[[start]][[solver]]$xbest
}

saveRDS(start_point_r3, 'start_points_r3.rds')

result_list_r3 <- list()

#iterate other the start points
for(i in 1:5){
  #for(i in 1:length(start_points)){
  start_parameters <- start_point_r3[[i]]
  
  #=========================
  #PROBLEM SPECIFICATIONS
  #=========================
  problem<-list(f="evaluation_function_meigo",
                x_0 = start_parameters,
                x_L = x_L,
                x_U = x_U,
                c_L = c_L, 
                c_U = c_U)
  
  result_list_r3[[i]] <- list()
  
  #iterate over the solver
  for(solver in local_solver){
    
    #evaluation_function_meigo(start_parameters, modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, SeasonList = SeasonList)$g
    
    #specify options for the solver
    opts<-list(#maxeval = 1000,
      maxtime = 60 * 2, 
      local_solver = solver, 
      log_var = 7:8,
      local_bestx = 1,
      inter_save = 0)
    
    Results_meigo<-MEIGO(problem,
                         opts,
                         algorithm="ESS", 
                         modelfn = custom_PhenoFlex_GDHwrapper,
                         bloomJDays = bloomJDays,
                         SeasonList = SeasonList)
    
    result_list_r3[[i]][[solver]] <- Results_meigo
    
  }

}

result_list_r3_6to10 <- readRDS("result_run3_6to10.rds")
result_list_r3_combined <- c(result_list_r3, result_list_r3_6to10[6:10])


fbest_r3 <- purrr::map(result_list_r3_combined, function(x){
  purrr::map(x, function(y){
    data.frame(
      f = c(y$f, y$fbest),
      neval = c(y$neval, y$numeval)
    )
  }) %>% 
    bind_rows(.id = 'solver')
  
}) %>% 
  bind_rows(.id = 'start_point')

#save start point to other column
fbest_r3$start_point <- factor(fbest_r3$start_point, levels = 1:10)
fbest_r3$start_point_r3 <- fbest_r3$start_point
fbest_r3$start_point_r2 <- fbest_r3$start_point

levels(fbest_r3$start_point) <- top_r2$start_point
levels(fbest_r3$start_point_r2) <- top_r2$start_point_r2

top_r2$start_point_r3 <- 1:10
test$start_point_r3 <- NA

fbest_combined <- fbest_r3 %>% 
  merge(top_r2[,c('start_point', 'start_point_r2', 'start_point_r3', 'neval')], 
      by = c('start_point', 'start_point_r2', 'start_point_r3')) %>% 
  mutate(neval = neval.x + neval.y,
         id = paste0(solver, start_point, start_point_r2, start_point_r3),
         run = 'run_3') %>% 
  select(start_point, start_point_r2, start_point_r3, solver, f, neval, id, run) %>% 
  rbind(test) %>% 
  mutate(start_point = factor(start_point, levels = 1:10),
         start_point_r2 = factor(start_point_r2, levels = 1:10),
         start_point_r3 = factor(start_point_r3, levels = 1:10))

ggplot(fbest_r3, aes(x = neval/1000, y = f, col = solver)) + 
  geom_point() +
  geom_line() +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  #geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~start_point_r3) +
  coord_cartesian(ylim = c(0, 2000)) +
  theme_bw()

ggsave('moderl_perofrmance_r3.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')


#combine fbest2 and fbest, also

#  rbind(start_r2) %>% 
#  ggplot(aes(x = neval/1000, y = f, col = solver)) + 
ggplot(fbest_combined, aes(x = neval/1000, y = f)) + 
  geom_point(aes(col = run)) +
  geom_line(aes(col = run, group = id)) +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  #geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~as.factor(start_point), scales = 'free_x') +
  #geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  coord_cartesian(ylim = c(0, 2000)) +
  theme_bw()

ggsave('moderl_perofrmance_combined.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')

#combine it with the other data.frame



top_r3 <- fbest_combined %>% 
  filter(run == 'run_3') %>% 
  group_by(solver, start_point, start_point_r2, start_point_r3) %>% 
  summarise(f = min(f),
            neval = max(neval)) %>%
  ungroup() %>% 
  distinct(f, .keep_all = TRUE) %>% 
  slice_min(order_by = f, n = 10)


#highlight
ggplot(fbest_combined, aes(x = neval/1000, y = f)) + 
  geom_point(aes(col = run)) +
  geom_line(aes(col = run, group = id)) +
  ylab('Residual Sum of Squares (RSS)')+
  xlab('Number of evaluations * 1000')+
  #geom_point(data = top_results, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  facet_grid(~as.factor(start_point), scales = 'free_x') +
  geom_point(data = top_r3, aes(x = neval/1000, y = f), col = 'red', shape = 1, size = 4, stroke = 2)+
  coord_cartesian(ylim = c(0, 2000)) +
  theme_bw()

ggsave('moderl_perofrmance_combined_highlight.jpeg', height = 10, width = 20, units = 'cm', device = 'jpeg')



for(i in seq_along(top_r3$f)){
  start <- as.numeric(top_r3$start_point_r3[i])
  solver <- top_r3$solver[i]
  par <- result_list_r3_combined[[start]][[solver]]$xbest
  p <- get_temp_response_plot(result_list_combined[[start]][[top_results$solver[i]]]$xbest, 
                              temp_values = seq(-10,40, by = 0.1), log_A = TRUE, hourtemps = hourtemps[hourtemps$Year %in% years,])
  fname <- paste0('figures/run3_start_point_',i, 'temp_response.jpeg')
  ggsave(plot = p, filename = fname, height = 10, width = 15, units = 'cm', device = 'jpeg')
  
  
  print(i)
  #training
  pred_days <- evaluation_function_meigo(par, modelfn = custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, SeasonList = SeasonList, return_bloom_days = TRUE)
  print(RMSEP(pred_days, bloomJDays))
  #evaluation
  pred_eval_days <- evaluation_function_meigo(par, modelfn = custom_PhenoFlex_GDHwrapper, bloomJDays = eval_bloomJDays, SeasonList = eval_SeasonList, return_bloom_days = TRUE)
  print(RMSEP(pred_eval_days, eval_bloomJDays))
  print('------')
}















#have a look at the temperature response curves
#and have a look at the rmse / rpiq for the parameters

purrr::map(results_list, function(x){
  purrr::map_dbl(x, function(y){
    if(length(y) == 0){
      return(NA)
    }
    pred <- evaluation_function_meigo(y$xbest, modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, 
                              SeasonList = SeasonList, return_bloom_days = TRUE)
    
    RMSEP(predicted = pred, observed = bloomJDays)
    
  })
}) 



purrr::map(results_list, function(x){
  purrr::map_dbl(x, function(y){
    if(length(y) == 0){
      return(NA)
    }
    pred <- evaluation_function_meigo(y$xbest, modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, 
                                      SeasonList = SeasonList, return_bloom_days = TRUE)
    
    RPIQ(predicted = pred, observed = bloomJDays)
    
  })
}) 

#see if the optimizer yield similar temperature response functions
source('code/99_helper_functions.R')

get_temp_response_plot(results_list[[7]][['DHC']]$xbest, temp_values = seq(-10,40, by = 0.1), log_A = TRUE)
get_temp_response_plot(results_list[[7]][['SA']]$xbest, temp_values = seq(-10,40, by = 0.1), log_A = TRUE)
#outcome is quite variable (and far from being close to acceptable)










#--> make a more thourough search with these 3 solvers
results_list_sa <- list()
#use the sofar fitted values as a start point

#iterate other the start points
for(i in 1:8){
  results_list_sa[[i]] <- list()
  
  #iterate over the solver
  for(solver in c('SA')){
    
    start_parameters <- results_list[[i]][[solver]]$xbest
    
    #=========================
    #PROBLEM SPECIFICATIONS
    #=========================
    problem<-list(f="evaluation_function_meigo",
                  x_0 = start_parameters,
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
      inter_save = 1)
    
    Results_meigo<-MEIGO(problem,
                         opts,
                         algorithm="ESS", 
                         modelfn = custom_PhenoFlex_GDHwrapper,
                         bloomJDays = bloomJDays,
                         SeasonList = SeasonList)
    
    results_list_sa[[i]][[solver]] <- Results_meigo
    
    # pred_days <- evaluation_function_meigo(x =  Results_meigo$xbest,
    #                           modelfn = custom_PhenoFlex_GDHwrapper,
    #                           bloomJDays = bloomJDays,
    #                           SeasonList = SeasonList, 
    #                           return_bloom_days = TRUE)
    # 
    # chillR::RPIQ(predicted = pred_days, observed = bloomJDays)
    # chillR::RMSEP(predicted = pred_days, observed = bloomJDays)
    
    
  }
  
  
  
}


fbest2 <- purrr::map(results_list_2, function(x){
  unlist(purrr::map(x, 'fbest'))
}) 

fbest2 <- do.call('rbind', fbest2[1:8])

fbest2_long <- reshape2::melt(fbest2)
colnames(fbest2_long) <- c('start_point', 'solver', 'f')

ggplot(fbest2_long, aes(x = solver, y = f)) +
  geom_point() +
  facet_grid(~start_point)

#get rank per run
fbest2_long %>% 
  group_by(start_point) %>% 
  mutate(rank = rank(f)) %>% 
  ungroup() %>% 
  group_by(solver) %>% 
  summarise(mean_rank = mean(rank),
            median_rank = median(rank),
            sd_rank = sd(rank))

#--> looks like DHC has best rank and simular sd
#SA and SOLNP also look decent




#have a look at the temperature response curves
#and have a look at the rmse / rpiq for the parameters

purrr::map(results_list_2, function(x){
  purrr::map_dbl(x, function(y){
    if(length(y) == 0){
      return(NA)
    }
    pred <- evaluation_function_meigo(y$xbest, modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, 
                                      SeasonList = SeasonList, return_bloom_days = TRUE)
    
    RMSEP(predicted = pred, observed = bloomJDays)
    
  })
}) 



purrr::map(results_list_2, function(x){
  purrr::map_dbl(x, function(y){
    if(length(y) == 0){
      return(NA)
    }
    pred <- evaluation_function_meigo(y$xbest, modelfn =  custom_PhenoFlex_GDHwrapper, bloomJDays = bloomJDays, 
                                      SeasonList = SeasonList, return_bloom_days = TRUE)
    
    RPIQ(predicted = pred, observed = bloomJDays)
    
  })
}) 



get_temp_response_plot(results_list_2[[6]][['DHC']]$xbest, temp_values = seq(-10,40, by = 0.1), log_A = TRUE)
results_list_2[[6]][['DHC']]$xbest

get_temp_response_plot(results_list_2[[6]][['SA']]$xbest, temp_values = seq(-10,40, by = 0.1), log_A = TRUE)

##-> best starting point does not necissarily lead to best solution, but they are all pretty bad though



#maybe follow up using DHC only

results_list_2[[7]][['DHC']]$f





df_performance <- purrr::map(1:8, function(x){
  purrr::map(results_list[[x]], function(y){
    
    if(length(y) == 0){
      return(NULL)
    } else{
      return(    data.frame(start_point = x,
                            f = y$f,
                            neval = y$neval))
    }

  }) %>% 
    bind_rows(.id = 'solver')
  
  
}) %>% 
  bind_rows()


#now add the information of results_list_2
df_performance_2 <- purrr::map(1:8, function(x){
  purrr::map(results_list_2[[x]], function(y){
    return(    data.frame(start_point = x,
                          f = y$f,
                          neval = y$neval))
    
  }) %>% 
    bind_rows(.id = 'solver')
}) %>% 
  bind_rows()

#now add the total number of runs to the df_performance_2

n_runs <- df_performance %>% 
  group_by(solver, start_point) %>% 
  summarise(max = max(neval))

df_performance <- merge.data.frame(df_performance_2, n_runs, by = c('solver', 'start_point')) %>% 
  mutate(neval = neval + max) %>% 
  select(solver, start_point, f, neval) %>% 
  rbind(df_performance)


ggplot(df_performance, aes(x = neval, y = f, col = as.factor(start_point))) + 
  geom_point() +
  geom_line() +
  facet_grid(~solver) +
  coord_cartesian(ylim = c(0, 5000))




#maybe give the optimizers more time?


#--> make a more thourough search with these 3 solvers
results_list_3 <- list()
#use the sofar fitted values as a start point

#iterate other the start points
for(i in 1:8){
  results_list_3[[i]] <- list()
  
  #iterate over the solver
  for(solver in c('DHC', 'SA', 'SOLNP')){
    
    start_parameters <- results_list_2[[i]][[solver]]$xbest
    
    #=========================
    #PROBLEM SPECIFICATIONS
    #=========================
    problem<-list(f="evaluation_function_meigo",
                  x_0 = start_parameters,
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
      inter_save = 1)
    
    Results_meigo<-MEIGO(problem,
                         opts,
                         algorithm="ESS", 
                         modelfn = custom_PhenoFlex_GDHwrapper,
                         bloomJDays = bloomJDays,
                         SeasonList = SeasonList)
    
    results_list_3[[i]][[solver]] <- Results_meigo
    
    # pred_days <- evaluation_function_meigo(x =  Results_meigo$xbest,
    #                           modelfn = custom_PhenoFlex_GDHwrapper,
    #                           bloomJDays = bloomJDays,
    #                           SeasonList = SeasonList, 
    #                           return_bloom_days = TRUE)
    # 
    # chillR::RPIQ(predicted = pred_days, observed = bloomJDays)
    # chillR::RMSEP(predicted = pred_days, observed = bloomJDays)
    
    
  }
  
  
  
}








#--> make a more thourough search with these 3 solvers
results_list_3 <- list()
#use the sofar fitted values as a start point

#iterate other the start points
for(i in 1:8){
  print(i)
  results_list_3[[i]] <- list()
  
  #iterate over the solver
  for(solver in c('DHC')){
    
    # start_parameters <- results_list_2[[i]][[solver]]$xbest
    start_parameters <- Results_meigo$xbest
    
    results_list_2[[i]][[solver]]$fbest
    
    #contract x_L and x_U
    #if the start point is then outside the point, then only contract the other (upper or lower)
    #if both then exclude the point, then don't contract the bounds at all
    contract_percent <- 0.1
    #zeros cannot be contracted with this method, maybe add 1% of upper limit to the lower limit?
    
    #nimporovement at n_contraction: 3, 
    
    #n_contraction=6
    
    x_L_old <- x_L
    x_U_old <- x_U
    
    x_L_new <- x_L_old * (1 + contract_percent)
    x_L_new <- ifelse(x_L_new == 0, yes = x_L_new + (contract_percent * x_U_old), no = x_L_new)
    
    x_U_new <- x_U_old * (1 - contract_percent)
    
    #prevent the start point to be outside the bounds
    x_L_new <- ifelse(x_L_new > start_parameters, yes = x_L_old, no = x_L_new) 
    x_U_new <- ifelse(x_U_new < start_parameters, yes = x_U_old, no = x_U_new) 

    
    #=========================
    #PROBLEM SPECIFICATIONS
    #=========================
    problem<-list(f="evaluation_function_meigo",
                  x_0 = start_parameters,
                  x_L = x_L_new,
                  x_U = x_U_new,
                  c_L = c_L, 
                  c_U = c_U)
    
    #specify options for the solver
    opts<-list(#maxeval = 1000,
      maxtime = 60 * 5, 
      ndiverse = 500,
      local_solver = solver, 
      log_var = 7:8,
      local_bestx = 1,
      inter_save = 0,
      local_n2 = 10,
      local_balance = 0.1
      )
    #increase the frequency the local solver is called
    
    Results_meigo<-MEIGO(problem,
                         opts,
                         algorithm="ESS", 
                         modelfn = custom_PhenoFlex_GDHwrapper,
                         bloomJDays = bloomJDays,
                         SeasonList = SeasonList)
    
    results_list_3[[i]] <- Results_meigo
    
    # pred_days <- evaluation_function_meigo(x =  Results_meigo$xbest,
    #                           modelfn = custom_PhenoFlex_GDHwrapper,
    #                           bloomJDays = bloomJDays,
    #                           SeasonList = SeasonList, 
    #                           return_bloom_days = TRUE)
    # 
    # chillR::RPIQ(predicted = pred_days, observed = bloomJDays)
    # chillR::RMSEP(predicted = pred_days, observed = bloomJDays)
    
    
  }
  
  
  
}

results_list_3


pred_days <- evaluation_function_meigo(x =  start_parameters,
                                       modelfn = custom_PhenoFlex_GDHwrapper,
                                       bloomJDays = bloomJDays,
                                       SeasonList = SeasonList, 
                                       return_bloom_days = TRUE)

chillR::RPIQ(predicted = pred_days, observed = bloomJDays)
#baseline 0.32
chillR::RMSEP(predicted = pred_days, observed = bloomJDays)












length(results_list[[1]][[8]])
#plot all of the different starting points together with their model performance?
results_list



pred_days <- evaluation_function_meigo(x =  results_list[[6]]$xbest,
                                       modelfn = custom_PhenoFlex_GDHwrapper,
                                       bloomJDays = bloomJDays,
                                       SeasonList = SeasonList, 
                                       return_bloom_days = TRUE)

chillR::RPIQ(predicted = pred_days, observed = bloomJDays)
#baseline 0.32
chillR::RMSEP(predicted = pred_days, observed = bloomJDays)
#baseline 29.03


plot(results_list[[1]]$f ~ results_list[[1]]$neval,ylim = c(0,20000))
lines(results_list[[1]]$f ~ results_list[[1]]$neval, lty = 1)
points(results_list[[2]]$f ~ results_list[[2]]$neval, pch = 2)
lines(results_list[[2]]$f ~ results_list[[2]]$neval, lty = 2)
points(results_list[[3]]$f ~ results_list[[3]]$neval, pch = 3)
lines(results_list[[3]]$f ~ results_list[[3]]$neval, lty = 3)
points(results_list[[4]]$f ~ results_list[[4]]$neval, pch = 4)
lines(results_list[[4]]$f ~ results_list[[4]]$neval, lty = 4)
points(results_list[[5]]$f ~ results_list[[5]]$neval, pch = 5)
lines(results_list[[5]]$f ~ results_list[[5]]$neval, lty = 5)
points(results_list[[6]]$f ~ results_list[[6]]$neval, pch = 6)
lines(results_list[[6]]$f ~ results_list[[6]]$neval, lty = 6)
points(results_list[[7]]$f ~ results_list[[7]]$neval, pch = 7)
lines(results_list[[7]]$f ~ results_list[[7]]$neval, lty = 7)

results_list[[5]]$xbest
results_list[[6]]$xbest




#run optimization with the new and the standard starting point








#container for the multistart output
# multistart_list <- list()
# 
# 
# local_solver <- c('NM', 'BFGS', 'CG', 'LBFGSB', 'SA', 'SOLNP', 'DHC', 'NL2SOL')
# for(solver in local_solver){
#   
#   #change solver in the options
#   opts<-list(maxeval=1000, maxtime = 120, local_solver = solver, local_n1=2, local_n2=3, local_bestx = 1, ndiverse = 100)
#   
#   #========================
#   #END OF PROBLEM SPECIFICATIONS
#   # ========================
#   
#   #run the multistart optimization
#   Results_multistart<-MEIGO(problem,
#                             opts,
#                             algorithm="MULTISTART", 
#                             modelfn = custom_PhenoFlex_GDHwrapper,
#                             bloomJDays = bloomJDays,
#                             SeasonList = SeasonList)
#   
#   #load multistart data
#   load('eSSR_multistart.RData')
#   
#   #index of best performing start points
#   best <- order(Results_multistart$f0)[1:10]
#   
#   #save to list
#   multistart_list[[solver]] <- list(x0 = Results_multistart$x0[best,],
#                                     f0 = Results_multistart$f0[best])
#   
# }
# 
# 
# f0 <- purrr::map(multistart_list, 'f0') %>% 
#   bind_cols() %>% 
#   reshape2::melt(variable.name = 'local_solver', value.name = 'f0')
# 
# ggplot(f0, aes(x = local_solver, y = f0)) +
#   geom_boxplot()
#   
# 
# #check the different parameters
# x0 <- purrr::map(multistart_list, 'x0') %>% 
#   do.call(what = 'rbind', args = .) %>% 
#   data.frame() 
# 
# colnames(x0) <-c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1', 'Tf', 'Tc', 'Tb', 'slope')
# 
# #Tb, Tu, Tc and zc from GDH model
# #s1 for transitioning between heat and chill
# #E0, E1, A0, A1, Tf, slope, yc dynamic model
# 
# param_df <- data.frame(parameter = c('Tb', 'Tu', 'Tc', 'zc', 's1', 'E0', 'E1', 'A0', 'A1', 'Tf', 'slope', 'yc'),
#                        model = c(rep('GDH', 4), 'PhenoFlex', rep('Dynamic model', 7)))
# 
# ranges_df <- data.frame(parameter = c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1', 'Tf', 'Tc', 'Tb', 'slope'),
#                         upper = x_U,
#                         lower = x_L)
# 
# multistart_df <- x0 %>% 
#   cbind(f0) %>% 
#   reshape2::melt(id.vars = c('local_solver', 'f0'), value.name = 'par_value', variable.name = 'parameter') %>% 
#   merge(param_df, by = 'parameter') %>% 
#   merge(ranges_df, by = 'parameter')
# 
# 
# 
# library(colorRamps)
# 
# multistart_df %>% 
#   ggplot(aes(x = local_solver, y = par_value)) +
#   geom_boxplot() +
#   geom_hline(aes(yintercept = upper), linetype = 'dashed') + 
#   geom_hline(aes(yintercept = lower), linetype = 'dashed') + 
#   geom_jitter(width = 0.2,aes(col = f0), )+
#   scale_color_gradientn(colours=alpha(matlab.like(15), alpha = .5),
#                        name="Model Performance")+
#   facet_wrap(model ~ parameter, scales = 'free_y')
# 
# #calculate the rmse and rpiq for each set of parameter and use that for plotting instead
# f0$rmse <- NA
# f0$rpiq <- NA
# for(i in 1:nrow(x0)){
#   pred_days <- evaluation_function_meigo(x = unlist(x0[i,]), 
#                             modelfn = custom_PhenoFlex_GDHwrapper,
#                             bloomJDays = bloomJDays,
#                             SeasonList = SeasonList, 
#                             return_bloom_days = TRUE)
#   
#   f0$rpiq[i] <- RPIQ(predicted = pred_days, observed = bloomJDays)
#   f0$rmse[i] <- RMSEP(predicted = pred_days, observed = bloomJDays)
#   
# }
# 
# f0_baseline <- evaluation_function_meigo(x = x_0, 
#                                          modelfn = custom_PhenoFlex_GDHwrapper,
#                                          bloomJDays = bloomJDays,
#                                          SeasonList = SeasonList)$F
# 
# x_0
# x0[which(f0$f0 <= f0_baseline),]
#multistart doesnt really help finding better starting points than the one we have
#-->only one point found which is a better start

