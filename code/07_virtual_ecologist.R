#prepare virtual ecologist

files <- list.files('data/hajar_parameters/', full.names = TRUE)
par_list <- lapply(files,read.csv) 

files <- list.files('data/training/', full.names = TRUE)
pheno_train <- lapply(files, read.csv)


#read temperature data
hourtemps <- read.csv('data/hourtemps.csv')

temp_list <- list()

#generate list of temperatures for each cultivar
for(i in 1:length(pheno_train)){
  temp_list[[i]] <- chillR::genSeasonList(hourtemps, years = pheno_train[[i]]$Year)
}


source('code/04_functions_meigo.R')

return_bloomdays <- function(x, 
         modelfn,
         SeasonList,
         na_penalty = 365){
  
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
  
  return(pred_bloom)
  
}

custom_PhenoFlex_GDHwrapper <- function (x, par){
  #x is one of the elements in season List
  #par are the parameters of the model
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

cultivars <- c('Blanquina', 'Clara', 'Collaos', 'Coloradona',
               'DelaRiega', 'Perezosa', 'Perico', 'Raxao',
               'Teorica', 'Verdialona', 'Xuanina')

pred_df <- data.frame()
for(i in 1:length(par_list)){
  
  par <- par_list[[i]]$Last_run

  pred_bloom <- return_bloomdays(x = par, 
                                          modelfn = custom_PhenoFlex_GDHwrapper,
                                          SeasonList = temp_list[[i]])
  
  pred_df <- rbind(pred_df,
                    data.frame(cultivar = cultivars[i],
                               year = pheno_train[[i]]$Year,
                               pheno = pred_bloom))
}


synthetic_pheno_data <- split(pred_df, f = pred_df$cultivar)




#      ('yc', 'zc', 's1', 'Tu', 'theta*', 'theta_c', 'Tau(thetha*)', 'pie_c',   'Tf', 'Tc', 'Tb',  'slope')
x_0 <- c(40,   190,   0.5,  25,   278.1,    285.1,   47.7,           30.6,        4,   36,     4,    1.60)
x_U <- c(80,   500,   1.0,  30,   281,      287,       48,             32,       10,   40,    10,    5.00)
x_L <- c(20,   100,   0.1,   0,   279,      284,       16,             24,        0,    0,     0,    0.05)

#limits for the inequality constraints
#         #gdh parameters   #q10 for E0 and E1
c_L <- c(  0,   0,   0,     1.5, 1.5)
c_U <- c(Inf, Inf, Inf,     3.5, 3.5)

source('code/05_functions_meigo_nonlin_equations.R')

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
      maxtime = 60 * 2, 
      local_solver = solver, 
      local_bestx = 1,
      inter_save = 0)
    
    res_list[[i]][[solver]]<-MEIGO(problem,
                                   opts,
                                   algorithm="ESS", 
                                   modelfn = custom_PhenoFlex_GDHwrapper,
                                   bloomJDays = synthetic_pheno_data[[i]]$pheno,
                                   SeasonList = temp_list[[i]])
    
  } 
}

purrr::map(res_list, function(x) x[['DHC']]$fbest)
purrr::map(res_list, function(x) x[['DHC']]$xbest)

convert_parameters <- function(par){
  params<-numeric(4)
  
  params[1] <- par[5]   #theta*
  params[2] <- par[6]    #theta_c
  params[3] <- par[7]    #Tau(thetha*)
  params[4] <- par[8]     #pi_c
  
  
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
  par[5:8] <- c(E0, E1, A0, A1)
  
  return(par)
}

ve_df <- data.frame()

for(i in 1:length(res_list)){
  
  ve_df <- rbind(ve_df ,
                 data.frame(cultivar = cultivars[i],
             hajar = par_list[[i]]$Last_run,
             meigo = convert_parameters(res_list[[i]]$DHC$xbest),
             parameter = c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1',   'Tf', 'Tc', 'Tb',  'slope'))
  )
  
}


ve_df %>% 
  ggplot(aes(x = cultivar)) + 
  geom_point(aes(y = hajar, col = 'hajar')) + 
  geom_point(aes(y = meigo, col = 'meigo')) + 
  facet_wrap(~parameter, scales = 'free_y')

#maybe there is a mistake in the conversion of the parameters and 