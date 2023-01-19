#this function takes the weather data (x) and the model parameters (par) and returns the expected day of bloom 
custom_PhenoFlex_GDHwrapper <- function (x, par){
  #x is one of the elements in season List
  #par are the parameters of the model
  


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

custom_PhenoFlex_GDHwrapper_log <- function (x, par){
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


custom_PhenoFlex_GDHwrapper_no_control <- function (x, par){
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

evaluation_function_meigo_vns <- function(x, 
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
    return(F)
  } else{
    return(pred_bloom)
  }
  
}


return_predicted_days <- function(par, 
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
  par 
  
  #calculate the predicted flower dates
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, par = par))
  
  #if the model returns no bloom day, then give penalty term instead
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
  
  return(pred_bloom)
  
}


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

#this function is used to calculate E0, E1, A0 and A1 based on certain values of theta*, theta_c, tau(theta*) and pi_c
solve_nle <- function(x,params){
  
  #Initialize the vector of residuals
  y<-numeric(2)
  
  
  #Introduce the values of the parameters actually optimized
  par1<-params[1]    #theta*
  par2<-params[2]    #theta_c
  par3<-params[3]    #Tau(thetha*)
  par4<-params[4]    #pi_c
  
  T1<-297
  T2<-279
  
  eta <- 1/3
  q=1/par1-1/par2
  
  
  #Equation 35 Fishman et al 87. In the paper it seems that the log is multiplying, but it is not.
  #y[1] <- (x[1]-x[2])/(exp((x[2]-x[1])*q)-1)/log(1-exp((x[1]-x[2])*q))-x[2]
  y[1] <- log((x[1]-x[2])/(exp((x[2]-x[1])*q)-1)/log(1-exp((x[1]-x[2])*q)))-log(x[2])
  
  
  #Calculate these terms needed for the next equation
  #From Eq 36, we get A1 depending only in E0 and E1. Here the log is indeed multiplying!!.
  A1 <- -exp(x[2]/par1)/par3*log(1-exp((x[1]-x[2])*q))
  k1T1<-A1*exp(-x[2]/T1)
  k1T2<-A1*exp(-x[2]/T2)
  
  # Equation 38. It only depends on E0 and E1. Also in A1 but we had its expresion in terms of E0 and E1
  lhs<-(exp((x[2]-x[1])/par2)-exp((x[2]-x[1])/T1))/(exp((x[2]-x[1])/T2)-exp((x[2]-x[1])/T1))
  rhs<-(1-exp(-k1T2*(1-eta)*par4))/(1-exp(-(k1T1*eta+k1T2*(1-eta))*par4))
  
  
  y[2] <- log(lhs)-log(rhs)   #Taking logs the problems is much easily solved
  
  return(y)
}


evaluation_function_meigo_nonlinear <- function(x, 
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
  
  params[1] <- x[5]   #theta*
  params[2] <- x[6]    #theta_c
  params[3] <- x[7]    #Tau(thetha*)
  params[4] <- x[8]     #pi_c
  
  
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
  par[5:8] <- c(E0, E1, A0, A1)
  
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
  g[1] <- par[4] - par[11]
  #Tx >= Tb
  g[2] <- par[10] - par[11]
  #Tc >= Tu
  g[3] <- par[10] - par[4]
  
  
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[4] <- exp((10 * par[5]) / (297 * 279))
  
  g[5] <- exp((10 * par[6]) / (297 * 279))
  
  
  if(return_bloom_days == FALSE){
    #output
    return(list(F=F, g=g))
  } else{
    return(pred_bloom)
  }
}