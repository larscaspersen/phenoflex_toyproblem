#modification made by Jose


#this function is needed to do the actual model calculations
#calculate accumulated heat and chill
#in the end give a day for bloom date

#modification JAE: make it run on logarithmic A0 and A1
PhenoFlex_GDHwrapper_JAE <- function (x, par) 
{
  
  if (par[4] <= par[11]) 
    return(NA)
  if (par[10] <= par[4]) 
    return(NA)
  
  
  #JAE transform A0 and A1 into their exponential since the optimization algorithm
  #works with their logarithms
  par[7]<-exp(par[7])
  par[8]<-exp(par[8])
  
  
  
  
  bloomindex <- PhenoFlex(temp = x$Temp, times = seq_along(x$Temp), 
                          yc = par[1], zc = par[2], s1 = par[3], Tu = par[4], 
                          E0 = par[5], E1 = par[6], A0 = par[7], A1 = par[8], 
                          Tf = par[9], Tc = par[10], Tb = par[11], slope = par[12], 
                          Imodel = 0L, basic_output = TRUE)$bloomindex
  if (bloomindex == 0) 
    return(NA)
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  n <- length(JDaylist)
  if (n == 1) 
    return(JDay)
  return(JDay + which(JDaylist == bloomindex)/n - 1/(n/ceiling(n/2)))
}


#this function returns the model performance with a given set of parameters

#modification JAE: added another constraint for E0 and E1 using the Q10 concept
chifull_JAE <- function (par, modelfn, bloomJDays, SeasonList, na_penalty = 365,
                         ...)
{
  sres <- predictBloomDays_JAE(par = par, SeasonList = SeasonList,
                               modelfn = modelfn, ...)
  s <- (sres - bloomJDays)
  nai <- which(is.na(sres))
  s[nai] <- na_penalty
  
  
  #JAE Calculation of the Q10 coefficients
  #Discard the solution (e.g. F(x)=1e20 if they are not in the range
  
  e0<-par[5]
  e1<-par[6]
  
  Q10_E0<-(exp(-e0/297)/exp(-e0/279))^(10/18)
  Q10_E1<-(exp(-e1/297)/exp(-e1/279))^(10/18)
  
  if (Q10_E0<1.5 | Q10_E0>3.5 | Q10_E1<1.5 | Q10_E1>3.5){
    return(1e20)
  }
  else{
    
    
    return(sum(s[!is.na(s)]^2))
  }
} 

#internal function, used to calculate the error
#this function calls the gdh wrapper and calculate the bloom days

#no changes by JAE
predictBloomDays_JAE <- function (par, SeasonList, modelfn, ...) 
{
  paravail <- requireNamespace("parallel")
  sres <- if (FALSE) {
    mc.cores <- as.numeric(Sys.getenv("OMP_NUM_THREADS"))
    if (is.na(mc.cores)) 
      mc.cores <- getOption("mc.cores", default = 1L)
    parallel::mclapply(X = SeasonList, FUN = modelfn, par = par, 
                       ..., mc.cores = mc.cores)
  }
  else {
    lapply(X = SeasonList, FUN = modelfn, par = par, ...)
  }
  return(invisible(simplify2array(x = sres)))
}


#this is the wrapper function for the whole fitting process
#jose only adjusted the called functions 

phenologyFitter_JAE <- function (par.guess = NULL, modelfn = PhenoFlex_GDHwrapper_JAE, 
                                 bloomJDays, SeasonList, control = list(smooth = FALSE, verbose = TRUE, 
                                                                        maxit = 1000, nb.stop.improvement = 250), lower, upper, 
                                 seed = 1235433, ...) 
{
  control$seed <- seed
  stopifnot(is.list(SeasonList))
  stopifnot(length(SeasonList) == length(bloomJDays))
  res <- phenologyFit()
  res$par.guess <- par.guess
  res$modelfn <- modelfn
  res$bloomJDays <- bloomJDays
  res$SeasonList <- SeasonList
  res$lower <- lower
  res$upper <- upper
  res$control <- control
  for (i in c(1:length(SeasonList))) {
    minJDay <- SeasonList[[i]]$JDay[1]
    maxJDay <- SeasonList[[i]]$JDay[length(SeasonList[[i]]$JDay)]
    if (maxJDay > minJDay) {
      stop(paste0("Season ", i, " is overlapping with the previous or following one. Aborting!"))
    }
    if (bloomJDays[i] > maxJDay && bloomJDays[i] < minJDay) {
      stop(paste0("In season ", i, " the bloomJDay is outside the provided JDay vector. Aborting!"))
    }
    dx <- diff(SeasonList[[i]]$JDay)
    mx <- min(dx)
    kmx <- which(dx == mx)
    SeasonList[[i]]$JDay[1:kmx] <- SeasonList[[i]]$JDay[1:kmx] + 
      mx - 1
    if (bloomJDays[i] > minJDay) 
      bloomJDays[i] <- bloomJDays[i] + mx - 1
    res$SeasonList[[i]]$JDayunwrapped <- SeasonList[[i]]$JDay
    res$bloomJDaysunwrapped <- bloomJDays
  }
  res$model_fit <- GenSA::GenSA(par = par.guess, fn = chifull_JAE, 
                                bloomJDays = bloomJDays, SeasonList = SeasonList, modelfn = modelfn, 
                                control = control, lower = lower, upper = upper)
  res$par <- res$model_fit$par
  res$pbloomJDays <- predictBloomDays_JAE(par = res$par, SeasonList = res$SeasonList, 
                                          modelfn = res$modelfn)
  return(res)
}