library(chillR)
library(ggplot2)


#read bloom data
pheno_data <- readxl::read_excel('data/Blanquina.xlsx')

#drop nas
pheno_data <- na.omit(pheno_data)

#have year as numeric
pheno_data$Year <- as.numeric(pheno_data$Year)

#read hourly temperature data
hourtemps <- read.csv('data/hourtemps.csv')



#JAE Set seed for random numbers to fix the validation years
set.seed(1235)


#split phenological data into training and calibration data
val_years <- sample(pheno_data$Year, 7)
val_data <- pheno_data[pheno_data$Year %in% val_years,]
cal_data <- pheno_data[!pheno_data$Year %in% val_years,]

#split temperature data into chill seasons
SeasonList <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data$Year)


######run 1#######

#JAE. Transform values of A0 and A1 into their logarithms and consider ranges
#from 1e3 to 1e15 for A0 and 1e13 to 1e19 for A1
#Ranges for E0 and E1 should be between 2000 and 6000 (E0) and 90000 and 12000 (E1)
#See Table 1 from https://doi.org/10.1016/S0022-5193(87)80237-0
#In my opinion, the upper value for the slope should be 5, not 50

#default parameters according to the vignette
#          yc,  zc,  s1, Tu,    E0,      E1,     A0,         A1,   Tf, Tc, Tb,  slope
par <-   c(40, 190, 0.5, 25, 3372.8,  9900.3, log(6319.5), log(5.939917e13),  4, 36,  4,  1.60)
upper <- c(80, 500, 1.0, 30, 6000.0, 12000.0, log(1e15),    log(1e19),       10, 40, 10, 5.00)
lower <- c(20, 100, 0.1, 0 , 2000.0, 9000.0, log(1e3),      log(1e13),        0,  0,  0,  0.05)


#JAE create a random initial parameter vector from lower and upper
par <-   runif(12)*(upper-lower)+lower

#JAE source auxliary functions that consider the logs of the A0 and A1 parameters
#and calculates the exponentials inside. Other changes commented inside the functions
source("code/03_PhenoFlex_GDHwrapper_JAE.R")
source("code/04_chifull_JAE.R")
source("code/05_predictBloomDays_JAE.R")
source("code/06_phenologyFitter_JAE.R")



#JAE try a set of 1000 initial solutions before applying the optimization algorithm
# I applied this filter below because I saw it in one of the intermediate functions
# but it seems to work without it, and causes some problems if we apply it.
# So, we skip it. The same filter is later applied during the optimization

# bloomJDays <- cal_data$pheno
# for (i in c(1:length(SeasonList))) {
#   minJDay <- SeasonList[[i]]$JDay[1]
#   maxJDay <- SeasonList[[i]]$JDay[length(SeasonList[[i]]$JDay)]
#   if (maxJDay > minJDay) {
#     stop(paste0("Season ", i, " is overlapping with the previous or following one. Aborting!"))
#   }
#   if (bloomJDays[i] > maxJDay && bloomJDays[i] < minJDay) {
#     stop(paste0("In season ", i, " the bloomJDay is outside the provided JDay vector. Aborting!"))
#   }
#   dx <- diff(SeasonList[[i]]$JDay)
#   mx <- min(dx)
#   kmx <- which(dx == mx)
#   SeasonList[[i]]$JDay[1:kmx] <- SeasonList[[i]]$JDay[1:kmx] + 
#     mx - 1
#   if (bloomJDays[i] > minJDay) 
#     bloomJDays[i] <- bloomJDays[i] + mx - 1
#   SeasonList[[i]]$JDayunwrapped <- SeasonList[[i]]$JDay
#   bloomJDaysunwrapped <- bloomJDays
# }



#JAE perform n_init initial simulations to find nice initial guesses for the parameters 
n_init <- 10000

#Initialize a matrix of initial solutions with zeros and add the first solution
#with our own initial guess
sol_x <- matrix(0, n_init, 12)
sol_x[1,]<- c(40, 190, 0.5, 25, 3372.8,  9900.3, log(6319.5), log(5.939917e13),  4, 36,  4,  1.60)
sol_f <- rep(Inf, n_init)

for (i in 1:n_init){
  sol_x[i,] <- runif(12)*(upper-lower)+lower
  sol_f[i] <- chifull_JAE(sol_x[i,], PhenoFlex_GDHwrapper_JAE, cal_data$pheno, SeasonList, na_penalty = 365)
}


#Retrieve the best solution among the n_init simulated as initial guess for the optimization
par<-sol_x[which.min(sol_f),]

#JAE Changed option VERBOSE to TRUE to check how the objective function evolves
#JAE now modelfn is PhenoFlex-GDHwrapper_JAE to undo the logarithmic transformation
Fit_res <- phenologyFitter_JAE(par.guess=par, 
                           modelfn = PhenoFlex_GDHwrapper_JAE,
                           bloomJDays = cal_data$pheno,
                           SeasonList = SeasonList,
                           lower = lower,
                           upper = upper,
                           seed = 001,
                           control = list(smooth=FALSE, verbose=TRUE, maxit=1000,
                                          nb.stop.improvement=250))


#check error 
RMSEP(Fit_res$pbloomJDays, cal_data$pheno)
RPIQ(Fit_res$pbloomJDays, cal_data$pheno)
mean(abs(cal_data$pheno - Fit_res$pbloomJDays))
mean(cal_data$pheno - Fit_res$pbloomJDays)



#read functions to easily create plots
source('code/99_helper_functions.R')

#values for plotting temperature response
temp_values = seq(-5, 40, 0.1)


#JAE Transform the final parameters log(A0) and log(A1)
#to A0 and A1 values
Fit_res$par[7]<-exp(Fit_res$par[7])
Fit_res$par[8]<-exp(Fit_res$par[8])



get_qqplot(cal_data, Fit_res$pbloomJDays)
get_hist_plot(cal_data$pheno - Fit_res$pbloomJDays)
get_temp_response_plot(Fit_res$par, temp_values)



# Not sure if this may be useful: 
# In the next fitting run, we re-define the parameters based on the results obtained here (i.e. Fit_res$par). Lower and upper
# bounds are then adjusted to keep the estimated par within the boundaries... We finally run the fitting function again,
# searching for lower RMSE and "plausible" response curves.


#run on validation data








