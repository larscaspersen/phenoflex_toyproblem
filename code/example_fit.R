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



#split phenological data into training and calibration data
val_years <- sample(pheno_data$Year, 7)
val_data <- pheno_data[pheno_data$Year %in% val_years,]
cal_data <- pheno_data[!pheno_data$Year %in% val_years,]

#split temperature data into chill seasons
SeasonList <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data$Year)


######run 1#######

#default parameters according to the vignette
#          yc,  zc,  s1, Tu,    E0,      E1,     A0,         A1,   Tf, Tc, Tb,  slope
par <-   c(40, 190, 0.5, 25, 3372.8,  9900.3, 6319.5, 5.939917e13,  4, 36,  4,  1.60)
upper <- c(80, 500, 1.0, 30, 4000.0, 10000.0, 7000.0,       6.e13, 10, 40, 10, 50.00)
lower <- c(20, 100, 0.1, 0 , 3000.0,  9000.0, 6000.0,       5.e13,  0,  0,  0,  0.05)


Fit_res <- phenologyFitter(par.guess=par, 
                           modelfn = PhenoFlex_GDHwrapper,
                           bloomJDays = cal_data$pheno,
                           SeasonList = SeasonList,
                           lower = lower,
                           upper = upper,
                           seed = 001,
                           control = list(smooth=FALSE, verbose=FALSE, maxit=1000,
                                          nb.stop.improvement=250))


#check error 
RMSEP(Fit_res$pbloomJDays, cal_data$pheno)
RPIQ(Fit_res$pbloomJDays, cal_data$pheno)
mean(abs(cal_data$pheno - Fit_res$pbloomJDays))
mean(cal_data$pheno - Fit_res$pbloomJDays)



#read functions to easily create plots
source('code/helper_functions.R')

#values for plotting temperature response
temp_values = seq(-5, 40, 0.1)


get_qqplot(cal_data, Fit_res$pbloomJDays)
get_hist_plot(cal_data$pheno - Fit_res$pbloomJDays)
get_temp_response_plot(Fit_res$par, temp_values)



# Not sure if this may be useful: 
# In the next fitting run, we re-define the parameters based on the results obtained here (i.e. Fit_res$par). Lower and upper
# bounds are then adjusted to keep the estimated par within the boundaries... We finally run the fitting function again,
# searching for lower RMSE and "plausible" response curves.






