#meigo in action v2

library(chillR)
library(MEIGOR)
library(tidyverse)

#functions needed for the fitting
source('code/05_functions_meigo_nonlin_equations.R')

#training data
pheno_train <- read.csv('data/training/calibration_data_Blanquina.csv')
bloomJDays <- pheno_train$pheno

#temperature data
hourtemps <- read.csv('data/hourtemps.csv')
SeasonList <- chillR::genSeasonList(hourtemps, years = pheno_train$Year)

#      ('yc', 'zc', 's1', 'Tu', 'theta*', 'theta_c', 'Tau(thetha*)', 'pie_c',   'Tf', 'Tc', 'Tb',  'slope')
x_0 <- c(40,   190,   0.5,  25,   278.1,    285.1,   47.7,           30.6,        4,   36,     4,    1.60)
x_U <- c(80,   500,   1.0,  30,   281,      287,       48,             32,       10,   40,    10,    5.00)
x_L <- c(20,   100,   0.1,   0,   279,      284,       16,             24,        0,    0,     0,    0.05)

#limits for the inequality constraints
#         #gdh parameters   #q10 for E0 and E1
c_L <- c(  0,   0,   0,     1.5, 1.5)
c_U <- c(Inf, Inf, Inf,     3.5, 3.5)

#=========================
#PROBLEM SPECIFICATIONS
#=========================
problem<-list(f="evaluation_function_meigo",
              x_0= x_0,
              x_L = x_L,
              x_U = x_U,
              c_L = c_L, 
              c_U = c_U)

#specify options for the solver
opts<-list(#maxeval = 1000,
  maxtime = 60 * 2, 
  #ndiverse = 10000,
  local_solver = 'SOLNP', 
  #local_bestx = 1,
  inter_save = 0)


Results_meigo<-MEIGO(problem = problem,
                     opts = opts,
                     algorithm="ESS", 
                     modelfn = custom_PhenoFlex_GDHwrapper,
                     bloomJDays = bloomJDays,
                     SeasonList = SeasonList)

source('code/99_helper_functions.R')

get_temp_response_plot_v2(par = Results_meigo$xbest, temp_values = seq(-10,40, by = 0.1), par_type = 'new', log_A = FALSE, hourtemps = hourtemps)
#somehow the old temperature response plot does not work with the new parameters, even after converting them

evaluation_function_meigo(x = Results_meigo$xbest,
                          modelfn = custom_PhenoFlex_GDHwrapper,
                          bloomJDays = bloomJDays,
                          SeasonList = SeasonList, return_bloom_days = TRUE)


par <- Results_meigo$xbest

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
  error('Could not find corresponding values of E0, E1, A0 and A1')
  
  #You would add here a flag to let your optimization procedure know
  #That this solution should be ignored by lack of convergence
  
} else {
  
  par[5] <- output$x[1]
  par[6] <- output$x[2]
  
  #A1 and A0 can be calculated through Equations 36 and 37
  
  q=1/params[1]-1/params[2]
  
  par[8] <- -exp(E1/params[1])/params[3]*log(1-exp((E0-E1)*q))
  par[7] <- A1*exp((E0-E1)/params[2])
}



res <- PhenoFlex(temp = SeasonList[[1]]$Temp, 
          times = 1:length(SeasonList[[1]]$Temp),
          yc = par[1],  zc = par[2], 
          s1 = par[3], Tu = par[4],
          E0 = par[5], E1 = par[6],
          A0 = par[7], A1 = par[8],
          Tf = par[9],  Tc = par[10],
          Tb = par[11],  slope = par[12],
          Imodel = 0L, basic_output=FALSE)

yc = par[1]
zc = par[2]


DBreakDay <- res$bloomindex
seasontemps<-SeasonList[[1]]
seasontemps[,"x"]<-res$x
seasontemps[,"y"]<-res$y
seasontemps[,"z"]<-res$z
seasontemps$Date <- lubridate::parse_date_time(x = paste(SeasonList[[1]]$Year, seasontemps$JDay), orders = "yj")


CR_full<-seasontemps$Date[which(seasontemps$y>=yc)[1]]
Bloom<-seasontemps$Date[which(seasontemps$z>=zc)[1]]

chillplot<-ggplot(data=seasontemps[1:DBreakDay,],aes(x=Date,y=y)) +
  geom_line(col="blue",lwd=1.5) +
  theme_bw(base_size=20) +
  geom_hline(yintercept=yc,lty=2,col="blue",lwd=1.2) +
  geom_vline(xintercept=CR_full,lty=3,col="blue",lwd=1.2) +
  ylab("Chill accumulation (y)") +
  labs(title="Chilling") +
  annotate("text",label="Chill req. (yc)", 
           x=ISOdate(1986,10,01),
           y=yc*1.1, col="blue",size=5)

heatplot<-ggplot(data=seasontemps[1:DBreakDay,],aes(x=Date,y=z)) +
  geom_line(col="red",lwd=1.5) +
  theme_bw(base_size=20) +
  scale_y_continuous(position = "right") +
  geom_hline(yintercept=zc,lty=2,col="red",lwd=1.2) +
  geom_vline(xintercept=CR_full,lty=3,col="blue",lwd=1.2) +
  geom_vline(xintercept=Bloom,lty=3,col="red",lwd=1.2) +
  ylab("Heat accumulation (z)") +
  labs(title="Forcing")  +
  annotate("text",label="Heat req. (zc)", 
           x=ISOdate(1986,10,01),
           y=zc*0.95, col="red",size=5)


library(patchwork)
chillplot + heatplot

