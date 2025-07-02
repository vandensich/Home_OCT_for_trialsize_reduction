try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

library(ggh4x)
library(dMod)
library(tidyverse)
library(tidyr)
library(tibble)
library(dplyr)
library(boot)
library(truncnorm)

# set the Trial conditions like dosing times and patient numbers

runin_schedule <- c(seq(28,84,56))
SoC_schedule <- c(seq(max(runin_schedule)+56, max(runin_schedule)+3*56, 56))
n_patients <- 1000

# generate patient parameters based on the Parameter specific distributions

Patient_pars <- c()

set.seed(42)

dummy <- data.frame()
for (i in 1:n_patients){
  dummy <- data.frame(
    ID = i,
    CSTDEL = log10(exp(rtruncnorm(n=1, b=6.9, mean=4.6, sd=sqrt(1.02)))),
    KIN = log10(exp(4.01)),
    EC50_Eylea = log10(exp(rnorm(1, mean=-1.67, sd=sqrt(2.82)))),
    EC50_Lucentis = log10(exp(rnorm(1, mean=-3.82, sd=sqrt(2.82)))),
    EC50_SoC = log10(exp(rnorm(1, mean=-3.82, sd=sqrt(2.82)))),
    CSTMIN = log10(exp(rtruncnorm(n=1, a=4.60517, mean=5.42, sd=sqrt(3.34e-2)))),
    CST0 = log10(exp(rnorm(1, mean=5.51, sd=sqrt(3.83e-2)))),
    PK0 = log10(exp(rtruncnorm(n=1, b=-0.6931, mean=-1.19, sd=sqrt(4.54)))),
    HILL = log10(exp(0.79)),
    sigma_CST_obs_abs = 5.49
  )
  Patient_pars <- rbind(Patient_pars, dummy)
  dummy <- data.frame()
}
Patient_pars <- as.data.frame(Patient_pars)
Patient_pars <- as.data.frame(Patient_pars[,-1],row.names=as.character(Patient_pars[,1]))


# Define the subset sizes you want (20, 30, 40, ... 200)
subset_sizes <- seq(100, 1000, by=100)

# Initialize an empty data frame to store means
Patient_means <- data.frame()

# Loop over each subset size
for (n in subset_sizes) {
  # Compute the mean of each parameter for the first n patients
  means_vec <- colMeans(Patient_pars[1:n, ])
  
  # Convert the result to a one-row data frame
  temp_df <- as.data.frame(t(means_vec))
  
  # Add a column indicating the number of patients included
  temp_df$n_patients <- n
  
  # Append to the Patient_means data frame
  Patient_means <- rbind(Patient_means, temp_df)
}

# Optionally rearrange columns so that n_patients is the first column
Patient_means <- Patient_means[, c(ncol(Patient_means), 1:(ncol(Patient_means)-1))]
Patient_means <- as.data.frame(Patient_means[,-1],row.names=as.character(Patient_means[,1]))

# Patient_means now contains mean parameter vectors for each subset size.
# Each row corresponds to the mean of the first 'n_patients' patients.


## Model Definition - Equations ------------------------------------------------

model_name <- "Notal"
flist <- NULL %>% 
  addReaction("PK_VH", "","(log(2)/T12) * PK_VH") %>%
  addReaction("CST", "","(KIN / (CSTMIN + CSTDEL))*(1 + (CSTDEL / CSTMIN)*((PK_VH + 1E-08) ^ HILL)/(EC50 ^ HILL+(PK_VH + 1E-08) ^ HILL)) * CST") %>%
  addReaction("0", "EC50", "0" ) %>%
  addReaction("0", "CSTMIN", "0" ) %>%
  addReaction("0", "T12", "0" ) %>%
  addReaction("", "CST","KIN") 


## Observables definition ---------------------

observables <- eqnvec( CST_obs = "CST")

g <- Y(observables, flist, compile=TRUE, modelname=paste0("g_",model_name))

# error model

obsError <- c(names(observables))
errors <- as.eqnvec(
  c(paste0("sigma_",obsError,"_abs")
  ), names = c(obsError)
)
errPars <- getSymbols(errors)
errPars <- errPars[which(grepl("sigma", errPars))]

myerr <- Y(errors,
           f = c(observables, as.eqnvec(flist)),
           states = names(errors),
           attach.input = FALSE,
           compile = TRUE,
           modelname = paste0("err_", model_name))

## event generation of the different dosings

myevents <- NULL 
myevents <-  addEvent(myevents, var = "EC50", time = runin_schedule[1], value = "EC50_Eylea", method = "replace") 
myevents <-  addEvent(myevents, var = "T12", time = runin_schedule[1], value = 9, method = "replace") 

for(i in runin_schedule){
  myevents <-  addEvent(myevents, var = "PK_VH", time = i, value = 2, method = "add") 
}




for(i in SoC_schedule){
  myevents <-  addEvent(myevents, var = "PK_VH", time = i, value = 2, method = "add") 
}
myevents <-  addEvent(myevents, var = "EC50", time = SoC_schedule[1], value = "EC50_Eylea", method = "replace") 
myevents <-  addEvent(myevents, var = "T12", time = SoC_schedule[1], value = 9, method = "replace") 

## Model Generation ---------------------

modelHBV <- odemodel(flist, 
                     events = myevents,
                     
                     modelname = paste0("odemodel_", model_name),
                     compile = TRUE)


constraints <- resolveRecurrence(c(
  PK_VH = "PK0",
  CST =  "CST0",
  EMAX = "CSTDEL / CSTMIN",
  KOUT = "KIN / (CSTMIN + CSTDEL)"
))

innerpars <- unique(c(getParameters(modelHBV), getSymbols(observables), errPars))
names(innerpars) <- innerpars

trafo <- replaceSymbols(names(observables), observables, innerpars)
trafo <- replaceSymbols(names(constraints), constraints, trafo)

trafoL <- branch(trafo, table=Patient_pars) %>%
  insert(x~10^y, x=c("EC50_Eylea"), y=EC50_Eylea) %>%

  insert(x~10^y, x=c("EC50"), y=EC50_Eylea) %>%
  insert(x~y, x=c("T12"), y=9) %>%
  
  insert(x~(10**(x)), x = .currentSymbols)




tolerances <- 1e-7

p0 <- x <- NULL
for (C in names(trafoL)) {
  p0 <- p0 + P(trafoL[[C]], condition = C)
  x <- x + Xs(modelHBV, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
              optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
              condition = C)
}

trafoL_mean <- branch(trafo, table=Patient_means) %>%
  insert(x~10^y, x=c("EC50_Eylea"), y=EC50_Eylea) %>%
  insert(x~10^y, x=c("EC50"), y=EC50_Eylea) %>%
  insert(x~y, x=c("T12"), y=9) %>%
  
  insert(x~(10**(x)), x = .currentSymbols)




tolerances <- 1e-7

p0_mean <- x_mean <- NULL
for (C in names(trafoL_mean)) {
  p0_mean <- p0_mean + P(trafoL_mean[[C]], condition = C)
  x_mean <- x_mean + Xs(modelHBV, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
              optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
              condition = C)
}

## Generate objective function and initial parameter set -------

times <- seq(0,max(SoC_schedule)+112,1)

# Set seed for reproducibility
set.seed(123)

# Initialize an empty vector to store the time points
time_points <- c()

# Loop through each week
for (week in 1:((max(SoC_schedule)+112)/7)) {
  # Generate a random number of days (5 or 6) for the current week
  num_days <- sample(5:6, 1)
  
  # Generate random days within the week (1 to 7)
  week_days <- sample(1:7, num_days)
  
  # Convert week and days to a single time point (assuming week starts at day 1)
  week_time_points <- (week - 1) * 7 + week_days
  
  # Append the generated time points to the main vector
  time_points <- c(time_points, week_time_points)
}

# Sort the time points (optional, depending on your use case)
times_OnSite_more <- seq(0,max(SoC_schedule)+112,14)

times_Notal <- unique(sort(c(0,runin_schedule,SoC_schedule,time_points, times_OnSite_more)))

times_OnSite <- seq(0,max(SoC_schedule)+112,28)

Pred_patients_means <- c()
rownames(Pred_patients_means)<- NULL
dummy <- c()
j <- 1
for (i in 1:nrow(Patient_means)){
  
  par <- unlist(Patient_means[i,], use.names = FALSE)
  names(par) <- names(Patient_means[1,])
  dummy <-  (g*x_mean*p0_mean)(times,par, conditions = rownames(Patient_means[i,]))[[1]]
  dummy <- as.data.frame(wide2long(dummy))
  dummy$Patients <- rownames(Patient_means[i,])
  Pred_patients_means <- rbind(Pred_patients_means, dummy)
  dummy <- c()
}
save(Pred_patients_means, file = "Pred_Patients_mean.RData")
save(Patient_means, file = "Patient_mean.RData")
Data_patients_Notal <- c()
rownames(Patient_pars)<- NULL
dummy <- c()
j <- 1
for (i in unique(row.names(Patient_pars))){

  par <- unlist(Patient_pars[i,], use.names = FALSE)
  names(par) <- names(Patient_pars[1,])
  dummy <-  (g*x*p0)(times_Notal,par, conditions = i)[[1]]
  dummy <- as.data.frame(wide2long(dummy))
  dummy$ID <- j
  j <- j+1
  dummy$noise <- par["sigma_CST_obs_abs"][[1]]* rnorm(length(dummy$value), 0 , sd = 1)
  dummy$value_clean <- dummy$value
  dummy$value <- dummy$value + dummy$noise
  Data_patients_Notal <- rbind(Data_patients_Notal, dummy)
  dummy <- c()
}

write.csv(Data_patients_Notal, file = "Data_patients_Notal_SoC.csv")

Data_patients_OnSite_more <- subset(Data_patients_Notal, time %in% times_OnSite_more)

write.csv(Data_patients_OnSite_more, file = "Data_patients_OnSite_more_SoC.csv")







Patient_pars_plot <- gather(Patient_pars, name, value)
Patient_pars_plot$value <- log(10**as.numeric(Patient_pars_plot$value)) 

P4 <- ggplot(data =Patient_pars_plot , aes(x= value))+
  facet_wrap(~name,scale="free" )+
  # geom_histogram( binwidth = .1)+
  geom_density(aes(y = after_stat(density), color = "Distribution of parameters"))+
  stat_theodensity(aes(y = after_stat(density), color = "Gaussian fit"), distri = "norm",geom = "line")+
  theme_dMod()+ 
  scale_color_dMod(name = "Color" )

ggsave("Distribution_all_parameters_SoC_drug.pdf",P4, width = 24, height= 20)


## NONMEM Datafile generation-------------------------------

data_onsite <- Data_patients_OnSite[c(-5)]
data_onsite <- subset(data_onsite, time >= min(SoC_schedule))
data_onsite$time <- data_onsite$time - min(SoC_schedule)
names(data_onsite)[names(data_onsite) == 'value'] <- 'DV'
names(data_onsite)[names(data_onsite) == 'time'] <- 'TIME'
data_onsite$CMT <- 2
data_onsite$AMT <- 0
data_onsite$DOSE <- 0
data_onsite$MDV <-0

data_onsite <- subset(data_onsite, name=="CST_obs")
data_onsite <- data_onsite[c(-2)]
data_onsite <- data_onsite[c(-4)]

dummy <- c()
j <- 1
for (i in unique(row.names(Patient_pars))){
  dummy <- data.frame(
    TIME = SoC_schedule - min(SoC_schedule),
    CMT = 1,
    AMT = 2,
    DV = 0,
    DOSE = 2,
    MDV = 1,
    ID = as.numeric(j) )
  j <- j+1
  data_onsite <- rbind(data_onsite, dummy)
  dummy <- c()
}
data_onsite$AGE <- 75
data_onsite$CENTER <- 3
data_onsite$DRUG <- 5
data_onsite$DIA <- 5
data_onsite$C <- 0
data_onsite$EYE <- 1
data_onsite$TAD <- 0

data_onsite <- arrange(data_onsite, ID, TIME, CMT)

write.csv(data_onsite, file = "Data_patients_OnSite_nonmem_SoC.csv", row.names = FALSE, quote = FALSE)

data_notal <- Data_patients_Notal[c(-5)]
data_notal <- subset(data_notal, time >= min(SoC_schedule))
data_notal$time <- data_notal$time - min(SoC_schedule)
names(data_notal)[names(data_notal) == 'value'] <- 'DV'
names(data_notal)[names(data_notal) == 'time'] <- 'TIME'
data_notal$CMT <- 2
data_notal$AMT <- 0
data_notal$DOSE <- 0
data_notal$MDV <-0

data_notal <- subset(data_notal, name=="CST_obs")
data_notal <- data_notal[c(-2)]
data_notal <- data_notal[c(-4)]


dummy <- c()
j <- 1
for (i in unique(row.names(Patient_pars))){
  dummy <- data.frame(
    TIME = SoC_schedule- min(SoC_schedule),
    CMT = 1,
    AMT = 2,
    DV = 0,
    DOSE = 2,
    MDV = 1,
    ID = as.numeric(j) )
  j <- j+1
  data_notal <- rbind(data_notal, dummy)
  dummy <- c()
}
data_notal$AGE <- 75
data_notal$CENTER <- 3
data_notal$DRUG <- 5
data_notal$DIA <- 5
data_notal$C <- 0
data_notal$EYE <- 1
data_notal$TAD <- 0

data_notal <- arrange(data_notal, ID, TIME, CMT)

write.csv(data_notal, file = "Data_patients_Notal_nonmem_SoC.csv", row.names = FALSE, quote = FALSE)

data_onsite_more <- Data_patients_OnSite_more[c(-5)]
data_onsite_more <- subset(data_onsite_more, time >= min(SoC_schedule))
data_onsite_more$time <- data_onsite_more$time - min(SoC_schedule)
names(data_onsite_more)[names(data_onsite_more) == 'value'] <- 'DV'
names(data_onsite_more)[names(data_onsite_more) == 'time'] <- 'TIME'
data_onsite_more$CMT <- 2
data_onsite_more$AMT <- 0
data_onsite_more$DOSE <- 0
data_onsite_more$MDV <-0

data_onsite_more <- subset(data_onsite_more, name=="CST_obs")
data_onsite_more <- data_onsite_more[c(-2)]
data_onsite_more <- data_onsite_more[c(-4)]

dummy <- c()
j <- 1
for (i in unique(row.names(Patient_pars))){
  dummy <- data.frame(
    TIME = SoC_schedule - min(SoC_schedule),
    CMT = 1,
    AMT = 2,
    DV = 0,
    DOSE = 2,
    MDV = 1,
    ID = as.numeric(j) )
  j <- j+1
  data_onsite_more <- rbind(data_onsite_more, dummy)
  dummy <- c()
}
data_onsite_more$AGE <- 75
data_onsite_more$CENTER <- 3
data_onsite_more$DRUG <- 5
data_onsite_more$DIA <- 5
data_onsite_more$C <- 0
data_onsite_more$EYE <- 1
data_onsite_more$TAD <- 0

data_onsite_more <- arrange(data_onsite_more, ID, TIME, CMT)

write.csv(data_onsite_more, file = "Data_patients_OnSite_more_nonmem_SoC.csv", row.names = FALSE, quote = FALSE)

