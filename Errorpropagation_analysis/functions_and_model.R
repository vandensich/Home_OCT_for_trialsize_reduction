try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

library(ggh4x)
library(dMod)
library(tidyverse)
library(tidyr)
library(tibble)
library(dplyr)
library(boot)
library(truncnorm)
library(MASS)  # For mvrnorm
library(Matrix)  # For nearPD
# library(ggbreak)

# set the Trial conditions like dosing times and patient numbers


trial_schedule <- c(seq(0, 2*56, 56))
n_patients <- 1

# generate patient parameters based on the Parameter specific distributions

Patient_pars <- c()

set.seed(42)

dummy <- c()

generate_dummy <- function(df) {
  # Initialize an empty list to store results
  dummy <- list()
  
  # Loop through the rows of the input data frame
  for (i in 1:nrow(df)) {
    # Extract the variable name and FE value
    variable <- unlist(as.character(df$Parameter[i]))  # Extract second part of the 'Name'
    fe_value <- df$Estimate[i]
    
    # Check if the variable is 'sigma_CST_obs_abs' and handle differently
    if (variable == "additive error CST") {
      dummy[[variable]] <- fe_value  # Use the raw value for sigma
    } else {
      # dummy[[variable]] <- log10(exp(fe_value))  # Apply the transformation for other variables
      dummy[[variable]] <- fe_value  # Apply the transformation for other variables
      
    }
  }
  
  # Convert the list to a data frame
  dummy <- as.data.frame(dummy)
  colnames(dummy)[ncol(dummy)] <- "sigma_CST_obs_abs"
  return(dummy)
}

## Model Definition - Equations ------------------------------------------------

model_name <- "Notal"
flist <- NULL %>% 
  addReaction("PK_VH", "","(log(2)/T12) * PK_VH") %>%
  addReaction("CST", "","(KIN / (CSTMIN + CSTDEL))*(1 + (CSTDEL / CSTMIN)*((PK_VH + 1E-08) ^ HILL)/(EC50 ^ HILL+(PK_VH + 1E-08) ^ HILL)) * CST") %>%
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




for(i in trial_schedule){
  myevents <-  addEvent(myevents, var = "PK_VH", time = i, value = 2, method = "add") 
}

## Model Generation ---------------------

modelHBV <- odemodel(flist, 
                     events = myevents,
                     
                     modelname = paste0("odemodel_", model_name),
                     compile = TRUE,  verbose = TRUE,deriv = TRUE, jacobian = c("full"))


constraints <- resolveRecurrence(c(
  PK_VH = "PK0",
  CST =  "CST0"
))

innerpars <- unique(c(getParameters(modelHBV), getSymbols(observables), errPars))
names(innerpars) <- innerpars

trafo <- replaceSymbols(names(observables), observables, innerpars)
trafo <- replaceSymbols(names(constraints), constraints, trafo)

## Generate objective function and initial parameter set -------

times <- seq(0,max(trial_schedule)+56,1)


create_jacobian_matrices <- function(df, jacobian) {
  # Extract time points
  time_points <- df$time
  jacobian_of_p0 <- jacobian
  # Remove the time column for matrix creation
  df_data <- df[, -1]
  
  # Split each entry name into parts using the period (.)
  split_column_names <- strsplit(names(df_data), "\\.")
  
  # Extract unique row and column names
  row_names <- unique(sapply(split_column_names, `[`, 1))
  col_names <- unique(sapply(split_column_names, `[`, 2))
  
  # Initialize a list to store matrices
  timepoint_matrices <- list()
  
  # Loop over each time point to create corresponding matrix
  for (i in 1:nrow(df_data)) {
    # Create an empty matrix
    time_matrix <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                          dimnames = list(row_names, col_names))
    
    # Populate the matrix with the respective values
    for (j in 1:ncol(df_data)) {
      row_name <- split_column_names[[j]][1]
      col_name <- split_column_names[[j]][2]
      time_matrix[row_name, col_name] <- as.numeric(df_data[i, j])
    }
    
    # divide each row of the deriv by the jacobian of the parameter transformation to obtain the sensitivity equation
    # for (f in 1:nrow(time_matrix)){
    #   time_matrix[f,]<- t(jacobian_of_p0**-1) * time_matrix[f,]
    # }
    
    timepoint_matrices[[i]] <- time_matrix
  }
  
  return(timepoint_matrices)
}


extract_theta_block <- function(filename) {
  # Read the lines from the file
  lines <- readLines(filename)
  
  # Header processing
  header_line <- unlist(strsplit(lines[2], split = "\\s+"))
  theta_indices <- which(grepl("THETA", header_line))
  
  # Extracting the THETA block
  theta_block <- matrix(0, nrow = length(theta_indices), ncol = length(theta_indices))
  
  for (i in 1:length(theta_indices)) {
    line_elements <- unlist(strsplit(lines[2 + i], split = "\\s+"))
    row_values <- as.numeric(line_elements[theta_indices])
    theta_block[i, ] <- row_values
  }
  
  return(theta_block)
}



get_errors_from_file <- function(filename) {
  # Extract number of patients and setting from the filename
  file_parts <- unlist(strsplit(filename, "_"))
  n_patients <- as.numeric(file_parts[1])
  setting <- ifelse(grepl("OS", filename), "OnSite", "Notal")
  drug <- ifelse(grepl("Trialdrug", filename), "Trialdrug", "Eylea")
  
  # Extract patient parameters
  Patient_pars <- dummy <- c()
  table <- read.csv(paste0(file_parts[1], "_", file_parts[2],"_",file_parts[3],"_",file_parts[4],"_",file_parts[5], "_results.csv"))
  dummy <- generate_dummy(table)
  Patient_pars <- as.data.frame(dummy)

  trafoL <- branch(trafo, table=Patient_pars) %>%
    insert(x~y, x= c("CSTMIN"), y= c("CSTMIN_Trialdrug")) %>%


    insert(x~y, x=c("T12"), y=9) %>%

    # insert(x~(10**(x)), x = .currentSymbols)
    insert(x~(exp(x)), x = .currentSymbols)


  tolerances <- 1e-6

  p0 <- x <- NULL
  for (C in names(trafoL)) {
    p0 <- p0 + P(trafoL[[C]], condition = C)
    x <- x + Xs(modelHBV, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
                optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                condition = C)
  }

  # Model definition and error calculation
  par <- unlist(Patient_pars[1,], use.names = FALSE)
  names(par) <- names(Patient_pars[1,])
  dummy <- (g * x * p0)(times, par)
  error <- getDerivs(dummy)

  # Extract the Jacobian of the parameter transformation function
  jacobian_of_p0 <- attr(p0(par)[[1]], "deriv")[which(attr(p0(par)[[1]], "deriv") != 0)]
  names(jacobian_of_p0) <- names(par)

  # Generate the sensitivity matrices
  timepoint_matrices <- create_jacobian_matrices(as.data.frame(error[[1]]), jacobian_of_p0)

  # Extract the theta block and transform covariance matrix
  theta_block <- read.csv(paste0(file_parts[1], "_", file_parts[2],"_",file_parts[3],"_",file_parts[4],"_",file_parts[5], "_covariance_matrix.csv"))
  # Transform covariance matrix from log-space to linear-space DELTA METHOD 1/ comes from deriv of ln() = Jacobian of parameter transformation
  # theta_block <- diag(10^(Patient_pars)) %*% theta_block %*% t(diag(10^(Patient_pars)))
  # theta_block <- diag(1 / (log(10^(Patient_pars)))) %*% theta_block %*% t(diag(1 / (log(10^(Patient_pars)))))
  # theta_block <- diag(exp(Patient_pars)) %*% theta_block %*% t(diag(exp(Patient_pars)))

  colnames(theta_block) <- names(Patient_pars)
  rownames(theta_block) <- names(Patient_pars)
  # print(filename)
  # print(timepoint_matrices[[155]])
  # print(theta_block)
  # print(length(timepoint_matrices))
  # Calculate absolute and relative errors
  abs_error <- sqrt(diag(timepoint_matrices[[169]] %*% as.matrix(theta_block) %*% t(timepoint_matrices[[169]])))
  rel_error <- abs_error * 100 / dummy[[1]][169, 2:4]
  res <- data.frame(
    Absolute_Error = abs_error,
    Relative_Error = rel_error,
    Number_of_Patients = n_patients,
    Setting = setting,
    Drug = drug
  )
  res$name <- rownames(res)
  rownames(res) <- NULL
  # Return results
  return(res)
}


get_errors_via_MC <- function(filename, N = 3000) {
  
  
  # Extract number of patients and setting from the filename
  file_parts <- unlist(strsplit(filename, "_"))
  n_patients <- as.numeric(file_parts[1])
  setting <- ifelse(grepl("OS", filename), "OnSite", "Notal")
  drug <- ifelse(grepl("Trialdrug", filename), "Trialdrug", "Eylea")
  
  # Extract patient parameters (log-space)
  table <- read.csv(paste0(file_parts[1], "_", file_parts[2],"_",file_parts[3],"_",file_parts[4],"_",file_parts[5], "_results.csv"))
  dummy <- generate_dummy(table)
  Patient_pars <- as.data.frame(dummy)
  
  # Read the covariance matrix (log-space)
  theta_block <- as.matrix(read.csv(paste0(file_parts[1], "_", file_parts[2],"_",file_parts[3],"_",file_parts[4],"_",file_parts[5], "_covariance_matrix.csv")))
  colnames(theta_block) <- names(Patient_pars)
  rownames(theta_block) <- names(Patient_pars)
  
  # For each simulation
  CST_values <- numeric(N)
  times <- seq(0,max(trial_schedule)+56,1)
  
  for (i in 1:N) {
    # Sample parameters from multivariate normal in log-space
    sampled_pars_log <- c(mvrnorm(1, mu = as.numeric(Patient_pars[1, ]), Sigma = theta_block))
    names(sampled_pars_log) <- names(Patient_pars)
    
    
    sampled_pars <- as.data.frame(as.list(sampled_pars_log))
    # Create the parameter transformation
    trafoL <- branch(trafo, table=sampled_pars) %>%
      insert(x~y, x= c("CSTMIN"), y= c("CSTMIN_Trialdrug")) %>%
      insert(x~y, x=c("T12"), y=9) %>%
      insert(x~exp(x), x = .currentSymbols)
    
    tolerances <- 1e-6
    
    p0 <- x <- NULL
    for (C in names(trafoL)) {
      p0 <- p0 + P(trafoL[[C]], condition = C)
      x <- x + Xs(modelHBV, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
                  optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                  condition = C)
    }
    
    # Run the model
    par <- sampled_pars  # Parameters are in linear scale
    result <- tryCatch({
      dummy <- (g * x * p0)(times, par)
      # Extract CST at desired timepoint, e.g., time = 154
      CST_values[i] <-dummy[[1]][169,2]
    }, error = function(e) {
      CST_values[i] <- NA
    })
  }
  print(paste0(filename, " done"))
  # Remove NA values from CST_values
  CST_values <- CST_values[!is.na(CST_values)]
  
  # Compute statistics
  mean_CST <- mean(CST_values)
  sd_CST <- sd(CST_values)
  
  # Bootstrap the SE itself
  
  boot_se_samples <- replicate(1000, {
    resampled_CST <- sample(CST_values, replace = TRUE)
    sd(resampled_CST)
  })
  ci_level <- 0.95
  
  lower_ci <- quantile(boot_se_samples, (1 - ci_level) / 2)
  upper_ci <- quantile(boot_se_samples, 1 - (1 - ci_level) / 2)
  
  
  median_CST <- median(CST_values)
  median_absolute_deviation <- mad(CST_values)
  rel_error <- sd_CST * 100 / mean_CST
  
  res <- data.frame(
    Mean_CST = mean_CST,
    Absolute_Error = sd_CST,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    Relative_Error = rel_error,
    Median_CST = median_CST,
    Median_Absolute_Deviation = median_absolute_deviation,
    Number_of_Patients = n_patients,
    Setting = setting,
    Drug = drug
  )
  
  return(res)
}

