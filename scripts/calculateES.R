#### Calculate Effect Sizes ####

# This script was written to calculate effect sizes in MetaLAb compatible spreadsheets on statistical learning studies (segmenting artificial mini-languages based on statistics). 
# Author: Christina Bergmann
# chbergma'at'gmail.com
# last modified: Feb 1, 2017

#### Load libraries ####

#None needed right now

#### Get data ####

source("scripts/ReadIn.R")

#### calculate paired ES based on available data ####

# The formulae here are based on those in MetaLab
# Link https://github.com/langcog/metalab/blob/master/scripts/compute_es.R

# We don't have any correlations in the dataset, so we will for now randomly impute them base on the mean and SD of correlations found in InWordDB. 

median_corr = .635
var_corr = .1 
#I am reducing the variance here, the actual SD is .177 

set.seed(111)
db$corr_imputed = rnorm(length(db$n_1), mean = median_corr, sd = var_corr)

# Initiate the additional columns for d and d_var

db$d_calc <- NA
db$d_var_calc <- NA
db$es_method <- "missing"

# Effectively, all correlations are now imputed. 
db$corr = ifelse(is.na(db$corr), db$corr_imputed, db$corr)

#Correlations are never perfect and cannot be higher than 1, so fixing some possible imputation issues first
db$corr = ifelse(db$corr>.99, .99, db$corr)

#Now compute effect sizes based on the data we have, row by row. MetaLAb solves this a bit differently, here I opt for readability.

for(line in 1:length(db$n_1)){
  if(db[line,]$participant_design == "within_two"){
    if (complete.cases(db[line,]$x_1, db[line,]$x_2, db[line,]$SD_1, db[line,]$SD_2)) {
      # Lipsey & Wilson, 3.14
      pooled_SD <- sqrt((db[line,]$SD_1 ^ 2 + db[line,]$SD_2 ^ 2) / 2)
      db[line,]$d_calc <- (db[line,]$x_1 - db[line,]$x_2) / pooled_SD
      db[line,]$es_method  <- "group_means_two"
    } else if (complete.cases(db[line,]$t)) {
      #Dunlap et al., 1996, p.171
      wc <- sqrt(2 * (1 - db[line,]$corr))
      db[line,]$d_calc <- (db[line,]$t / sqrt(db[line,]$n_1)) * wc
      db[line,]$es_method  <- "t_two"
    } else if (complete.cases(db[line,]$F)) {
      #No case here, reduces variance from varied ES calculations
      wc <- sqrt(2 * (1 - db[line,]$corr))
      db[line,]$d_calc <- sqrt(db[line,]$F / db[line,]$n_1) * wc
      db[line,]$es_method  <- "f_two"
    }
    #Next step: effect size variance (needed for weighting the effect sizes)
    #Lipsey & Wilson (2001) 
    if (complete.cases(db[line,]$n_1, db[line,]$d_calc)) {
      db[line,]$d_var_calc <- ((1 / db[line,]$n_1) + (db[line,]$d_calc ^ 2 / (2 * db[line,]$n_1))) * 2 * (1 - db[line,]$corr)
    } 
  }else if(db[line,]$participant_design == "between"){
    #we only have one case here so far which has means and SD, so we just enter the approriate formula
    # Lipsey & Wilson (2001), 3.23
    pooled_SD <- sqrt(((db[line,]$n_1 - 1) * db[line,]$SD_1 ^ 2 + (db[line,]$n_2 - 1) * db[line,]$SD_2 ^ 2) / (db[line,]$n_1 + db[line,]$n_2 - 2))
    db[line,]$d_calc <- (db[line,]$x_1 - db[line,]$x_2) / pooled_SD 
    db[line,]$d_var_calc <- ((db[line,]$n_1 + db[line,]$n_2) / (db[line,]$n_1 * db[line,]$n_2)) + (db[line,]$d_calc ^ 2 / (2 * (db[line,]$n_1 + db[line,]$n_2)))
    db[line,]$es_method  <- "between"
  }
}

#### Remove missing values ####

db = db[!is.na(db$d_calc),]


#### Compute Hedge's g based on Cohen's d ####

# Morris, 2010, p. 21
J <- 1 - 3 / (4 * (db$n_1 - 1 - 1))
db$g_calc <- db$d_calc * J
db$g_var_calc <- J ^ 2 * db$d_var_calc

#### Add weights ####

db$weights_g <-  1/(db$g_var_calc)^2

db$weights_d <-  1/(db$d_var_calc)^2


#### Outlier removal ####

# Following InWordD/standard practice, we will remove effect sizes more than 3 SD away from the median effect (in both positive and negative directions)

db$nooutlier = ifelse(db$g_calc > mean(db$g_calc, na.rm = TRUE) + 3*sd(db$g_calc, na.rm = TRUE) 
                         | db$g_calc < mean(db$g_calc, na.rm = TRUE) - 3*sd(db$g_calc, na.rm = TRUE),FALSE, TRUE)

db = db[db$nooutlier,]
#Turns out there are no outlier in the dataset

#### Add weights ####

db$weights_g <-  1/(db$g_var_calc)^2

db$weights_d <-  1/(db$d_var_calc)^2