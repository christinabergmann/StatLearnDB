#### Bias and Replicability Analyses ####

# This script was written to analyze bias and replicability in statistical learning studies (segmenting artificial mini-languages based on statistics). 
# Author: Christina Bergmann
# chbergma'at'gmail.com
# last modified: January 20, 2017

#### Load libraries ####

library(metafor)
library(pwr)

#### get data ####

source("scripts/MetaAnalysis.R")

#### Publication Bias ####

# Publication bias can become visible in funnel plot asymmetry. 
# Metafor comes with several options to check this, I chose ranktest and regtest

ranktest(rma_all)
ranktest(rma_DR)

#No longer possible with switch from rma to rma.mv
#regtest(rma_all)
#regtest(rma_DR)

# Both datasets seem ok. 



#### Export Data for p-curve ####
#The very pedestrian way because reasons. 


for(line in 1:length(db$n_1)){
  if(!is.na(db[line,]$t)){
    newline = paste("t(", db[line,]$n_1-1, ")=",db[line,]$t, sep = "")
    #write(newline, file = textfile, append = TRUE)
    print(newline)
    #write(newline, file = textfile, append = TRUE)    
}
}

#### Power ####

# What is the average power to detect an effect (to be compared with p-curve estimate)
d = rma.mv(g_calc, g_var_calc, data = db_DR, random = ~ 1 | short_cite)$b[,1]
# Let's try to extract the ES automagically later, for now I enter the value I had when running fir model by hand

pwr.t.test(n = median(db$n_1, na.rm=TRUE), d = .228, sig.level = .05, type = "paired", alternative = "two.sided")

pwr.t.test( d = .228, sig.level = .05, type = "paired", power = .8, alternative = "two.sided")


#### Original effect sizes vs MA effect sizes ####

original_effects = db[db$study_ID=="SaffranAslinNewport1996",]


#### What was the power of the original study? #####

pwr.t.test(n = median(db$n_1, na.rm=TRUE), d = median(original_effects$g_calc), sig.level = .05, type = "paired", alternative = "two.sided")
pwr.t.test(d = median(original_effects$g_calc), power = .8, sig.level = .05, type = "paired", alternative = "two.sided")
