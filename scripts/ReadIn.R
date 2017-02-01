#### Read in data and clean up ####

# This script was written to pre-process MetaLab compatible spreadsheets on statistical learning studies (segmenting artificial mini-languages based on statistics). 
# Author: Christina Bergmann
# chbergma'at'gmail.com
# last modified: Feb 1, 2017

#### Data read in ####

db = read.csv("data/StatLEarnDB.csv")

#### Some cleanup ####

# Rename columns for the dependent variables. We expect a novelty prefrence, so x_1 will correspond to the novel items, and x_2 to the familiar (statistical words with high TP)

db$x_1 = db$LTP_x_2
db$x_2 = db$HTP_x_1

# All other cleanup takes place in the .Rmd for maximum transparency