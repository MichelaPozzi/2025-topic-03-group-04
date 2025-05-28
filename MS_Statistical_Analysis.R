#### LICENSE OF USE ####
# BSD 3-Clause License
# Copyright (c) 2019, Maiwen Caudron-Herger
# All rights reserved.
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#  1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#################################################
#### Mass spectrometry analysis Step-by-Step ####
#################################################


# Remark:
# A lot of help for the development of this analysis was found within the stack overflow community (https://stackoverflow.com/)

# Save the script under the name filename.r of your choice
# Execute the analysis with the command: source("filename.r")
# It will take about 20-30 min for a dataset with 5000 proteins or shorter, depending on the computer
# It will take about 5 to 10 min for the provided sample dataset


#### This script is divided into several parts as listed below: ####
# PART 1 : read raw data; normalization between replicates; sliding window; normalization to sum = 100
# PART 2 : Finding the fit parameters from average curves (every local maximum above a value of 2)
# PART 3 : Gaussian fit on each curves (three control and three rnase replicates for each protein)
# PART 4 : quality control of the fit and adjustment if possible
# PART 5 : Assessment of p_values and FDR-adjusted p-values for all maxima
# PART 6 : Evaluation of the shifts - use of selection criteria
# PART 7 - Computing of the last table and export as CSV file
# GRAPHICS 1 - average curves for one protein of interest - function curves_plot  
# GRAPHICS 2 - all curves for one protein of interest - function curves_all


#*************************************************************************************************************
#### PART 1 : read raw data; normalization between replicates; sliding window; normalization to sum = 100 ####
#*************************************************************************************************************

# Set working directory - the directory where you store your raw data
# setwd("path to the directory")

#**************************************************************************************************
# Read the raw data from CSV file
# There is a header and the name of the proteins are stored in the first column. The values are separated by a ";"
table <- read.table("Mass_Spec_RawData_Sample.csv", header=TRUE, row.names=1, sep = ";")

# Run head(table) to see the first lines of the object table.
# In our case, we have a table with 150 columns (3x 25 control fractions and 3x 25 rnase fractions)
																					

#**************************************************************************************************
# Define the information found in the table
# Use factors to define the treatment in the samples as control (CTRL) or rnase (RNASE) - there are as many treatments as columns in the table. Here, we have 150 treatments
treatment <- factor(rep(c("CTRL", "RNASE"),75))
# Use factors to define the condition in the samples. Here we have all six replicates for 1 fraction after each other and then the next fractions
condition <- factor(rep(c("ctrl1","rnase1","ctrl2","rnase2","ctrl3","rnase3"),25))
# Use factors to define the fractions from fraction1 to fraction25 through the 150 columns. Here, we have first all six fraction1, then all six fraction2 etc ... until all six fraction25
fraction <- as.vector(matrix(rep(paste("fraction",1:25,sep=""),6), nrow = 6, ncol=25, byrow = TRUE))


# Store the number of rows in a variable
n_row <- nrow(table)
# Store the row names (here the name sof the proteins) in a variable
row_names <- rownames(table)


# Make dataframe with the information (treatment, condition, fraction) for each column - here 150 columns
data <- data.frame(row.names = colnames(table), treatment = treatment, condition = condition, fraction = fraction)


#**********************************************************************************************************************
#### Normalization step between the replicates in each fraction and for each treatment using the mean value method ####

# Define subtables for each fraction (with three CTRL and three RNASE replicates)
table.frxn1 <- table[,data$fraction =="fraction1"]
table.frxn2 <- table[,data$fraction =="fraction2"]
table.frxn3 <- table[,data$fraction =="fraction3"]
table.frxn4 <- table[,data$fraction =="fraction4"]
table.frxn5 <- table[,data$fraction =="fraction5"]
table.frxn6 <- table[,data$fraction =="fraction6"]
table.frxn7 <- table[,data$fraction =="fraction7"]
table.frxn8 <- table[,data$fraction =="fraction8"]
table.frxn9 <- table[,data$fraction =="fraction9"]
table.frxn10 <- table[,data$fraction =="fraction10"]
table.frxn11 <- table[,data$fraction =="fraction11"]
table.frxn12 <- table[,data$fraction =="fraction12"]
table.frxn13 <- table[,data$fraction =="fraction13"]
table.frxn14 <- table[,data$fraction =="fraction14"]
table.frxn15 <- table[,data$fraction =="fraction15"]
table.frxn16 <- table[,data$fraction =="fraction16"]
table.frxn17 <- table[,data$fraction =="fraction17"]
table.frxn18 <- table[,data$fraction =="fraction18"]
table.frxn19 <- table[,data$fraction =="fraction19"]
table.frxn20 <- table[,data$fraction =="fraction20"]
table.frxn21 <- table[,data$fraction =="fraction21"]
table.frxn22 <- table[,data$fraction =="fraction22"]
table.frxn23 <- table[,data$fraction =="fraction23"]
table.frxn24 <- table[,data$fraction =="fraction24"]
table.frxn25 <- table[,data$fraction =="fraction25"]


# Re-define the information for the fraction subtables (dataframe named "data3") and create one table per treatment and per fraction. Here there are 50 subtables in total
data3 <- data.frame(row.names=colnames(table.frxn1), treatment = rep(c("CTRL", "RNASE"),3), condition = c("ctrl1","rnase1","ctrl2","rnase2","ctrl3","rnase3"))
table.frxn1.CTRL <- table.frxn1[,data3$treatment == "CTRL"] 
table.frxn1.RNASE <- table.frxn1[,data3$treatment == "RNASE"]
table.frxn2.CTRL <- table.frxn2[,data3$treatment == "CTRL"] 
table.frxn2.RNASE <- table.frxn2[,data3$treatment == "RNASE"]
table.frxn3.CTRL <- table.frxn3[,data3$treatment == "CTRL"] 
table.frxn3.RNASE <- table.frxn3[,data3$treatment == "RNASE"]
table.frxn4.CTRL <- table.frxn4[,data3$treatment == "CTRL"] 
table.frxn4.RNASE <- table.frxn4[,data3$treatment == "RNASE"]
table.frxn5.CTRL <- table.frxn5[,data3$treatment == "CTRL"] 
table.frxn5.RNASE <- table.frxn5[,data3$treatment == "RNASE"]
table.frxn6.CTRL <- table.frxn6[,data3$treatment == "CTRL"] 
table.frxn6.RNASE <- table.frxn6[,data3$treatment == "RNASE"]
table.frxn7.CTRL <- table.frxn7[,data3$treatment == "CTRL"] 
table.frxn7.RNASE <- table.frxn7[,data3$treatment == "RNASE"]
table.frxn8.CTRL <- table.frxn8[,data3$treatment == "CTRL"] 
table.frxn8.RNASE <- table.frxn8[,data3$treatment == "RNASE"]
table.frxn9.CTRL <- table.frxn9[,data3$treatment == "CTRL"] 
table.frxn9.RNASE <- table.frxn9[,data3$treatment == "RNASE"]
table.frxn10.CTRL <- table.frxn10[,data3$treatment == "CTRL"] 
table.frxn10.RNASE <- table.frxn10[,data3$treatment == "RNASE"]
table.frxn11.CTRL <- table.frxn11[,data3$treatment == "CTRL"] 
table.frxn11.RNASE <- table.frxn11[,data3$treatment == "RNASE"]
table.frxn12.CTRL <- table.frxn12[,data3$treatment == "CTRL"] 
table.frxn12.RNASE <- table.frxn12[,data3$treatment == "RNASE"]
table.frxn13.CTRL <- table.frxn13[,data3$treatment == "CTRL"] 
table.frxn13.RNASE <- table.frxn13[,data3$treatment == "RNASE"]
table.frxn14.CTRL <- table.frxn14[,data3$treatment == "CTRL"] 
table.frxn14.RNASE <- table.frxn14[,data3$treatment == "RNASE"]
table.frxn15.CTRL <- table.frxn15[,data3$treatment == "CTRL"] 
table.frxn15.RNASE <- table.frxn15[,data3$treatment == "RNASE"]
table.frxn16.CTRL <- table.frxn16[,data3$treatment == "CTRL"] 
table.frxn16.RNASE <- table.frxn16[,data3$treatment == "RNASE"]
table.frxn17.CTRL <- table.frxn17[,data3$treatment == "CTRL"] 
table.frxn17.RNASE <- table.frxn17[,data3$treatment == "RNASE"]
table.frxn18.CTRL <- table.frxn18[,data3$treatment == "CTRL"] 
table.frxn18.RNASE <- table.frxn18[,data3$treatment == "RNASE"]
table.frxn19.CTRL <- table.frxn19[,data3$treatment == "CTRL"] 
table.frxn19.RNASE <- table.frxn19[,data3$treatment == "RNASE"]
table.frxn20.CTRL <- table.frxn20[,data3$treatment == "CTRL"] 
table.frxn20.RNASE <- table.frxn20[,data3$treatment == "RNASE"]
table.frxn21.CTRL <- table.frxn21[,data3$treatment == "CTRL"] 
table.frxn21.RNASE <- table.frxn21[,data3$treatment == "RNASE"]
table.frxn22.CTRL <- table.frxn22[,data3$treatment == "CTRL"] 
table.frxn22.RNASE <- table.frxn22[,data3$treatment == "RNASE"]
table.frxn23.CTRL <- table.frxn23[,data3$treatment == "CTRL"] 
table.frxn23.RNASE <- table.frxn23[,data3$treatment == "RNASE"]
table.frxn24.CTRL <- table.frxn24[,data3$treatment == "CTRL"] 
table.frxn24.RNASE <- table.frxn24[,data3$treatment == "RNASE"]
table.frxn25.CTRL <- table.frxn25[,data3$treatment == "CTRL"] 
table.frxn25.RNASE <- table.frxn25[,data3$treatment == "RNASE"]


# Determine average values for each treatment, each fraction and each replicate. Here, the values for the three replicates are stored in the same table
avg.table.frxn1.CTRL <- sapply(table.frxn1.CTRL, function(x) mean(x))
avg.table.frxn2.CTRL <- sapply(table.frxn2.CTRL, function(x) mean(x))
avg.table.frxn3.CTRL <- sapply(table.frxn3.CTRL, function(x) mean(x))
avg.table.frxn4.CTRL <- sapply(table.frxn4.CTRL, function(x) mean(x))
avg.table.frxn5.CTRL <- sapply(table.frxn5.CTRL, function(x) mean(x))
avg.table.frxn6.CTRL <- sapply(table.frxn6.CTRL, function(x) mean(x))
avg.table.frxn7.CTRL <- sapply(table.frxn7.CTRL, function(x) mean(x))
avg.table.frxn8.CTRL <- sapply(table.frxn8.CTRL, function(x) mean(x))
avg.table.frxn9.CTRL <- sapply(table.frxn9.CTRL, function(x) mean(x))
avg.table.frxn10.CTRL <- sapply(table.frxn10.CTRL, function(x) mean(x))
avg.table.frxn11.CTRL <- sapply(table.frxn11.CTRL, function(x) mean(x))
avg.table.frxn12.CTRL <- sapply(table.frxn12.CTRL, function(x) mean(x))
avg.table.frxn13.CTRL <- sapply(table.frxn13.CTRL, function(x) mean(x))
avg.table.frxn14.CTRL <- sapply(table.frxn14.CTRL, function(x) mean(x))
avg.table.frxn15.CTRL <- sapply(table.frxn15.CTRL, function(x) mean(x))
avg.table.frxn16.CTRL <- sapply(table.frxn16.CTRL, function(x) mean(x))
avg.table.frxn17.CTRL <- sapply(table.frxn17.CTRL, function(x) mean(x))
avg.table.frxn18.CTRL <- sapply(table.frxn18.CTRL, function(x) mean(x))
avg.table.frxn19.CTRL <- sapply(table.frxn19.CTRL, function(x) mean(x))
avg.table.frxn20.CTRL <- sapply(table.frxn20.CTRL, function(x) mean(x))
avg.table.frxn21.CTRL <- sapply(table.frxn21.CTRL, function(x) mean(x))
avg.table.frxn22.CTRL <- sapply(table.frxn22.CTRL, function(x) mean(x))
avg.table.frxn23.CTRL <- sapply(table.frxn23.CTRL, function(x) mean(x))
avg.table.frxn24.CTRL <- sapply(table.frxn24.CTRL, function(x) mean(x))
avg.table.frxn25.CTRL <- sapply(table.frxn25.CTRL, function(x) mean(x))

avg.table.frxn1.RNASE <- sapply(table.frxn1.RNASE, function(x) mean(x))
avg.table.frxn2.RNASE <- sapply(table.frxn2.RNASE, function(x) mean(x))
avg.table.frxn3.RNASE <- sapply(table.frxn3.RNASE, function(x) mean(x))
avg.table.frxn4.RNASE <- sapply(table.frxn4.RNASE, function(x) mean(x))
avg.table.frxn5.RNASE <- sapply(table.frxn5.RNASE, function(x) mean(x))
avg.table.frxn6.RNASE <- sapply(table.frxn6.RNASE, function(x) mean(x))
avg.table.frxn7.RNASE <- sapply(table.frxn7.RNASE, function(x) mean(x))
avg.table.frxn8.RNASE <- sapply(table.frxn8.RNASE, function(x) mean(x))
avg.table.frxn9.RNASE <- sapply(table.frxn9.RNASE, function(x) mean(x))
avg.table.frxn10.RNASE <- sapply(table.frxn10.RNASE, function(x) mean(x))
avg.table.frxn11.RNASE <- sapply(table.frxn11.RNASE, function(x) mean(x))
avg.table.frxn12.RNASE <- sapply(table.frxn12.RNASE, function(x) mean(x))
avg.table.frxn13.RNASE <- sapply(table.frxn13.RNASE, function(x) mean(x))
avg.table.frxn14.RNASE <- sapply(table.frxn14.RNASE, function(x) mean(x))
avg.table.frxn15.RNASE <- sapply(table.frxn15.RNASE, function(x) mean(x))
avg.table.frxn16.RNASE <- sapply(table.frxn16.RNASE, function(x) mean(x))
avg.table.frxn17.RNASE <- sapply(table.frxn17.RNASE, function(x) mean(x))
avg.table.frxn18.RNASE <- sapply(table.frxn18.RNASE, function(x) mean(x))
avg.table.frxn19.RNASE <- sapply(table.frxn19.RNASE, function(x) mean(x))
avg.table.frxn20.RNASE <- sapply(table.frxn20.RNASE, function(x) mean(x))
avg.table.frxn21.RNASE <- sapply(table.frxn21.RNASE, function(x) mean(x))
avg.table.frxn22.RNASE <- sapply(table.frxn22.RNASE, function(x) mean(x))
avg.table.frxn23.RNASE <- sapply(table.frxn23.RNASE, function(x) mean(x))
avg.table.frxn24.RNASE <- sapply(table.frxn24.RNASE, function(x) mean(x))
avg.table.frxn25.RNASE <- sapply(table.frxn25.RNASE, function(x) mean(x))


# Determine normalization factor for each condition, as the mean of the 2 most similar replicates
# Create a function norm_fact for this step:
norm_fact <- function(x) {
				if( (abs(x[1]-x[2])<abs(x[1]-x[3])) && (abs(x[1]-x[2])<abs(x[2]-x[3])) ) 
					{mean(c(x[1],x[2]))} else if( (abs(x[1]-x[3])<abs(x[1]-x[2])) && (abs(x[1]-x[3])<abs(x[2]-x[3])) )
												  {mean(c(x[1],x[3]))} else {mean(c(x[2],x[3]))} 
                        }


# Determine the normalization factor for the CTRL fractions and for the RNASE fractions
norm_mean_frxn1_CTRL <-  norm_fact(avg.table.frxn1.CTRL)/avg.table.frxn1.CTRL
norm_mean_frxn2_CTRL <-  norm_fact(avg.table.frxn2.CTRL)/avg.table.frxn2.CTRL
norm_mean_frxn3_CTRL <-  norm_fact(avg.table.frxn3.CTRL)/avg.table.frxn3.CTRL
norm_mean_frxn4_CTRL <-  norm_fact(avg.table.frxn4.CTRL)/avg.table.frxn4.CTRL
norm_mean_frxn5_CTRL <-  norm_fact(avg.table.frxn5.CTRL)/avg.table.frxn5.CTRL
norm_mean_frxn6_CTRL <-  norm_fact(avg.table.frxn6.CTRL)/avg.table.frxn6.CTRL
norm_mean_frxn7_CTRL <-  norm_fact(avg.table.frxn7.CTRL)/avg.table.frxn7.CTRL
norm_mean_frxn8_CTRL <-  norm_fact(avg.table.frxn8.CTRL)/avg.table.frxn8.CTRL
norm_mean_frxn9_CTRL <-  norm_fact(avg.table.frxn9.CTRL)/avg.table.frxn9.CTRL
norm_mean_frxn10_CTRL <-  norm_fact(avg.table.frxn10.CTRL)/avg.table.frxn10.CTRL
norm_mean_frxn11_CTRL <-  norm_fact(avg.table.frxn11.CTRL)/avg.table.frxn11.CTRL
norm_mean_frxn12_CTRL <-  norm_fact(avg.table.frxn12.CTRL)/avg.table.frxn12.CTRL
norm_mean_frxn13_CTRL <-  norm_fact(avg.table.frxn13.CTRL)/avg.table.frxn13.CTRL
norm_mean_frxn14_CTRL <-  norm_fact(avg.table.frxn14.CTRL)/avg.table.frxn14.CTRL
norm_mean_frxn15_CTRL <-  norm_fact(avg.table.frxn15.CTRL)/avg.table.frxn15.CTRL
norm_mean_frxn16_CTRL <-  norm_fact(avg.table.frxn16.CTRL)/avg.table.frxn16.CTRL
norm_mean_frxn17_CTRL <-  norm_fact(avg.table.frxn17.CTRL)/avg.table.frxn17.CTRL
norm_mean_frxn18_CTRL <-  norm_fact(avg.table.frxn18.CTRL)/avg.table.frxn18.CTRL
norm_mean_frxn19_CTRL <-  norm_fact(avg.table.frxn19.CTRL)/avg.table.frxn19.CTRL
norm_mean_frxn20_CTRL <-  norm_fact(avg.table.frxn20.CTRL)/avg.table.frxn20.CTRL
norm_mean_frxn21_CTRL <-  norm_fact(avg.table.frxn21.CTRL)/avg.table.frxn21.CTRL
norm_mean_frxn22_CTRL <-  norm_fact(avg.table.frxn22.CTRL)/avg.table.frxn22.CTRL
norm_mean_frxn23_CTRL <-  norm_fact(avg.table.frxn23.CTRL)/avg.table.frxn23.CTRL
norm_mean_frxn24_CTRL <-  norm_fact(avg.table.frxn24.CTRL)/avg.table.frxn24.CTRL
norm_mean_frxn25_CTRL <-  norm_fact(avg.table.frxn25.CTRL)/avg.table.frxn25.CTRL

norm_mean_frxn1_RNASE <-  norm_fact(avg.table.frxn1.RNASE)/avg.table.frxn1.RNASE
norm_mean_frxn2_RNASE <-  norm_fact(avg.table.frxn2.RNASE)/avg.table.frxn2.RNASE
norm_mean_frxn3_RNASE <-  norm_fact(avg.table.frxn3.RNASE)/avg.table.frxn3.RNASE
norm_mean_frxn4_RNASE <-  norm_fact(avg.table.frxn4.RNASE)/avg.table.frxn4.RNASE
norm_mean_frxn5_RNASE <-  norm_fact(avg.table.frxn5.RNASE)/avg.table.frxn5.RNASE
norm_mean_frxn6_RNASE <-  norm_fact(avg.table.frxn6.RNASE)/avg.table.frxn6.RNASE
norm_mean_frxn7_RNASE <-  norm_fact(avg.table.frxn7.RNASE)/avg.table.frxn7.RNASE
norm_mean_frxn8_RNASE <-  norm_fact(avg.table.frxn8.RNASE)/avg.table.frxn8.RNASE
norm_mean_frxn9_RNASE <-  norm_fact(avg.table.frxn9.RNASE)/avg.table.frxn9.RNASE
norm_mean_frxn10_RNASE <-  norm_fact(avg.table.frxn10.RNASE)/avg.table.frxn10.RNASE
norm_mean_frxn11_RNASE <-  norm_fact(avg.table.frxn11.RNASE)/avg.table.frxn11.RNASE
norm_mean_frxn12_RNASE <-  norm_fact(avg.table.frxn12.RNASE)/avg.table.frxn12.RNASE
norm_mean_frxn13_RNASE <-  norm_fact(avg.table.frxn13.RNASE)/avg.table.frxn13.RNASE
norm_mean_frxn14_RNASE <-  norm_fact(avg.table.frxn14.RNASE)/avg.table.frxn14.RNASE
norm_mean_frxn15_RNASE <-  norm_fact(avg.table.frxn15.RNASE)/avg.table.frxn15.RNASE
norm_mean_frxn16_RNASE <-  norm_fact(avg.table.frxn16.RNASE)/avg.table.frxn16.RNASE
norm_mean_frxn17_RNASE <-  norm_fact(avg.table.frxn17.RNASE)/avg.table.frxn17.RNASE
norm_mean_frxn18_RNASE <-  norm_fact(avg.table.frxn18.RNASE)/avg.table.frxn18.RNASE
norm_mean_frxn19_RNASE <-  norm_fact(avg.table.frxn19.RNASE)/avg.table.frxn19.RNASE
norm_mean_frxn20_RNASE <-  norm_fact(avg.table.frxn20.RNASE)/avg.table.frxn20.RNASE
norm_mean_frxn21_RNASE <-  norm_fact(avg.table.frxn21.RNASE)/avg.table.frxn21.RNASE
norm_mean_frxn22_RNASE <-  norm_fact(avg.table.frxn22.RNASE)/avg.table.frxn22.RNASE
norm_mean_frxn23_RNASE <-  norm_fact(avg.table.frxn23.RNASE)/avg.table.frxn23.RNASE
norm_mean_frxn24_RNASE <-  norm_fact(avg.table.frxn24.RNASE)/avg.table.frxn24.RNASE
norm_mean_frxn25_RNASE <-  norm_fact(avg.table.frxn25.RNASE)/avg.table.frxn25.RNASE


# Compute vectors containing the normalization factors for the respective treatment and replicates and fractions 1 to 25
norm.ctrl1 <- c(norm_mean_frxn1_CTRL[1],norm_mean_frxn2_CTRL[1],norm_mean_frxn3_CTRL[1],norm_mean_frxn4_CTRL[1],norm_mean_frxn5_CTRL[1],norm_mean_frxn6_CTRL[1],norm_mean_frxn7_CTRL[1],norm_mean_frxn8_CTRL[1],norm_mean_frxn9_CTRL[1],norm_mean_frxn10_CTRL[1],norm_mean_frxn11_CTRL[1],norm_mean_frxn12_CTRL[1],norm_mean_frxn13_CTRL[1],norm_mean_frxn14_CTRL[1],norm_mean_frxn15_CTRL[1],norm_mean_frxn16_CTRL[1],norm_mean_frxn17_CTRL[1],norm_mean_frxn18_CTRL[1],norm_mean_frxn19_CTRL[1],norm_mean_frxn20_CTRL[1],norm_mean_frxn21_CTRL[1],norm_mean_frxn22_CTRL[1],norm_mean_frxn23_CTRL[1],norm_mean_frxn24_CTRL[1],norm_mean_frxn25_CTRL[1])
norm.ctrl2 <- c(norm_mean_frxn1_CTRL[2],norm_mean_frxn2_CTRL[2],norm_mean_frxn3_CTRL[2],norm_mean_frxn4_CTRL[2],norm_mean_frxn5_CTRL[2],norm_mean_frxn6_CTRL[2],norm_mean_frxn7_CTRL[2],norm_mean_frxn8_CTRL[2],norm_mean_frxn9_CTRL[2],norm_mean_frxn10_CTRL[2],norm_mean_frxn11_CTRL[2],norm_mean_frxn12_CTRL[2],norm_mean_frxn13_CTRL[2],norm_mean_frxn14_CTRL[2],norm_mean_frxn15_CTRL[2],norm_mean_frxn16_CTRL[2],norm_mean_frxn17_CTRL[2],norm_mean_frxn18_CTRL[2],norm_mean_frxn19_CTRL[2],norm_mean_frxn20_CTRL[2],norm_mean_frxn21_CTRL[2],norm_mean_frxn22_CTRL[2],norm_mean_frxn23_CTRL[2],norm_mean_frxn24_CTRL[2],norm_mean_frxn25_CTRL[2])
norm.ctrl3 <- c(norm_mean_frxn1_CTRL[3],norm_mean_frxn2_CTRL[3],norm_mean_frxn3_CTRL[3],norm_mean_frxn4_CTRL[3],norm_mean_frxn5_CTRL[3],norm_mean_frxn6_CTRL[3],norm_mean_frxn7_CTRL[3],norm_mean_frxn8_CTRL[3],norm_mean_frxn9_CTRL[3],norm_mean_frxn10_CTRL[3],norm_mean_frxn11_CTRL[3],norm_mean_frxn12_CTRL[3],norm_mean_frxn13_CTRL[3],norm_mean_frxn14_CTRL[3],norm_mean_frxn15_CTRL[3],norm_mean_frxn16_CTRL[3],norm_mean_frxn17_CTRL[3],norm_mean_frxn18_CTRL[3],norm_mean_frxn19_CTRL[3],norm_mean_frxn20_CTRL[3],norm_mean_frxn21_CTRL[3],norm_mean_frxn22_CTRL[3],norm_mean_frxn23_CTRL[3],norm_mean_frxn24_CTRL[3],norm_mean_frxn25_CTRL[3])

norm.rnase1 <- c(norm_mean_frxn1_RNASE[1],norm_mean_frxn2_RNASE[1],norm_mean_frxn3_RNASE[1],norm_mean_frxn4_RNASE[1],norm_mean_frxn5_RNASE[1],norm_mean_frxn6_RNASE[1],norm_mean_frxn7_RNASE[1],norm_mean_frxn8_RNASE[1],norm_mean_frxn9_RNASE[1],norm_mean_frxn10_RNASE[1],norm_mean_frxn11_RNASE[1],norm_mean_frxn12_RNASE[1],norm_mean_frxn13_RNASE[1],norm_mean_frxn14_RNASE[1],norm_mean_frxn15_RNASE[1],norm_mean_frxn16_RNASE[1],norm_mean_frxn17_RNASE[1],norm_mean_frxn18_RNASE[1],norm_mean_frxn19_RNASE[1],norm_mean_frxn20_RNASE[1],norm_mean_frxn21_RNASE[1],norm_mean_frxn22_RNASE[1],norm_mean_frxn23_RNASE[1],norm_mean_frxn24_RNASE[1],norm_mean_frxn25_RNASE[1])
norm.rnase2 <- c(norm_mean_frxn1_RNASE[2],norm_mean_frxn2_RNASE[2],norm_mean_frxn3_RNASE[2],norm_mean_frxn4_RNASE[2],norm_mean_frxn5_RNASE[2],norm_mean_frxn6_RNASE[2],norm_mean_frxn7_RNASE[2],norm_mean_frxn8_RNASE[2],norm_mean_frxn9_RNASE[2],norm_mean_frxn10_RNASE[2],norm_mean_frxn11_RNASE[2],norm_mean_frxn12_RNASE[2],norm_mean_frxn13_RNASE[2],norm_mean_frxn14_RNASE[2],norm_mean_frxn15_RNASE[2],norm_mean_frxn16_RNASE[2],norm_mean_frxn17_RNASE[2],norm_mean_frxn18_RNASE[2],norm_mean_frxn19_RNASE[2],norm_mean_frxn20_RNASE[2],norm_mean_frxn21_RNASE[2],norm_mean_frxn22_RNASE[2],norm_mean_frxn23_RNASE[2],norm_mean_frxn24_RNASE[2],norm_mean_frxn25_RNASE[2])
norm.rnase3 <- c(norm_mean_frxn1_RNASE[3],norm_mean_frxn2_RNASE[3],norm_mean_frxn3_RNASE[3],norm_mean_frxn4_RNASE[3],norm_mean_frxn5_RNASE[3],norm_mean_frxn6_RNASE[3],norm_mean_frxn7_RNASE[3],norm_mean_frxn8_RNASE[3],norm_mean_frxn9_RNASE[3],norm_mean_frxn10_RNASE[3],norm_mean_frxn11_RNASE[3],norm_mean_frxn12_RNASE[3],norm_mean_frxn13_RNASE[3],norm_mean_frxn14_RNASE[3],norm_mean_frxn15_RNASE[3],norm_mean_frxn16_RNASE[3],norm_mean_frxn17_RNASE[3],norm_mean_frxn18_RNASE[3],norm_mean_frxn19_RNASE[3],norm_mean_frxn20_RNASE[3],norm_mean_frxn21_RNASE[3],norm_mean_frxn22_RNASE[3],norm_mean_frxn23_RNASE[3],norm_mean_frxn24_RNASE[3],norm_mean_frxn25_RNASE[3])


# Define subtables for each treatment and each replicate 
data.ctrl1 <- data$condition =="ctrl1"
data.ctrl2 <- data$condition =="ctrl2"
data.ctrl3 <- data$condition =="ctrl3"
data.rnase1 <- data$condition =="rnase1"
data.rnase2 <- data$condition =="rnase2"
data.rnase3 <- data$condition =="rnase3"


#**************************************************************************************************
# Normalization step, fraction-wise
table.ctrl1 <- data.frame(mapply('*', table[,data.ctrl1], norm.ctrl1, SIMPLIFY=FALSE))
table.ctrl2 <- data.frame(mapply('*', table[,data.ctrl2], norm.ctrl2, SIMPLIFY=FALSE))
table.ctrl3 <- data.frame(mapply('*', table[,data.ctrl3], norm.ctrl3, SIMPLIFY=FALSE))

table.rnase1 <- data.frame(mapply('*', table[,data.rnase1], norm.rnase1, SIMPLIFY=FALSE))
table.rnase2 <- data.frame(mapply('*', table[,data.rnase2], norm.rnase2, SIMPLIFY=FALSE))
table.rnase3 <- data.frame(mapply('*', table[,data.rnase3], norm.rnase3, SIMPLIFY=FALSE))


# Get the proper rownames for the tables
rownames(table.ctrl1) <- row_names
rownames(table.ctrl2) <- row_names
rownames(table.ctrl3) <- row_names
rownames(table.rnase1) <- row_names
rownames(table.rnase2) <- row_names
rownames(table.rnase3) <- row_names


#**************************************************************************************************
#### Apply a sliding window/moving average of 3 points to the data ####
table.ctrl1.SW <- data.frame(table.ctrl1[1],(table.ctrl1[1:23]+table.ctrl1[2:24]+table.ctrl1[3:25])/3,table.ctrl1[25])
table.ctrl2.SW <- data.frame(table.ctrl2[1],(table.ctrl2[1:23]+table.ctrl2[2:24]+table.ctrl2[3:25])/3,table.ctrl2[25])
table.ctrl3.SW <- data.frame(table.ctrl3[1],(table.ctrl3[1:23]+table.ctrl3[2:24]+table.ctrl3[3:25])/3,table.ctrl3[25])
table.rnase1.SW <- data.frame(table.rnase1[1],(table.rnase1[1:23]+table.rnase1[2:24]+table.rnase1[3:25])/3,table.rnase1[25])
table.rnase2.SW <- data.frame(table.rnase2[1],(table.rnase2[1:23]+table.rnase2[2:24]+table.rnase2[3:25])/3,table.rnase2[25])
table.rnase3.SW <- data.frame(table.rnase3[1],(table.rnase3[1:23]+table.rnase3[2:24]+table.rnase3[3:25])/3,table.rnase3[25])


# Get the proper rownames for the tables
colnames(table.ctrl1.SW) <- colnames(table.ctrl1)
colnames(table.ctrl2.SW) <- colnames(table.ctrl2)
colnames(table.ctrl3.SW) <- colnames(table.ctrl3)
colnames(table.rnase1.SW) <- colnames(table.rnase1)
colnames(table.rnase2.SW) <- colnames(table.rnase2)
colnames(table.rnase3.SW) <- colnames(table.rnase3)


#**************************************************************************************************
#### Normaliztion over the amount to sum = 100 (%) ####

table.ctrl1.SW.norm <- table.ctrl1.SW * 100 / rowSums(table.ctrl1.SW)
table.ctrl2.SW.norm <- table.ctrl2.SW * 100 / rowSums(table.ctrl2.SW)
table.ctrl3.SW.norm <- table.ctrl3.SW * 100 / rowSums(table.ctrl3.SW)

table.rnase1.SW.norm <- table.rnase1.SW * 100 / rowSums(table.rnase1.SW)
table.rnase2.SW.norm <- table.rnase2.SW * 100 / rowSums(table.rnase2.SW)
table.rnase3.SW.norm <- table.rnase3.SW * 100 / rowSums(table.rnase3.SW)


#**************************************************************************************************
# Replace NA, NaN with 0 if any
table.ctrl1.SW.norm <- rapply(table.ctrl1.SW.norm, f=function(x) ifelse(is.na(x),0,x), how="replace") 
table.ctrl2.SW.norm <- rapply(table.ctrl2.SW.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.ctrl3.SW.norm <- rapply(table.ctrl3.SW.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase1.SW.norm <- rapply(table.rnase1.SW.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase2.SW.norm <- rapply(table.rnase2.SW.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase3.SW.norm <- rapply(table.rnase3.SW.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )

table.ctrl1.SW.norm <- rapply(table.ctrl1.SW.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace") 
table.ctrl2.SW.norm <- rapply(table.ctrl2.SW.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.ctrl3.SW.norm <- rapply(table.ctrl3.SW.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase1.SW.norm <- rapply(table.rnase1.SW.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase2.SW.norm <- rapply(table.rnase2.SW.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase3.SW.norm <- rapply(table.rnase3.SW.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )


#**************************************************************************************************
# Average value for normalization of CTRL and RNASE samples to 100
my.list.ctrl.norm <- list(table.ctrl1.SW, table.ctrl2.SW, table.ctrl3.SW)
my.list.rnase.norm <- list(table.rnase1.SW, table.rnase2.SW, table.rnase3.SW)

ctrl_norm_mean <- Reduce("+", my.list.ctrl.norm)/length(my.list.ctrl.norm)
rnase_norm_mean <- Reduce("+", my.list.rnase.norm)/length(my.list.rnase.norm)


# Change names of the columns: from "fraction1" to "fraction25"
col_fractions <- paste("fraction",1:25,sep="")
colnames(ctrl_norm_mean) <- col_fractions
colnames(rnase_norm_mean) <- col_fractions


# Normalization to 100 for each protein (= the sum of the amount of protein from fraction1 to fraction25 is 100)
ctrl_norm_mean <- ctrl_norm_mean*100/rowSums(ctrl_norm_mean)
rnase_norm_mean <- rnase_norm_mean*100/rowSums(rnase_norm_mean)


# Replace the NaN and NA values by 0 if any
ctrl_norm_mean <- rapply(ctrl_norm_mean, f=function(x) ifelse(is.nan(x),0,x), how="replace")
rnase_norm_mean <- rapply(rnase_norm_mean, f=function(x) ifelse(is.nan(x),0,x), how="replace")
ctrl_norm_mean <- rapply(ctrl_norm_mean, f=function(x) ifelse(is.na(x),0,x), how="replace")
rnase_norm_mean <- rapply(rnase_norm_mean, f=function(x) ifelse(is.na(x),0,x), how="replace")


# If one of the curve is 0 over all the fractions of one condition, the other condition will be set to 0 to take it out of the analysis
ctrl_norm_mean[rowSums(rnase_norm_mean[1:25])==0,] <- 0
rnase_norm_mean[rowSums(ctrl_norm_mean[1:25])==0,] <- 0


#*********************************************************************************************************
#### PART 2 : Finding the fit parameters from average curves (every local maximum above a value of 2) ####
#*********************************************************************************************************

#*********************************************************************************************************
#### Find local maxima in each protein mean profile (peak calling step - function find_peaks) ####
# Restriction to values above an absolut threshold of 2 (to get rid of the noise)
# Identify also shoulders not found by peak calling
#*********************************************************************************************************

# Define the function find_peaks
# A "peak" is defined as a local maxima with m points either side of it being smaller than it.
# Hence, the bigger the parameter m, the more stringent the peak finding procedure
# m set to 2

find_peaks <- function (x, m = 2){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
     })
   pks <- unlist(pks)
     
   # part added to deal with plateaux
     n = length(pks)
     if (n>1) {
     			rem <- numeric(0)
     			for (i in 1:(n-1)) { if  ((pks[i]+1) == pks[(i+1)]) {rem <- c(rem, pks[i+1]) } }
     			for (i in 1:(n-1)) { if  ((pks[i]+2) == pks[(i+1)]) {rem <- c(rem, pks[i+1]) } }
     			pks <- pks[! pks %in% rem]
     		  }		   
     
   # part added to deal with the 1st and 25th values, in case they are max
     if ( sum(x[1]>x[2:(m+1)]) == 2 ) {pks <- c(pks, 1)} 
     if ( sum(x[25]>x[24:(25-m)]) == 2 ) {pks <- c(pks, 25)}
     pks <- unlist(pks)
     pks
}


#*********************************************************************************************************
# Apply the function to the data set and retrieve the values
# restriction to values above an absolut threshold of 2%
# New column: "maxima" (column 26)

ctrl_norm_mean$maxima <- apply(ctrl_norm_mean, 1, function(x) {
															   list <- find_peaks(x)
															   list <- list[x[list] > 2] 
															   list <- unlist(list)
															   list
															   })
															   
rnase_norm_mean$maxima <- apply(rnase_norm_mean, 1, function(x) {
															   list <- find_peaks(x)
															   list <- list[x[list] > 2] 
															   list <- unlist(list)
															   list
															   })


#**************************************************************************************************
# Step to retrieve the so-called "shoulder" - Group of fraction where there is a larger amount of protein that was not identified as "maxima" using the function "find_peaks"


# Get fractions where (RNASE fraction value > 2%). Dataframe with 0/1 and RNASE maxima
# Define new tables called ctrl_3 and rnase_3 that will be used temporary
rnase_3 <- rnase_norm_mean
rnase_3[1:25] <- (rnase_norm_mean[1:25] > 2)*1


# Get fractions where (CTRL fraction value > 2%). Dataframe with 0/1 and CTRL maxima
ctrl_3 <- ctrl_norm_mean
ctrl_3[1:25] <- (ctrl_norm_mean[1:25] > 2)*1


# Remove regions around maxima (3 fractions) in order not to identify those regions as "shoulder"
th_max_reg <- function(x) {
							ls <- as.numeric(unlist(x$maxima))
							n <- length(ls)
							dt <- as.data.frame(matrix(rep(1,25),1,25))
							if (n == 0) {dt[1:25] = 1} else {
							for (i in 1:n) { 
											if (ls[i] < 1) {dt[1:25] = 1}
											else if (ls[i] == 1) {dt[1:4] = 0}
											else if (ls[i] == 2) {dt[1:5] = 0}
											else if (ls[i] == 3) {dt[1:6] = 0}
											else if (ls[i] == 4) {dt[1:7] = 0}
											else if (ls[i] == 5) {dt[1:8] = 0}
											
											else if (ls[i]>=6  && ls[i]<=20) {dt[(ls[i]-3):(ls[i]+3)] = 0}
											
											else if (ls[i] == 21) {dt[18:24] = 0}
											else if (ls[i] == 22) {dt[19:25] = 0}
											else if (ls[i] == 23) {dt[20:25] = 0}
											else if (ls[i] == 24) {dt[21:25] = 0}
											else {dt[22:25] = 0}
											}
														   }	
							x[1:25] <- x[1:25]*dt
							x <- unlist(x)
							x[1:25]
						  }

rnase_3[1:25] <- t(apply( rnase_3, 1, function(x) {th_max_reg(x)} ))
ctrl_3[1:25] <- t(apply( ctrl_3, 1, function(x) {th_max_reg(x)} ))

# Get the regions (marked with 1, at least 4 consecutive fractions) and select the middle of it
# Use the function rle -> $lengths and $values
# Apply on the whole dataframe

peaks_reg <- function(x) {
							res <- rle(as.numeric(x[1:25]))
							pos_res <- which(res$lengths>=4 & res$values==1)
							if (length(pos_res) == 0) {
														vct <- numeric(0)
														}
							else {							
								vct <- numeric(0)
								for (i in 1:length(pos_res)) { vct <- c( vct,   sum(res$lengths[1:pos_res[i]]) - as.integer(res$lengths[pos_res[1]]/2) ) }
								  }
								  
							vct <- unlist(vct)
							vct
						 }
						 
						 
# get the "shoulders" in the curves
# New column: "peaks" (column 27)

rnase_3$peaks <- apply(rnase_3, 1, function(x) { peaks_reg(x) } )
ctrl_3$peaks <- apply(ctrl_3, 1, function(x) { peaks_reg(x) } )


#**************************************************************************************************
# Get the list of all maxima (maxima and peaks) for RNASE and CTRL
# New column: "ctrl_max" and "rnase_max" (column 28) respectively

						   						
rnase_3$rnase_max <- apply(rnase_3, 1, function(x) {
												 ls_max <- as.numeric(unlist(x$maxima))
												 ls_peaks <- as.numeric(unlist(x$peaks))
												 
												 rnase_max <- c(ls_max, ls_peaks)
												 rnase_max <- unlist(rnase_max)
												 rnase_max <- rnase_max[rnase_max!=0]
												 rnase_max <- sort(rnase_max, decreasing = FALSE)
												 if (length(rnase_max) == 0) {0} else {rnase_max}
												 })

ctrl_3$ctrl_max <- apply(ctrl_3, 1, function(x) {
												 ls_max <- as.numeric(unlist(x$maxima))
												 ls_peaks <- as.numeric(unlist(x$peaks))
												 
												 ctrl_max <- c(ls_max, ls_peaks)
												 ctrl_max <- unlist(ctrl_max)
												 ctrl_max <- ctrl_max[ctrl_max!=0]
												 ctrl_max <- sort(ctrl_max, decreasing = FALSE)
												 if (length(ctrl_max) == 0) {0} else {ctrl_max}
												 })												 


#**************************************************************************************************
# Get the number of maxima for RNASE and CTRL
# New column: "nb_max" (column 29)

rnase_3$nb_max <- apply(rnase_3, 1, function(x) {
												 ls_max <- as.numeric(unlist(x$rnase_max))
												 n = length(ls_max)
												 if (sum(ls_max) == 0) {0} else {n}
												 })
												 
ctrl_3$nb_max <- apply(ctrl_3, 1, function(x) {
												 ls_max <- as.numeric(unlist(x$ctrl_max))
												 n = length(ls_max)
												 if (sum(ls_max) == 0) {0} else {n}
												 })
										 

#**************************************************************************************************
#### Gaussian fit of the mean curves (mean of the three replicates per condition)              ####
#************************************************************************************************** 

# Addition of new columns in the table ctrl_3 and rnase_3
# Use of the information stored in the table ctrl_norm_mean and rnase_norm_mean
# Function to fit a gaussian on the raw mean data to obtain: fit_c, fit_mean, fit_sigma
# Function to get the residue value: fit_res
# Function to check if the fit was successful: fitted (TRUE/FALSE)

# Number of rows of the table - is the total number of proteins
lg <- dim(rnase_3)[1]
# Definition of a matrix with numbers ranging from 1 to lg (matrix with 1 row and lg columns)
vect <- matrix(c(1:lg), 1, lg)

# Functions to optimize for Gaussian fitting - Gaussians with 1 to 6 peaks
# Each parameter C, mean and sigma will be optimized for each CTRL and RNASE profile
f1 <- function(data, q) {
					C1 <- q[1]
					mean1 <- q[2]
					sigma1 <- q[3]
					res <- (C1 * exp(-(data$y-mean1)**2/(2 * sigma1**2))) - data$df.y
					sum(res * res)
				 }

f2 <- function(data, q) {
					C1 <- q[1]
					mean1 <- q[2]
					sigma1 <- q[3]
					C2 <- q[4]
					mean2 <- q[5]
					sigma2 <- q[6]
					res <- ( C1 * exp(-(data$y-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(data$y-mean2)**2/(2 * sigma2**2)) ) - data$df.y
					sum(res * res)
				 }

f3 <- function(data, q) {
					C1 <- q[1]
					mean1 <- q[2]
					sigma1 <- q[3]
					C2 <- q[4]
					mean2 <- q[5]
					sigma2 <- q[6]
					C3 <- q[7]
					mean3 <- q[8]
					sigma3 <- q[9]
					res <- ( C1 * exp(-(data$y-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(data$y-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(data$y-mean3)**2/(2 * sigma3**2)) ) - data$df.y
					sum(res * res)
				 }

f4 <- function(data, q) {
					C1 <- q[1]
					mean1 <- q[2]
					sigma1 <- q[3]
					C2 <- q[4]
					mean2 <- q[5]
					sigma2 <- q[6]
					C3 <- q[7]
					mean3 <- q[8]
					sigma3 <- q[9]
					C4 <- q[10]
					mean4 <- q[11]
					sigma4 <- q[12]
					res <- ( C1 * exp(-(data$y-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(data$y-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(data$y-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(data$y-mean4)**2/(2 * sigma4**2)) ) - data$df.y
					sum(res * res)
				 }

f5 <- function(data, q) {
					C1 <- q[1]
					mean1 <- q[2]
					sigma1 <- q[3]
					C2 <- q[4]
					mean2 <- q[5]
					sigma2 <- q[6]
					C3 <- q[7]
					mean3 <- q[8]
					sigma3 <- q[9]
					C4 <- q[10]
					mean4 <- q[11]
					sigma4 <- q[12]
					C5 <- q[13]
					mean5 <- q[14]
					sigma5 <- q[15]
					res <- ( C1 * exp(-(data$y-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(data$y-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(data$y-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(data$y-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(data$y-mean5)**2/(2 * sigma5**2)) ) - data$df.y
					sum(res * res)
				 }

f6 <- function(data, q) {
					C1 <- q[1]
					mean1 <- q[2]
					sigma1 <- q[3]
					C2 <- q[4]
					mean2 <- q[5]
					sigma2 <- q[6]
					C3 <- q[7]
					mean3 <- q[8]
					sigma3 <- q[9]
					C4 <- q[10]
					mean4 <- q[11]
					sigma4 <- q[12]
					C5 <- q[13]
					mean5 <- q[14]
					sigma5 <- q[15]
					C6 <- q[16]
					mean6 <- q[17]
					sigma6 <- q[18]
					res <- ( C1 * exp(-(data$y-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(data$y-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(data$y-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(data$y-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(data$y-mean5)**2/(2 * sigma5**2)) + C6 * exp(-(data$y-mean6)**2/(2 * sigma6**2)) ) - data$df.y
					sum(res * res)
				 }												  

#**************************************************************************************************
# Fit the RNASE mean curves and get the amplitude of the fit at each peak - new column: fit_c
rnase_3$fit_c <- apply(vect, 2, function(t) {
                          df.y <- c(as.numeric(rnase_norm_mean[t,1:25]))
                          y <- c(1:25)
                          data <- data.frame(y = y, df.y = df.y)
                          Mean_rnase_list <- as.numeric(unlist(rnase_3[t,"rnase_max"]))
                          n <- length(Mean_rnase_list)
                          gauss <- numeric(0)
                          # Before optimization of the parameters, check the number of maxima for the profile of the protein
												  if (n==0) { gauss <- 0 } else {
												  				  								if (n==1) {
												  		  	         									  vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]], Mean_rnase_list[1], 1 )
												  		  	         									  gauss <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$par[1],digits = 1 ),error=function(e) numeric(0))
												  				           								  } else {
												  				  												  if (n==2) {
												  				                 											 vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1 )
												  				                 											 gauss <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4)],digits = 1 ),error=function(e) numeric(0))
												  	  		                     											 } else {
												  				  																	if (n==3) {
												  				  																			   vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1 )
												  				  																			   gauss <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4,7)],digits = 1 ),error=function(e) numeric(0))
												  				 																			   } else {
												  				  																					  if (n==4) {
												  				  																								 vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1 )
												  				  																								 gauss <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4,7,10)],digits = 1 ),error=function(e) numeric(0))
												  				  																								 } else {
												                  																										 if (n==5) {
												                  																													vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1,rnase_norm_mean[t,Mean_rnase_list[5]],Mean_rnase_list[5],1 )
												                  																													gauss <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4,7,10,13)],digits = 1 ),error=function(e) numeric(0))
												                																													} else {
												                   																														    vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1,rnase_norm_mean[t,Mean_rnase_list[5]],Mean_rnase_list[5],1,rnase_norm_mean[t,Mean_rnase_list[6]],Mean_rnase_list[6],1 )
												                   																														    gauss <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4,7,10,13,16)],digits = 1 ),error=function(e) numeric(0))
												                  																															}
												  				  																										 }								
												  				  																						}
												  				  																		}
												  				  												   }
												  				  								}
												  gauss <- unlist(gauss)
												  if (length(gauss) == n) {as.numeric(gauss)} else {0}
												  })

#**************************************************************************************************
# Fit the RNASE mean curves and get the position of the peak(s) - new column: fit_mean
rnase_3$fit_mean <- apply(vect, 2, function(t) {
												  df.y <- c(as.numeric(rnase_norm_mean[t,1:25]))
												  y <- c(1:25)
												  data <- data.frame(y = y, df.y = df.y)
												  Mean_rnase_list <- as.numeric(unlist(rnase_3[t,"rnase_max"]))
												  n <- length(Mean_rnase_list)
												  gauss <- numeric(0)
										      if (n==0) { gauss <- 0 } else {
												  				  								if (n==1) {
												  		  	         									  vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]], Mean_rnase_list[1], 2 )
												  		  	         									  gauss <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$par[2],digits = 1 ),error=function(e) numeric(0))
												  				           								  } else {
												  				  												  if (n==2) {
												  				                 											 vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1 )
												  				                 											 gauss <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5)],digits = 1 ),error=function(e) numeric(0))
												  	  		                     											 } else {
												  				  																	if (n==3) {
												  				  																			   vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1 )
												  				  																			   gauss <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5,8)],digits = 1 ),error=function(e) numeric(0))
												  				 																			   } else {
												  				  																					  if (n==4) {
												  				  																								 vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1 )
												  				  																								 gauss <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5,8,11)],digits = 1 ),error=function(e) numeric(0))
												  				  																								 } else {
												                  																										 if (n==5) {
												                  																													vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1,rnase_norm_mean[t,Mean_rnase_list[5]],Mean_rnase_list[5],1 )
												                  																													gauss <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5,8,11,14)],digits = 1 ),error=function(e) numeric(0))
												                																													} else {
												                   																														    vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1,rnase_norm_mean[t,Mean_rnase_list[5]],Mean_rnase_list[5],1,rnase_norm_mean[t,Mean_rnase_list[6]],Mean_rnase_list[6],1 )
												                   																														    gauss <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5,8,11,14,17)],digits = 1 ),error=function(e) numeric(0))
												                  																															}
												  				  																										 }								
												  				  																						}
												  				  																		}
												  				  												   }
												  				  								}
												  gauss <- unlist(gauss)
												  if (length(gauss) == n) {as.numeric(gauss)} else {0}
												  })

#**************************************************************************************************
# Fit the RNASE mean curves and get the covariance (sigma) of the peak(s) - new column: fit_sigma
rnase_3$fit_sigma <- apply(vect, 2, function(t) {
												  df.y <- c(as.numeric(rnase_norm_mean[t,1:25]))
												  y <- c(1:25)
												  data <- data.frame(y = y, df.y = df.y)
												  Mean_rnase_list <- as.numeric(unlist(rnase_3[t,"rnase_max"]))
												  n <- length(Mean_rnase_list)
												  gauss <- numeric(0)
												  if (n==0) { gauss <- 0 } else {
												  				  								if (n==1) {
												  		  	         									  vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]], Mean_rnase_list[1], 2 )
												  		  	         									  gauss <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$par[3],digits = 1 ),error=function(e) numeric(0))
												  				           								  } else {
												  				  												  if (n==2) {
												  				                 											 vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1 )
												  				                 											 gauss <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6)],digits = 1 ),error=function(e) numeric(0))
												  	  		                     											 } else {
												  				  																	if (n==3) {
												  				  																			   vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1 )
												  				  																			   gauss <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6,9)],digits = 1 ),error=function(e) numeric(0))
												  				 																			   } else {
												  				  																					  if (n==4) {
												  				  																								 vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1 )
												  				  																								 gauss <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6,9,12)],digits = 1 ),error=function(e) numeric(0))
												  				  																								 } else {
												                  																										 if (n==5) {
												                  																													vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1,rnase_norm_mean[t,Mean_rnase_list[5]],Mean_rnase_list[5],1 )
												                  																													gauss <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6,9,12,15)],digits = 1 ),error=function(e) numeric(0))
												                																													} else {
												                   																														    vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1,rnase_norm_mean[t,Mean_rnase_list[5]],Mean_rnase_list[5],1,rnase_norm_mean[t,Mean_rnase_list[6]],Mean_rnase_list[6],1 )
												                   																														    gauss <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6,9,12,15,18)],digits = 1 ),error=function(e) numeric(0))
												                  																															}
												  				  																										 }								
												  				  																						}
												  				  																		}
												  				  												   }
												  				  								}
												  gauss <- unlist(gauss)
												  if (length(gauss) == n) {as.numeric(gauss)} else {0}
												  })

#**************************************************************************************************
# Fit the RNASE mean curves and get the residual sum of the squares - new column: fit_res
rnase_3$fit_res <- apply(vect, 2, function(t) {
												  df.y <- c(as.numeric(rnase_norm_mean[t,1:25]))
												  y <- c(1:25)
												  data <- data.frame(y = y, df.y = df.y)
												  Mean_rnase_list <- as.numeric(unlist(rnase_3[t,"rnase_max"]))
												  n <- length(Mean_rnase_list)
												  gauss <- numeric(0)
                          if (n==0) { gauss <- 0 } else {
												  				  								if (n==1) {
												  		  	         									  vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]], Mean_rnase_list[1], 1 )
												  		  	         									  gauss <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												  				           								  } else {
												  				  												  if (n==2) {
												  				                 											 vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1 )
												  				                 											 gauss <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												  	  		                     											 } else {
												  				  																	if (n==3) {
												  				  																			   vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1 )
												  				  																			   gauss <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												  				 																			   } else {
												  				  																					  if (n==4) {
												  				  																								 vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1 )
												  				  																								 gauss <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												  				  																								 } else {
												                  																										 if (n==5) {
												                  																													vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1,rnase_norm_mean[t,Mean_rnase_list[5]],Mean_rnase_list[5],1 )
												                  																													gauss <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												                																													} else {
												                   																														    vector <- c( rnase_norm_mean[t,Mean_rnase_list[1]],Mean_rnase_list[1],1,rnase_norm_mean[t,Mean_rnase_list[2]],Mean_rnase_list[2],1,rnase_norm_mean[t,Mean_rnase_list[3]],Mean_rnase_list[3],1,rnase_norm_mean[t,Mean_rnase_list[4]],Mean_rnase_list[4],1,rnase_norm_mean[t,Mean_rnase_list[5]],Mean_rnase_list[5],1,rnase_norm_mean[t,Mean_rnase_list[6]],Mean_rnase_list[6],1 )
												                   																														    gauss <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												                  																															}
												  				  																										 }								
												  				  																						}
												  				  																		}
												  				  												   }
												  				  								}
												  gauss <- unlist(gauss)
												  if (length(gauss) == 1) {as.numeric(gauss)} else {0}
												  })

#**************************************************************************************************
# Check whether a fit curve could be found or not for the protein - new column: fitted
rnase_3$fitted <- apply(rnase_3, 1, function(x) { 
																				 list_c <- as.numeric(unlist(x$fit_c))
																				 list_m <- as.numeric(unlist(x$fit_mean))
																				 list_s <- as.numeric(unlist(x$fit_sigma))
																				 list <- c(list_c, list_m, list_s)
																				 
																				 if (sum(list) == 0) {FALSE} else {TRUE}
																				 })												  


#**************************************************************************************************
# Repeat the same with the CTRL condition

#**************************************************************************************************
# Fit the CTRL mean curves and get the amplitude of the fit at each peak - new column: fit_c
ctrl_3$fit_c <- apply(vect, 2, function(t) {
  
  df.y <- c(as.numeric(rnase_norm_mean[t,1:25]))
                          df.y <- c(as.numeric(ctrl_norm_mean[t,1:25]))
												  y <- c(1:25)
												  data <- data.frame(y = y, df.y = df.y)
												  Mean_ctrl_list <- as.numeric(unlist(ctrl_3[t,"ctrl_max"]))
												  n <- length(Mean_ctrl_list)
												  gauss <- numeric(0)
												  if (n==0) { gauss <- 0 } else {
												  				  								if (n==1) {
												  		  	         									  vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]], Mean_ctrl_list[1], 1 )
												  		  	         									  gauss <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$par[1],digits = 1 ),error=function(e) numeric(0))
												  				           								  } else {
												  				  												  if (n==2) {
												  				                 											 vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1 )
												  				                 											 gauss <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4)],digits = 1 ),error=function(e) numeric(0))
												  	  		                     											 } else {
												  				  																	if (n==3) {
												  				  																			   vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1 )
												  				  																			   gauss <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4,7)],digits = 1 ),error=function(e) numeric(0))
												  				 																			   } else {
												  				  																					  if (n==4) {
												  				  																								 vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1 )
												  				  																								 gauss <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4,7,10)],digits = 1 ),error=function(e) numeric(0))
												  				  																								 } else {
												                  																										 if (n==5) {
												                  																													vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1,ctrl_norm_mean[t,Mean_ctrl_list[5]],Mean_ctrl_list[5],1 )
												                  																													gauss <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4,7,10,13)],digits = 1 ),error=function(e) numeric(0))
												                																													} else {
												                   																														    vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1,ctrl_norm_mean[t,Mean_ctrl_list[5]],Mean_ctrl_list[5],1,ctrl_norm_mean[t,Mean_ctrl_list[6]],Mean_ctrl_list[6],1 )
												                   																														    gauss <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(1,4,7,10,13,16)],digits = 1 ),error=function(e) numeric(0))
												                  																															}
												  				  																										 }								
												  				  																						}
												  				  																		}
												  				  												   }
												  				  								}
												  gauss <- unlist(gauss)
												  if (length(gauss) == n) {as.numeric(gauss)} else {0}
												  })

#**************************************************************************************************
# Fit the CTRL mean curves and get the position of the peak(s) - new column: fit_mean
ctrl_3$fit_mean <- apply(vect, 2, function(t) {
												  df.y <- c(as.numeric(ctrl_norm_mean[t,1:25]))
												  y <- c(1:25)
												  data <- data.frame(y = y, df.y = df.y)
												  Mean_ctrl_list <- as.numeric(unlist(ctrl_3[t,"ctrl_max"]))
												  n <- length(Mean_ctrl_list)
												  gauss <- numeric(0)
												  if (n==0) { gauss <- 0 } else {
												  				  								if (n==1) {
												  		  	         									  vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]], Mean_ctrl_list[1], 1 )
												  		  	         									  gauss <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$par[2],digits = 1 ),error=function(e) numeric(0))
												  				           								  } else {
												  				  												  if (n==2) {
												  				                 											 vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1 )
												  				                 											 gauss <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5)],digits = 1 ),error=function(e) numeric(0))
												  	  		                     											 } else {
												  				  																	if (n==3) {
												  				  																			   vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1 )
												  				  																			   gauss <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5,8)],digits = 1 ),error=function(e) numeric(0))
												  				 																			   } else {
												  				  																					  if (n==4) {
												  				  																								 vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1 )
												  				  																								 gauss <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5,8,11)],digits = 1 ),error=function(e) numeric(0))
												  				  																								 } else {
												                  																										 if (n==5) {
												                  																													vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1,ctrl_norm_mean[t,Mean_ctrl_list[5]],Mean_ctrl_list[5],1 )
												                  																													gauss <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5,8,11,14)],digits = 1 ),error=function(e) numeric(0))
												                																													} else {
												                   																														    vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1,ctrl_norm_mean[t,Mean_ctrl_list[5]],Mean_ctrl_list[5],1,ctrl_norm_mean[t,Mean_ctrl_list[6]],Mean_ctrl_list[6],1 )
												                   																														    gauss <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(2,5,8,11,14,17)],digits = 1 ),error=function(e) numeric(0))
												                  																															}
												  				  																										 }								
												  				  																						}
												  				  																		}
												  				  												   }
												  				  								}
												  gauss <- unlist(gauss)
												  if (length(gauss) == n) {as.numeric(gauss)} else {0}
												  })

#**************************************************************************************************
# Fit the CTRL mean curves and get the covariance (sigma) of the peak(s) - new column: fit_sigma
ctrl_3$fit_sigma <- apply(vect, 2, function(t) {
												  df.y <- c(as.numeric(ctrl_norm_mean[t,1:25]))
												  y <- c(1:25)
												  data <- data.frame(y = y, df.y = df.y)
												  Mean_ctrl_list <- as.numeric(unlist(ctrl_3[t,"ctrl_max"]))
												  n <- length(Mean_ctrl_list)
												  gauss <- numeric(0)
												  if (n==0) { gauss <- 0 } else {
												  				  								if (n==1) {
												  		  	         									  vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]], Mean_ctrl_list[1], 1 )
												  		  	         									  gauss <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$par[3],digits = 1 ),error=function(e) numeric(0))
												  				           								  } else {
												  				  												  if (n==2) {
												  				                 											 vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1 )
												  				                 											 gauss <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6)],digits = 1 ),error=function(e) numeric(0))
												  	  		                     											 } else {
												  				  																	if (n==3) {
												  				  																			   vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1 )
												  				  																			   gauss <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6,9)],digits = 1 ),error=function(e) numeric(0))
												  				 																			   } else {
												  				  																					  if (n==4) {
												  				  																								 vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1 )
												  				  																								 gauss <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6,9,12)],digits = 1 ),error=function(e) numeric(0))
												  				  																								 } else {
												                  																										 if (n==5) {
												                  																													vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1,ctrl_norm_mean[t,Mean_ctrl_list[5]],Mean_ctrl_list[5],1 )
												                  																													gauss <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6,9,12,15)],digits = 1 ),error=function(e) numeric(0))
												                																													} else {
												                   																														    vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1,ctrl_norm_mean[t,Mean_ctrl_list[5]],Mean_ctrl_list[5],1,ctrl_norm_mean[t,Mean_ctrl_list[6]],Mean_ctrl_list[6],1 )
												                   																														    gauss <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$par[c(3,6,9,12,15,18)],digits = 1 ),error=function(e) numeric(0))
												                  																															}
												  				  																										 }								
												  				  																						}
												  				  																		}
												  				  												   }
												  				  								}
												  gauss <- unlist(gauss)
												  if (length(gauss) == n) {as.numeric(gauss)} else {0}
												  })

#**************************************************************************************************
# Fit the CTRL mean curves and get the residual sum of the squares - new column: fit_res
ctrl_3$fit_res <- apply(vect, 2, function(t) {
												  df.y <- c(as.numeric(ctrl_norm_mean[t,1:25]))
												  y <- c(1:25)
												  data <- data.frame(y = y, df.y = df.y)
												  Mean_ctrl_list <- as.numeric(unlist(ctrl_3[t,"ctrl_max"]))
												  n <- length(Mean_ctrl_list)
												  gauss <- numeric(0)
												  if (n==0) { gauss <- 0 } else {
												  				  								if (n==1) {
												  		  	         									  vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]], Mean_ctrl_list[1], 1 )
												  		  	         									  gauss <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												  				           								  } else {
												  				  												  if (n==2) {
												  				                 											 vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1 )
												  				                 											 gauss <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												  	  		                     											 } else {
												  				  																	if (n==3) {
												  				  																			   vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1 )
												  				  																			   gauss <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												  				 																			   } else {
												  				  																					  if (n==4) {
												  				  																								 vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1 )
												  				  																								 gauss <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												  				  																								 } else {
												                  																										 if (n==5) {
												                  																													vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1,ctrl_norm_mean[t,Mean_ctrl_list[5]],Mean_ctrl_list[5],1 )
												                  																													gauss <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												                																													} else {
												                   																														    vector <- c( ctrl_norm_mean[t,Mean_ctrl_list[1]],Mean_ctrl_list[1],1,ctrl_norm_mean[t,Mean_ctrl_list[2]],Mean_ctrl_list[2],1,ctrl_norm_mean[t,Mean_ctrl_list[3]],Mean_ctrl_list[3],1,ctrl_norm_mean[t,Mean_ctrl_list[4]],Mean_ctrl_list[4],1,ctrl_norm_mean[t,Mean_ctrl_list[5]],Mean_ctrl_list[5],1,ctrl_norm_mean[t,Mean_ctrl_list[6]],Mean_ctrl_list[6],1 )
												                   																														    gauss <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) numeric(0))
												                  																															}
												  				  																										 }								
												  				  																						}
												  				  																		}
												  				  												   }
												  				  								}
												  gauss <- unlist(gauss)
												  if (length(gauss) == 1) {as.numeric(gauss)} else {0}
												  })												  

#**************************************************************************************************
# Check whether a fit curve could be found or not for the protein - new column: fitted
ctrl_3$fitted <- apply(ctrl_3, 1, function(x) { 
																				 list_c <- as.numeric(unlist(x$fit_c))
																				 list_m <- as.numeric(unlist(x$fit_mean))
																				 list_s <- as.numeric(unlist(x$fit_sigma))
																				 list <- c(list_c, list_m, list_s)
																				 
																				 if (sum(list) == 0) {FALSE} else {TRUE}
																				 })	


#**************************************************************************************************
# Control that the mean values are sorted.
# Fit_c and fit_sigma need to be sorted accordingly.
# Create a sorted parameters vector - new column: fit_param

ctrl_3$fit_param <- apply(vect, 2, function(z) {
												# Get the parameters after fitting
                        list_c <- as.numeric(unlist(ctrl_3[z,"fit_c"]))
												list_mean <- as.numeric(unlist(ctrl_3[z,"fit_mean"]))
												list_sigma <- as.numeric(unlist(ctrl_3[z,"fit_sigma"]))
												n <- length(list_c)
												vector <- numeric(0)
												
												# Sort the parameters
												if ((n==1) && c(list_c, list_mean, list_sigma) == c(0,0,0)) {fit_param <- c(0,0,0)} else { fit_param <- as.data.frame(matrix( c(list_c, list_mean, list_sigma),n,3,byrow="FALSE" ))
																																 fit_param <- fit_param[order(fit_param$V2),]
																																 fit_param <- unlist(fit_param)
																																 }
												fit_param
												})


rnase_3$fit_param <- apply(vect, 2, function(z) {
                        # Get the parameters after fitting
												list_c <- as.numeric(unlist(rnase_3[z,"fit_c"]))
												list_mean <- as.numeric(unlist(rnase_3[z,"fit_mean"]))
												list_sigma <- as.numeric(unlist(rnase_3[z,"fit_sigma"]))
												n <- length(list_c)
												vector <- numeric(0)
												
												# Sort the parameters
												if ((n==1) && c(list_c, list_mean, list_sigma) == c(0,0,0)) {fit_param <- c(0,0,0)} else { fit_param <- as.data.frame(matrix( c(list_c, list_mean, list_sigma),n,3,byrow="FALSE" ))
																																 fit_param <- fit_param[order(fit_param$V2),]
																																 fit_param <- unlist(fit_param)
																																 }
												fit_param
												})


#**************************************************************************************************
# Take care of the edges (for fitted mean <1 or >25).
# Set the fraction to min 1 and max 25.
# Calculate accordingly the amplitude of the fit at these positions.

ctrl_3$fit_mean_fxn <- apply(vect, 2, function(z) {
												   fit_param <- as.numeric(unlist(ctrl_3[z,"fit_param"]))
												   n <- length(fit_param)/3
												   fit_param <- as.data.frame(matrix( fit_param,n,3,byrow="FALSE" ))
												   list_c <- unlist(fit_param[,1])
												   list_mean <- unlist(fit_param[,2])
												   list_sigma <- unlist(fit_param[,3])
												   
												   if (list_mean[n] > 25) {list_mean[n] <- 25}
												   if (list_mean[1] < 1) {list_mean[1] <- 1}
												   unlist(list_mean)
												   if (ctrl_3[z,"fitted"] == "FALSE") {0} else {list_mean}
												   })
												   

ctrl_3$fit_c_fxn <- apply(vect, 2, function(z) {
												   fit_param <- as.numeric(unlist(ctrl_3[z,"fit_param"]))
												   n <- length(fit_param)/3
												   fit_param <- as.data.frame(matrix( fit_param,n,3,byrow="FALSE" ))
												   list_c <- unlist(fit_param[,1])
												   list_mean <- unlist(fit_param[,2])
												   list_sigma <- unlist(fit_param[,3])
												   
												   if (list_mean[n] > 25) {
												   						   q <- c(list_c[n],list_mean[n],list_sigma[n])
															 			   list_c[n] <- (q[1] * exp(-(25-q[2])**2/(2 * q[3]**2)))
												   						   }
												   list_mean[n] <- 25
												   if (list_mean[1] < 1) {
												   						  q <- c(list_c[1],list_mean[1],list_sigma[1])
															 			  list_c[1] <- (q[1] * exp(-(1-q[2])**2/(2 * q[3]**2)))
												   						  }
												   unlist(list_c)
												   list_c
												   })												   
												   

rnase_3$fit_mean_fxn <- apply(vect, 2, function(z) {
												   fit_param <- as.numeric(unlist(rnase_3[z,"fit_param"]))
												   n <- length(fit_param)/3
												   fit_param <- as.data.frame(matrix( fit_param,n,3,byrow="FALSE" ))
												   list_c <- unlist(fit_param[,1])
												   list_mean <- unlist(fit_param[,2])
												   list_sigma <- unlist(fit_param[,3])
												   
												   if (list_mean[n] > 25) {list_mean[n] <- 25}
												   if (list_mean[1] < 1) {list_mean[1] <- 1}
												   unlist(list_mean)
												   if (rnase_3[z,"fitted"] == "FALSE") {0} else {list_mean}
												   })
												   

rnase_3$fit_c_fxn <- apply(vect, 2, function(z) {
												   fit_param <- as.numeric(unlist(rnase_3[z,"fit_param"]))
												   n <- length(fit_param)/3
												   fit_param <- as.data.frame(matrix( fit_param,n,3,byrow="FALSE" ))
												   list_c <- unlist(fit_param[,1])
												   list_mean <- unlist(fit_param[,2])
												   list_sigma <- unlist(fit_param[,3])
												   
												   if (list_mean[n] > 25) {
												   						   q <- c(list_c[n],list_mean[n],list_sigma[n])
															 			   list_c[n] <- (q[1] * exp(-(25-q[2])**2/(2 * q[3]**2)))
												   						   }
												   list_mean[n] <- 25
												   if (list_mean[1] < 1) {
												   						  q <- c(list_c[1],list_mean[1],list_sigma[1])
															 			  list_c[1] <- (q[1] * exp(-(1-q[2])**2/(2 * q[3]**2)))
												   						  }
												   unlist(list_c)
												   list_c
												   })	


#****************************************************************
#### PART 3 : Gaussian fit on each curves                    ####
#****************************************************************

# Gaussian fit for 3x CTRL and 3x RNASE
# parameters to be optimized: C, mean, sigma -> fit_c, fit_mean, fit_sigma (from the fit of the mean curves)
# Function to get the sum of residual squares -> fit_res
# Number of peaks -> nb_peaks

#raw data:
#table.ctrl1.SW.norm 
#table.ctrl2.SW.norm 
#table.ctrl3.SW.norm 
#table.rnase1.SW.norm 
#table.rnase2.SW.norm 
#table.rnase3.SW.norm

# reference mean data with additional columns (temporary tables as defined earlier in the script):
# ctrl_3
# rnase_3

#define a function to calculate the gaussian fit for each curve - get the fit curve
#**************************************************************************************************
gauss_fit <- function(pos, table_mean, table_rep) {
  df.y <- c(as.numeric(table_rep[pos,1:25]))
  y <- c(1:25)
  data <- data.frame(y = y, df.y = df.y)
  list_c <- as.numeric(unlist(table_mean[pos,"fit_c"]))
  list_mean <- as.numeric(unlist(table_mean[pos,"fit_mean"]))
  list_sigma <- as.numeric(unlist(table_mean[pos,"fit_sigma"]))
  n <- length(list_c)
  vector <- numeric(0)
  dffit <- data.frame(y=seq(1, 25, 1))
  z <- seq(1, 25)
  if ((n==1) && (list_c == 100)) {dffit$df.y <- rep(0,25)} else {
    if (n==1) {
      vector <- c(list_c[1],list_mean[1],list_sigma[1])
      q <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) numeric(0))
      C1 <- q[1]
      mean1 <- q[2]
      sigma1 <- q[3]
      dffit$df.y <- C1 * exp(-(z-mean1)**2/(2 * sigma1**2))
    } else {
      if (n==2) {
        vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2])
        q <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) numeric(0))
        C1 <- q[1]
        mean1 <- q[2]
        sigma1 <- q[3]
        C2 <- q[4]
        mean2 <- q[5]
        sigma2 <- q[6]
        dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) ) 
      } else {
        if (n==3) {
          vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3])
          q <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) numeric(0))
          C1 <- q[1]
          mean1 <- q[2]
          sigma1 <- q[3]
          C2 <- q[4]
          mean2 <- q[5]
          sigma2 <- q[6]
          C3 <- q[7]
          mean3 <- q[8]
          sigma3 <- q[9]
          dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) ) 
        } else {
          if (n==4) {
            vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4])
            q <- tryCatch(round( (fit <- optim(vector, f4, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) numeric(0))
            C1 <- q[1]
            mean1 <- q[2]
            sigma1 <- q[3]
            C2 <- q[4]
            mean2 <- q[5]
            sigma2 <- q[6]
            C3 <- q[7]
            mean3 <- q[8]
            sigma3 <- q[9]
            C4 <- q[10]
            mean4 <- q[11]
            sigma4 <- q[12]
            dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) )
          } else {
            if (n==5) {
              vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5])
              q <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) numeric(0))
              C1 <- q[1]
              mean1 <- q[2]
              sigma1 <- q[3]
              C2 <- q[4]
              mean2 <- q[5]
              sigma2 <- q[6]
              C3 <- q[7]
              mean3 <- q[8]
              sigma3 <- q[9]
              C4 <- q[10]
              mean4 <- q[11]
              sigma4 <- q[12]
              C5 <- q[13]
              mean5 <- q[14]
              sigma5 <- q[15]
              dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) )
            } else {
              vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5],list_c[6],list_mean[6],list_sigma[6])
              q <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) numeric(0))
              C1 <- q[1]
              mean1 <- q[2]
              sigma1 <- q[3]
              C2 <- q[4]
              mean2 <- q[5]
              sigma2 <- q[6]
              C3 <- q[7]
              mean3 <- q[8]
              sigma3 <- q[9]
              C4 <- q[10]
              mean4 <- q[11]
              sigma4 <- q[12]
              C5 <- q[13]
              mean5 <- q[14]
              sigma5 <- q[15]
              C6 <- q[16]
              mean6 <- q[17]
              sigma6 <- q[18]
              dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) + C6 * exp(-(z-mean6)**2/(2 * sigma6**2)) ) 
            }
          }
        }
      }
    }
  }
  q <- unlist(q)
  if (length(q) == (3*n)) {dffit$df.y} else {rep(0,25)}
}

#**************************************************************************************************
# Make new tables with the fitted values for each curve

table.ctrl1.fit <- t(apply(vect, 2, function(d) { gauss_fit(d, ctrl_3, table.ctrl1.SW.norm)  }))
rownames(table.ctrl1.fit) <- rownames(table.ctrl1.SW.norm)
colnames(table.ctrl1.fit) <- colnames(table.ctrl1.SW.norm[1:25])
table.ctrl1.fit <- as.data.frame(table.ctrl1.fit)

table.ctrl2.fit <- t(apply(vect, 2, function(d) { gauss_fit(d, ctrl_3, table.ctrl2.SW.norm)  }))
rownames(table.ctrl2.fit) <- rownames(table.ctrl2.SW.norm)
colnames(table.ctrl2.fit) <- colnames(table.ctrl2.SW.norm[1:25])
table.ctrl2.fit <- as.data.frame(table.ctrl2.fit)

table.ctrl3.fit <- t(apply(vect, 2, function(d) { gauss_fit(d, ctrl_3, table.ctrl3.SW.norm)  }))
rownames(table.ctrl3.fit) <- rownames(table.ctrl3.SW.norm)
colnames(table.ctrl3.fit) <- colnames(table.ctrl3.SW.norm[1:25])
table.ctrl3.fit <- as.data.frame(table.ctrl3.fit)

table.rnase1.fit <- t(apply(vect, 2, function(d) { gauss_fit(d, rnase_3, table.rnase1.SW.norm)  }))
rownames(table.rnase1.fit) <- rownames(table.rnase1.SW.norm)
colnames(table.rnase1.fit) <- colnames(table.rnase1.SW.norm[1:25])
table.rnase1.fit <- as.data.frame(table.rnase1.fit)

table.rnase2.fit <- t(apply(vect, 2, function(d) { gauss_fit(d, rnase_3, table.rnase2.SW.norm)  }))
rownames(table.rnase2.fit) <- rownames(table.rnase2.SW.norm)
colnames(table.rnase2.fit) <- colnames(table.rnase2.SW.norm[1:25])
table.rnase2.fit <- as.data.frame(table.rnase2.fit)

table.rnase3.fit <- t(apply(vect, 2, function(d) { gauss_fit(d, rnase_3, table.rnase3.SW.norm)  }))
rownames(table.rnase3.fit) <- rownames(table.rnase3.SW.norm)
colnames(table.rnase3.fit) <- colnames(table.rnase2.SW.norm[1:25])
table.rnase3.fit <- as.data.frame(table.rnase3.fit)


#**************************************************************************************************
# Replace NA, NaN with 0 if any
table.ctrl1.fit <- rapply(table.ctrl1.fit, f=function(x) ifelse(is.na(x),0,x), how="replace") 
table.ctrl2.fit <- rapply(table.ctrl2.fit, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.ctrl3.fit <- rapply(table.ctrl3.fit, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase1.fit <- rapply(table.rnase1.fit, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase2.fit <- rapply(table.rnase2.fit, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase3.fit <- rapply(table.rnase3.fit, f=function(x) ifelse(is.na(x),0,x), how="replace" )

table.ctrl1.fit <- rapply(table.ctrl1.fit, f=function(x) ifelse(is.nan(x),0,x), how="replace") 
table.ctrl2.fit <- rapply(table.ctrl2.fit, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.ctrl3.fit <- rapply(table.ctrl3.fit, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase1.fit <- rapply(table.rnase1.fit, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase2.fit <- rapply(table.rnase2.fit, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase3.fit <- rapply(table.rnase3.fit, f=function(x) ifelse(is.nan(x),0,x), how="replace" )


#**************************************************************************************************
# check the sum of the fitted curves. Should be close to 100 (%)

table.ctrl1.fit$check_sum <- apply(vect, 2, function(x) { check <- round(sum(as.numeric(table.ctrl1.fit[x,1:25])), digits=1)
														check
														})
table.ctrl2.fit$check_sum <- apply(vect, 2, function(x) { check <- round(sum(as.numeric(table.ctrl2.fit[x,1:25])), digits=1)
														check
														})
table.ctrl3.fit$check_sum <- apply(vect, 2, function(x) { check <- round(sum(as.numeric(table.ctrl3.fit[x,1:25])), digits=1)
														check
														})
table.rnase1.fit$check_sum <- apply(vect, 2, function(x) { check <- round(sum(as.numeric(table.rnase1.fit[x,1:25])), digits=1)
														check
														})
table.rnase2.fit$check_sum <- apply(vect, 2, function(x) { check <- round(sum(as.numeric(table.rnase2.fit[x,1:25])), digits=1)
														check
														})
table.rnase3.fit$check_sum <- apply(vect, 2, function(x) { check <- round(sum(as.numeric(table.rnase3.fit[x,1:25])), digits=1)
														check
														})

#**************************************************************************************************
# get the fit parameters
# all at once in a vector (in the order: C1, mean1, sigma1, C2, mean2, sigma2, ..)

gauss_fit_param <- function(pos, table_mean, table_rep) {
								df.y <- c(as.numeric(table_rep[pos,1:25]))
								y <- c(1:25)
								data <- data.frame(y = y, df.y = df.y)
								list_c <- as.numeric(unlist(table_mean[pos,"fit_c"]))
								list_mean <- as.numeric(unlist(table_mean[pos,"fit_mean"]))
								list_sigma <- as.numeric(unlist(table_mean[pos,"fit_sigma"]))
								list <- c(list_c, list_mean, list_sigma)
								n <- length(list_c)
								vector <- numeric(0)
									  if ((n==1) && (list == c(0,0,0))) {q <- rep(0,3)} else {
												  											if (n==1) {
												  														vector <- c(list_c[1],list_mean[1],list_sigma[1])
												  														q <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) c(0,0,0))
												  														
												  														} else {
												  																if (n==2) {
												  																			vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2])
																															q <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) c(0,0,0))
																															
												  																			} else {
												  																					if (n==3) {
												  																								vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3])
												  																								q <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) c(0,0,0))
												  				
												  																								} else {
												  																										if (n==4) {
												  																													vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4])
												  																													q <- tryCatch(round( (fit <- optim(vector, f4,data = data,  method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) c(0,0,0))
												  																											
												  																													} else {
												  																															if (n==5) {
												  																																		vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5])
												  																																		q <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) c(0,0,0))
												  																																	
												  																																		} else {
												  																																				vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5],list_c[6],list_mean[6],list_sigma[6])
												  																																				q <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) c(0,0,0))
												  																																			
												  																																				}
												  																															}
												  																										}
												  																					}
												  																}
												  											}

																#vector mit parameters
																q
								}


table.ctrl1.fit$fit_param <- apply(vect, 2, function(x) { gauss_fit_param(x, ctrl_3, table.ctrl1.SW.norm) })
table.ctrl2.fit$fit_param <- apply(vect, 2, function(x) { gauss_fit_param(x, ctrl_3, table.ctrl2.SW.norm) })
table.ctrl3.fit$fit_param <- apply(vect, 2, function(x) { gauss_fit_param(x, ctrl_3, table.ctrl3.SW.norm) })
table.rnase1.fit$fit_param <- apply(vect, 2, function(x) { gauss_fit_param(x, rnase_3, table.rnase1.SW.norm) })
table.rnase2.fit$fit_param <- apply(vect, 2, function(x) { gauss_fit_param(x, rnase_3, table.rnase2.SW.norm) })
table.rnase3.fit$fit_param <- apply(vect, 2, function(x) { gauss_fit_param(x, rnase_3, table.rnase3.SW.norm) })


# define a function to calculate the gaussian fit for each curve - get the residue
#**************************************************************************************************

gauss_fit_res <- function(pos, table_mean, table_rep) {
  df.y <- c(as.numeric(table_rep[pos,1:25]))
  y <- c(1:25)
  data <- data.frame(y = y, df.y = df.y)
  list_c <- as.numeric(unlist(table_mean[pos,"fit_c"]))
  list_mean <- as.numeric(unlist(table_mean[pos,"fit_mean"]))
  list_sigma <- as.numeric(unlist(table_mean[pos,"fit_sigma"]))
  list <- c(list_c, list_mean, list_sigma)
  n <- length(list_c)
  vector <- numeric(0)
  if ((n==1) && (list == c(0,0,0))) {q <- rep(0,3)} else {
    if (n==1) {
      vector <- c(list_c[1],list_mean[1],list_sigma[1])
      q <- tryCatch(round( (fit <- optim(vector, f1, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) c(0,0,0))
      
    } else {
      if (n==2) {
        vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2])
        q <- tryCatch(round( (fit <- optim(vector, f2, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) c(0,0,0))
        
      } else {
        if (n==3) {
          vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3])
          q <- tryCatch(round( (fit <- optim(vector, f3, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) c(0,0,0))
          
        } else {
          if (n==4) {
            vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4])
            q <- tryCatch(round( (fit <- optim(vector, f4,data = data,  method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) c(0,0,0))
            
          } else {
            if (n==5) {
              vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5])
              q <- tryCatch(round( (fit <- optim(vector, f5, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) c(0,0,0))
              
            } else {
              vector <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5],list_c[6],list_mean[6],list_sigma[6])
              q <- tryCatch(round( (fit <- optim(vector, f6, data = data, method="BFGS", control=list(reltol=1e-9)))$value,digits = 1 ),error=function(e) c(0,0,0))
              
            }
          }
        }
      }
    }
  }
  
  #vector mit parameters
  q <- unlist(q)
  if (length(q) == 1) {
    if (q > 1000) {0} else {
      as.numeric(q)}
  } else {0}
}


table.ctrl1.fit$fit_res <- apply(vect, 2, function(x) { gauss_fit_res(x, ctrl_3, table.ctrl1.SW.norm) })
table.ctrl2.fit$fit_res <- apply(vect, 2, function(x) { gauss_fit_res(x, ctrl_3, table.ctrl2.SW.norm) })
table.ctrl3.fit$fit_res <- apply(vect, 2, function(x) { gauss_fit_res(x, ctrl_3, table.ctrl3.SW.norm) })
												  
table.rnase1.fit$fit_res <- apply(vect, 2, function(x) { gauss_fit_res(x, rnase_3, table.rnase1.SW.norm) })
table.rnase2.fit$fit_res <- apply(vect, 2, function(x) { gauss_fit_res(x, rnase_3, table.rnase2.SW.norm) })
table.rnase3.fit$fit_res <- apply(vect, 2, function(x) { gauss_fit_res(x, rnase_3, table.rnase3.SW.norm) })


#**************************************************************************************************
# Normalize the fit to 100

table.ctrl1.fit.norm <- table.ctrl1.fit[1:25] * 100 / rowSums(table.ctrl1.fit[1:25])
table.ctrl2.fit.norm <- table.ctrl2.fit[1:25] * 100 / rowSums(table.ctrl2.fit[1:25])
table.ctrl3.fit.norm <- table.ctrl3.fit[1:25] * 100 / rowSums(table.ctrl3.fit[1:25])

table.rnase1.fit.norm <- table.rnase1.fit[1:25] * 100 / rowSums(table.rnase1.fit[1:25])
table.rnase2.fit.norm <- table.rnase2.fit[1:25] * 100 / rowSums(table.rnase2.fit[1:25])
table.rnase3.fit.norm <- table.rnase3.fit[1:25] * 100 / rowSums(table.rnase3.fit[1:25])

# Replace NA, NaN with 0
table.ctrl1.fit.norm <- rapply(table.ctrl1.fit.norm, f=function(x) ifelse(is.na(x),0,x), how="replace") 
table.ctrl2.fit.norm <- rapply(table.ctrl2.fit.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.ctrl3.fit.norm <- rapply(table.ctrl3.fit.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase1.fit.norm <- rapply(table.rnase1.fit.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase2.fit.norm <- rapply(table.rnase2.fit.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
table.rnase3.fit.norm <- rapply(table.rnase3.fit.norm, f=function(x) ifelse(is.na(x),0,x), how="replace" )

table.ctrl1.fit.norm <- rapply(table.ctrl1.fit.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace") 
table.ctrl2.fit.norm <- rapply(table.ctrl2.fit.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.ctrl3.fit.norm <- rapply(table.ctrl3.fit.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase1.fit.norm <- rapply(table.rnase1.fit.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase2.fit.norm <- rapply(table.rnase2.fit.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
table.rnase3.fit.norm <- rapply(table.rnase3.fit.norm, f=function(x) ifelse(is.nan(x),0,x), how="replace" )


#***************************************************************************************************
# Add columns to the tables with the fitted values:
# nb_max
# fit_param
# fit_res
# check_sum

table.ctrl1.fit.norm$nb_max <- ctrl_3$nb_max
table.ctrl2.fit.norm$nb_max <- ctrl_3$nb_max
table.ctrl3.fit.norm$nb_max <- ctrl_3$nb_max

table.rnase1.fit.norm$nb_max <- rnase_3$nb_max
table.rnase2.fit.norm$nb_max <- rnase_3$nb_max
table.rnase3.fit.norm$nb_max <- rnase_3$nb_max


table.ctrl1.fit.norm$fit_param <- table.ctrl1.fit$fit_param
table.ctrl2.fit.norm$fit_param <- table.ctrl2.fit$fit_param
table.ctrl3.fit.norm$fit_param <- table.ctrl3.fit$fit_param

table.rnase1.fit.norm$fit_param <- table.rnase1.fit$fit_param
table.rnase2.fit.norm$fit_param <- table.rnase2.fit$fit_param
table.rnase3.fit.norm$fit_param <- table.rnase3.fit$fit_param


table.ctrl1.fit.norm$fit_res <- table.ctrl1.fit$fit_res
table.ctrl2.fit.norm$fit_res <- table.ctrl2.fit$fit_res
table.ctrl3.fit.norm$fit_res <- table.ctrl3.fit$fit_res

table.rnase1.fit.norm$fit_res <- table.rnase1.fit$fit_res
table.rnase2.fit.norm$fit_res <- table.rnase2.fit$fit_res
table.rnase3.fit.norm$fit_res <- table.rnase3.fit$fit_res


table.ctrl1.fit.norm$check_sum <- table.ctrl1.fit$check_sum
table.ctrl2.fit.norm$check_sum <- table.ctrl2.fit$check_sum
table.ctrl3.fit.norm$check_sum <- table.ctrl3.fit$check_sum

table.rnase1.fit.norm$check_sum <- table.rnase1.fit$check_sum
table.rnase2.fit.norm$check_sum <- table.rnase2.fit$check_sum
table.rnase3.fit.norm$check_sum <- table.rnase3.fit$check_sum

												  
#***************************************************************************************************
# create tables with fit values, defined to 0.1 increments - will help to later calculate the p_values
# normalized to 100 (area under the curve)

# general function to create the tables
fine_fit <- function(fit)			{ 
												  dffit <- data.frame(y=seq(0, 25, 0.1))
												  z <- round(seq(0, 25, 0.1), digits=1)
												  n <- length(fit)/3
												  if (n==0) {dffit$df.y <- rep(0,251)} else {
												  											if (n==1) {
												  														q <- fit
												  														C1 <- q[1]
																										mean1 <- q[2]
																										sigma1 <- q[3]
																										dffit$df.y <- C1 * exp(-(z-mean1)**2/(2 * sigma1**2))
												  														} else {
												  																if (n==2) {
												  																			q <- fit
												  																			C1 <- q[1]
																															mean1 <- q[2]
																															sigma1 <- q[3]
																															C2 <- q[4]
																															mean2 <- q[5]
																															sigma2 <- q[6]
																															dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) ) 
												  																			} else {
												  																					if (n==3) {
												  																								q <- fit
												  																								C1 <- q[1]
																																				mean1 <- q[2]
																																				sigma1 <- q[3]
																																				C2 <- q[4]
																																				mean2 <- q[5]
																																				sigma2 <- q[6]
																																				C3 <- q[7]
																																				mean3 <- q[8]
																																				sigma3 <- q[9]
																																				dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) ) 
												  																								} else {
												  																										if (n==4) {
												  																													q <- fit
												  																													C1 <- q[1]
																																									mean1 <- q[2]
																																									sigma1 <- q[3]
																																									C2 <- q[4]
																																									mean2 <- q[5]
																																									sigma2 <- q[6]
																																									C3 <- q[7]
																																									mean3 <- q[8]
																																									sigma3 <- q[9]
																																									C4 <- q[10]
																																									mean4 <- q[11]
																																									sigma4 <- q[12]
																																									dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) )
												  																													} else {
												  																															if (n==5) {
												  																																		q <- fit
												  																																		C1 <- q[1]
																																														mean1 <- q[2]
																																														sigma1 <- q[3]
																																														C2 <- q[4]
																																														mean2 <- q[5]
																																														sigma2 <- q[6]
																																														C3 <- q[7]
																																														mean3 <- q[8]
																																														sigma3 <- q[9]
																																														C4 <- q[10]
																																														mean4 <- q[11]
																																														sigma4 <- q[12]
																																														C5 <- q[13]
																																														mean5 <- q[14]
																																														sigma5 <- q[15]
																																														dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) )
												  																																		} else {
												  																																				q <- fit
												  																																				C1 <- q[1]
																																																mean1 <- q[2]
																																																sigma1 <- q[3]
																																																C2 <- q[4]
																																																mean2 <- q[5]
																																																sigma2 <- q[6]
																																																C3 <- q[7]
																																																mean3 <- q[8]
																																																sigma3 <- q[9]
																																																C4 <- q[10]
																																																mean4 <- q[11]
																																																sigma4 <- q[12]
																																																C5 <- q[13]
																																																mean5 <- q[14]
																																																sigma5 <- q[15]
																																																C6 <- q[16]
																																																mean6 <- q[17]
																																																sigma6 <- q[18]
																																																dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) + C6 * exp(-(z-mean6)**2/(2 * sigma6**2)) ) 
												  																																				}
												  																															}
												  																										}
												  																					}
												  																}
												  											}
											# *0.1 to scale the resulting curve at the level of the normalized raw data. The sum here is at the end 1000 and not 100, but the values fits to the raw data.
												  dffit$df.y*100/(sum(dffit$df.y)*0.1)
											}							


table.ctrl1.fit.fine <- t(apply(vect, 2, function(t) {
													fit_param <- as.numeric(unlist(table.ctrl1.fit[t,"fit_param"]))
													fine_fit(fit_param)
													}))
rownames(table.ctrl1.fit.fine) <- rownames(table.ctrl1.fit)
colnames(table.ctrl1.fit.fine) <- lapply(seq(0,25,0.1), as.character)

table.ctrl2.fit.fine <- t(apply(vect, 2, function(t) {
													fit_param <- as.numeric(unlist(table.ctrl2.fit[t,"fit_param"]))
													fine_fit(fit_param)
													}))
rownames(table.ctrl2.fit.fine) <- rownames(table.ctrl2.fit)
colnames(table.ctrl2.fit.fine) <- lapply(seq(0,25,0.1), as.character)

table.ctrl3.fit.fine <- t(apply(vect, 2, function(t) {
													fit_param <- as.numeric(unlist(table.ctrl3.fit[t,"fit_param"]))
													fine_fit(fit_param)
													}))
rownames(table.ctrl3.fit.fine) <- rownames(table.ctrl3.fit)
colnames(table.ctrl3.fit.fine) <- lapply(seq(0,25,0.1), as.character)


table.rnase1.fit.fine <- t(apply(vect, 2, function(t) {
													fit_param <- as.numeric(unlist(table.rnase1.fit[t,"fit_param"]))
													fine_fit(fit_param)
													}))
rownames(table.rnase1.fit.fine) <- rownames(table.rnase1.fit)
colnames(table.rnase1.fit.fine) <- lapply(seq(0,25,0.1), as.character)

table.rnase2.fit.fine <- t(apply(vect, 2, function(t) {
													fit_param <- as.numeric(unlist(table.rnase2.fit[t,"fit_param"]))
													fine_fit(fit_param)
													}))
rownames(table.rnase2.fit.fine) <- rownames(table.rnase2.fit)
colnames(table.rnase2.fit.fine) <- lapply(seq(0,25,0.1), as.character)

table.rnase3.fit.fine <- t(apply(vect, 2, function(t) {
													fit_param <- as.numeric(unlist(table.rnase3.fit[t,"fit_param"]))
													fine_fit(fit_param)
													}))
rownames(table.rnase3.fit.fine) <- rownames(table.rnase3.fit)
colnames(table.rnase3.fit.fine) <- lapply(seq(0,25,0.1), as.character)


#***************************************************************************************************
# Function to calculate the fit value at a given peak (with 0.1 instead of 1 fraction precision)
# The function will be used in the next step to evaluate the p-value at each maximum
# Return the amplitude of the fit curve at a given position

fine_fit_corr <- function(pos,vector,table_rep)	 { 
  df.y <- c(as.numeric(table_rep[pos,1:25]))
  y <- c(1:25)
  q <- vector
  n <- length(vector)/3
  gauss <- numeric(0)
  value <- numeric(0)
  
  fun_f1 <- function(q) {
    C1 <- q[1]
    mean1 <- q[2]
    sigma1 <- q[3]
    res <- (C1 * exp(-(y-mean1)**2/(2 * sigma1**2))) - df.y
    sum(res * res)
  }
  
  if (n==0) {df.z <- rep(0,251)}
  else {		
    gauss <- tryCatch(round( (fit <- optim(vector, fun_f1, method="BFGS", control=list(reltol=1e-9)))$par,digits = 1 ),error=function(e) numeric(0))
    z <- round(seq(0, 25, 0.1),digits=1)
    if (length(gauss) == 0) {df.z <- rep(0,251)} else { df.z <- (gauss[1] * exp(-(z-gauss[2])**2/(2 * gauss[3]**2))) }
    }
  
  df.z <- df.z*100/sum(df.z)
  value <- df.z[q[2]/0.1+1]
  value
}


#***********************************************************************
#### PART 4 : quality control of the fit and adjustment if possible ####
#***********************************************************************

# data used:
# ctrl_3 and ctrl_norm_mean
# rnase_3 and rnase_norm_mean

#************************************************************************
# Create summary tables

ctrl_mean <- ctrl_3
ctrl_mean[1:25] <- ctrl_norm_mean[1:25]

rnase_mean <- rnase_3
rnase_mean[1:25] <- rnase_norm_mean[1:25]

#***************************************************************************************************
# Area under the curve: here, area under the respective maxima

fun_gauss_area <- function(x) {
						 C <- x[1]
						 Mean <- x[2]
						 Sigma <- x[3]
						 y <- round(seq(1,25,0.2), digits=1)
						 if (Sigma == 0) {res <- 0} else {
						 								  res <- C * exp(-(y-Mean)**2/(2 * Sigma**2))
						 								  }
						 res <- round(sum(res)*0.2, digits=1)
						 res
						 }

# add a column with areas under each maximum and sum of the areas
# column "fit_area"
# column "sum_area" - the sum of the areas is expected to be close to 100 (due to the previous normalization step)
ctrl_mean$fit_area <- apply(vect, 2, function(z) { 
												  list_c <- as.numeric(unlist(ctrl_mean[z,"fit_c"]))
												  list_mean <- as.numeric(unlist(ctrl_mean[z,"fit_mean"]))
												  list_sigma <- as.numeric(unlist(ctrl_mean[z,"fit_sigma"]))
												  n <- as.numeric(unlist(ctrl_mean[z,"nb_max"]))
												  area <- numeric(0)
												  vector <- numeric(0)
												  area_list <- numeric(0)
												  
												  if (n == 0) {area_list <- 0} else {
												  						for (i in 1:n) {
												  										vector <- c(list_c[i],list_mean[i],list_sigma[i])
												  										area <- fun_gauss_area(vector)
												  										area_list <- c(area_list, area)
												  										}
												  						}
												  area_list
												  }
							)

ctrl_mean$sum_area <- apply(vect, 2, function(z) {
												 list_area <- as.numeric(unlist(ctrl_mean[z,"fit_area"]))
												 sum_area <- sum(list_area)
												 sum_area
												})

rnase_mean$fit_area <- apply(vect, 2, function(z) { 
												  list_c <- as.numeric(unlist(rnase_mean[z,"fit_c"]))
												  list_mean <- as.numeric(unlist(rnase_mean[z,"fit_mean"]))
												  list_sigma <- as.numeric(unlist(rnase_mean[z,"fit_sigma"]))
												  n <- as.numeric(unlist(rnase_mean[z,"nb_max"]))
												  area <- numeric(0)
												  vector <- numeric(0)
												  area_list <- numeric(0)
												  
												  if (n == 0) {area_list <- 0} else {
												  						for (i in 1:n) {
												  										vector <- c(list_c[i],list_mean[i],list_sigma[i])
												  										area <- fun_gauss_area(vector)
												  										area_list <- c(area_list, area)
												  										}
												  						}
												  area_list
												  }
							)

rnase_mean$sum_area <- apply(vect, 2, function(z) {
												 list_area <- as.numeric(unlist(rnase_mean[z,"fit_area"]))
												 sum_area <- sum(list_area)
												 sum_area
												})

#****************************************************************************************************************************
# Label peaks which have either a residual value > 260 (no good fitting) or negative area (area under the curve < 0 after fitting)
# add column "pb_fit"

ctrl_mean$pb_fit <- apply(vect, 2, function(x) {
												res <- as.numeric(unlist(ctrl_mean[x,"fit_res"]))
												sigma <- as.numeric(unlist(ctrl_mean[x,"fit_sigma"]))
												list_area <- unlist(ctrl_mean[x, "fit_area"])
												neg <- any(list_area < 0)
												nl <- any(sigma == 0) && (res != 0)
												if (res > 260 || neg == "TRUE" || nl == "TRUE") {TRUE} else FALSE
												})
												
rnase_mean$pb_fit <- apply(vect, 2, function(x) {
												res <- as.numeric(unlist(rnase_mean[x,"fit_res"]))
												sigma <- as.numeric(unlist(rnase_mean[x,"fit_sigma"]))
												list_area <- unlist(rnase_mean[x, "fit_area"])
												neg <- any(list_area < 0)
												nl <- any(sigma == 0) && (res != 0)
												if (res > 260 || neg == "TRUE" || nl == "TRUE") {TRUE} else FALSE
												})												


#***************************************************************************************************
#### Make correction for peaks labeled as "problematic" peaks in the column "pb_fit"            ####

# Create a subtable with the problematic peaks from the control table
ctrl_mean_corr <- ctrl_mean[ctrl_mean$pb_fit == "TRUE",]
ctrl_mean_corr <- ctrl_mean_corr[, c(1:25,28)]

# Define a counter (vect3) adapted to the size of the table ctrl_mean_corr
# It will be used in the function "apply" below.
vect3 <- matrix(c(1:dim(ctrl_mean_corr)[1]), 1, dim(ctrl_mean_corr)[1])

ctrl_mean_corr$fit_param_corr <- apply(vect3, 2, function(x) {
															  list_mean <- as.numeric(unlist(ctrl_mean_corr[x,"ctrl_max"]))
															  
															  list_c <- ctrl_mean_corr[x,as.numeric(unlist(ctrl_mean_corr[x,]$ctrl_max))]
															  
															  n <- length(list_c)
															  
															  fit_param <- numeric(0)
															  
															  if (n == 0) {fit_param <- 0} else {
															  		          for (i in 1:n) {
															  				  		    # use raw data to fit
															  				  		    z <- c(1:25)
															  				  		    df.z <- ctrl_mean_corr[x,1:25]
															  				  
															  				  		    # values after maxima
															  				  	      if (list_mean[i] == 25) {df.z <- df.z} else {
															  				  		                            p <- list_mean[i]+1
															  				  		                            df.z[p:25] <- (df.z[p:25]-df.z[(list_mean[i]):24]<=0)*df.z[p:25]
															  				  		                            m <- which(df.z[p:25] == 0)
															  				  		                            if (length(m) == 0) {df.z[p:25] <- df.z[p:25]} else {df.z[(p+min(m)-1):25] <- 0}
															  				  		                            }
															  				  
															  				  		    # values before maxima
															  				 		      p <- list_mean[i]-1
															  				  		    df.z[1:p] <- (df.z[2:list_mean[i]]-df.z[1:p]>=0)*df.z[1:p]
															  				  		    m <- which(df.z[1:p] == 0)      
															  				  		    if (length(m) == 0) {df.z[1:p] <- df.z[1:p]} else {df.z[1:max(m)] <- 0}
															  				  
															  				  		    # function to be minimized
															  				  		    fun1 <- function(q) {
															  				  					                 C1 <- q[1]
															  				  					                 mean1 <- q[2]
															  				  					                 sigma1 <- q[3]
															  				  					                 res <- (C1 * exp(-(z-mean1)**2/(2 * sigma1**2))) - df.z
															  				 				 	                   sum(res * res)
															  				  					                 }
																			  
																			  		      # find fit parameters
																			  		      vector <- c(list_c[i], list_mean[i],2)
															  				  		    fit <- optim(vector, fun1, method="BFGS", control=list(reltol=1e-9))
															  				  		    fit_param <- c(fit_param, fit$par)
															  				  		    fit_param <- round(as.numeric(fit_param), digits=1)
															  				  		    } # end for loop
															  					
															  					    fit_param <- unique(matrix(fit_param,n, 3, byrow=TRUE))
															  					    fit_param <- as.vector(t(fit_param))
															  					    fit_param <- as.data.frame(matrix(fit_param,1, length(fit_param), byrow=TRUE))
															  					    }
															                fit_param
															  })

# Distribute the corrected fit parameters in respective columns for later usage.
ctrl_mean_corr$fit_c_corr <- apply(vect3, 2, function(x) {
															 list_param <- unlist(ctrl_mean_corr[x, "fit_param_corr"])
															 n <- length(list_param)/3
															 fit_c <- unlist(matrix(list_param,n, 3, byrow=TRUE)[,1])
															 fit_c <- as.data.frame(matrix(fit_c,1, length(fit_c), byrow=TRUE))
															 fit_c
															 })

ctrl_mean_corr$fit_mean_corr <- apply(vect3, 2, function(x) {
															 list_param <- unlist(ctrl_mean_corr[x, "fit_param_corr"])
															 n <- length(list_param)/3
															 fit_mean <- unlist(matrix(list_param,n, 3, byrow=TRUE)[,2])
															 fit_mean <- as.data.frame(matrix(fit_mean,1, length(fit_mean), byrow=TRUE))
															 fit_mean
															 })

ctrl_mean_corr$fit_sigma_corr <- apply(vect3, 2, function(x) {
															 list_param <- unlist(ctrl_mean_corr[x, "fit_param_corr"])
															 n <- length(list_param)/3
															 fit_s <- unlist(matrix(list_param,n, 3, byrow=TRUE)[,3])
															 fit_s <- as.data.frame(matrix(fit_s,1, length(fit_s), byrow=TRUE))
															 fit_s
															 })

# Function to calculate the area under the curve at a specific peak
# as defined by the parameters C, Mean and Sigma
gauss_area <- function(x) {
  C <- x[1]
  Mean <- x[2]
  Sigma <- x[3]
  y <- round(seq(1,25,0.2), digits=1)
  
  res <- C * exp(-(y-Mean)**2/(2 * Sigma**2))
  
  res <- round(sum(res)*0.2, digits=1)
  res
}

ctrl_mean_corr$fit_area <- apply(vect3, 2, function(t) { 
												  list <- as.numeric(unlist(ctrl_mean_corr[t,"fit_param_corr"]))
												  n <- length(list)/3
												  vector <- numeric(0)
												  area_list <- numeric(0)
												  
												  for (i in 1:n) {
												  				  vector <- list[(i*3-2):(i*3)]
												  				  area <- gauss_area(vector)
												  				  area_list <- c(area_list, area)
												  				  }
												  unlist(area_list)
												  area_list <- as.data.frame(matrix(area_list,1, length(area_list), byrow=TRUE))
												  area_list
												  })

												  
#***************************************************************************************************												  
# Same for the rnase table

rnase_mean_corr <- rnase_mean[rnase_mean$pb_fit == "TRUE",]
rnase_mean_corr <- rnase_mean_corr[, c(1:25,28)]

vect4 <- matrix(c(1:dim(rnase_mean_corr)[1]), 1, dim(rnase_mean_corr)[1])

rnase_mean_corr$fit_param_corr <- apply(vect4, 2, function(x) {
															  list_mean <- as.numeric(unlist(rnase_mean_corr[x,"rnase_max"]))
															  
															  list_c <- rnase_mean_corr[x,as.numeric(unlist(rnase_mean_corr[x,]$rnase_max))]
															  
															  n <- length(list_c)
															  
															  fit_param <- numeric(0)
															  
															  if (n == 0) {fit_param <- 0} else {
															  		for (i in 1:n) {
															  				  		# raw data to fit
															  				  		z <- c(1:25)
															  				  		df.z <- rnase_mean_corr[x,1:25]
															  				  
															  				  		
															  				  		# values after maxima
															  				  		if (list_mean[i] == 25) {df.z <- df.z} else {
															  				  		p <- list_mean[i]+1
															  				  		df.z[p:25] <- (df.z[p:25]-df.z[(list_mean[i]):24]<=0)*df.z[p:25]
															  				  		m <- which(df.z[p:25] == 0)
															  				  		if (length(m) == 0) {df.z[p:25] <- df.z[p:25]} else {df.z[(p+min(m)-1):25] <- 0}}
															  				  
															  				  		# values before maxima
															  				 		 p <- list_mean[i]-1
															  				  		df.z[1:p] <- (df.z[2:list_mean[i]]-df.z[1:p]>=0)*df.z[1:p]
															  				  		m <- which(df.z[1:p] == 0)
															  				  		if (length(m) == 0) {df.z[1:p] <- df.z[1:p]} else {df.z[1:max(m)] <- 0}
															  				  
															  				  		# function to be minimized
															  				  		f1 <- function(q) {
															  				  					C1 <- q[1]
															  				  					mean1 <- q[2]
															  				  					sigma1 <- q[3]
															  				  					res <- (C1 * exp(-(z-mean1)**2/(2 * sigma1**2))) - df.z
															  				 				 	sum(res * res)
															  				  					}
																			  
																			  		# find fit parameters
																			  		vector <- c(list_c[i], list_mean[i],2)
															  				  		fit <- optim(vector, f1, method="BFGS", control=list(reltol=1e-9))
															  				  		fit_param <- c(fit_param, fit$par)
															  				  		fit_param <- round(as.numeric(fit_param), digits=1)
															  				  		}
															  					
															            fit_param <- unique(matrix(fit_param,n, 3, byrow=TRUE))
															            fit_param <- (as.vector(t(fit_param)))
															            fit_param <- as.data.frame(matrix(fit_param,1, length(fit_param), byrow=TRUE)) 
															  					}
															  fit_param
															  })

rnase_mean_corr$fit_c_corr <- apply(vect4, 2, function(x) {
															 list_param <- unlist(rnase_mean_corr[x, "fit_param_corr"])
															 n <- length(list_param)/3
															 fit_c <- unlist(matrix(list_param,n, 3, byrow=TRUE)[,1])
															 fit_c <- as.data.frame(matrix(fit_c,1, length(fit_c), byrow=TRUE))
															 fit_c
															 })

rnase_mean_corr$fit_mean_corr <- apply(vect4, 2, function(x) {
															 list_param <- unlist(rnase_mean_corr[x, "fit_param_corr"])
															 n <- length(list_param)/3
															 fit_mean <- unlist(matrix(list_param,n, 3, byrow=TRUE)[,2])
															 fit_mean <- as.data.frame(matrix(fit_mean,1, length(fit_mean), byrow=TRUE))
															 fit_mean
															 })

rnase_mean_corr$fit_sigma_corr <- apply(vect4, 2, function(x) {
															 list_param <- unlist(rnase_mean_corr[x, "fit_param_corr"])
															 n <- length(list_param)/3
															 fit_s <- unlist(matrix(list_param,n, 3, byrow=TRUE)[,3])
															 fit_s <- as.data.frame(matrix(fit_s,1, length(fit_s), byrow=TRUE))
															 fit_s
															 })

rnase_mean_corr$fit_area <- apply(vect4, 2, function(t) { 
												  list <- as.numeric(unlist(rnase_mean_corr[t,"fit_param_corr"]))
												  n <- length(list)/3
												  vector <- numeric(0)
												  area_list <- numeric(0)
												  
												  for (i in 1:n) {
												  				  vector <- list[(i*3-2):(i*3)]
												  				  area <- gauss_area(vector)
												  				  area_list <- c(area_list, area)
												  				  }
												  unlist(area_list)
												  area_list <- as.data.frame(matrix(area_list,1, length(area_list), byrow=TRUE))
												  area_list
												  })


#**************************************************************************************************
# Take care of the edge (for mean values <1 or >25)

ctrl_mean_corr$fit_mean_fxn <- apply(vect3, 2, function(z) {
												   list_c <- as.numeric(unlist(ctrl_mean_corr[z,"fit_c_corr"]))
												   list_mean <- as.numeric(unlist(ctrl_mean_corr[z,"fit_mean_corr"]))
												   list_sigma <- as.numeric(unlist(ctrl_mean_corr[z,"fit_sigma_corr"]))
												   n <- length(list_c)
												   if (list_mean[n] > 25) {list_mean[n] <- 25}
												   if (list_mean[1] < 1) {list_mean[1] <- 1}
												   list_mean <- as.data.frame(matrix(list_mean,1, length(list_mean), byrow=TRUE))
												   list_mean
												   })
												   
ctrl_mean_corr$fit_c_fxn <- apply(vect3, 2, function(z) {
												   list_c <- as.numeric(unlist(ctrl_mean_corr[z,"fit_c_corr"]))
												   list_mean <- as.numeric(unlist(ctrl_mean_corr[z,"fit_mean_corr"]))
												   list_sigma <- as.numeric(unlist(ctrl_mean_corr[z,"fit_sigma_corr"]))
												   n <- length(list_c)
												   if (list_mean[n] > 25) {
												   						   q <- c(list_c[n],list_mean[n],list_sigma[n])
															 			   list_c[n] <- (q[1] * exp(-(25-q[2])**2/(2 * q[3]**2)))
												   						   }
												   if (list_mean[1] < 1) {
												   						  q <- c(list_c[1],list_mean[1],list_sigma[1])
															 			  list_c[1] <- (q[1] * exp(-(1-q[2])**2/(2 * q[3]**2)))
												   						  }
												   list_c <- as.data.frame(matrix(list_c,1, length(list_c), byrow=TRUE))
												   list_c
												   })												   

rnase_mean_corr$fit_mean_fxn <- apply(vect4, 2, function(z) {
												   list_c <- as.numeric(unlist(rnase_mean_corr[z,"fit_c_corr"]))
												   list_mean <- as.numeric(unlist(rnase_mean_corr[z,"fit_mean_corr"]))
												   list_sigma <- as.numeric(unlist(rnase_mean_corr[z,"fit_sigma_corr"]))
												   n <- length(list_c)
												   if (list_mean[n] > 25) {list_mean[n] <- 25}
												   if (list_mean[1] < 1) {list_mean[1] <- 1}
												   list_mean <- as.data.frame(matrix(list_mean,1, length(list_mean), byrow=TRUE))
												   list_mean
												   })
												   
rnase_mean_corr$fit_c_fxn <- apply(vect4, 2, function(z) {
												   list_c <- as.numeric(unlist(rnase_mean_corr[z,"fit_c_corr"]))
												   list_mean <- as.numeric(unlist(rnase_mean_corr[z,"fit_mean_corr"]))
												   list_sigma <- as.numeric(unlist(rnase_mean_corr[z,"fit_sigma_corr"]))
												   n <- length(list_c)
												   if (list_mean[n] > 25) {
												   						   q <- c(list_c[n],list_mean[n],list_sigma[n])
															 			   list_c[n] <- (q[1] * exp(-(25-q[2])**2/(2 * q[3]**2)))
															 			   }
												   if (list_mean[1] < 1) {
												   						  q <- c(list_c[1],list_mean[1],list_sigma[1])
															 			  list_c[1] <- (q[1] * exp(-(1-q[2])**2/(2 * q[3]**2)))
												   						  }
												   list_c <- as.data.frame(matrix(list_c,1, length(list_c), byrow=TRUE))
												   list_c
												   })


#***************************************************************************************************
# Assemble now all the information

# Parts from the original table
table_ctrl_mean <- ctrl_mean[1:25]
table_ctrl_mean$ctrl_max <- ctrl_mean$ctrl_max
table_ctrl_mean$nb_max <- ctrl_mean$nb_max
table_ctrl_mean$fit_c <- ctrl_mean$fit_c
table_ctrl_mean$fit_c_fxn <- ctrl_mean$fit_c_fxn
table_ctrl_mean$fit_mean <- ctrl_mean$fit_mean
table_ctrl_mean$fit_mean_fxn <- ctrl_mean$fit_mean_fxn
table_ctrl_mean$fit_sigma <- ctrl_mean$fit_sigma
table_ctrl_mean$fit_res <- ctrl_mean$fit_res
table_ctrl_mean$fit_area <- ctrl_mean$fit_area
table_ctrl_mean$sum_area <- ctrl_mean$sum_area
table_ctrl_mean$fitted <- ctrl_mean$fitted
table_ctrl_mean$pb_fit <- ctrl_mean$pb_fit

# Parts from the adjusted table
table_ctrl_mean$ctrl_max_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_ctrl_mean[x,])
															 max <- 1
															 if (dim(ctrl_mean_corr)[1] == 0){
															 max <- numeric()
															 } else {
															        max <- as.numeric(unlist(ctrl_mean_corr[pos,]$ctrl_max)*(table_ctrl_mean[pos,"pb_fit"] == "TRUE"))
															        }
															        max
															 })

table_ctrl_mean$fit_c_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_ctrl_mean[x,])
															 fit_c <- as.numeric(unlist(ctrl_mean_corr[pos,]$fit_c_corr)*(table_ctrl_mean[pos,"pb_fit"] == "TRUE"))
															 fit_c
															 })

table_ctrl_mean$fit_mean_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_ctrl_mean[x,])
															 fit_mean <- as.numeric(unlist(ctrl_mean_corr[pos,]$fit_mean_corr)*(table_ctrl_mean[pos,"pb_fit"] == "TRUE"))
															 fit_mean
															 })
															 
table_ctrl_mean$fit_sigma_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_ctrl_mean[x,])
															 fit_sigma <- as.numeric(unlist(ctrl_mean_corr[pos,]$fit_sigma_corr)*(table_ctrl_mean[pos,"pb_fit"] == "TRUE"))
															 fit_sigma
															 })
															 
table_ctrl_mean$fit_area_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_ctrl_mean[x,])
															 fit_area <- as.numeric(unlist(ctrl_mean_corr[pos,]$fit_area)*(table_ctrl_mean[pos,"pb_fit"] == "TRUE"))
															 fit_area
															 })															 

table_ctrl_mean$fit_mean_fxn_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_ctrl_mean[x,])
															 fit_mean_fxn <- as.numeric(unlist(ctrl_mean_corr[pos,]$fit_mean_fxn)*(table_ctrl_mean[pos,"pb_fit"] == "TRUE"))
															 fit_mean_fxn
															 })

table_ctrl_mean$fit_c_fxn_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_ctrl_mean[x,])
															 fit_c_fxn <- as.numeric(unlist(ctrl_mean_corr[pos,]$fit_c_fxn)*(table_ctrl_mean[pos,"pb_fit"] == "TRUE"))
															 fit_c_fxn
															 })

table_ctrl_mean$nb_max <- apply(vect, 2, function(x) {
												 		if (table_ctrl_mean[x,"pb_fit"] == "FALSE") {ls_max <- as.numeric(unlist(table_ctrl_mean[x,"fit_mean_fxn"])) }
												 		else {ls_max <- as.numeric(unlist(table_ctrl_mean[x,"fit_mean_fxn_corr"]))}
												 		
												 		n = length(ls_max)
												 		if (n == 0) {0} else {
												 							  if ((n == 1) && (ls_max == 0)) {0} else {n}
												 							  }
												 		})

# Same procedure for the RNase table
# Parts from the original table
table_rnase_mean <- rnase_mean[1:25]
table_rnase_mean$rnase_max <- rnase_mean$rnase_max
table_rnase_mean$nb_max <- rnase_mean$nb_max
table_rnase_mean$fit_c <- rnase_mean$fit_c
table_rnase_mean$fit_c_fxn <- rnase_mean$fit_c_fxn
table_rnase_mean$fit_mean <- rnase_mean$fit_mean
table_rnase_mean$fit_mean_fxn <- rnase_mean$fit_mean_fxn
table_rnase_mean$fit_sigma <- rnase_mean$fit_sigma
table_rnase_mean$fit_res <- rnase_mean$fit_res
table_rnase_mean$fit_area <- rnase_mean$fit_area
table_rnase_mean$sum_area <- rnase_mean$sum_area
table_rnase_mean$fitted <- rnase_mean$fitted
table_rnase_mean$pb_fit <- rnase_mean$pb_fit

# Parts from the adjusted table
table_rnase_mean$rnase_max_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_rnase_mean[x,])
															 max <- as.numeric(unlist(rnase_mean_corr[pos,]$rnase_max)*(table_rnase_mean[pos,"pb_fit"] == "TRUE"))
															 max
															 })

table_rnase_mean$fit_c_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_rnase_mean[x,])
															 fit_c <- as.numeric(unlist(rnase_mean_corr[pos,]$fit_c_corr)*(table_rnase_mean[pos,"pb_fit"] == "TRUE"))
															 fit_c
															 })

table_rnase_mean$fit_mean_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_rnase_mean[x,])
															 fit_mean <- as.numeric(unlist(rnase_mean_corr[pos,]$fit_mean_corr)*(table_rnase_mean[pos,"pb_fit"] == "TRUE"))
															 fit_mean
															 })
															 
table_rnase_mean$fit_sigma_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_rnase_mean[x,])
															 fit_sigma <- as.numeric(unlist(rnase_mean_corr[pos,]$fit_sigma_corr)*(table_rnase_mean[pos,"pb_fit"] == "TRUE"))
															 fit_sigma
															 })
															 
table_rnase_mean$fit_area_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_rnase_mean[x,])
															 fit_area <- as.numeric(unlist(rnase_mean_corr[pos,]$fit_area)*(table_rnase_mean[pos,"pb_fit"] == "TRUE"))
															 fit_area
															 })

table_rnase_mean$fit_mean_fxn_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_rnase_mean[x,])
															 fit_mean_fxn <- as.numeric(unlist(rnase_mean_corr[pos,]$fit_mean_fxn)*(table_rnase_mean[pos,"pb_fit"] == "TRUE"))
															 fit_mean_fxn
															 })

table_rnase_mean$fit_c_fxn_corr <- apply(vect, 2, function(x) {
															 pos <- rownames(table_rnase_mean[x,])
															 fit_c_fxn <- as.numeric(unlist(rnase_mean_corr[pos,]$fit_c_fxn)*(table_rnase_mean[pos,"pb_fit"] == "TRUE"))
															 fit_c_fxn
															 })
															 
table_rnase_mean$nb_max <- apply(vect, 2, function(x) {
												 		if (table_rnase_mean[x,"pb_fit"] == "FALSE") {ls_max <- as.numeric(unlist(table_rnase_mean[x,"fit_mean_fxn"])) }
												 		else {ls_max <- as.numeric(unlist(table_rnase_mean[x,"fit_mean_fxn_corr"]))}
												 		
												 		n = length(ls_max)
												 		if (n == 0) {0} else {
												 							  if ((n == 1) && (ls_max == 0)) {0} else {n}
												 							  }
												 		})


#****************************************************************
#### PART 5 : Obtain p_values for all maxima                 ####
#****************************************************************

# calculate the p_values for all maxima
# including those with high res values,
# negative areas and null sigma (all previously corrected in PART 4).

# Add a new column with p-values for the maxima in the CTRL table "table_ctrl_mean".
# The p-values are calculated protein per protein
# Reminder: "vect" is a matrix of dimension 1 and 1000

table_ctrl_mean$p_values <- apply(vect, 2, function(z) {
						# get the list of maxima and the corrected maxima if they exist (pb_fit = TRUE)
						list_max <- numeric(0)
						if (table_ctrl_mean[z,"pb_fit"] == "TRUE") { list_max <- as.numeric(unlist(table_ctrl_mean[z,"fit_mean_fxn_corr"]))	} else {list_max <- as.numeric(unlist(table_ctrl_mean[z,"fit_mean_fxn"]))}
												  											  
						n <- (length(list_max))
												  
						# get the residue and sigma values for ctrl maxima
						res <- as.numeric(unlist(table_ctrl_mean[z, "fit_res"]))
						sigma <- as.numeric(unlist(table_ctrl_mean[z,"fit_sigma"]))
						nl <- any(sigma == 0) && (res != 0)
						high_res <- (res > 260)
						p_value <- numeric(0)
												  
						if (sum(list_max) == 0) {p_value <- 1} else {
												  						                  # control the quality of the fit in the CTRL sample (mean curve)
						                                            if (high_res == "TRUE" || nl == "TRUE") {
												  																	for (i in 1:n) {
												  																	fit_mean <- list_max[i]
												  																	if (fit_mean > 25) {p_value <- c(p_value,1)} else {
												  																	                           # control the quality of the fit in the RNase sample (mean curve)
												  																														 res_R <- as.numeric(unlist(table_rnase_mean[z, "fit_res"]))
												  																														 sigma_R <- as.numeric(unlist(table_rnase_mean[z,"fit_sigma"]))
												  																														 nl_R <- any(sigma_R == 0) && (res_R != 0)
												  																														 high_res_R <- (res_R > 260)
												  																		                         if (high_res_R == "TRUE" || nl_R == "TRUE") {
												  																																										fit_c_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_c_fxn_corr"]))
												  																																										fit_mean_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_mean_fxn_corr"]))
												  																																										fit_sigma_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_sigma_corr"]))
												  																																										fit_c <- fit_c_list[fit_mean_list == fit_mean]
												  																																										fit_sigma <- fit_sigma_list[fit_mean_list == fit_mean]
												  																																										vector <- c(fit_c, fit_mean, fit_sigma)
												  																																										# Call the function "fine_fit_corr" to access the corrected fit values
												  																																										vc1 <- fine_fit_corr(z,vector,table.ctrl1.SW.norm)
												  																																										vc2 <- fine_fit_corr(z,vector,table.ctrl2.SW.norm)
												  																																										vc3 <- fine_fit_corr(z,vector,table.ctrl3.SW.norm)
												  																																										vr1 <- fine_fit_corr(z,vector,table.rnase1.SW.norm)
												  																																										vr2 <- fine_fit_corr(z,vector,table.rnase2.SW.norm)
												  																																										vr3 <- fine_fit_corr(z,vector,table.rnase3.SW.norm)
												  																																										x <- c(vc1, vc2, vc3, vr1, vr2, vr3)
												  																		
												  																																										if ( all(abs(x[1:3]-x[1])<1e-13, na.rm=TRUE) ) {x[1:3] <- c(x[1], x[1]+0.0001, x[1])}
												  																																										if ( all(abs(x[4:6]-x[4])<1e-13, na.rm=TRUE) ) {x[4:6] <- c(x[4], x[4]+0.0001, x[4])}
												  																		
												  																																										ftest <- if ( sum(is.na(x))>1 )  {1} else {var.test(x[1:3],x[4:6],na.action=na.omit)$p.value}
												  																																										pv <- if ( sum(is.na(x))>1 )  {1} else if ( all(abs(x[1:6]-x[1])<1e-13, na.rm=TRUE) ) {1} else {t.test(x[1:3],x[4:6],paired=FALSE, na.action=na.omit, var.equal = (ftest >0.05) )$p.value}
												  																																										p_value <- c(p_value, pv)
												  																																										
												  																																										} else {
												  																																													fit_c_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_c_fxn_corr"]))
												  																																													fit_mean_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_mean_fxn_corr"]))
												  																																													fit_sigma_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_sigma_corr"]))
												  																																													fit_c <- fit_c_list[fit_mean_list == fit_mean]
												  																																													fit_sigma <- fit_sigma_list[fit_mean_list == fit_mean]
												  																																													vector <- c(fit_c, fit_mean, fit_sigma)
												  																																													vc1 <- fine_fit_corr(z,vector,table.ctrl1.SW.norm)
												  																																													vc2 <- fine_fit_corr(z,vector,table.ctrl2.SW.norm)
												  																																													vc3 <- fine_fit_corr(z,vector,table.ctrl3.SW.norm)
												  																																													vr1 <- table.rnase1.fit.fine[z,(list_max[i])/0.1+1]
												  																																													vr2 <- table.rnase2.fit.fine[z,(list_max[i])/0.1+1]
												  																																													vr3 <- table.rnase3.fit.fine[z,(list_max[i])/0.1+1]
												  																																													x <- c(vc1, vc2, vc3, vr1, vr2, vr3)
												  																		
												  																																													if ( all(abs(x[1:3]-x[1])<1e-13, na.rm=TRUE) ) {x[1:3] <- c(x[1], x[1]+0.0001, x[1])}
												  																																													if ( all(abs(x[4:6]-x[4])<1e-13, na.rm=TRUE) ) {x[4:6] <- c(x[4], x[4]+0.0001, x[4])}
												  																		
												  																																													ftest <- if ( sum(is.na(x))>1 )  {1} else {var.test(x[1:3],x[4:6],na.action=na.omit)$p.value}
												  																																													pv <- if ( sum(is.na(x))>1 )  {1} else if ( all(abs(x[1:6]-x[1])<1e-13, na.rm=TRUE) ) {1} else {t.test(x[1:3],x[4:6],paired=FALSE, na.action=na.omit, var.equal = (ftest >0.05) )$p.value}
												  																																													p_value <- c(p_value, pv)
																										 																															} #end if else loop
																										 																} #end if else loop									  						
												  														    					} #end for loop
												  												 				} else {
												  										
												  																			for (i in 1:n) {
												  																								fit_mean <- list_max[i]
												  																								if (fit_mean > 25) {p_value <- c(p_value,1)} else {
												  																																						res_R <- as.numeric(unlist(table_rnase_mean[z, "fit_res"]))
												  																																						sigma_R <- as.numeric(unlist(table_rnase_mean[z,"fit_sigma"]))
												  																																						nl_R <- any(sigma_R == 0) && (res_R != 0)
												  																																						high_res_R <- (res_R > 260)
												  																		
												  																																						if (high_res_R == "TRUE" || nl_R == "TRUE") {
												  																																																		fit_c_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_c_fxn"]))
												  																																																		fit_mean_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_mean_fxn"]))
												  																																																		fit_sigma_list <- as.numeric(unlist(table_ctrl_mean[z,"fit_sigma"]))
												  																																																		fit_c <- fit_c_list[fit_mean_list == fit_mean]
												  																																																		fit_sigma <- fit_sigma_list[fit_mean_list == fit_mean]
												  																																																		vector <- c(fit_c, fit_mean, fit_sigma)
												  																																																		vc1 <- table.ctrl1.fit.fine[z,(list_max[i])/0.1+1]
												  																																																		vc2 <- table.ctrl2.fit.fine[z,(list_max[i])/0.1+1]
												  																																																		vc3 <- table.ctrl3.fit.fine[z,(list_max[i])/0.1+1]
												  																																																		vr1 <- fine_fit_corr(z,vector,table.rnase1.SW.norm)
												  																																																		vr2 <- fine_fit_corr(z,vector,table.rnase2.SW.norm)
												  																																																		vr3 <- fine_fit_corr(z,vector,table.rnase3.SW.norm)
												  																																																		x <- c(vc1, vc2, vc3, vr1, vr2, vr3)
												  																		
												  																																																		if ( all(abs(x[1:3]-x[1])<1e-13, na.rm=TRUE) ) {x[1:3] <- c(x[1], x[1]+0.0001, x[1])}
												  																																																		if ( all(abs(x[4:6]-x[4])<1e-13, na.rm=TRUE) ) {x[4:6] <- c(x[4], x[4]+0.0001, x[4])}
												  																		
												  																																																		ftest <- if ( sum(is.na(x))>1 )  {1} else {var.test(x[1:3],x[4:6],na.action=na.omit)$p.value}
												  																																																		pv <- if ( sum(is.na(x))>1 )  {1} else if ( all(abs(x[1:6]-x[1])<1e-13, na.rm=TRUE) ) {1} else {t.test(x[1:3],x[4:6],paired=FALSE, na.action=na.omit, var.equal = (ftest >0.05) )$p.value}
												  																																																		p_value <- c(p_value, pv)
												  																																																		
												  																																																		} else {
												  																																																					vc1 <- table.ctrl1.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vc2 <- table.ctrl2.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vc3 <- table.ctrl3.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vr1 <- table.rnase1.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vr2 <- table.rnase2.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vr3 <- table.rnase3.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					x <- c(vc1, vc2, vc3, vr1, vr2, vr3)
												  																		
												  																																																					if ( all(abs(x[1:3]-x[1])<1e-13, na.rm=TRUE) ) {x[1:3] <- c(x[1], x[1]+0.0001, x[1])}
												  																																																					if ( all(abs(x[4:6]-x[4])<1e-13, na.rm=TRUE) ) {x[4:6] <- c(x[4], x[4]+0.0001, x[4])}
												  																		
												  																																																					ftest <- if ( sum(is.na(x))>1 )  {1} else {var.test(x[1:3],x[4:6],na.action=na.omit)$p.value}
												  																																																					pv <- if ( sum(is.na(x))>1 )  {1} else if ( all(abs(x[1:6]-x[1])<1e-13, na.rm=TRUE) ) {1} else {t.test(x[1:3],x[4:6],paired=FALSE, na.action=na.omit, var.equal = (ftest >0.05) )$p.value}
												  																																																					p_value <- c(p_value, pv)
												  													    																																						} #end if else loop
												  													    																							} #end if else loop
												  													     									} #end for loop
												  											 	 						} #end if else loop
												  					} #end if else loop						 	 						
												  p_value
												  })
												  
																						  
# for RNASE matrix
table_rnase_mean$p_values <- apply(vect, 2, function(z) {
						# get the list of maxima
						list_max <- numeric(0)
						if (table_rnase_mean[z,"pb_fit"] == "TRUE") { list_max <- as.numeric(unlist(table_rnase_mean[z,"fit_mean_fxn_corr"]))	} else {list_max <- as.numeric(unlist(table_rnase_mean[z,"fit_mean_fxn"]))}
												  											  
						n <- (length(list_max))
												  
						# get the residue and sigma values for rnase maxima
						res <- as.numeric(unlist(table_rnase_mean[z, "fit_res"]))
						sigma <- as.numeric(unlist(table_rnase_mean[z,"fit_sigma"]))
						nl <- any(sigma == 0) && (res != 0)
						high_res <- (res > 260)
						p_value <- numeric(0)
												  
						if (sum(list_max) == 0) {p_value <- 1} else {
						                                            # control the quality of the fit in the RNase sample (mean curve)					  						
						                                            if (high_res == "TRUE" || nl == "TRUE") {
												  																	for (i in 1:n) {
												  																	fit_mean <- list_max[i]
												  																	if (fit_mean > 25) {p_value <- c(p_value,1)} else {
												  																														 res_C <- as.numeric(unlist(table_ctrl_mean[z, "fit_res"]))
												  																														 sigma_C <- as.numeric(unlist(table_ctrl_mean[z,"fit_sigma"]))
												  																														 nl_C <- any(sigma_C == 0) && (res_C != 0)
												  																														 high_res_C <- (res_C > 260)
												  																		
												  																														 # control the quality of the fit in the CTRL sample (mean curve)	
												  																														 if (high_res_C == "TRUE" || nl_C == "TRUE") {
												  																																										fit_c_list <- as.numeric(unlist(table_rnase_mean[z,"fit_c_fxn_corr"]))
												  																																										fit_mean_list <- as.numeric(unlist(table_rnase_mean[z,"fit_mean_fxn_corr"]))
												  																																										fit_sigma_list <- as.numeric(unlist(table_rnase_mean[z,"fit_sigma_corr"]))
												  																																										fit_c <- fit_c_list[fit_mean_list == fit_mean]
												  																																										fit_sigma <- fit_sigma_list[fit_mean_list == fit_mean]
												  																																										vector <- c(fit_c, fit_mean, fit_sigma)
												  																																										vc1 <- fine_fit_corr(z,vector,table.ctrl1.SW.norm)
												  																																										vc2 <- fine_fit_corr(z,vector,table.ctrl2.SW.norm)
												  																																										vc3 <- fine_fit_corr(z,vector,table.ctrl3.SW.norm)
												  																																										vr1 <- fine_fit_corr(z,vector,table.rnase1.SW.norm)
												  																																										vr2 <- fine_fit_corr(z,vector,table.rnase2.SW.norm)
												  																																										vr3 <- fine_fit_corr(z,vector,table.rnase3.SW.norm)
												  																																										x <- c(vc1, vc2, vc3, vr1, vr2, vr3)
												  																		
												  																																										if ( all(abs(x[1:3]-x[1])<1e-13, na.rm=TRUE) ) {x[1:3] <- c(x[1], x[1]+0.0001, x[1])}
												  																																										if ( all(abs(x[4:6]-x[4])<1e-13, na.rm=TRUE) ) {x[4:6] <- c(x[4], x[4]+0.0001, x[4])}
												  																		
												  																																										ftest <- if ( sum(is.na(x))>1 )  {1} else {var.test(x[1:3],x[4:6],na.action=na.omit)$p.value}
												  																																										pv <- if ( sum(is.na(x))>1 )  {1} else if ( all(abs(x[1:6]-x[1])<1e-13, na.rm=TRUE) ) {1} else {t.test(x[1:3],x[4:6],paired=FALSE, na.action=na.omit, var.equal = (ftest >0.05) )$p.value}
												  																																										p_value <- c(p_value, pv)
												  																																										
												  																																										} else {
												  																																													fit_c_list <- as.numeric(unlist(table_rnase_mean[z,"fit_c_fxn_corr"]))
												  																																													fit_mean_list <- as.numeric(unlist(table_rnase_mean[z,"fit_mean_fxn_corr"]))
												  																																													fit_sigma_list <- as.numeric(unlist(table_rnase_mean[z,"fit_sigma_corr"]))
												  																																													fit_c <- fit_c_list[fit_mean_list == fit_mean]
												  																																													fit_sigma <- fit_sigma_list[fit_mean_list == fit_mean]
												  																																													vector <- c(fit_c, fit_mean, fit_sigma)
												  																																													vr1 <- fine_fit_corr(z,vector,table.rnase1.SW.norm)
												  																																													vr2 <- fine_fit_corr(z,vector,table.rnase2.SW.norm)
												  																																													vr3 <- fine_fit_corr(z,vector,table.rnase3.SW.norm)
												  																																													vc1 <- table.ctrl1.fit.fine[z,(list_max[i])/0.1+1]
												  																																													vc2 <- table.ctrl2.fit.fine[z,(list_max[i])/0.1+1]
												  																																													vc3 <- table.ctrl3.fit.fine[z,(list_max[i])/0.1+1]
												  																																													x <- c(vc1, vc2, vc3, vr1, vr2, vr3)
												  																		
												  																																													if ( all(abs(x[1:3]-x[1])<1e-13, na.rm=TRUE) ) {x[1:3] <- c(x[1], x[1]+0.0001, x[1])}
												  																																													if ( all(abs(x[4:6]-x[4])<1e-13, na.rm=TRUE) ) {x[4:6] <- c(x[4], x[4]+0.0001, x[4])}
												  																		
												  																																													ftest <- if ( sum(is.na(x))>1 )  {1} else {var.test(x[1:3],x[4:6],na.action=na.omit)$p.value}
												  																																													pv <- if ( sum(is.na(x))>1 )  {1} else if ( all(abs(x[1:6]-x[1])<1e-13, na.rm=TRUE) ) {1} else {t.test(x[1:3],x[4:6],paired=FALSE, na.action=na.omit, var.equal = (ftest >0.05) )$p.value}
												  																																													p_value <- c(p_value, pv)
																										 																															} #end if else loop
																										 																} #end if else loop									  						
												  														    					} #end for loop
												  												 				} else {
												  										
												  																			for (i in 1:n) {
												  																								fit_mean <- list_max[i]
												  																								if (fit_mean > 25) {p_value <- c(p_value,1)} else {
												  																																						res_C <- as.numeric(unlist(table_ctrl_mean[z, "fit_res"]))
												  																																						sigma_C <- as.numeric(unlist(table_ctrl_mean[z,"fit_sigma"]))
												  																																						nl_C <- any(sigma_C == 0) && (res_C != 0)
												  																																						high_res_C <- (res_C > 260)
												  																		
												  																																						if (high_res_C == "TRUE" || nl_C == "TRUE") {
												  																																																		fit_c_list <- as.numeric(unlist(table_rnase_mean[z,"fit_c_fxn"]))
												  																																																		fit_mean_list <- as.numeric(unlist(table_rnase_mean[z,"fit_mean_fxn"]))
												  																																																		fit_sigma_list <- as.numeric(unlist(table_rnase_mean[z,"fit_sigma"]))
												  																																																		fit_c <- fit_c_list[fit_mean_list == fit_mean]
												  																																																		fit_sigma <- fit_sigma_list[fit_mean_list == fit_mean]
												  																																																		vector <- c(fit_c, fit_mean, fit_sigma)
												  																																																		vr1 <- table.rnase1.fit.fine[z,(list_max[i])/0.1+1]
												  																																																		vr2 <- table.rnase2.fit.fine[z,(list_max[i])/0.1+1]
												  																																																		vr3 <- table.rnase3.fit.fine[z,(list_max[i])/0.1+1]
												  																																																		vc1 <- fine_fit_corr(z,vector,table.ctrl1.SW.norm)
												  																																																		vc2 <- fine_fit_corr(z,vector,table.ctrl2.SW.norm)
												  																																																		vc3 <- fine_fit_corr(z,vector,table.ctrl3.SW.norm)
												  																																																		x <- c(vc1, vc2, vc3, vr1, vr2, vr3)
												  																		
												  																																																		if ( all(abs(x[1:3]-x[1])<1e-13, na.rm=TRUE) ) {x[1:3] <- c(x[1], x[1]+0.0001, x[1])}
												  																																																		if ( all(abs(x[4:6]-x[4])<1e-13, na.rm=TRUE) ) {x[4:6] <- c(x[4], x[4]+0.0001, x[4])}
												  																		
												  																																																		ftest <- if ( sum(is.na(x))>1 )  {1} else {var.test(x[1:3],x[4:6],na.action=na.omit)$p.value}
												  																																																		pv <- if ( sum(is.na(x))>1 )  {1} else if ( all(abs(x[1:6]-x[1])<1e-13, na.rm=TRUE) ) {1} else {t.test(x[1:3],x[4:6],paired=FALSE, na.action=na.omit, var.equal = (ftest >0.05) )$p.value}
												  																																																		p_value <- c(p_value, pv)
												  																																																		
												  																																																		} else {
												  																																																					vc1 <- table.ctrl1.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vc2 <- table.ctrl2.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vc3 <- table.ctrl3.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vr1 <- table.rnase1.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vr2 <- table.rnase2.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					vr3 <- table.rnase3.fit.fine[z,(list_max[i])/0.1+1]
												  																																																					x <- c(vc1, vc2, vc3, vr1, vr2, vr3)
												  																		
												  																																																					if ( all(abs(x[1:3]-x[1])<1e-13, na.rm=TRUE) ) {x[1:3] <- c(x[1], x[1]+0.0001, x[1])}
												  																																																					if ( all(abs(x[4:6]-x[4])<1e-13, na.rm=TRUE) ) {x[4:6] <- c(x[4], x[4]+0.0001, x[4])}
												  																		
												  																																																					ftest <- if ( sum(is.na(x))>1 )  {1} else {var.test(x[1:3],x[4:6],na.action=na.omit)$p.value}
												  																																																					pv <- if ( sum(is.na(x))>1 )  {1} else if ( all(abs(x[1:6]-x[1])<1e-13, na.rm=TRUE) ) {1} else {t.test(x[1:3],x[4:6],paired=FALSE, na.action=na.omit, var.equal = (ftest >0.05) )$p.value}
												  																																																					p_value <- c(p_value, pv)
												  													    																																						} #end if else loop
												  													    																							} #end if else loop
												  													     									} #end for loop
												  											 	 						} #end if else loop
												  					} #end if else loop						 	 						
												  p_value
												  })

#***************************************************************************************************
# adjusted p-values

table_ctrl_pv <- unlist(table_ctrl_mean$p_values)
adj_table_ctrl_pv <- p.adjust(table_ctrl_pv, method="fdr")

table_ctrl_mean$p_values_fdr <- apply(vect, 2, function(x) {
													  
													  n <- sum(table_ctrl_mean$nb_max[1:x]) + sum(table_ctrl_mean$nb_max[1:x] == 0)
													  m <- table_ctrl_mean[x, "nb_max"] + (table_ctrl_mean[x, "nb_max"] == 0)*1
													  list_adj_pv <- adj_table_ctrl_pv[1:n]
													  adj_pv <- list_adj_pv[(n-m+1):n]
													  adj_pv
													  })

table_rnase_pv <- unlist(table_rnase_mean$p_values)
adj_table_rnase_pv <- p.adjust(table_rnase_pv, method="fdr")

table_rnase_mean$p_values_fdr <- apply(vect, 2, function(x) {
													  n <- sum(table_rnase_mean$nb_max[1:x]) + sum(table_rnase_mean$nb_max[1:x] == 0)
													  m <- table_rnase_mean[x, "nb_max"] + (table_rnase_mean[x, "nb_max"] == 0)*1
													  list_adj_pv <- adj_table_rnase_pv[1:n]
													  adj_pv <- list_adj_pv[(n-m+1):n]
													  adj_pv
													  })



#**************************************************************************************************
#### PART 6 : Evaluation of the shifts and selection of significant shifts                     ####
#**************************************************************************************************

#****************************************************************
# Select peaks with amplitude >= 25% of amplitude of global maxima (max_th)

table_ctrl_mean$max_th <- apply(vect, 2, function(z) {
												if (table_ctrl_mean[z,"pb_fit"]) {
															 list_c <- as.numeric(unlist(table_ctrl_mean[z,"fit_c_fxn_corr"]))
															 list_mean <- as.numeric(unlist(table_ctrl_mean[z,"fit_mean_fxn_corr"]))
															 n <- length(list_c)
															 glob_max_pos <- which.max(list_c)
															 thres <- list_c[glob_max_pos]/4
															 list_mean <- list_mean[list_c>=thres]
															 list_mean <- list_mean[list_mean >1]
															 list_mean <- list_mean[list_mean <= 25]
															 if (n == 0) {0} else {list_mean}
															 } else {
												list_c <- as.numeric(unlist(table_ctrl_mean[z,"fit_c_fxn"]))
												list_mean <- as.numeric(unlist(table_ctrl_mean[z,"fit_mean_fxn"]))
												n <- length(list_c)
												glob_max_pos <- which.max(list_c)
												thres <- list_c[glob_max_pos]/4
												list_mean <- list_mean[list_c>=thres]
												n <- length(list_mean)
												
												# Get rid of peaks which are next to each other (ex 18.4 and 18.6 - can happen during fitting with Gaussian). At least 1 fraction in between. Remove one peak of the list.
												if (n>1) {
															to_keep <- c("TRUE", (abs(list_mean[2:n]-list_mean[1:(n-1)])>=1))
															list_mean <- list_mean[to_keep == "TRUE"]
														}	
												
												if (n == 0) {0} else {list_mean}
																	 }
												})


table_rnase_mean$max_th <- apply(vect, 2, function(z) {
												if (table_rnase_mean[z,"pb_fit"]) {
															 list_c <- as.numeric(unlist(table_rnase_mean[z,"fit_c_fxn_corr"]))
															 list_mean <- as.numeric(unlist(table_rnase_mean[z,"fit_mean_fxn_corr"]))
															 n <- length(list_c)
															 glob_max_pos <- which.max(list_c)
															 thres <- list_c[glob_max_pos]/4
															 list_mean <- list_mean[list_c>=thres]
															 list_mean <- list_mean[list_mean >1]
															 list_mean <- list_mean[list_mean <= 25]
															 if (n == 0) {0} else {list_mean}
															 } else {
												list_c <- as.numeric(unlist(table_rnase_mean[z,"fit_c_fxn"]))
												list_mean <- as.numeric(unlist(table_rnase_mean[z,"fit_mean_fxn"]))
												n <- length(list_c)
												glob_max_pos <- which.max(list_c)
												thres <- list_c[glob_max_pos]/4
												list_mean <- list_mean[list_c>=thres]
												n <- length(list_mean)
												
												# Get rid of peaks which are next to each other (ex 18.4 and 18.6 - can happen during fitting with Gaussian). At least 1 fraction in between. Remove one peak of the list.
												if (n>1) {
															to_keep <- c("TRUE", (abs(list_mean[2:n]-list_mean[1:(n-1)])>=1))
															list_mean <- list_mean[to_keep == "TRUE"]
														}
												
												if (n == 0) {0} else {list_mean}
																	 }
												})

# nb of peaks >= threshold

table_ctrl_mean$nb_max_th <- apply(vect, 2, function(x) {
												 		ls_max <- as.numeric(unlist(table_ctrl_mean[x,"max_th"]))
												 		n = length(ls_max)
												 		if (n == 0) {0} else {
												 							  if ((n == 1) && (ls_max == 0)) {0} else {n}
												 							  }
												 		})

table_rnase_mean$nb_max_th <- apply(vect, 2, function(x) {
												 		ls_max <- as.numeric(unlist(table_rnase_mean[x,"max_th"]))
												 		n = length(ls_max)
												 		if (n == 0) {0} else {
												 							  if ((n == 1) && (ls_max == 0)) {0} else {n}
												 							  }
												 		})


#**************************************************************************************************
# function to obtain a fine fit (0.1 precision) of the ctrl_mean curve for each protein.

fine_fit_fun_ctrl <- function(t) { 
												  dffit <- data.frame(y=round(seq(1, 25, 0.1),digits=1))
												  list_c <- as.numeric(unlist(ctrl_mean[t,"fit_c"]))
												  list_mean <- as.numeric(unlist(ctrl_mean[t,"fit_mean"]))
												  list_sigma <- as.numeric(unlist(ctrl_mean[t,"fit_sigma"]))
												  n <- length(list_c)
												  z <- round(seq(1, 25, 0.1),digits=1)
												  
												  if (rnase_mean[t,"fitted"] == "FALSE") {dffit$df.y <- rep(0,241)} else {
												  											if (n==1) {
												  														q <- c(list_c[1],list_mean[1],list_sigma[1])
												  														C1 <- q[1]
																										mean1 <- q[2]
																										sigma1 <- q[3]
																										dffit$df.y <- C1 * exp(-(z-mean1)**2/(2 * sigma1**2))
												  														} else {
												  																if (n==2) {
												  																			q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2])
												  																			C1 <- q[1]
																															mean1 <- q[2]
																															sigma1 <- q[3]
																															C2 <- q[4]
																															mean2 <- q[5]
																															sigma2 <- q[6]
																															dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) ) 
												  																			} else {
												  																					if (n==3) {
												  																								q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3])
												  																								C1 <- q[1]
																																				mean1 <- q[2]
																																				sigma1 <- q[3]
																																				C2 <- q[4]
																																				mean2 <- q[5]
																																				sigma2 <- q[6]
																																				C3 <- q[7]
																																				mean3 <- q[8]
																																				sigma3 <- q[9]
																																				dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) ) 
												  																								} else {
												  																										if (n==4) {
												  																													q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4])
												  																													C1 <- q[1]
																																									mean1 <- q[2]
																																									sigma1 <- q[3]
																																									C2 <- q[4]
																																									mean2 <- q[5]
																																									sigma2 <- q[6]
																																									C3 <- q[7]
																																									mean3 <- q[8]
																																									sigma3 <- q[9]
																																									C4 <- q[10]
																																									mean4 <- q[11]
																																									sigma4 <- q[12]
																																									dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) )
												  																													} else {
												  																															if (n==5) {
												  																																		q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5])
												  																																		C1 <- q[1]
																																														mean1 <- q[2]
																																														sigma1 <- q[3]
																																														C2 <- q[4]
																																														mean2 <- q[5]
																																														sigma2 <- q[6]
																																														C3 <- q[7]
																																														mean3 <- q[8]
																																														sigma3 <- q[9]
																																														C4 <- q[10]
																																														mean4 <- q[11]
																																														sigma4 <- q[12]
																																														C5 <- q[13]
																																														mean5 <- q[14]
																																														sigma5 <- q[15]
																																														dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) )
												  																																		} else {
												  																																				q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5],list_c[6],list_mean[6],list_sigma[6])
												  																																				C1 <- q[1]
																																																mean1 <- q[2]
																																																sigma1 <- q[3]
																																																C2 <- q[4]
																																																mean2 <- q[5]
																																																sigma2 <- q[6]
																																																C3 <- q[7]
																																																mean3 <- q[8]
																																																sigma3 <- q[9]
																																																C4 <- q[10]
																																																mean4 <- q[11]
																																																sigma4 <- q[12]
																																																C5 <- q[13]
																																																mean5 <- q[14]
																																																sigma5 <- q[15]
																																																C6 <- q[16]
																																																mean6 <- q[17]
																																																sigma6 <- q[18]
																																																dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) + C6 * exp(-(z-mean6)**2/(2 * sigma6**2)) ) 
												  																																				}
												  																															}
												  																										}
												  																					}
												  																}
												  											}
											dffit
											}

#**************************************************************************************************
# function to obtain a fine fit (0.1 precision) of the rnase_mean curve for each protein.

fine_fit_fun_rnase <- function(t) { 
												  dffit <- data.frame(y=round(seq(1, 25, 0.1),digits=1))
												  list_c <- as.numeric(unlist(rnase_mean[t,"fit_c"]))
												  list_mean <- as.numeric(unlist(rnase_mean[t,"fit_mean"]))
												  list_sigma <- as.numeric(unlist(rnase_mean[t,"fit_sigma"]))
												  n <- length(list_c)
												  z <- round(seq(1, 25, 0.1),digits=1)
												  
												  if (rnase_mean[t,"fitted"] == "FALSE") {dffit$df.y <- rep(0,241)} else {
												  											if (n==1) {
												  														q <- c(list_c[1],list_mean[1],list_sigma[1])
												  														C1 <- q[1]
																										mean1 <- q[2]
																										sigma1 <- q[3]
																										dffit$df.y <- C1 * exp(-(z-mean1)**2/(2 * sigma1**2))
												  														} else {
												  																if (n==2) {
												  																			q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2])
												  																			C1 <- q[1]
																															mean1 <- q[2]
																															sigma1 <- q[3]
																															C2 <- q[4]
																															mean2 <- q[5]
																															sigma2 <- q[6]
																															dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) ) 
												  																			} else {
												  																					if (n==3) {
												  																								q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3])
												  																								C1 <- q[1]
																																				mean1 <- q[2]
																																				sigma1 <- q[3]
																																				C2 <- q[4]
																																				mean2 <- q[5]
																																				sigma2 <- q[6]
																																				C3 <- q[7]
																																				mean3 <- q[8]
																																				sigma3 <- q[9]
																																				dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) ) 
												  																								} else {
												  																										if (n==4) {
												  																													q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4])
												  																													C1 <- q[1]
																																									mean1 <- q[2]
																																									sigma1 <- q[3]
																																									C2 <- q[4]
																																									mean2 <- q[5]
																																									sigma2 <- q[6]
																																									C3 <- q[7]
																																									mean3 <- q[8]
																																									sigma3 <- q[9]
																																									C4 <- q[10]
																																									mean4 <- q[11]
																																									sigma4 <- q[12]
																																									dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) )
												  																													} else {
												  																															if (n==5) {
												  																																		q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5])
												  																																		C1 <- q[1]
																																														mean1 <- q[2]
																																														sigma1 <- q[3]
																																														C2 <- q[4]
																																														mean2 <- q[5]
																																														sigma2 <- q[6]
																																														C3 <- q[7]
																																														mean3 <- q[8]
																																														sigma3 <- q[9]
																																														C4 <- q[10]
																																														mean4 <- q[11]
																																														sigma4 <- q[12]
																																														C5 <- q[13]
																																														mean5 <- q[14]
																																														sigma5 <- q[15]
																																														dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) )
												  																																		} else {
												  																																				q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5],list_c[6],list_mean[6],list_sigma[6])
												  																																				C1 <- q[1]
																																																mean1 <- q[2]
																																																sigma1 <- q[3]
																																																C2 <- q[4]
																																																mean2 <- q[5]
																																																sigma2 <- q[6]
																																																C3 <- q[7]
																																																mean3 <- q[8]
																																																sigma3 <- q[9]
																																																C4 <- q[10]
																																																mean4 <- q[11]
																																																sigma4 <- q[12]
																																																C5 <- q[13]
																																																mean5 <- q[14]
																																																sigma5 <- q[15]
																																																C6 <- q[16]
																																																mean6 <- q[17]
																																																sigma6 <- q[18]
																																																dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) + C6 * exp(-(z-mean6)**2/(2 * sigma6**2)) ) 
												  																																				}
												  																															}
												  																										}
												  																					}
												  																}
												  											}
											dffit
											}


# Function to obtain the amplitude of the CTRL and RNASE fit curves at each maximum for each protein (1 position, 1 value in each curve, 0.1 precision)
ampl_max <- function(max,pos) {
								value <- numeric(0)
								
								res_C <- as.numeric(unlist(table_ctrl_mean[pos, "fit_res"]))
								sigma_C <- as.numeric(unlist(table_ctrl_mean[pos,"fit_sigma"]))
								nl_C <- any(sigma_C == 0) && (res_C != 0)
								
								res_R <- as.numeric(unlist(table_rnase_mean[pos, "fit_res"]))
								sigma_R <- as.numeric(unlist(table_rnase_mean[pos,"fit_sigma"]))
								nl_R <- any(sigma_R == 0) && (res_R != 0)
								
								if ((table_ctrl_mean[pos, "fit_res"] > 260 || nl_C) && (table_rnase_mean[pos, "fit_res"] <= 260) && nl_R == "FALSE") {
							  				fit_c_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_c_corr"]))
							  				fit_mean_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_mean_corr"]))
							  				fit_sigma_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_sigma_corr"]))
							  				n <- length(fit_c_ctrl)
							  				Fractions <- round(seq(1,25,0.1),digits=1)
							  				dfpb_ctrl <- as.data.frame(Fractions)
							  				dfpb <- as.data.frame(Fractions)
							  				dfpb$Am <- rep(0,241)
							  				dfpb_ctrl$Amount <- rep(0,241)
							  				
							  				for (i in 1:n) {
							  				dfpb$Am <- ( fit_c_ctrl[i] * exp(-(Fractions-fit_mean_ctrl[i])**2/(2 * fit_sigma_ctrl[i]**2)) )
							  				dfpb$Amount <- dfpb_ctrl$Amount
							  				dfpb_ctrl$Amount <- apply(dfpb,1, function(x) {max(x[2:3])})
							  								}
							  				
							  				dffit_ctrl <- dfpb_ctrl
							  				dffit_rnase <- fine_fit_fun_rnase(pos) 
											  colnames(dffit_rnase) <- c("Fractions","Amount")
							  				
							  				} else {
							  						if ((table_ctrl_mean[pos, "fit_res"] <= 260 && nl_C == "FALSE") && (table_rnase_mean[pos, "fit_res"] > 260 || nl_R)) {
							  									fit_c_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_c_corr"]))
							  									fit_mean_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_mean_corr"]))
							  									fit_sigma_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_sigma_corr"]))
							  									n <- length(fit_c_rnase)
							  									Fractions <- round(seq(1,25,0.1),digits=1)
							  									dfpb_rnase <- as.data.frame(Fractions)
							  									dfpb <- as.data.frame(Fractions)
							  									dfpb$Am <- rep(0,241)
							  									dfpb_rnase$Amount <- rep(0,241)
							  				
							  									for (i in 1:n) {
							  									dfpb$Am <- ( fit_c_rnase[i] * exp(-(Fractions-fit_mean_rnase[i])**2/(2 * fit_sigma_rnase[i]**2)) )
							  									dfpb$Amount <- dfpb_rnase$Amount
							  									dfpb_rnase$Amount <- apply(dfpb,1, function(x) {max(x[2:3])})
							  													}
							  				
							  									dffit_rnase <- dfpb_rnase
							  									dffit_ctrl <- fine_fit_fun_ctrl(pos) 
																  colnames(dffit_ctrl) <- c("Fractions","Amount")
																
																} else {
																		if ((table_ctrl_mean[pos, "fit_res"] > 260 || nl_C) && (table_rnase_mean[pos, "fit_res"] > 260) || nl_R) {
																						  fit_c_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_c_corr"]))
							  															fit_mean_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_mean_corr"]))
							  															fit_sigma_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_sigma_corr"]))
							  															n <- length(fit_c_ctrl)
							  															Fractions <- round(seq(1,25,0.1),digits=1)
							  															dfpb_ctrl <- as.data.frame(Fractions)
							  															dfpb <- as.data.frame(Fractions)
							  															dfpb$Am <- rep(0,241)
							  															dfpb_ctrl$Amount <- rep(0,241)
							  				
							  															for (i in 1:n) {
							  															dfpb$Am <- ( fit_c_ctrl[i] * exp(-(Fractions-fit_mean_ctrl[i])**2/(2 * fit_sigma_ctrl[i]**2)) )
							  															dfpb$Amount <- dfpb_ctrl$Amount
							  															dfpb_ctrl$Amount <- apply(dfpb,1, function(x) {max(x[2:3])})
							  																			}
							  															
																						  fit_c_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_c_corr"]))
							  															fit_mean_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_mean_corr"]))
							  															fit_sigma_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_sigma_corr"]))
							  															n <- length(fit_c_rnase)
							  															Fractions <- round(seq(1,25,0.1),digits=1)
							  															dfpb_rnase <- as.data.frame(Fractions)
							  															dfpb <- as.data.frame(Fractions)
							  															dfpb$Am <- rep(0,241)
							  															dfpb_rnase$Amount <- rep(0,241)
							  				
							  															for (i in 1:n) {
							  															dfpb$Am <- ( fit_c_rnase[i] * exp(-(Fractions-fit_mean_rnase[i])**2/(2 * fit_sigma_rnase[i]**2)) )
							  															dfpb$Amount <- dfpb_rnase$Amount
							  															dfpb_rnase$Amount <- apply(dfpb,1, function(x) {max(x[2:3])})
							  																			}
																						
																						dffit_ctrl <- dfpb_ctrl
																						dffit_rnase <- dfpb_rnase
																						
																						} else {
																								dffit_ctrl <- fine_fit_fun_ctrl(pos) 
																								colnames(dffit_ctrl) <- c("Fractions","Amount")
																								dffit_rnase <- fine_fit_fun_rnase(pos) 
																								colnames(dffit_rnase) <- c("Fractions","Amount")
																								}
																		
																		}
												}
												
								value <- c(dffit_ctrl[dffit_ctrl$Fractions == max,"Amount"],dffit_rnase[dffit_rnase$Fractions == max,"Amount"])				
								value
								
								# end of the function
								}


# Evaluation of the parameters used to define a shift.
evaluate_shifts <- function(pos) 		{
																	   shift <- numeric(0)
																	   
																	   # get the maxima (< threshold)
																	   list_rnase_max_th <- as.numeric(unlist(table_rnase_mean[pos,"max_th"]))
																	   list_rnase_p_values <- as.numeric(unlist(table_rnase_mean[pos,"p_values_fdr"]))
																	   
																	   if (table_rnase_mean[pos,"pb_fit"]) {
																	   										list_rnase_max <- as.numeric(unlist(table_rnase_mean[pos,"fit_mean_fxn_corr"]))
																	   										list_rnase_area <- as.numeric(unlist(table_rnase_mean[pos,"fit_area_corr"]))
																	   										} else {
																	   												list_rnase_max <- as.numeric(unlist(table_rnase_mean[pos,"fit_mean_fxn"]))
																	   												list_rnase_area <- as.numeric(unlist(table_rnase_mean[pos,"fit_area"]))
																	   												}
																	   
																	   list_ctrl_max_th <- as.numeric(unlist(table_ctrl_mean[pos,"max_th"]))
																	   list_ctrl_p_values <- as.numeric(unlist(table_ctrl_mean[pos,"p_values_fdr"]))
																	   
																	   if (table_ctrl_mean[pos,"pb_fit"]) {
																	   										list_ctrl_max <- as.numeric(unlist(table_ctrl_mean[pos,"fit_mean_fxn_corr"]))
																	   										list_ctrl_area <- as.numeric(unlist(table_ctrl_mean[pos,"fit_area_corr"]))
																	   										} else {
																	   												list_ctrl_max <- as.numeric(unlist(table_ctrl_mean[pos,"fit_mean_fxn"]))
																	   												list_ctrl_area <- as.numeric(unlist(table_ctrl_mean[pos,"fit_area"]))
																	   												}
																	   
																	   
																	   
																	   if (sum(list_rnase_max_th) == 0) {shift <- 0} else {
																	   
																	   R <- length(list_rnase_max_th)
																	   C <- length(list_ctrl_max_th)
																	   
																	   
																	   #**************************************************************************************************
																	   # calculate the distance between the peaks
																	   dist <- numeric(0)
																	   for (i in 1:R) {
																	   				   dist <- c(dist, (list_rnase_max_th[i]-list_ctrl_max_th))
																	   				   }				   
																	   dist <- unlist(dist)
																	   #dist
																	   
																	   
																	   #**************************************************************************************************
																	   # get the corresponding list of amounts
																	   list_rnase_amount <- numeric(0)
																	   list_ctrl_amount <- numeric(0)
																	   list_rnase_amount <- list_rnase_area[list_rnase_max %in% list_rnase_max_th]
																	   list_ctrl_amount <- list_ctrl_area[list_ctrl_max %in% list_ctrl_max_th]
																	   
																	   amount_C <- list_ctrl_amount
																	   amount_R <- list_rnase_amount
																	   amount_C <- round(amount_C, digits=1)
																	   amount_R <- round(amount_R, digits=1)


																	   #**************************************************************************************************
																	   # calculate loss of amplitude for each ctrl maximum
																	   # calculate gain of amplitude for each rnase maximum
																	   
																	   loss_C <- numeric(0)
																	   for (i in 1:C) {
																	   				   loss <- (ampl_max(list_ctrl_max_th[i], pos)[1] - ampl_max(list_ctrl_max_th[i], pos)[2])*100/ampl_max(list_ctrl_max_th[i], pos)[1]
																	   				   loss_C <- c(loss_C, loss)
																	   				   }				   
																	   loss_C <- unlist(loss_C)
																	   loss_C <- round(loss_C, digits=4)
																	   
																	   gain_R <- numeric(0)
																	   for (i in 1:R) {
																	   				   gain <- (ampl_max(list_rnase_max_th[i], pos)[2] - ampl_max(list_rnase_max_th[i], pos)[1])*100/ampl_max(list_rnase_max_th[i], pos)[2]
																	   				   gain_R <- c(gain_R, gain)
																	   				   }				   
																	   gain_R <- unlist(gain_R)
																	   gain_R <- round(gain_R, digits=4)
																	   
																	
																	   #**************************************************************************************************
																	   p_value_C <- list_ctrl_p_values[list_ctrl_max %in% list_ctrl_max_th]
																	   p_value_R <- list_rnase_p_values[list_rnase_max %in% list_rnase_max_th]
																	   
																	   
																	   #**************************************************************************************************
																	   # compute shift vector with all parameters (dist,amount_C,amount_R,p_value_C,p_value_R)
																	   shift <- c(dist, rep(list_ctrl_max_th,R),rep(amount_C,R),rep(loss_C,R),rep(p_value_C,R))
																	   
																	   for (i in 1:R) {
																	   				   shift <- c(shift, rep(list_rnase_max_th[i], C))
																	   				   }
																	   
																	   for (i in 1:R) {
																	   				   shift <- c(shift, rep(amount_R[i], C))
																	   				   }
																	  
																	   for (i in 1:R) {
																	   				   shift <- c(shift, rep(gain_R[i], C))
																	   				   }
																	   
																	   for (i in 1:R) {
																	   				   shift <- c(shift, rep(p_value_R[i], C))
																	   				   }
																	   
																	   shift <- unlist(shift)
																													  		}
										shift
										} 


#**************************************************************************************************
# Select significant shift:
# dist > 1 AND one p_value < 0.05
# do not select shift if gain/loss are both negative!
shift_select <- function(shift) {
									value <- numeric(0)
									n <- length(shift)/9
									if (n<1) {value <- 0} else {
																shift_mat <- matrix(shift,n,9,byrow=FALSE)
																shift_mat <- as.data.frame(shift_mat)
																
																# select dist > +/- 1 fraction
																shift_mat <- shift_mat[abs(shift_mat$V1) > 1,]
																
																# select one of the p_value < 0.05
																shift_high_pv <- shift_mat[(shift_mat$V5 >= 0.05),]
																shift_high_pv <- shift_high_pv[(shift_high_pv$V9 >= 0.05),]
																
																value <- shift_mat[! rownames(shift_mat) %in% rownames(shift_high_pv),]
																}
								value <- as.numeric(unlist(value))
								if (length(value) == 0) {value <- 0} else {
								value <- value}
								value
								}


#**************************************************************************************************
# Apply functions evaluate_shifts and shift_select to all proteins

table_ctrl_mean$selected_shift <- apply(vect, 2, function(pos) {
															  shifts <- evaluate_shifts(pos)
															  shifts <- unlist(shifts)
															  shifts
															  selected_shift <- shift_select(shifts)
															  selected_shift <- unlist(selected_shift)
															  selected_shift
															  })

table_rnase_mean$selected_shift <- apply(vect, 2, function(pos) {
															  shifts <- evaluate_shifts(pos)
															  selected_shift <- shift_select(shifts)
															  selected_shift <- unlist(selected_shift)
															  selected_shift
															  })


#**************************************************************************************************
# Add some columns to help to select for interesting proteins

# Number of shift per protein
table_ctrl_mean$nb_shift <- apply(vect, 2, function(pos) {
														  shift <- as.numeric(unlist(table_ctrl_mean[pos, "selected_shift"]))
														  n <- length(shift)/9
														  
														  shift_mat <- matrix(shift,n,9,byrow=FALSE)
														  shift_mat <- as.data.frame(shift_mat)
														  
														  # Compute real number of shifts
														  shift_mat <- shift_mat[shift_mat[,2]>0,]
														  shift_mat <- shift_mat[shift_mat[,4]>0,]
														  shift_mat <- shift_mat[shift_mat[,6]>0,]
														  shift_mat <- shift_mat[shift_mat[,8]>0,]
														  n <- dim(shift_mat)[1]
														  
														  if (n<1) {0} else {n}
														  
														  })


table_rnase_mean$nb_shift <- apply(vect, 2, function(pos) {
														  shift <- as.numeric(unlist(table_rnase_mean[pos, "selected_shift"]))
														  n <- length(shift)/9
														  
														  shift_mat <- matrix(shift,n,9,byrow=FALSE)
														  shift_mat <- as.data.frame(shift_mat)
														  
														  # Compute real number of shifts
														  shift_mat <- shift_mat[shift_mat[,2]>0,]
														  shift_mat <- shift_mat[shift_mat[,4]>0,]
														  shift_mat <- shift_mat[shift_mat[,6]>0,]
														  shift_mat <- shift_mat[shift_mat[,8]>0,]
														  n <- dim(shift_mat)[1]
														  
														  if (n<1) {0} else {n}
														  
														  })


#**************************************************************************************************
#### PART 7 - Computing of the last table and exprot as CSV file                               ####
#**************************************************************************************************

# function to add lines with shifts from a protein "pos" to an existing dataframe "tab"
shift_table <- function(tab,pos) {
							  		prot_name <- as.vector(tab$protein_name)
							  		
							  		shift <- as.numeric(unlist(table_ctrl_mean[pos, "selected_shift"]))
									n <- length(shift)/9
														  
									if (n<1) {shift_mat <- rep(0,9)
											  } else {
													  shift_mat <- matrix(shift,n,9,byrow=FALSE)
													  }					  					 		
									shift_mat <- as.data.frame(shift_mat)
								
									if (n<1) { 
												tab[nrow(tab)+1,2:26] <- as.numeric(unlist(table_ctrl_mean[pos,1:25]))
												tab[nrow(tab),]$max_ctrl <- table_ctrl_mean[pos,"max_th"]
												tab[nrow(tab),]$nb_max_ctrl <- table_ctrl_mean[pos,"nb_max_th"]
												
												tab[nrow(tab),29:53] <- as.numeric(unlist(table_rnase_mean[pos,1:25]))
												tab[nrow(tab),]$max_rnase <- table_rnase_mean[pos,"max_th"]
												tab[nrow(tab),]$nb_max_rnase <- table_rnase_mean[pos,"nb_max_th"]
												
												tab[nrow(tab),]$nb_shift <- table_rnase_mean[pos,"nb_shift"]
												
												tab[nrow(tab),57:65] <- as.numeric(unlist(shift_mat))
												
												prot_name <- c(prot_name, pos)
												
												} else {
														for (i in 1:n) {
																		tab[nrow(tab)+1,2:26] <- as.numeric(unlist(table_ctrl_mean[pos,1:25]))
																		tab[nrow(tab),]$max_ctrl <- table_ctrl_mean[pos,"max_th"]
																		tab[nrow(tab),]$nb_max_ctrl <- table_ctrl_mean[pos,"nb_max_th"]
												
																		tab[nrow(tab),29:53] <- as.numeric(unlist(table_rnase_mean[pos,1:25]))
																		tab[nrow(tab),]$max_rnase <- table_rnase_mean[pos,"max_th"]
																		tab[nrow(tab),]$nb_max_rnase <- table_rnase_mean[pos,"nb_max_th"]
												
																		tab[nrow(tab),]$nb_shift <- table_rnase_mean[pos,"nb_shift"]
												
																		tab[nrow(tab),57:65] <- as.numeric(unlist(shift_mat[i,]))

																		prot_name <- c(prot_name, pos)
																		}
														}
							 		tab$protein_name <- prot_name
							 		tab
							 		
							  		}
							  		
							  		
# Creation of a summary table

# First row of the table:
rownames <- rownames(table)
lgth <- length(rownames)
sum_tab <- numeric(0)
sum_tab <- matrix(as.numeric(table_ctrl_mean[rownames[1],1:25]), 1, 25, byrow = "TRUE")
sum_tab <- as.data.frame(sum_tab)
pos <- rownames[1]
colnames(sum_tab)[1:25] <- c("C_mean_fxn1","C_mean_fxn2","C_mean_fxn3","C_mean_fxn4","C_mean_fxn5","C_mean_fxn6","C_mean_fxn7","C_mean_fxn8","C_mean_fxn9","C_mean_fxn10","C_mean_fxn11","C_mean_fxn12","C_mean_fxn13","C_mean_fxn14","C_mean_fxn15","C_mean_fxn16","C_mean_fxn17","C_mean_fxn18","C_mean_fxn19","C_mean_fxn20","C_mean_fxn21","C_mean_fxn22","C_mean_fxn23","C_mean_fxn24","C_mean_fxn25")
sum_tab$max_ctrl <- table_ctrl_mean[pos,"max_th"]
sum_tab$nb_max_ctrl <- table_ctrl_mean[pos,"nb_max_th"]
sum_tab[,28:52] <- as.numeric(table_rnase_mean[pos,1:25])
colnames(sum_tab)[28:52] <- c("R_mean_fxn1","R_mean_fxn2","R_mean_fxn3","R_mean_fxn4","R_mean_fxn5","R_mean_fxn6","R_mean_fxn7","R_mean_fxn8","R_mean_fxn9","R_mean_fxn10","R_mean_fxn11","R_mean_fxn12","R_mean_fxn13","R_mean_fxn14","R_mean_fxn15","R_mean_fxn16","R_mean_fxn17","R_mean_fxn18","R_mean_fxn19","R_mean_fxn20","R_mean_fxn21","R_mean_fxn22","R_mean_fxn23","R_mean_fxn24","R_mean_fxn25")

sum_tab$max_rnase <- table_rnase_mean[pos,"max_th"]
sum_tab$nb_max_rnase <- table_rnase_mean[pos,"nb_max_th"]
							  		
sum_tab$nb_shift <- as.numeric(table_rnase_mean[pos,"nb_shift"])
							  		
shift <- as.numeric(unlist(table_ctrl_mean[pos, "selected_shift"]))
									n <- length(shift)/9
									shift_mat <- matrix(shift,n,9,byrow=FALSE)				
									shift_mat <- as.data.frame(shift_mat)

sum_tab <- cbind(sum_tab, shift_mat)

colnames(sum_tab)[56:64] <- c("dist","ctrl_peak","ctrl_peak_amount","ctrl_peak_amount_loss","ctrl_peak_p_value","rnase_peak","rnase_peak_amount","rnase_peak_amount_gain","rnase_peak_p_value")


sum_tab <- cbind(as.data.frame(pos),sum_tab)
colnames(sum_tab)[1] <- "protein_name"
sum_tab <- as.data.frame(sum_tab)


# addition of the other rows:
for (i in 2:lgth) {
				   pos <- rownames[i]
				   sum_tab <- shift_table(sum_tab, pos)
				   }


# any shift towards larger fractions? (positive shift)
lgth <- dim(sum_tab)[1]
vect6 <- matrix(c(1:lgth), 1, lgth)

sum_tab$right_shift <- apply(vect6, 2, function(pos) {
														  n <- sum_tab[pos,"nb_shift"]
														  dist <- sum_tab[pos,"dist"]
														 
														  if (n<1) {FALSE} else {
																			 	if (dist>0) {
																			 			  	TRUE
																			 				} else {FALSE}
														  					    }
														  })
														  

#**************************************************************************************************
# Add one column with the information about difficult fitting of the raw data. All data with fitted = "FALSE" and pb_fit = "TRUE"
sum_tab$pb_fit <- apply(sum_tab, 1, function(x) {
												 pos <- x$protein_name
												 pb_fit <- any(!table_ctrl_mean[pos,"fitted"],table_ctrl_mean[pos,"pb_fit"],!table_rnase_mean[pos,"fitted"],table_rnase_mean[pos,"pb_fit"])
												 pb_fit
												 })

#**************************************************************************************************

# Collapse the two list ctrl_max and rnase_max to export the table to a file.
sum_tab_exp <- sum_tab

sum_tab_exp$max_ctrl <- apply(vect6, 2, function(x) {
																	list <- unlist(sum_tab_exp[x,"max_ctrl"])
																	list <- paste(list, sep=", ", collapse=", ")
																	list
																	}) 

sum_tab_exp$max_rnase <- apply(vect6, 2, function(x) {
																	list <- unlist(sum_tab_exp[x,"max_rnase"])
																	list <- paste(list, sep=", ", collapse=", ")
																	list
																	})


#**************************************************************************************************
# FILTERS on the last table (sum_tab -> contains information about the 4968 starting proteins)

# remove gain/loss double negative (I verified that there is a corresponding counterpart gain/loss double positive)
# replace shift information for those proteins with 0.

# Remove ctrl peak info for negative loss. No shift. 
sum_tab_exp[sum_tab_exp$ctrl_peak_amount_loss <= 0,][, c("nb_shift","dist","ctrl_peak","ctrl_peak_amount","ctrl_peak_amount_loss","ctrl_peak_p_value")] <- 0
# Remove right_shift info for negative loss. 
sum_tab_exp[sum_tab_exp$ctrl_peak_amount_loss <= 0,][, c("right_shift")] <- FALSE
# Remove rnase peak info for non_significant gain -> not the case in the sample file: error message
# Error in value[[jvseq[[jjj]]]] : subscript out of bounds
sum_tab_exp[((sum_tab_exp$ctrl_peak_amount_loss <= 0)+(sum_tab_exp$rnase_peak_p_value > 0.05)) == 2,][, c("rnase_peak","rnase_peak_amount","rnase_peak_amount_gain","rnase_peak_p_value")] <- 0

# Remove rnase peak info for negative gain. No shift.
sum_tab_exp[sum_tab_exp$rnase_peak_amount_gain <= 0,][, c("nb_shift","dist","rnase_peak","rnase_peak_amount","rnase_peak_amount_gain","rnase_peak_p_value")] <- 0

# Remove ctrl peak info for non_significant loss -> not the case in the sample file: error message
# Error in value[[jvseq[[jjj]]]] : subscript out of bounds
# sum_tab_exp[((sum_tab_exp$rnase_peak_amount_gain <= 0)+(sum_tab_exp$ctrl_peak_p_value > 0.05)) == 2,][, c("ctrl_peak","ctrl_peak_amount","ctrl_peak_amount_loss","ctrl_peak_p_value")] <- 0

# Remove right_shift info for negative gain. 
sum_tab_exp[sum_tab_exp$rnase_peak_amount_gain <= 0,][, c("right_shift")] <- FALSE

# Remove right_shift info for non_significant gain. Just in case -> apparently there are no cases like this - therefore error message:
# Error in `[<-.data.frame`(`*tmp*`, , c("right_shift"), value = FALSE) : replacement has 1 row, data has 0
# sum_tab_exp[((sum_tab_exp$ctrl_peak_amount_loss <= 0)+(sum_tab_exp$rnase_peak_p_value > 0.05)) == 2,][, c("right_shift")] <- FALSE


#**************************************************************************************************
# Remove Gain or Loss which appear in duplicate
# EXPORT TABLES
prot_name <- unique(sum_tab_exp$protein_name)
prot_lg <- length(prot_name)

i<-1
pos <- prot_name[i]
					temp <- sum_tab_exp[sum_tab_exp$protein_name == pos,]
					temp <- unique(temp)
					ctrl_peak <- temp$ctrl_peak
					rnase_peak <- temp$rnase_peak
					n <- length(ctrl_peak)

					if (n==1) {write.table(temp, file="MS_Analysis_Shifts.csv", row.names = FALSE, col.names = TRUE, append = FALSE, sep = ";")} else {
								   																		  temp <- temp[!((temp$ctrl_peak==0)+(temp$rnase_peak==0))==2,]
								   																		  
								   																		  ctrl_peak <- temp$ctrl_peak
								   																		  rnase_peak <- temp$rnase_peak
								   																		  n <- length(ctrl_peak)
								   																		  
								   																		  C_peak <- temp[((ctrl_peak>0)+(rnase_peak>0))==2,]$ctrl_peak
								   																		  R_peak <- temp[((ctrl_peak>0)+(rnase_peak>0))==2,]$rnase_peak
								   																		  temp <- temp[!((temp$rnase_peak %in% R_peak)+(temp$ctrl_peak %in% C_peak))==1,]
								   																		  
								   																		  write.table(temp, file="MS_Analysis_Shifts.csv", row.names = FALSE, col.names = TRUE, append = FALSE, sep = ";")
														   															  }
						
for (i in 2:(prot_lg)) {
 # for (i in 2:(100)) {
					pos <- prot_name[i]
					temp <- sum_tab_exp[sum_tab_exp$protein_name == pos,]
					temp <- unique(temp)
					ctrl_peak <- temp$ctrl_peak
					rnase_peak <- temp$rnase_peak
					n <- length(ctrl_peak)

					if (n==1) {write.table(temp, file="MS_Analysis_Shifts.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ";")} else {
								   																		  temp <- temp[!((temp$ctrl_peak==0)+(temp$rnase_peak==0))==2,]
								   																		  
								   																		  ctrl_peak <- temp$ctrl_peak
								   																		  rnase_peak <- temp$rnase_peak
								   																		  n <- length(ctrl_peak)
								   																		  
								   																		  C_peak <- temp[((ctrl_peak>0)+(rnase_peak>0))==2,]$ctrl_peak
								   																		  R_peak <- temp[((ctrl_peak>0)+(rnase_peak>0))==2,]$rnase_peak
								   																		  temp <- temp[!((temp$rnase_peak %in% R_peak)+(temp$ctrl_peak %in% C_peak))==1,]
								   																		  
								   																		  write.table(temp, file="MS_Analysis_Shifts.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ";")
														   															  }
						}

# Read the final table from the file
#**************************************************************************************************					
MS_Analysis_table <- read.table("MS_Analysis_Shifts.csv", header = TRUE, sep = ";")


#**************************************************************************************************
#### GRAPHICS 1 - average curves for one protein of interest                                   ####
#**************************************************************************************************

# CURVES - AMOUNTS (MEAN)
#**************************************************************************************************
# Load packages
library(sp)
library(raster)
library(grid)
library(ggplot2)
library(gridExtra)
library(lattice)

# Curves for a particular protein
# pos = "CTCF_HUMAN"


# define a function to calculate gaussian fit of the curves
#***************************************************************************************************
gauss_fit_fun_ctrl <- function(t) { 
												  dffit <- data.frame(y=round(seq(1, 25, 0.5),digits=1))
												  list_c <- as.numeric(unlist(ctrl_mean[t,"fit_c"]))
												  list_mean <- as.numeric(unlist(ctrl_mean[t,"fit_mean"]))
												  list_sigma <- as.numeric(unlist(ctrl_mean[t,"fit_sigma"]))
												  n <- length(list_c)
												  z <- round(seq(1, 25, 0.5),digits=1)
												  
												  if (rnase_mean[t,"fitted"] == "FALSE") {dffit$df.y <- rep(0,49)} else {
												  											if (n==1) {
												  														q <- c(list_c[1],list_mean[1],list_sigma[1])
												  														C1 <- q[1]
																										mean1 <- q[2]
																										sigma1 <- q[3]
																										dffit$df.y <- C1 * exp(-(z-mean1)**2/(2 * sigma1**2))
												  														} else {
												  																if (n==2) {
												  																			q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2])
												  																			C1 <- q[1]
																															mean1 <- q[2]
																															sigma1 <- q[3]
																															C2 <- q[4]
																															mean2 <- q[5]
																															sigma2 <- q[6]
																															dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) ) 
												  																			} else {
												  																					if (n==3) {
												  																								q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3])
												  																								C1 <- q[1]
																																				mean1 <- q[2]
																																				sigma1 <- q[3]
																																				C2 <- q[4]
																																				mean2 <- q[5]
																																				sigma2 <- q[6]
																																				C3 <- q[7]
																																				mean3 <- q[8]
																																				sigma3 <- q[9]
																																				dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) ) 
												  																								} else {
												  																										if (n==4) {
												  																													q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4])
												  																													C1 <- q[1]
																																									mean1 <- q[2]
																																									sigma1 <- q[3]
																																									C2 <- q[4]
																																									mean2 <- q[5]
																																									sigma2 <- q[6]
																																									C3 <- q[7]
																																									mean3 <- q[8]
																																									sigma3 <- q[9]
																																									C4 <- q[10]
																																									mean4 <- q[11]
																																									sigma4 <- q[12]
																																									dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) )
												  																													} else {
												  																															if (n==5) {
												  																																		q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5])
												  																																		C1 <- q[1]
																																														mean1 <- q[2]
																																														sigma1 <- q[3]
																																														C2 <- q[4]
																																														mean2 <- q[5]
																																														sigma2 <- q[6]
																																														C3 <- q[7]
																																														mean3 <- q[8]
																																														sigma3 <- q[9]
																																														C4 <- q[10]
																																														mean4 <- q[11]
																																														sigma4 <- q[12]
																																														C5 <- q[13]
																																														mean5 <- q[14]
																																														sigma5 <- q[15]
																																														dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) )
												  																																		} else {
												  																																				q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5],list_c[6],list_mean[6],list_sigma[6])
												  																																				C1 <- q[1]
																																																mean1 <- q[2]
																																																sigma1 <- q[3]
																																																C2 <- q[4]
																																																mean2 <- q[5]
																																																sigma2 <- q[6]
																																																C3 <- q[7]
																																																mean3 <- q[8]
																																																sigma3 <- q[9]
																																																C4 <- q[10]
																																																mean4 <- q[11]
																																																sigma4 <- q[12]
																																																C5 <- q[13]
																																																mean5 <- q[14]
																																																sigma5 <- q[15]
																																																C6 <- q[16]
																																																mean6 <- q[17]
																																																sigma6 <- q[18]
																																																dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) + C6 * exp(-(z-mean6)**2/(2 * sigma6**2)) ) 
												  																																				}
												  																															}
												  																										}
												  																					}
												  																}
												  											}
											dffit
											}							
	
#***************************************************************************************************	
gauss_fit_fun_rnase <- function(t) { 
												  dffit <- data.frame(y=round(seq(1, 25, 0.5),digits=1))
												  list_c <- as.numeric(unlist(rnase_mean[t,"fit_c"]))
												  list_mean <- as.numeric(unlist(rnase_mean[t,"fit_mean"]))
												  list_sigma <- as.numeric(unlist(rnase_mean[t,"fit_sigma"]))
												  n <- length(list_c)
												  z <- round(seq(1, 25, 0.5),digits=1)
												  
												  if (rnase_mean[t,"fitted"] == "FALSE") {dffit$df.y <- rep(0,49)} else {
												  											if (n==1) {
												  														q <- c(list_c[1],list_mean[1],list_sigma[1])
												  														C1 <- q[1]
																										mean1 <- q[2]
																										sigma1 <- q[3]
																										dffit$df.y <- C1 * exp(-(z-mean1)**2/(2 * sigma1**2))
												  														} else {
												  																if (n==2) {
												  																			q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2])
												  																			C1 <- q[1]
																															mean1 <- q[2]
																															sigma1 <- q[3]
																															C2 <- q[4]
																															mean2 <- q[5]
																															sigma2 <- q[6]
																															dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) ) 
												  																			} else {
												  																					if (n==3) {
												  																								q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3])
												  																								C1 <- q[1]
																																				mean1 <- q[2]
																																				sigma1 <- q[3]
																																				C2 <- q[4]
																																				mean2 <- q[5]
																																				sigma2 <- q[6]
																																				C3 <- q[7]
																																				mean3 <- q[8]
																																				sigma3 <- q[9]
																																				dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) ) 
												  																								} else {
												  																										if (n==4) {
												  																													q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4])
												  																													C1 <- q[1]
																																									mean1 <- q[2]
																																									sigma1 <- q[3]
																																									C2 <- q[4]
																																									mean2 <- q[5]
																																									sigma2 <- q[6]
																																									C3 <- q[7]
																																									mean3 <- q[8]
																																									sigma3 <- q[9]
																																									C4 <- q[10]
																																									mean4 <- q[11]
																																									sigma4 <- q[12]
																																									dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) )
												  																													} else {
												  																															if (n==5) {
												  																																		q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5])
												  																																		C1 <- q[1]
																																														mean1 <- q[2]
																																														sigma1 <- q[3]
																																														C2 <- q[4]
																																														mean2 <- q[5]
																																														sigma2 <- q[6]
																																														C3 <- q[7]
																																														mean3 <- q[8]
																																														sigma3 <- q[9]
																																														C4 <- q[10]
																																														mean4 <- q[11]
																																														sigma4 <- q[12]
																																														C5 <- q[13]
																																														mean5 <- q[14]
																																														sigma5 <- q[15]
																																														dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) )
												  																																		} else {
												  																																				q <- c(list_c[1],list_mean[1],list_sigma[1],list_c[2],list_mean[2],list_sigma[2],list_c[3],list_mean[3],list_sigma[3],list_c[4],list_mean[4],list_sigma[4],list_c[5],list_mean[5],list_sigma[5],list_c[6],list_mean[6],list_sigma[6])
												  																																				C1 <- q[1]
																																																mean1 <- q[2]
																																																sigma1 <- q[3]
																																																C2 <- q[4]
																																																mean2 <- q[5]
																																																sigma2 <- q[6]
																																																C3 <- q[7]
																																																mean3 <- q[8]
																																																sigma3 <- q[9]
																																																C4 <- q[10]
																																																mean4 <- q[11]
																																																sigma4 <- q[12]
																																																C5 <- q[13]
																																																mean5 <- q[14]
																																																sigma5 <- q[15]
																																																C6 <- q[16]
																																																mean6 <- q[17]
																																																sigma6 <- q[18]
																																																dffit$df.y <- ( C1 * exp(-(z-mean1)**2/(2 * sigma1**2)) + C2 * exp(-(z-mean2)**2/(2 * sigma2**2)) + C3 * exp(-(z-mean3)**2/(2 * sigma3**2)) + C4 * exp(-(z-mean4)**2/(2 * sigma4**2)) + C5 * exp(-(z-mean5)**2/(2 * sigma5**2)) + C6 * exp(-(z-mean6)**2/(2 * sigma6**2)) ) 
												  																																				}
												  																															}
												  																										}
												  																					}
												  																}
												  											}
											dffit
											}							

												  
#***************************************************************************************************
# Standard deviation for each mean curve for graphics
ctrl_norm_sd <- t(apply(vect, 2, function(x) { 
											list <- rep(0,25)
											if (table_ctrl_mean[x, "fitted"] == "FALSE") {list <- list} else {
																										for (i in 1:25) {
																									    				list[i] <- sd(c(table.ctrl1.SW.norm[x,i],table.ctrl2.SW.norm[x,i],table.ctrl3.SW.norm[x,i]),na.rm=TRUE)
															 															}
															 										  }					
											list
											}
					  ))
rownames(ctrl_norm_sd) <- rownames(table.ctrl1.SW.norm)
colnames(ctrl_norm_sd) <- col_fractions
ctrl_norm_sd <- as.data.frame(ctrl_norm_sd)


rnase_norm_sd <- t(apply(vect, 2, function(x) { 
											list <- rep(0,25)
											if (table_rnase_mean[x, "fitted"] == "FALSE") {list <- list} else {
																										for (i in 1:25) {
																									    				list[i] <- sd(c(table.ctrl1.SW.norm[x,i],table.ctrl2.SW.norm[x,i],table.ctrl3.SW.norm[x,i]),na.rm=TRUE)
															 															}
															 										  }					
											list
											}
					  ))
rownames(rnase_norm_sd) <- rownames(table.rnase1.SW.norm)
colnames(rnase_norm_sd) <- col_fractions
rnase_norm_sd <- as.data.frame(rnase_norm_sd)


# define a function to plot the mean curves (CTRL and RNASE) (with single peaks for proteins with pb_fit)
# usage: curves_plot (pos) where "pos" is the Entry name (UniProt) of the protein - and only Entry name!
# Example: curves_plot("CTCF_HUMAN")
#********************************************************************************************************
curves_plot <- function(pos) {
  
  # check if protein (pos) is listed in the table
  if (pos %in% rownames(table)) {
    
    # single peak fit for all proteins with pb_fit = TRUE
    if ((table_ctrl_mean[pos, "pb_fit"] == "TRUE") && (table_rnase_mean[pos, "pb_fit"] == "FALSE")) {
      fit_c_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_c_corr"]))
      fit_mean_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_mean_corr"]))
      fit_sigma_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_sigma_corr"]))
      n <- length(fit_c_ctrl)
      Fractions <- round(seq(1,25,0.5),digits=1)
      dfpb_ctrl <- as.data.frame(Fractions)
      plot_ctrl <- ggplot()
      
      for (i in 1:n) {
        dfpb_ctrl$Amount <- ( fit_c_ctrl[i] * exp(-(Fractions-fit_mean_ctrl[i])**2/(2 * fit_sigma_ctrl[i]**2)) )
        plot_ctrl <- plot_ctrl + geom_line(data=dfpb_ctrl, aes(x=Fractions, y=Amount), color="darkgreen", size=1)
      }
      
      dffit_rnase <- gauss_fit_fun_rnase(pos) 
      colnames(dffit_rnase) <- c("Fractions","Amount")
      
    } else {
      if ((table_ctrl_mean[pos, "pb_fit"] == "FALSE") && (table_rnase_mean[pos, "pb_fit"] == "TRUE")) {
        fit_c_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_c_corr"]))
        fit_mean_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_mean_corr"]))
        fit_sigma_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_sigma_corr"]))
        n <- length(fit_c_rnase)
        Fractions <- round(seq(1,25,0.5),digits=1)
        dfpb_rnase <- as.data.frame(Fractions)
        plot_rnase <- ggplot()
        
        for (i in 1:n) {
          dfpb_rnase$Amount <- ( fit_c_rnase[i] * exp(-(Fractions-fit_mean_rnase[i])**2/(2 * fit_sigma_rnase[i]**2)) )
          plot_rnase <- plot_rnase + geom_line(data=dfpb_rnase, aes(x=Fractions, y=Amount), color="darkred", size=1)
        }
        
        dffit_ctrl <- gauss_fit_fun_ctrl(pos) 
        colnames(dffit_ctrl) <- c("Fractions","Amount")
        
      } else {
        if ((table_ctrl_mean[pos, "pb_fit"] == "TRUE") && (table_rnase_mean[pos, "pb_fit"] == "TRUE")) {
          fit_c_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_c_corr"]))
          fit_mean_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_mean_corr"]))
          fit_sigma_ctrl <- as.numeric(unlist(table_ctrl_mean[pos,"fit_sigma_corr"]))
          n <- length(fit_c_ctrl)
          Fractions <- round(seq(1,25,0.5),digits=1)
          dfpb_ctrl <- as.data.frame(Fractions)
          plot_CR <- ggplot()
          
          for (i in 1:n) {
            dfpb_ctrl$Amount <- ( fit_c_ctrl[i] * exp(-(Fractions-fit_mean_ctrl[i])**2/(2 * fit_sigma_ctrl[i]**2)) )
            plot_CR <- plot_CR + geom_line(data=dfpb_ctrl, aes(x=Fractions, y=Amount), color="darkgreen", size=1)
          }
          
          fit_c_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_c_corr"]))
          fit_mean_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_mean_corr"]))
          fit_sigma_rnase <- as.numeric(unlist(table_rnase_mean[pos,"fit_sigma_corr"]))
          n <- length(fit_c_rnase)
          Fractions <- round(seq(1,25,0.5),digits=1)
          dfpb_rnase <- as.data.frame(Fractions)
          
          for (i in 1:n) {
            dfpb_rnase$Amount <- ( fit_c_rnase[i] * exp(-(Fractions-fit_mean_rnase[i])**2/(2 * fit_sigma_rnase[i]**2)) )
            plot_CR <- plot_CR + geom_line(data=dfpb_rnase, aes(x=Fractions, y=Amount), color="darkred", size=1)
          }
          
        } else {
          dffit_ctrl <- gauss_fit_fun_ctrl(pos) 
          colnames(dffit_ctrl) <- c("Fractions","Amount")
          dffit_rnase <- gauss_fit_fun_rnase(pos) 
          colnames(dffit_rnase) <- c("Fractions","Amount")
        }
        
      }
    }					
    
    # raw data - Amounts 
    x <- c(1:25)
    ctrl <- as.numeric(ctrl_norm_mean[pos,1:25])
    title <- pos
    ctrl <- as.numeric(ctrl)
    df_ctrl <- data.frame(x, ctrl)
    colnames(df_ctrl) <- c("Fractions","Amount")
    rnase <- as.numeric(rnase_norm_mean[pos,1:25])
    df_rnase <- data.frame(x, rnase)
    colnames(df_rnase) <- c("Fractions","Amount")
    
    # standard deviation of mean values
    df_sd_ctrl <- df_ctrl
    df_sd_ctrl$upper <- ctrl + as.numeric(ctrl_norm_sd[pos,1:25])
    df_sd_ctrl$lower <- ctrl - as.numeric(ctrl_norm_sd[pos,1:25])
    
    df_sd_rnase <- df_rnase
    df_sd_rnase$upper <- rnase + as.numeric(rnase_norm_sd[pos,1:25])
    df_sd_rnase$lower <- rnase - as.numeric(rnase_norm_sd[pos,1:25])
    
    # scatterPlot - depending on the "fitted", "fit_res", "pb_fit" values	
    
    #***************************************************************************************************
    if ((table_ctrl_mean[pos, "fitted"] == "FALSE") || (table_rnase_mean[pos, "fitted"] == "FALSE") ) {
      
      scatterPlot <- ggplot() + 
        
        # curve CTRL
        geom_line(data=df_ctrl, aes(x=Fractions, y=Amount, color="green"), size=.5) +
        geom_point(data=df_ctrl, mapping=aes(x=Fractions, y=Amount), color="green", shape=1, size=1) +
        
        # shadow CTRL
        geom_area(data=df_ctrl, aes(x=Fractions, y=Amount), fill="green", alpha=.1, show.legend=F) +
        
        # shadded errors
        geom_ribbon(data=df_sd_ctrl,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
        
        # curve RNASE
        geom_line(data=df_rnase, aes(x=Fractions, y=Amount, color="red"), size=.5) +
        geom_point(data=df_rnase, mapping=aes(x=Fractions, y=Amount), color="red", shape=2, size=1) +
        
        # shadded errors
        geom_ribbon(data=df_sd_rnase,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
        
        # shadow RNASE
        geom_area(data=df_rnase, aes(x=Fractions, y=Amount), fill="red", alpha=.1, show.legend=F) +
        
        # fit CTRL
        geom_line(data=dffit_ctrl, aes(x=Fractions, y=Amount), color="darkgreen", size=1) +
        
        # fit RNASE
        geom_line(data=dffit_rnase, aes(x=Fractions, y=Amount), color="darkred", size=1) +
        
        # New axis label
        ylab("normalized protein amount") +
        
        # title
        title_pos +
        
        # Legend
        scale_color_manual(values = c("green","red"), labels=c("ctrl","RNases")) +
        theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank(), axis.title.x=element_blank(), plot.title = element_text(size = rel(1.5), colour = "darkred"), legend.background = element_rect(colour = "black"),legend.key = element_rect(fill = "white")) +
        
        # panel options
        theme(panel.background = element_rect(fill = "white", color = "darkred", size = 3), panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey")) +
        
        # x-axis scale and ticks
        scale_x_discrete(name ="Fractions", limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))
      
      
      
    } else {
      
      #***************************************************************************************************
      if ((table_ctrl_mean[pos, "pb_fit"] == "TRUE") && (table_rnase_mean[pos, "pb_fit"] == "FALSE")) {
        
        scatterPlot <- plot_ctrl + 
          
          # curve CTRL
          geom_line(data=df_ctrl, aes(x=Fractions, y=Amount, color="green"), size=.5) +
          geom_point(data=df_ctrl, mapping=aes(x=Fractions, y=Amount), color="green", shape=1, size=1) +
          
          # shadow CTRL
          geom_area(data=df_ctrl, aes(x=Fractions, y=Amount), fill="green", alpha=.1, show.legend=F) +
          
          # shadded errors
          geom_ribbon(data=df_sd_ctrl,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
          
          # curve RNASE
          geom_line(data=df_rnase, aes(x=Fractions, y=Amount, color="red"), size=.5) +
          geom_point(data=df_rnase, mapping=aes(x=Fractions, y=Amount), color="red", shape=2, size=1) +
          
          # shadded errors
          geom_ribbon(data=df_sd_rnase,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
          
          # shadow RNASE
          geom_area(data=df_rnase, aes(x=Fractions, y=Amount), fill="red", alpha=.1, show.legend=F) +
          
          # fit CTRL
          #plot_ctrl +
          
          # fit RNASE
          geom_line(data=dffit_rnase, aes(x=Fractions, y=Amount), color="darkred", size=1) +
          
          # New axis label
          ylab("normalized protein amount") +
          
          # title
          labs(title = title) +
          
          # Legend
          scale_color_manual(values = c("green","red"), labels=c("ctrl","RNases")) +
          theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank(), axis.title.x=element_blank(), plot.title = element_text(size = rel(1.5), colour = "darkred"), legend.background = element_rect(colour = "black"),legend.key = element_rect(fill = "white")) +
          
          # panel options
          theme(panel.background = element_rect(fill = "white", color = "darkred", size = 3), panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey")) +
          
          # x-axis scale and ticks
          scale_x_discrete(name ="Fractions", limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))
        
        
      } else {
        
        #***************************************************************************************************
        if ((table_ctrl_mean[pos, "pb_fit"] == "FALSE") && (table_rnase_mean[pos, "pb_fit"] == "TRUE")) {
          
          scatterPlot <- plot_rnase + 
            
            # curve CTRL
            geom_line(data=df_ctrl, aes(x=Fractions, y=Amount, color="green"), size=.5) +
            geom_point(data=df_ctrl, mapping=aes(x=Fractions, y=Amount), color="green", shape=1, size=1) +
            
            # shadow CTRL
            geom_area(data=df_ctrl, aes(x=Fractions, y=Amount), fill="green", alpha=.1, show.legend=F) +
            
            # shadded errors
            geom_ribbon(data=df_sd_ctrl,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
            
            # curve RNASE
            geom_line(data=df_rnase, aes(x=Fractions, y=Amount, color="red"), size=.5) +
            geom_point(data=df_rnase, mapping=aes(x=Fractions, y=Amount), color="red", shape=2, size=1) +
            
            # shadded errors
            geom_ribbon(data=df_sd_rnase,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
            
            # shadow RNASE
            geom_area(data=df_rnase, aes(x=Fractions, y=Amount), fill="red", alpha=.1, show.legend=F) +
            
            # fit CTRL
            geom_line(data=dffit_ctrl, aes(x=Fractions, y=Amount), color="darkgreen", size=1) +
            
            # fit RNASE
            # plot_rnase +
            
            # New axis label
            ylab("normalized protein amount") +
            
            # title
            labs(title = title) +
            
            # Legend
            scale_color_manual(values = c("green","red"), labels=c("ctrl","RNases")) +
            theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank(), axis.title.x=element_blank(), plot.title = element_text(size = rel(1.5), colour = "darkred"), legend.background = element_rect(colour = "black"),legend.key = element_rect(fill = "white")) +
            
            # panel options
            theme(panel.background = element_rect(fill = "white", color = "darkred", size = 3), panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey")) +
            
            # x-axis scale and ticks
            scale_x_discrete(name ="Fractions", limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))
          
          
        } else {
          
          #***************************************************************************************************
          if ((table_ctrl_mean[pos, "pb_fit"] == "TRUE") && (table_rnase_mean[pos, "pb_fit"] == "TRUE")) {
            
            scatterPlot <- plot_CR + 
              
              # curve CTRL
              geom_line(data=df_ctrl, aes(x=Fractions, y=Amount, color="green"), size=.5) +
              geom_point(data=df_ctrl, mapping=aes(x=Fractions, y=Amount), color="green", shape=1, size=1) +
              
              # shadow CTRL
              geom_area(data=df_ctrl, aes(x=Fractions, y=Amount), fill="green", alpha=.1, show.legend=F) +
              
              # shadded errors
              geom_ribbon(data=df_sd_ctrl,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
              
              # curve RNASE
              geom_line(data=df_rnase, aes(x=Fractions, y=Amount, color="red"), size=.5) +
              geom_point(data=df_rnase, mapping=aes(x=Fractions, y=Amount), color="red", shape=2, size=1) +
              
              # shadded errors
              geom_ribbon(data=df_sd_rnase,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
              
              # shadow RNASE
              geom_area(data=df_rnase, aes(x=Fractions, y=Amount), fill="red", alpha=.1, show.legend=F) +
              
              # fit CTRL
              # plot_CR +
              
              # fit RNASE
              # plot_CR +
              
              # New axis label
              ylab("normalized protein amount") +
              
              # title
              labs(title = title) +
              
              # Legend
              scale_color_manual(values = c("green","red"), labels=c("ctrl","RNases")) +
              theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank(), axis.title.x=element_blank(), plot.title = element_text(size = rel(1.5), colour = "darkred"), legend.background = element_rect(colour = "black"),legend.key = element_rect(fill = "white")) +
              
              # panel options
              theme(panel.background = element_rect(fill = "white", color = "darkred", size = 3), panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey")) +
              
              # x-axis scale and ticks
              scale_x_discrete(name ="Fractions", limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))
            
            
          } else {
            
            
            #***************************************************************************************************
            # if everything is OK ;-)
            
            scatterPlot <- ggplot() +
              
              # curve CTRL
              geom_line(data=df_ctrl, aes(x=Fractions, y=Amount, color="green"), size=.5) +
              geom_point(data=df_ctrl, mapping=aes(x=Fractions, y=Amount), color="green", shape=1, size=1) +
              
              # shadow CTRL
              geom_area(data=df_ctrl, aes(x=Fractions, y=Amount), fill="green", alpha=.1, show.legend=F) +
              
              # shadded errors
              geom_ribbon(data=df_sd_ctrl,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
              
              # curve RNASE
              geom_line(data=df_rnase, aes(x=Fractions, y=Amount, color="red"), size=.5) +
              geom_point(data=df_rnase, mapping=aes(x=Fractions, y=Amount), color="red", shape=2, size=1) +
              
              # shadded errors
              geom_ribbon(data=df_sd_rnase,aes(x=Fractions,ymin=lower,ymax=upper),alpha=0.3) +
              
              # shadow RNASE
              geom_area(data=df_rnase, aes(x=Fractions, y=Amount), fill="red", alpha=.1, show.legend=F) +
              
              # fit CTRL
              geom_line(data=dffit_ctrl, aes(x=Fractions, y=Amount), color="darkgreen", size=1) +
              
              # fit RNASE
              geom_line(data=dffit_rnase, aes(x=Fractions, y=Amount), color="darkred", size=1) +
              
              # New axis label
              ylab("normalized protein amount") +
              
              # title
              labs(title = title) +
              
              # Legend
              scale_color_manual(values = c("#09b957","#f42730"), labels=c("Control","RNase")) +
              theme(legend.position=c(0.02,0.98),
                    legend.direction="horizontal",
                    legend.title.align = '0',
                    legend.justification=c(0,1),
                    legend.title=element_blank(),
                    axis.title.x=element_blank(),
                    plot.title = element_text(size = rel(1.5),
                                              colour = "black",
                                              hjust = 0.5, face = "bold"),
                    legend.background = element_rect(colour = "black"),
                    legend.key = element_rect(fill = "white")) +
              
              # panel options
              theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey"), axis.title.y = element_text(color="black", face="bold"), plot.margin = margin(t=1, r=0, b=1, l=0, "cm"))  +
              
              # x-axis scale and ticks
              scale_x_discrete(name ="Fractions", limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))
            
            
            
          }		
        }
      }		
    }				
    
    
    # plot a raster with ggplot
    dataf1 <- data.frame(ctrl_norm_mean[pos,1:25],rnase_norm_mean[pos,1:25])
    r1 <- raster(xmn=0, xmx = 25, ymn = 0, ymx = 1, nrows = 1, ncols = 25)
    r2 <- raster(xmn=0, xmx = 25, ymn = 0, ymx = 1, nrows = 1, ncols = 25)
    r1[] <- as.numeric(dataf1[1,1:25])
    r2[] <- as.numeric(dataf1[1,26:50])
    r.spdf1 <- as(r1, "SpatialPixelsDataFrame")
    r1.df <- as.data.frame(r.spdf1)
    colnames(r1.df) <- c("Amount","Fraction","ctrl")
    r.spdf2 <- as(r2, "SpatialPixelsDataFrame")
    r2.df <- as.data.frame(r.spdf2)
    colnames(r2.df) <- c("Amount","Fractions","RNases")
    
    # now plot the whole
    plot1 <- ggplot(r1.df, aes(x=Fraction+0.5, y=ctrl)) + geom_tile(aes(fill = Amount)) + coord_equal(2) + theme(legend.position="none") + scale_fill_gradient(low="white", high="green") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank()) + scale_x_continuous(breaks=c(seq(1,25,by=1)))
    plot2 <- ggplot(r2.df, aes(x=Fractions+0.5, y=RNases)) + geom_tile(aes(fill = Amount)) + coord_equal(2) + theme(legend.position="none") + scale_fill_gradient(low="white", high="red") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank()) + scale_x_continuous(breaks=c(seq(1,25,by=1)))
    
    # combine all graphics in one panel
    grid.arrange(scatterPlot, plot1, plot2, ncol=1, nrow=3, heights=c(1.5,.5,.5), bottom = "Fractions")
    
    
    
    # end from "if (pos %in% rownames(table))"
    # if the protein is not listed in the data 
  } else {
    scatterPlot <- ggplot() + geom_text(aes(x=5, y=5, label=paste(pos, "not listed")), color="darkred", size=10) +
      theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
    scatterPlot
  }
  
  # end of the function "curves_plot()" - curves_plot("CTCF_HUMAN")
}


#**************************************************************************************************
#### GRAPHICS 2 - all curves for one protein of interest                                       ####
#**************************************************************************************************

# Load packages
library(sp)
library(raster)
library(grid)
library(ggplot2)
library(gridExtra)
library(lattice)

# Function to plot all profiles (with sliding window and normalized) of one protein -> curves_all(pos)
#***************************************************************************************************
curves_all <- function(pos) {

# check if protein (pos) is listed in the table
if (pos %in% rownames(table)) {
  
datafa <- data.frame(table.ctrl1.SW.norm[pos,1:25],table.rnase1.SW.norm[pos,1:25])
datafb <- data.frame(table.ctrl2.SW.norm[pos,1:25],table.rnase2.SW.norm[pos,1:25]) 
datafc <- data.frame(table.ctrl3.SW.norm[pos,1:25],table.rnase3.SW.norm[pos,1:25])

title <- rownames(datafa)

#**************************************************************************************************
# shadded curves 
x <- c(1:25, 1:25)
y1 <- as.numeric(datafa[])
y2 <- as.numeric(datafb[])
y3 <- as.numeric(datafc[])

group <- as.factor(rep(c(1,2), each=25))
df1 <- data.frame(x, y1, group)
df2 <- data.frame(x, y2, group)
df3 <- data.frame(x, y3, group)

colnames(df1) <- c("Fraction","amount","Group")
colnames(df2) <- c("Fraction","amount","Group")
colnames(df3) <- c("Fraction","amount","Group")

scatterPlot1 <- ggplot(df1,aes(Fraction, amount, color=Group, shape=Group)) +   geom_area(data=df1[c(1:25),], fill="green", alpha=.1, show.legend=F) + geom_area(data=df1[c(26:50),], fill="red", alpha=.1, show.legend=F) +
 scale_color_manual(values = c("green","red"), labels=c("ctrl","RNases")) + geom_point() + scale_shape_discrete(labels=c("ctrl","RNases")) +
 theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank(),axis.title.x=element_blank()) + 
  # x-axis scale and ticks
  scale_x_discrete(name ="Fractions", limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))

scatterPlot2 <- ggplot(df2,aes(Fraction, amount, color=Group, shape=Group)) +   geom_area(data=df2[c(1:25),], fill="green", alpha=.1, show.legend=F) + geom_area(data=df2[c(26:50),], fill="red", alpha=.1, show.legend=F) +
 scale_color_manual(values = c("green","red"), labels=c("ctrl","RNases")) + geom_point() + scale_shape_discrete(labels=c("ctrl","RNases")) +
 theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank(),axis.title.x=element_blank()) +
  # x-axis scale and ticks
  scale_x_discrete(name ="Fractions", limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))

scatterPlot3 <- ggplot(df3,aes(Fraction, amount, color=Group, shape=Group)) +   geom_area(data=df3[c(1:25),], fill="green", alpha=.1, show.legend=F) + geom_area(data=df3[c(26:50),], fill="red", alpha=.1, show.legend=F) +
 scale_color_manual(values = c("green","red"), labels=c("ctrl","RNases")) + geom_point() + scale_shape_discrete(labels=c("ctrl","RNases")) +
 theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank(),axis.title.x=element_blank()) +
  # x-axis scale and ticks
  scale_x_discrete(name ="Fractions", limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))

# combine all graphics in one panel
grid.arrange(scatterPlot1, scatterPlot2,scatterPlot3, ncol=1, nrow=3, heights=c(1,1,1), top = title, bottom = "fractions")							

# if the protein is not listed in the data
} else {
scatterPlot <- ggplot() + geom_text(aes(x=5, y=5, label=paste(pos, "not listed")), color="darkred", size=10) +
theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
scatterPlot
}										
								}		#end of the function




