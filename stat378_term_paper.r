# install.packages("reshape2")
library(reshape2)
library(car)
library(qpcR)
data = read.table("C:\\Users\\adhig\\OneDrive\\Desktop\\University Courses\\STAT 378\\Term_Paper_-_DATA.txt", header = TRUE) # nolint
attach(data)
View(data) # nolint
# plot(data, main="Scatterplot of variables in Data set") # nolint

#######################################################################################################################################################################################################
# SPLITTING DATA USING DUPLEX ALGORITHM
s1 = (length(ID) - 1) * var(ID)
s2 = (length(LA) - 1) * var(LA)
s3 = (length(TP) - 1) * var(TP)
s4 = (length(PPCC) - 1) * var(PPCC)
s5 = (length(PPA) - 1) * var(PPA)
s6 = (length(NAP) - 1) * var(NAP)
s7 = (length(NHB) - 1) * var(NHB)
s8 = (length(PHSG) - 1) * var(PHSG)
s9 = (length(CLF) - 1) * var(CLF)
s10 = (length(TPI) - 1) * var(TPI)
s11 = (length(TSC) - 1) * var(TSC)
s12 = (length(GR) - 1) * var(GR)
Z=cbind((ID-mean(ID))/sqrt(s1), (LA-mean(LA))/sqrt(s2), (TP-mean(TP))/sqrt(s3), (PPCC-mean(PPCC))/sqrt(s4), (PPA-mean(PPA))/sqrt(s5), (NAP-mean(NAP))/sqrt(s6), (NHB-mean(NHB))/sqrt(s7), (PHSG-mean(PHSG))/sqrt(s8), (CLF-mean(CLF))/sqrt(s9), (TPI-mean(TPI))/sqrt(s10), (TSC-mean(TSC))/sqrt(s11), (GR-mean(GR))/sqrt(s12)) # nolint
#t(Z)%*%Z
T=chol(t(Z)%*%Z)
#t(T)%*%T
#solve(T)
W=Z%*%solve(T)
distances = dist(W)
# Creating a dataframe with all distances and observation pairs using reshape2 package
distances_df <- melt(as.matrix(distances), varnames = c("row_num", "col_num"))
# Removing all observations whose distance is 0
distances_df <- subset(distances_df, value != 0)
# Creating two lists to contain the observation number for prediction and estimation points
estimation_index <- c()
prediction_index <- c()
# Assign 24 pairs of observations equally into each data set
for (i in 1:24){
    # max_index is the index for the observations with max distance
    max_index = which.max(distances_df$value)
    # obs1 and obs2 are the observations with max_distance
    obs1 = distances_df[max_index, "row_num"] # row_num
    obs2 = distances_df[max_index, "col_num"] # col_num
    if ((i %% 2) == 1){
        # if i is odd, add observation to estimation data
        estimation_index <- append(estimation_index, obs1)
        estimation_index <- append(estimation_index, obs2)
    } else {
        # if i is even, add observation to prediction data
        prediction_index <- append(prediction_index, obs1)
        prediction_index <- append(prediction_index, obs2)
    }
    # Update the dataframe by removing all observations with obs1 or obs2
    distances_df <- subset(distances_df, row_num != obs1 & row_num != obs2)
    distances_df <- subset(distances_df, col_num != obs1 & col_num != obs2)
}
# Adding one observation each of the next pair to each data set
max_index = which.max(distances_df$value)
obs1 = distances_df[max_index, "row_num"] # row_num
obs2 = distances_df[max_index, "col_num"] # col_num
estimation_index <- append(estimation_index, obs1)
prediction_index <- append(prediction_index, obs2)
distances_df <- subset(distances_df, row_num != obs1 & row_num != obs2)
distances_df <- subset(distances_df, col_num != obs1 & col_num != obs2)

# Creating data frames for prediction and estimation data sets
pred_data <- data[prediction_index,]
est_data <- data[-prediction_index,]

#########################################################################################################################################################################################################

attach(est_data)
full_model = lm(TSC~ID+LA+TP+PPCC+PPA+PHSG+TPI+NAP+NHB+CLF+GR, data = est_data)

