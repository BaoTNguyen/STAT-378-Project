library(olsrr)
library(car)
library(qpcR)
library(reshape2)
source('load_data.R')
data <- load_data()
attach(data)
options(max.print = .Machine$integer.max)

colnames(data) <- cbind('ID',
                        'LA',
                        'TP',
                        'PPCC',
                        'PPA',
                        'NAP',
                        'NHB',
                        'PHSG',
                        'CLF',
                        'TPI',
                        'TSC',
                        'GR')

# Creating three indicator variables to represent GR variable consisting of 4 categories
# We choose 4==W as our reference category for Geographical Region
NE <- ifelse(data$GR == 1, 1, 0)
NC <- ifelse(data$GR == 2, 1, 0)
S <- ifelse(data$GR == 3, 1, 0)
data <- cbind(data, NE, NC, S)

# Correlation matrix
cor(data)

# DUPLEX Algorithm
s1 = (length(data$S) - 1) * var(data$S)
s2 = (length(data$LA) - 1) * var(data$LA)
s3 = (length(data$TP) - 1) * var(data$TP)
s4 = (length(data$PPCC) - 1) * var(data$PPCC)
s5 = (length(data$PPA) - 1) * var(data$PPA)
s6 = (length(data$NAP) - 1) * var(data$NAP)
s7 = (length(data$NHB) - 1) * var(data$NHB)
s8 = (length(data$PHSG) - 1) * var(data$PHSG)
s9 = (length(data$CLF) - 1) * var(data$CLF)
s10 = (length(data$TPI) - 1) * var(data$TPI)
s11 = (length(data$TSC) - 1) * var(data$TSC)
s12 = (length(data$NE) - 1) * var(data$NE)
s13 = (length(data$NC) - 1) * var(data$NC)

Z = cbind((data$S-mean(data$S))/sqrt(s1), (data$LA-mean(data$LA))/sqrt(s2), (data$TP-mean(data$TP))/sqrt(s3), (data$PPCC-mean(data$PPCC))/sqrt(s4), (data$PPA-mean(data$PPA))/sqrt(s5), (data$NAP-mean(data$NAP))/sqrt(s6), (data$NHB-mean(data$NHB))/sqrt(s7), (data$PHSG-mean(data$PHSG))/sqrt(s8), (data$CLF-mean(data$CLF))/sqrt(s9), (data$TPI-mean(data$TPI))/sqrt(s10), (data$TSC-mean(data$TSC))/sqrt(s11), (data$NE-mean(data$NE))/sqrt(s12), (data$NC-mean(data$NC))/sqrt(s13))
T=chol(t(Z)%*%Z)
W=Z%*%solve(T)
distances = dist(W)

# Creating a dataframe with all distances and observation pairs using reshape2 package
distances_df <- melt(as.matrix(distances), varnames = c("row_num", "col_num"))

# Removing all observations whose distance is 0
distances_df <- subset(distances_df, value != 0)

# Creating two lists to contain the observation number for prediction and estimation points
estimation_index <- c()
prediction_index <- c()
temp_df <- distances_df

# Assign 1 pair of observations into each data set
for (i in 1:2){
  # max_index is the index for the observations with max distance
  max_index = which.max(temp_df$value)
  
  # obs1 and obs2 are the observations with max_distance
  obs1 = temp_df[max_index, "row_num"] # row_num
  obs2 = temp_df[max_index, "col_num"] # col_num
  
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
  temp_df <- subset(temp_df, row_num != obs1 & row_num != obs2)
  temp_df <- subset(temp_df, col_num != obs1 & col_num != obs2)
}
rm(temp_df)

# Add 23 more observations to each data set
for (i in 1:(2*23)){
  if ((i %% 2) == 1){
    temp_df <- subset(distances_df, row_num %in% estimation_index)
    temp_df <- subset(temp_df, !(col_num %in% estimation_index))
    temp_df <- subset(temp_df, !(col_num %in% prediction_index))
    max_index = which.max(temp_df$value)
    obs = temp_df[max_index, "col_num"]
    estimation_index <- append(estimation_index, obs)
  } else {
    temp_df <- subset(distances_df, row_num %in% prediction_index)
    temp_df <- subset(temp_df, !(col_num %in% prediction_index))
    temp_df <- subset(temp_df, !(col_num %in% estimation_index))
    max_index = which.max(temp_df$value)
    obs = temp_df[max_index, "col_num"]
    prediction_index <- append(prediction_index, obs)
  }
  rm(temp_df)
}

# Now we create the data sets for training group and holdout group
holdout_data <- data[prediction_index,]

# All remaining observations will be used for training data
training_data <- data[-prediction_index,]

# FULL MODEL
full_model = lm(TSC ~ LA+TP+PPCC+PPA+NAP+NHB+PHSG+CLF+TPI+NE+NC+S, data = training_data)

# Checking model adequacy
# QQ Plot
e = full_model$resid
qqnorm(e, col = "blue")
qqline(e, col = 2)

# Residual versus predicted response
t = rstudent(full_model)
y_hat = full_model$fit
plot(y_hat, t, col="blue", title("Plot of R-student residuals versus predicted response"))
abline(h=0, col = "red")

# Residual versus regressor
par(mfrow = c(2,2))
#plot(training_data$ID, t, col="blue", xlab="ID", title("Plot of R-student residuals versus ID"))
#abline(h=0, col="red")
plot(training_data$LA, t, col="blue", xlab="LA", title("Plot of R-student residuals versus LA"))
abline(h=0, col="red")
plot(training_data$TP, t, col="blue", xlab="TP", title("Plot of R-student residuals versus TP"))
abline(h=0, col="red")
plot(training_data$PPCC, t, col="blue", xlab="PPCC", title("Plot of R-student residuals versus PPCC"))
abline(h=0, col="red")
par(mfrow = c(2,2))
plot(training_data$PPA, t, col="blue", xlab="PPA", title("Plot of R-student residuals versus PPA"))
abline(h=0, col="red")
plot(training_data$NAP, t, col="blue", xlab="NAP", title("Plot of R-student residuals versus NAP"))
abline(h=0, col="red")
plot(training_data$NHB, t, col="blue", xlab="NHB", title("Plot of R-student residuals versus NHB"))
abline(h=0, col="red")
plot(training_data$PHSG, t, col="blue", xlab="PHSG", title("Plot of R-student residuals versus PHSG"))
abline(h=0, col="red")
par(mfrow = c(3,2))
plot(training_data$CLF, t, col="blue", xlab="CLF", title("Plot of R-student residuals versus CLF"))
abline(h=0, col="red")
plot(training_data$TPI, t, col="blue", xlab="TPI", title("Plot of R-student residuals versus TPI"))
abline(h=0, col="red")
plot(training_data$NE, t, col="blue", xlab="GR=NE", title("Plot of R-student residuals versus GR=NE"))
abline(h=0, col="red")
plot(training_data$NC, t, col="blue", xlab="GR=NC", title("Plot of R-student residuals versus GR=NC"))
abline(h=0, col="red")
plot(training_data$S, t, col="blue", xlab="GR=S", title("Plot of R-student residuals versus GR=S"))
abline(h=0, col="red")

# Partial Regression plots
avPlots(full_model)

# Partial Residual plots
crPlots(full_model)

# PLOTS OF VARIABLES
plot(training_data$LA, training_data$TSC, title="TSC vs LA")
plot(training_data$TP, training_data$TSC, title="TSC vs TP")
plot(training_data$PPCC, training_data$TSC, title="TSC vs PPCC")
plot(training_data$PPA, training_data$TSC, title="TSC vs PPA")
plot(training_data$NAP, training_data$TSC, title="TSC vs NAP")
plot(training_data$NHB, training_data$TSC, title="TSC vs NHB")
plot(training_data$PHSG, training_data$TSC, title="TSC vs PHSG")
plot(training_data$CLF, training_data$TSC, title="TSC vs CLF")
plot(training_data$TPI, training_data$TSC, title="TSC vs TPI")

test_data_predictions = predict(full_model, holdout_data)
(cor(holdout_data$TSC, test_data_predictions))^2

TSC_t = unlist(log(training_data[11]))
LA_t = unlist(sqrt(training_data[2]))
TP_t = unlist(log(training_data[3]))
PPCC_t = unlist(training_data[4])
PPA_t = unlist(training_data[5])
NAP_t = unlist(log(training_data[6]))
NHB_t = unlist(log(training_data[7]))
PHSG_t = unlist(training_data[8])
CLF_t = unlist(log(training_data[9]))
TPI_t = unlist(log(training_data[10]))
NE_t = unlist(training_data[12])
NC_t = unlist(training_data[13])
S_t = unlist(training_data[14])

transformed_model <- lm(log(TSC) ~ sqrt(LA)+log(TP)+PPCC+PPA+log(NAP)+log(NHB)+PHSG+log(CLF)+log(TPI)+NE+NC+S, data = training_data)
e = transformed_model$resid
qqnorm(e, col = "blue")
qqline(e, col = 2)
summary(transformed_model)
predicted_test_data = predict(transformed_model, holdout_data)
(cor(holdout_data$TSC, predicted_test_data))^2

# Residual versus predicted response
t = rstudent(transformed_model)
y_hat = transformed_model$fit
plot(y_hat, t, col="blue", title("Plot of R-student residuals versus predicted response"))
abline(h=0, col = "red")

# Residual versus regressor
par(mfrow = c(2,2))
plot(LA_t, t, col="blue", xlab="LA", title("Plot of R-student residuals versus LA"))
abline(h=0, col="red")
plot(TP_t, t, col="blue", xlab="TP", title("Plot of R-student residuals versus TP"))
abline(h=0, col="red")
plot(PPCC_t, t, col="blue", xlab="PPCC", title("Plot of R-student residuals versus PPCC"))
abline(h=0, col="red")
plot(PPA_t, t, col="blue", xlab="PPA", title("Plot of R-student residuals versus PPA"))
abline(h=0, col="red")
par(mfrow = c(2,2))
plot(NAP_t, t, col="blue", xlab="NAP", title("Plot of R-student residuals versus NAP"))
abline(h=0, col="red")
plot(NHB_t, t, col="blue", xlab="NHB", title("Plot of R-student residuals versus NHB"))
abline(h=0, col="red")
plot(PHSG_t, t, col="blue", xlab="PHSG", title("Plot of R-student residuals versus PHSG"))
abline(h=0, col="red")
plot(CLF_t, t, col="blue", xlab="CLF", title("Plot of R-student residuals versus CLF"))
abline(h=0, col="red")
par(mfrow = c(2,2))
plot(TPI_t, t, col="blue", xlab="TPI", title("Plot of R-student residuals versus TPI"))
abline(h=0, col="red")
plot(training_data$NE, t, col="blue", xlab="GR=NE", title("Plot of R-student residuals versus GR=NE"))
abline(h=0, col="red")
plot(training_data$NC, t, col="blue", xlab="GR=NC", title("Plot of R-student residuals versus GR=NC"))
abline(h=0, col="red")
plot(training_data$S, t, col="blue", xlab="GR=S", title("Plot of R-student residuals versus GR=S"))
abline(h=0, col="red")

# Check Interaction Terms
# PHSG * TPI
transformed_model_PHSG_TPI <- lm(log(TSC) ~ sqrt(LA)+log(TP)+PPCC+PPA+log(NAP)+log(NHB)+PHSG+log(CLF)+log(TPI)+NE+NC+S + PHSG*log(TPI), data = training_data)
summary(transformed_model_PHSG_TPI)
anova(transformed_model, transformed_model_PHSG_TPI)

# PPCC * PPA
transformed_model_PPCC_PPA <- lm(log(TSC) ~ sqrt(LA)+log(TP)+PPCC+PPA+log(NAP)+log(NHB)+PHSG+log(CLF)+log(TPI)+NE+NC+S + PPCC*PPA, data = training_data)
summary(transformed_model_PPCC_PPA)
anova(transformed_model, transformed_model_PPCC_PPA)

# Final Full Model
final_model <- lm(log(TSC) ~ sqrt(LA)+log(TP)+PPCC+PPA+log(NAP)+log(NHB)+PHSG+log(CLF)+log(TPI)+NE+NC+S + PPCC*PPA + PHSG*log(TPI), data = training_data)
summary(final_model)

# Build Model
# Forward
null_model <- lm(log(TSC) ~ 1, data = training_data)
forward <- step(null_model, scope = formula(final_model), direction = "forward", trace = 0, k = qchisq(0.1, 1, lower = FALSE))
summary(forward)

# Backward
backward <- step(final_model, scope = formula(null_model), direction = "backward", trace = 0, k = qchisq(0.1, 1, lower = FALSE))
summary(backward)

# Stepwise
stepwise <- step(null_model, scope = formula(final_model), direction = "both", trace = 0, k = qchisq(0.1, 1, lower = FALSE))
summary(stepwise)

# Criteria
# AIC
AIC(final_model, transformed_model, forward, backward, stepwise)

# BIC
BIC(final_model, transformed_model, forward, backward, stepwise)

# Adjusted R^2
adjR2 = function(model) {
  1 - (1 - summary(model)$r.squared) * ((nrow(training_data) - 1)/(nrow(training_data) - length(model$coefficients) - 1))
}
adjR2(final_model)
adjR2(transformed_model)
adjR2(forward)
adjR2(backward)
adjR2(stepwise)

# Mallows Cp
ols_mallows_cp(transformed_model, final_model)
ols_mallows_cp(forward, final_model)
ols_mallows_cp(backward, final_model)
ols_mallows_cp(stepwise, final_model)

# PRESS
ols_press(transformed_model)
ols_press(final_model)
ols_press(forward)
ols_press(backward)
ols_press(stepwise)

# Check Collinearity
coll_transformed_model <- ols_coll_diag(transformed_model)
coll_final_model <- ols_coll_diag(final_model)
coll_forward <- ols_coll_diag(forward)
coll_backward <- ols_coll_diag(backward)
coll_stepwise <- ols_coll_diag(stepwise)

# Check Alternative Model
ID_t <- training_data$ID
LA_t <- training_data$LA
TP_t <- training_data$TP
PPCC_t <- training_data$PPCC
PPA_t <- training_data$PPA
NAP_t <- training_data$NAP
NHB_t <- training_data$NHB
PHSG_t <- training_data$PHSG
CLF_t <- training_data$CLF
TPI_t <- training_data$TPI
TSC_t <- training_data$TSC
GR_t <- training_data$GR

center_PPCC <- PPCC_t - mean(PPCC_t)
center_PPA <- PPA_t - mean(PPA_t)
center_PHSG <- PHSG_t - mean(PHSG_t)
center_NAP <- log(NAP_t) - mean(log(NAP_t))
center_TPI <- log(TPI_t) - mean(log(TPI_t))
center_PPCC_PPA <- center_PPCC * center_PPA
center_PHSG_TPI <- center_PHSG * center_TPI
center_NAP_TPI <- center_NAP * center_TPI

backward_alt <- lm(log(TSC_t) ~ PPCC_t+PPA_t+log(NAP_t)+PHSG_t+log(NHB_t)+center_PPCC_PPA + center_PHSG_TPI + center_NAP_TPI)
anova(backward_alt)
ols_coll_diag(backward_alt)
summary(backward_alt)
qqline(residuals(backward_alt), col = "red")
crPlots(backward_alt)
backward_alt_r2_pred <- 1 - ols_press(backward_alt)/sum(anova(backward_alt)$`Sum Sq`)

# Predict
transformed_model_pred <- predict(transformed_model, newdata = holdout_data)
final_model_pred = predict(final_model, newdata = holdout_data)
forward_pred <- predict(forward, newdata = holdout_data)
backward_pred <- predict(backward, newdata = holdout_data)
stepwise_pred <- predict(stepwise, newdata = holdout_data)

# Correlation
transformed_model_est_pred <- summary(transformed_model)$r.squared - cor(transformed_model_pred, holdout_data$TSC)^2
final_model_est_pred <- summary(final_model)$r.squared - cor(final_model_pred, holdout_data$TSC)^2
forward_est_pred <- summary(forward)$r.squared - cor(forward_pred, holdout_data$TSC)^2
backward_est_pred <- summary(backward)$r.squared - cor(backward_pred, holdout_data$TSC)^2
stepwise_est_pred <- summary(stepwise)$r.squared - cor(stepwise_pred, holdout_data$TSC)^2

# Correlation Matrix
cor(transformed_model_pred, holdout_data)
cor(final_model_pred, holdout_data)
cor(forward_pred, holdout_data)
cor(backward_pred, holdout_data)
cor(stepwise_pred, holdout_data)

# Fit All Data
transformed_model_all <- lm(transformed_model, data = data)
final_model_all <- lm(final_model, data = data)
forward_all <- lm(forward, data = data)
backward_all <- lm(backward, data = data)
stepwise_all <- lm(stepwise, data = data)

# Confidence Intervals
confint(backward, level = 0.95)

# R^2 Predicted
transformed_r2_pred <- 1 - ols_press(transformed_model)/sum(anova(transformed_model)$"Sum Sq")
final_r2_pred <- 1 - ols_press(final_model)/sum(anova(final_model)$"Sum Sq")
backward_r2_pred <- 1 - ols_press(backward)/sum(anova(backward)$"Sum Sq")
forward_r2_pred <- 1 - ols_press(forward)/sum(anova(forward)$"Sum Sq")
stepwise_r2_pred <- 1 - ols_press(stepwise)/sum(anova(stepwise)$"Sum Sq")

# Cooks Distance
cooks.distance(transformed_model_all)
cooks.distance(final_model_all)
cooks.distance(forward_all)
cooks.distance(backward_all)
cooks.distance(stepwise_all)

# DFFITS
dffits(transformed_model_all)
dffits(final_model_all)
dffits(forward_all)
dffits(backward_all)
dffits(stepwise_all)

# DFBETAS
dfbetas(transformed_model_all)
dfbetas(final_model_all)
dfbetas(forward_all)
dfbetas(backward_all)
dfbetas(stepwise_all)

# Covariance Ratio
covratio(transformed_model_all)
covratio(final_model_all)
covratio(forward_all)
covratio(backward_all)
covratio(stepwise_all)

# Influence Measures
# Training
influence.measures(transformed_model)
influence.measures(final_model)
influence.measures(forward)
influence.measures(backward)
influence.measures(stepwise)

# All
influence.measures(transformed_model_all)
influence.measures(final_model_all)
influence.measures(forward_all)
influence.measures(backward_all)
influence.measures(stepwise_all)

# Outliers
resid(transformed_model_all)/summary(transformed_model_all)$sigma
resid(final_model_all)/summary(final_model_all)$sigma
resid(forward_all)/summary(forward_all)$sigma
resid(backward_all)/summary(backward_all)$sigma
resid(stepwise_all)/summary(stepwise_all)$sigma
