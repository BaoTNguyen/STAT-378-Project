library(tidyverse)
library(ggsignif)
library(dplyr)
library(readxl)
library(chisq.posthoc.test)
library(reshape)
library(zoo)
library(olsrr)

data <- read.table("C:\\Users\\baotn\\Uni Stuff\\University Subjects\\STAT\\STAT 378\\Project\\Term Paper -  Fall 2023 - Data.txt")

factors <- c('id',
               'area',
               'total pop',
               'central pop',
               '65 or older pop',
               'active physicians',
               'hospital beds',
               '% graduates',
               'labor force',
               'income',
               'crimes',
               'region')

colnames(data) <- factors

y <- data$id
x1 <- data$area
x2 <- data$`total pop`
x3 <- data$`central pop`
x4 <- data$`65 or older pop`
x5 <- data$`active physicians`
x6 <- data$`hospital beds`
x7 <- data$`labor force`
x8 <- data$income
x9 <- data$crimes
x10 <- data$region
x11 <- data$`% graduates`

full_model <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11)
ols_step_all_possible(full_model)
anova(full_model)
# check for independence
# check for similar variance
# check for linearity
# check for normality
