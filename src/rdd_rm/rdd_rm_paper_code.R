
# Installing necessary packages and libraries -----------------------------

#install.packages("dplyr")
#install.packages("random")
#install.packages("stargazer")
#install.packages("reshape2")
#install.packages("readxl")
library(rddtools)
library(rdd)
library(boot)
library(random)
library(dplyr)
library(ggplot2)
library(reshape2)
library(readxl)


# Define Global Variables -----------------------------------------------

# beta is vec of constant, treatment effect and rest are coefficients that 
# are used in data generating process with the increasing degree of the 
# forcing variable.
beta <- c(0.42, 0.56 , 0.84, -3.00, 7.99, -9.01, 3.56, 6.11)

# n is the total sample of the population
n <- 100000

# k is the maximum degree of polynomial that is used in the simulations
k <-  6
set.seed(101)


# Define functions --------------------------------------------------------


# function 1: function to calculate the weights for x depending on --------

#' weights 
#' 
#' @param X : vector of forcing variable 
#' @param k : (global variable=6) degree of polynomial
#'
#' @return weight : weight matrix of dimensions k and lenght of forcing 
#'                  variable

set.seed(101)

weights <- function(X_pos, k){
  
  # number of observations above between [c,h].c is threshold and h is 
  # bandwidth value
  len <- length(X_pos)
  
  # sum Vector
  sum_vec <- matrix(rep(1,len),len,1)
  
  # vector where all elements other than the first one are 0 and the 
  # first one is 1
  e_vec <- matrix(c(1,rep(0,k)),k+1,len) 
  lm <- matrix(NA, k, len)
  
  for (i in 1:k){
    lm[i,] <-  X_pos^i
  }
  
  # matrix of K*n dimension, where K is the order of polynomial and n is number 
  # of observations in X_pos. first row equal to one and subsequent rows are 
  # the previous rows product by the respective unit of X_pos until K order of 
  # X_pos is reached 
  X_vec <-  rbind(rep(1,len),lm)  

  # obtain n (number of observations in X_pos) K*K matrix. in each matrix, the 
  # first columns first element is 1 and subsequent elements are the product of 
  # previous element by the respective X_pos unit until K order of X_pos is 
  # reached. the subsequent columns thereby are the product of the previous 
  # column by respective X_pos unit until K order of X_pos is reached.
  array1 <-  array(NA,dim=c(k+1,k+1,len))
  for(m in 1:len){
    for (i in 1:(k+1)){
      array1[,,m][i,]=X_pos[m]^(i-1)
      for(j in 1:k){
        array1[,,m][,j+1]= X_pos[m]*array1[,,m][,j]
      }
    }
  }
  
  # matrix to sum all matrices obtained above 
  sum_all <-  apply(array1,1:2,sum)
  
  # calculates the weight vector for order K
  weight <-  t(sum_vec)%*%t(e_vec)%*%solve(sum_all)%*%X_vec 
  return(weight)
}


# function 2: function for creating rdd_data ------------------------------

#' rdd_dgf
#'
#' @param scale1 : value for alpha for beta distribution
#' @param scale2: value for beta for beta distribution
#' @param mu_e : (default = 0) mean of error term for normal distribution
#' @param sd_e : standard deviation of error term for normal distribution
#' @param poly_x : degree of forcing variable
#' @param n_sample : number of observations, default is 100000
#' @param rdd : TRUE(default)- makes the data discontinuous (rdd setting)
#'              FALSE- allows for no discontinuity
#'
#' @return rdd_data : (default) a regression discontinuity dataframe with 
#'                    true Y, X up to the degree of polynomial specified 
#'                    by user, and the error term
#'         poly_data : (rdd=FALSE) a polynomial dataframe with true Y, X 
#'                     up to the degree of polynomial specified by user, 
#'                     and the error term.
rdd_dgf <- function(scale1, scale2, mu_e = 0, sd_e, poly_x, 
                    n_sample = n, rdd=TRUE){
  set.seed(101)
  X <- 2*(rbeta(n_sample, scale1, scale2)) - 1
  err <- rnorm(n_sample, mu_e, sd_e)
  
  # x_matrix stores the values of x^order_of_poly * beta(of order of 
  # polynomial)
  x_matrix <- matrix(NA, n_sample, poly_x)
  # x_poly just stores the x^order_of_poly
  x_poly <- matrix(NA, n_sample, poly_x)
  
  for (i in 1:poly_x) {
    x_matrix[,i] <- beta[i+2] * I(X^i)
    x_poly[,i] <- I(X^i)
  }
  
  D <- ifelse(X >= 0,1,0)
  if (rdd==TRUE) {
    # if rdd is TRUE, Y is discontinuous at the threshold
    Y <- beta[1] + beta[2]*D + rowSums(data.frame(x_matrix)) + err
    rdd_data <- data.frame(Y, D, x_poly, err)
    return(rdd_data)
  } else{
    # if not then the function renders continuous Y as a funciton of X
    Y <- beta[1] + rowSums(data.frame(x_matrix)) + err
    poly_data <- data.frame(Y, x_poly, err)
    return(poly_data)
  }
}


# function 3: simulation study for estimates and standard errors ----------

#' est_sd_sim
#'
#' @param scale_x1 : value for alpha for beta distribution
#' @param scale_x2: value for beta for beta distribution
#' @param e_mean : (default = 0)mean of e that is given to the 
#'                 function rdd_dgf
#' @param sd_e_ : standard deviation of error term that is given to the
#'                function rdd_dgf
#' @param rep : number of replications for simulation study
#' @param n_sample : number of samples in each simulation study
#' @param x_poly : degree of polynomial of forcing variable that is used 
#'                 in data generating function.
#'
#' @return ss : data frame of average of estimates and standard errors 
#'              of global and local polynomials.
est_sd_sim <- function(scale_x1, scale_x2, e_mean = 0, sd_e_, rep, 
                       n_sample, x_poly){
  # to store estimates and standard errors of global and local polynomials 
  est <- matrix(NA,k,rep)
  se <- matrix(NA,k,rep)
  est_local <- matrix(NA,k,rep)
  se_local <- matrix(NA,k,rep)
  
  # original data created by rdd_dgf described in above function
  data1 <- rdd_dgf(scale1 = scale_x1, scale2 = scale_x2, mu_e = e_mean, 
                   sd_e = sd_e_, poly_x = x_poly)
  
  for(j in 1:rep){
    # slice the data for each iteration from original data of n sample size
    data_sim <- slice_sample(data1, n = n_sample, replace = TRUE) 
    
    #calculate bandwidth using Imbens and Kalyanaraman Bandwidth selector
    bandwidth <- IKbandwidth(data_sim[,3],data_sim[,1],cutpoint=0, 
                             kernel = "triangular")
    
    # trim the sampled data based on bandwidth for local linear and 
    # quadratic regressions
    data2 <-  filter(data_sim,0-bandwidth<=X1 & X1<=0+bandwidth)
    
    # for up to six orders of global polynomial
    for (i in 1:k){
      # define the linear regression taking rdd in account and store the
      # estimates and standard errors for each iteration and degree of
      # global polynomial.
      lm_model <- lm(data_sim[,1]~data_sim[,2]*poly(data_sim[,3],i,raw = T), 
                     data = data_sim)
      est[i,j] <- lm_model$coefficients[2]
      se[i,j] <-sqrt(diag(vcov(lm_model)))[2]
    } 
    
    # for up to two orders of local polynomial
    for (i in 1:2){
      # define the linear regression taking rdd in account and store the
      # estimates and standard errors for each iteration and degree of
      # local polynomial.
      lm_model <- lm(data2[,1]~data2[,2]*poly(data2[,3],i,raw =T),
                     data=data2)
      est_local[i,j] <-lm_model$coefficients[2]
      se_local[i,j] <-sqrt(diag(vcov(lm_model)))[2]
    }
  }
  
  # store the averages estimates and standard errors of the treatment 
  # effect of global polynomials in sensitivity data frame.
  sensitivity <- data.frame("Order_of_Polynomial"=c(1:k),
                           "estimate"=round(rowMeans(est),3),
                           "std_error"=round(rowMeans(se),3))
  
  # store the averages estimates and standard errors of the treatment 
  # effect of local polynomials in sensitivity data frame.
  sensitivity_local <- data.frame("estimate"=round(rowMeans(est_local),3),
                                 "std_error" = round(rowMeans(se_local),3))
  
  # create one data frame to return both results
  ss <- bind_cols(sensitivity, sensitivity_local)
  colnames(ss) <- c("Order_of_Polynomial","est_global","se_global",
                    "est_local","se_local")
  
  return(list(ss = ss))
}


# function 4: Coverage Probability and rejection rate simulation ----------

#' rej_rate_sim
#'
#' @param n_threshold : number of replications for the simulations study
#' @param n_sample : number of samples in each simulation
#' @param poly : degree of polynomial of forcing variable used to to 
#'               generate polynomial_data using rdd_dgf
#' @param scale_x1 : value for scale 1 for beta distribution
#' @param scale_x2: value for scale 2 for beta distribution
#' @param sd_e_ : standard deviation of error term used in rdd_dgf
#'
#' @return median_se : median of standard errors for global and 
#'                     local polynomials
#'         rejectrate : rejection rate for the global and local 
#'                      polynomials
#'         estimate_sim : estimates of the treatment effect for every
#'                        iteration for the global and local polynomials. 
rej_rate_sim <- function(n_threshold, n_sample, poly, scale_x1, scale_x2, 
                         sd_e_){
  # create continuous data using rdd_dgf as rdd is FALSE
  data_poly <- rdd_dgf(scale1 = scale_x1, scale2 = scale_x2, sd_e = sd_e_, 
                       poly_x = poly, rdd = FALSE)
  
  # choose pseudo-threshold between 0.25 adn 0.75 percentile of forcing
  # variable
  q.25 <- quantile(data_poly[,2])[2]
  q.75 <- quantile(data_poly[,2])[4]
  thresholdVec <- runif(n_threshold, q.25, q.75)
  
  # containers to store x_normalized, dummy variable, estimates and standard
  # errors of treatment effect, upper confidence intervals, lower confidence 
  # intervals and rejections for global and local polynomails. 
  X_normalised_sim <- matrix(NA, n_sample, n_threshold)
  D_sim <- matrix(NA, nrow = n_sample, ncol = n_threshold)
  estimate_sim <- matrix(NA, k+2,n_threshold)
  std_error_sim <- matrix(NA,k+2,n_threshold)
  CI_lower <- matrix(NA,k,n_threshold)
  CI_upper <- matrix(NA,k,n_threshold)
  reject <-  matrix(NA,k,n_threshold)
  CI1_local_upper <-  matrix(NA,2,n_threshold)
  CI1_local_lower <-  matrix(NA,2,n_threshold)
  reject1 <-  matrix(NA,2,n_threshold)
  
  for (i in 1:n_threshold){
    # get a sample of the data for each iteration
    data <- slice_sample(data_poly, n = n_sample)
    
    # normalize the forcing variable by subtracting pseudo-threshold
    X_normalised_sim[,i] <- data[,2] - thresholdVec[i]
    
    # create and store the dummy variable for the same
    D_sim[,i] <- ifelse(data[,2] >= thresholdVec[i],1,0)
    
    # create a data frame for each iteration to pass in lm_model()
    df <- data.frame(y = data[,1], d = D_sim[,i], x_norm = 
                       X_normalised_sim[,i], x = data[,2])
    
    # for up to six orders of global polynomial
    for (j in 1:k) {
      # define the linear regression and store the estimates, standard errors,
      # upper and lower confidence intervals for each iteration and degree of
      # global polynomial
      lm_model <- lm(y~d*poly(x_norm,j,raw =T),data=df)
      estimate_sim[j,i] <- lm_model$coefficients[2]
      std_error_sim[j,i] <- sqrt(diag(vcov(lm_model)))[2]
      CI_lower[j,i] <- confint(lm_model)[2,1]
      CI_upper[j,i] <- confint(lm_model)[2,2]
      
      # check if 0 actually lies in the corresponsing confidence interval
      reject[j,i] = ifelse((0<CI_lower[j,i] | 0>CI_upper[j,i]),1,0)
    }
    
    #calculate bandwidth using Imbens and Kalyanaraman Bandwidth selector
    bandwidth1 <- IKbandwidth(data[,2],data[,1],cutpoint = thresholdVec[i],
                                 kernel = "triangular")
    
    # trim the sampled data based on bandwidth for local linear and 
    # quadratic regressions
    data2 <-  filter(df,thresholdVec[i]-bandwidth1<=df[,4] & 
                       df[,4]<=thresholdVec[i]+bandwidth1)
    
    # for up to two orders of local polynomial
    for (j in 1:2){
      # define the linear regression and store the estimates, standard errors,
      # upper and lower confidence intervals for each iteration and degree of
      # local polynomial
      lm2 <- lm(data2[,1] ~ data2[,2]*poly(data2[,3],j,raw=T),data=data2)
      CI1_local_lower[j,i] <- confint(lm2)[2,1]
      CI1_local_upper[j,i] <- confint(lm2)[2,2]
      estimate_sim[j+6,i] <- lm2$coefficients[2]
      std_error_sim[j+6,i] <- sqrt(diag(vcov(lm2)))[2]
      
      # check if 0 actually lies in the corresponsing confidence interval
      reject1[j,i] <- ifelse((0<CI1_local_lower[j,i] | 
                                0>CI1_local_upper[j,i]),1,0)
    }
  }
  
  # create a matrix to store the rejection rates of global and local
  # polynomials
  rejectrate <- matrix(NA,8,1)
  for (j in 1:2){
    rejectrate[j+6,] = sum(reject1[j,])/n_threshold
  }
  for(j in 1:6){
    rejectrate[j,] = sum(reject[j,])/n_threshold
  }
  
  # store the median standard errors of global and local polynomials
  sd_sim <- data.frame(std_error_sim)
  median_se <- rowMeans(sd_sim)
  
  # return the median standard errors, rejection rate and estimates of
  # treatment effect for each iteration.
  return(list(median_se = median_se, rejectrate = rejectrate, 
         estimate_sim=estimate_sim))
}


# function 5 --------------------------------------------------------------

# this function is exactly same as function 4, what is changed is that it
# does not generates data using rdd_dgf function instead takes a defined 
# dataset with first column as dependent variable and second column as 
# forcing variable. And then performs all the same steps as function 4 does.

#' rej_rate_emp
#'
#' @param n_threshold : number of replications for the simulations study
#' @param n_sample : number of samples in each simulation
#' @param data : takes a dataframe with first column as dependent variable
#'               and second column as forcing variable
#'
#' @return median_se : median of standard errors for global and 
#'                     local polynomials
#'         rejectrate : rejection rate for the global and local 
#'                      polynomials
#'         estimate_sim : estimates of the treatment effect for every
#'                        iteration for the global and local polynomials. 
rej_rate_emp <- function(n_threshold, n_sample, data){
  
  # takes the data and store it in a vaiable
  data_poly <- data
  
  # choose pseudo-threshold between 0.25 adn 0.75 percentile of forcing
  # variable
  q.25 <- quantile(data_poly[,2])[2]
  q.75 <- quantile(data_poly[,2])[4]
  thresholdVec <- runif(n_threshold, q.25, q.75)
  
  # containers to store x_normalized, dummy variable, estimates and standard
  # errors of treatment effect, upper confidence intervals, lower confidence 
  # intervals and rejections for global and local polynomails. 
  X_normalised_sim <- matrix(NA, n_sample, n_threshold)
  D_sim <- matrix(NA, nrow = n_sample, ncol = n_threshold)
  estimate_sim <- matrix(NA, k+2,n_threshold)
  std_error_sim <- matrix(NA,k+2,n_threshold)
  CI_lower <- matrix(NA,k,n_threshold)
  CI_upper <- matrix(NA,k,n_threshold)
  reject <-  matrix(NA,k,n_threshold)
  CI1_local_upper <-  matrix(NA,2,n_threshold)
  CI1_local_lower <-  matrix(NA,2,n_threshold)
  reject1 <-  matrix(NA,2,n_threshold)
  
  for (i in 1:n_threshold){
    # get a sample of the data for each iteration
    data <- slice_sample(data_poly, n = n_sample)
    
    # normalize the forcing variable by subtracting pseudo-threshold
    X_normalised_sim[,i] <- data[,2] - thresholdVec[i]
    
    # create and store the dummy variable for the same
    D_sim[,i] <- ifelse(data[,2] >= thresholdVec[i],1,0)
    
    # create a data frame for each iteration to pass in lm_model()
    df <- data.frame(y = data[,1], d = D_sim[,i], x_norm = 
                       X_normalised_sim[,i], x = data[,2])
    
    # for up to six orders of global polynomial
    for (j in 1:k) {
      # define the linear regression and store the estimates, standard errors,
      # upper and lower confidence intervals for each iteration and degree of
      # global polynomial
      lm_model <- lm(y~d*poly(x_norm,j,raw =T),data=df)
      estimate_sim[j,i] <- lm_model$coefficients[2]
      std_error_sim[j,i] <- sqrt(diag(vcov(lm_model)))[2]
      CI_lower[j,i] <- confint(lm_model)[2,1]
      CI_upper[j,i] <- confint(lm_model)[2,2]
      
      # check if 0 actually lies in the corresponding confidence interval
      reject[j,i] = ifelse((0<CI_lower[j,i] | 0>CI_upper[j,i]),1,0)
    }
    
    #calculate bandwidth using Imbens and Kalyanaraman Bandwidth selector
    bandwidth1 <- IKbandwidth(data[,2],data[,1],cutpoint = thresholdVec[i],
                              kernel = "triangular")
    
    # trim the sampled data based on bandwidth for local linear and 
    # quadratic regressions
    data2 <-  filter(df,thresholdVec[i]-bandwidth1<=df[,4] & 
                       df[,4]<=thresholdVec[i]+bandwidth1)
    
    # for up to two orders of local polynomial
    for (j in 1:2){
      # define the linear regression and store the estimates, standard errors,
      # upper and lower confidence intervals for each iteration and degree of
      # local polynomial
      lm2 <- lm(data2[,1] ~ data2[,2]*poly(data2[,3],j,raw=T),data=data2)
      CI1_local_lower[j,i] <- confint(lm2)[2,1]
      CI1_local_upper[j,i] <- confint(lm2)[2,2]
      estimate_sim[j+6,i] <- lm2$coefficients[2]
      std_error_sim[j+6,i] <- sqrt(diag(vcov(lm2)))[2]
      
      # check if 0 actually lies in the corresponsing confidence interval
      reject1[j,i] <- ifelse((0<CI1_local_lower[j,i] | 
                                0>CI1_local_upper[j,i]),1,0)
    }
  }
  
  # create a matrix to store the rejection rates of global and local
  # polynomials
  rejectrate <- matrix(NA,8,1)
  for (j in 1:2){
    rejectrate[j+6,] = sum(reject1[j,])/n_threshold
  }
  for(j in 1:6){
    rejectrate[j,] = sum(reject[j,])/n_threshold
  }
  
  # store the median standard errors of global and local polynomials
  sd_sim <- data.frame(std_error_sim)
  median_se <- rowMeans(sd_sim)
  
  # return the median standard errors, rejection rate and estimates of
  # treatment effect for each iteration.
  return(list(median_se = median_se, rejectrate = rejectrate, 
              estimate_sim=estimate_sim))
}


# Argument 1: Noisy Weights -----------------------------------------------

# SIMULATION STUDY

# global polynomial rdd regression
# dataset
data <-  rdd_dgf(scale1 = 2, scale2 = 5, sd_e = 0.1295, poly_x = 6) 

# Observations between [c,h]. c=0 & h=infinity
X_pos = sort(subset(data[,3],data[,3]>=0))
mu1 = weights(X_pos,k=1)  #Weight for K = 1
mu2 = weights(X_pos,k=2)  #Weight for K = 2
mu3 = weights(X_pos,k=3)  #Weight for K = 3
mu4 = weights(X_pos,k=4)  #Weight for K = 4
mu5 = weights(X_pos,k=5)  #Weight for K = 5
mu6 = weights(X_pos,k=6)  #Weight for K = 6

# merge them into dataframe
weight.frame = data.frame(mu1,mu2,mu3,mu4,mu5,mu6,X_pos)

# save plot of weights
arg1_sim = ggplot(weight.frame,aes(x=X_pos))+
  geom_line(aes(y=mu1,colour = "1"))+
  geom_line(aes(y=mu2,colour = "2"))+
  geom_line(aes(y=mu3,colour = "3"))+
  geom_line(aes(y=mu4,colour = "4"))+
  geom_line(aes(y=mu5,colour = "5"))+
  geom_line(aes(y=mu6,colour = "6"))+
  xlab("Forcing Variable")+
  ylab("Weights")+
  scale_color_manual(name = "Degree of Polynomial",
                     values =c("1"= "green","2"= "blue",
                     "3"="red","4"= "magenta","5" ="yellow","6"= "orange"))+
  theme_classic()+
  theme(legend.background=element_rect(colour="black"),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

# save the graph in current working directory
tiff("arg1_sim", units="in", width=6, height=5, res=500)
arg1_sim
dev.off()

# local linear rdd regression
# bandwidth Vector
bandwidth = matrix(NA,1,2) 

# bandwidth using triangular kernel
bandwidth[,1] =  IKbandwidth(data[,3], data[,1], cutpoint = 0, 
                             kernel="triangular") 

# bandwidth using rectangular kernel
bandwidth[,2] =  IKbandwidth(data[,3], data[,1], cutpoint = 0, 
                             kernel = "rectangular") 

# observations between [c,h] for triangular kernel
X_pos1 = sort(subset(data[,3], 0<=data[,3] & data[,3]<=bandwidth[,1])) 

# observations between [c,h] for rectangular kernel
X_pos2 = sort(subset(data[,3], 0<=data[,3] & data[,3]<=bandwidth[,2])) 

# weights for K = 1, 2 for triangular kernel and K = 1 for rectangular
# kernel respectively and store them in a dataframe
mu7 = matrix(c(weights(X_pos1,k=1),rep(0,length(X_pos) - length(X_pos1)))) 
mu8 = matrix(c(weights(X_pos1,k=2),rep(0,length(X_pos) - length(X_pos1)))) 
mu9 = matrix(c(weights(X_pos2,k=1),rep(0,length(X_pos) - length(X_pos2))))
weight.local.frame = data.frame(mu7,mu8,mu9,X_pos)

# save the plot of weights
arg1_sim_local = ggplot(weight.local.frame,aes(x=X_pos))+
  geom_line(aes(y=mu7,colour = "Linear and Triangular"))+
  geom_line(aes(y=mu8,colour = "Quadratic and Triangular"))+
  geom_line(aes(y=mu9,colour = "Linear and Rectangular"))+
  geom_hline(yintercept = 0)+ xlab("Forcing Variable")+ylab("Weights")+
  theme_classic()+
  scale_color_manual(name = "Degree of Polynomial and Kernel Type",
                     values =c("Linear and Triangular"= "green",
                               "Quadratic and Triangular"= "blue",
                               "Linear and Rectangular"="red"))+
  theme_classic()+
  theme(legend.background=element_rect(colour="black"),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

# save the graph in current working directory
tiff("arg1_sim_local", units="in", width=8, height=5, res=500)
arg1_sim_local
dev.off()

# make a weights table for largest value of X_pos for global high-order 
# polynomial case
round(X_pos[which.max(X_pos)],3)

# resulting Table
table1 = matrix(rbind(c(1:6),round(c(mu1[which.max(X_pos)], 
                                     mu2[which.max(X_pos)], 
                                     mu3[which.max(X_pos)], 
                                     mu4[which.max(X_pos)], 
                                     mu5[which.max(X_pos)], 
                                     mu6[which.max(X_pos)]),
                                   3)),6,2,byrow = TRUE)
colnames(table1)=c("Order of Global Polynomial","Normalized Weight of X =0.88")
table1


# Argument 2: Sensitivity of Estimates ------------------------------------

# SIMULATION STUDY

# simulate a data using rdd_dgf, then order the data to plot and see the
# implied discontinuity in the data
data1 <- rdd_dgf(scale1 = 2, scale2 = 5, sd_e = 0.1295, poly_x = 6,
                 n_sample = 100000)
data <- data1[order(data1$X1),]
plot(data[,3], data[,1])

# create buckets to store values of estimates and standard errors of the 
# treatment effect of the simulation of global and local polynomials 
est <- matrix(NA,15,6)
est_local <- matrix(NA,15,2)
se <- matrix(NA,15,6)
se_local <- matrix(NA,15,2)
c = 1

# this simulation can take approximately 15 mins (or more depending on the 
# system specifications)

# simulation for increasing standard deviation in the error term
for (i in seq(0.1, 3, 0.2)) {
  
  # data stores the rejection rate for local and global polynomials 
  data <- est_sd_sim(scale_x1 = 2, scale_x2 = 5, sd_e = i, rep = 2000, 
                     n_sample = 5000, x_poly = 5)
  
  # storing the results with increaing standard deviations in error term
  est[c,] <- data$ss[,2]
  est_local[c,] <- data$ss[1:2,4]
  se[c,] <- data$ss[,3]
  se_local[c,] <- data$ss[1:2,5]
  c = c + 1
}

# matrix to store the average results of the above simulation
results <- matrix(NA, 8,2)
colnames(results) <- c("estimates", "std_err")
rownames(results) <- c("local 1", "local 2", "global 1", "global 2", 
                       "global 3", "global 4", "global 5", "global 6")
results[1:2,1] <- round(colMeans(est_local),2)
results[3:8,1] <- round(colMeans(est),2)
results[1:2,2] <- round(colMeans(se_local),2)
results[3:8,2] <- round(colMeans(se),2)

# print the results
print(results)


# Argument 3: Confidence Interval -----------------------------------------

# this simulation can take 10-15 mins (or more depending on the system 
# specifications)

# the following code returns the median standard errors, rejection rate for 
# and global and local polynomial and estimates for each iteration for global
# and local polynomial.
rej_rate <- rej_rate_sim(n_threshold = 20000, n_sample = 5000, poly = 6, 
                         scale_x1 = 2, scale_x2 =  5, sd_e = 0.1295)

# store the rejection rate
rejrate_matrix <- matrix(c(round(rej_rate$rejectrate,3), 
                           round(rej_rate$median_se,3)), nrow = 8, ncol = 2)
colnames(rejrate_matrix) <- c("Rejection Rate", "Median S.E.")
rownames(rejrate_matrix) <- c("global 1", "global 2", "global 3", 
                                "global 4", "global 5", "global 6", 
                                "local 1", "local 2")

# print the rejection rates for global and local polynomial
print(rejrate_matrix)

# extract the estimates of treatment effect for every iteration
data <- rej_rate$estimate_sim
data_tf <- data.frame(data)
row_index <- rep(1:8)
data_tf <- cbind(row_index, data_tf)

# creating x-labs as the rejection rate for the corresponding global and 
# local polynomials
data_tf$row_index <- factor(data_tf$row_index, 
                            levels = c(1, 2, 3, 4, 5, 6,7, 8), 
                            labels = c("0.976", "1.000", "0.998", "0.163", 
                                       "0.096", "0.051", "0.327", "0.053"))

# subset to plot all but first column which is just row_index
data_subset <- data_tf[,-1]

# do the plot of estimates of treatment effect
ggplot(melt(data_tf, id.vars="row_index"), 
       aes(x = row_index, 
           y = value, 
           group = row_index, 
           fill = factor(row_index))) +
  geom_boxplot() + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("Estimates of Treatment Effect for Global and Local Polynomials") +
  xlab("Rejection Rate") +
  ylab("Estimates of Treatment Effect") +
  scale_y_continuous(limits = c(-0.1, 0.4)) +
  scale_fill_manual(values =c("#ffa600","#ff7c43","#f95d6a","#a05195",
                              "#d45087","#665191","#81BECE","#378BA4"),
                    labels = c("Global1","Global2", "Global3","Global4",
                               "Global5","Global6","Local1", "Local2"), 
                    guide = guide_legend(title = "Polynomials"))


# Empirical Application ---------------------------------------------------

set.seed(101)

# read and extract the excel data, don't forget to change the PATH
PATH = "../RDD_Guide_Dataset_0.xls"
rdd_data <- data.frame(read_xls(PATH, sheet = "Data"))

# normalizing the forcing variable
rdd_data$pretest <- rdd_data$pretest - 215


# Argument 1 --------------------------------------------------------------

# observations between [c,h]
X_pos_emp = sort(subset(rdd_data$pretest,rdd_data$pretest>=0))

# weights for global polynomial with order K and store in a dataframe
mu1_emp = weights(X_pos_emp,k=1) 
mu2_emp = weights(X_pos_emp,k=2)
mu3_emp = weights(X_pos_emp,k=3)
mu4_emp = weights(X_pos_emp,k=4)
mu5_emp = weights(X_pos_emp,k=5)
weight.frame.emp = data.frame(mu1_emp,mu2_emp,mu3_emp,mu4_emp,mu5_emp,X_pos_emp)

# plot the weights
arg1_emp = ggplot(weight.frame.emp,aes(x=X_pos_emp))+
  geom_line(aes(y=mu1_emp,colour="1"))+
  geom_line(aes(y=mu2_emp,colour="2"))+
  geom_line(aes(y=mu3_emp,colour="3"))+
  geom_line(aes(y=mu4_emp,colour="4"))+
  geom_line(aes(y=mu5_emp,colour = "5"))+
  xlab("Forcing Variable")+ylab("Weights")+
  scale_color_manual(name = "Degree of Polynomial",values =c("1"= "green",
                     "2"= "blue","3"="red","4"= "magenta","5" ="yellow"))+
  theme_classic()+
  theme(legend.background=element_rect(colour="black"),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

# save the graph in current working directory
tiff("arg1_emp", units="in", width=6, height=5, res=500)
arg1_emp
dev.off()

# local linear rdd regression
# bandwidth value for local linear and local quadratic
bandwidth_emp = matrix(NA,1,2)

# using triangular kernel
bandwidth_emp[,1] =  IKbandwidth(rdd_data$pretest, rdd_data$posttest, 
                                 cutpoint = 0, kernel="triangular")

# using rectangular kernel
bandwidth_emp[,2] =  IKbandwidth(rdd_data$pretest, rdd_data$posttest, 
                                 cutpoint = 0, kernel="rectangular")

# observations between [c,h] for triangular kernel
X_pos1_emp = sort(subset(rdd_data$pretest, 0<=rdd_data$pretest & 
                           rdd_data$pretest<=bandwidth_emp[,1]))

# observations between [c,h] for rectangular kernel
X_pos2_emp = sort(subset(rdd_data$pretest, 0<=rdd_data$pretest & 
                           rdd_data$pretest<=bandwidth_emp[,2]))

# weights for K = 1, 2 for triangular kernel and K = 1 for rectangular
# kernel respectively and store them in a dataframe
mu7_emp = matrix(c(weights(X_pos1_emp,k=1),
                   rep(0,length(X_pos_emp) - length(X_pos1_emp))))
mu8_emp = matrix(c(weights(X_pos1_emp,k=2),
                   rep(0,length(X_pos_emp) - length(X_pos1_emp))))
mu9_emp = matrix(c(weights(X_pos2_emp,k=1),
                   rep(0,length(X_pos_emp) - length(X_pos2_emp))))
weight.local.frame.emp = data.frame(mu7_emp,mu8_emp,mu9_emp,X_pos_emp)

# generate plot for weights
arg1_emp_loc = ggplot(weight.local.frame.emp,aes(x=X_pos_emp))+
  geom_line(aes(y=mu7_emp,colour="Linear and Triangular"))+
  geom_line(aes(y=mu8_emp,colour="Quadratic and Triangular"))+
  geom_line(aes(y=mu9_emp,color="Linear and Rectangular"))+
  geom_hline(yintercept = 0)+ xlab("Forcing Variable")+ylab("Weights")+
  theme_classic()+
  scale_color_manual(name = "Degree of Polynomial and Kernel Type",
                     values =c("Linear and Triangular"= "green",
                               "Quadratic and Triangular"= "blue",
                               "Linear and Rectangular"="red"))+
  theme_classic()+
  theme(legend.background=element_rect(colour="black"),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

# save the graph in current working directory
tiff("arg1_emp_loc", units="in", width=8, height=5, res=500)
arg1_emp_loc
dev.off()

# make a weight table for largest value of X_pos for global high-order 
# polynomial case
round(X_pos_emp[which.max(X_pos_emp)],3) #Largest value

# resulting table
table2 = matrix(rbind(c(1:5),round(c(mu1_emp[which.max(X_pos_emp)],
                                     mu2_emp[which.max(X_pos_emp)],
                                     mu3_emp[which.max(X_pos_emp)],
                                     mu4_emp[which.max(X_pos_emp)],
                                     mu5_emp[which.max(X_pos_emp)]),
                                   3)),5,2,byrow = TRUE)
colnames(table2) = c("Order of Global Polynomial","Normalized Weight of X =52")
table2


# Argument 2 --------------------------------------------------------------

# buckets to store the empirical results
results_emp <- matrix(NA, k+2, 2)
colnames(results_emp) <- c("estimates", "std err")

# for up to six order of global polynomial
for (i in 1:k) {
  # define linear regression model taking rdd in account to calculate 
  # estimates and standard errors of the treatment effect for global
  # polynomial
  lm_model <- lm(rdd_data$posttest ~ rdd_data$treat*poly
                 (rdd_data$pretest,i,raw = TRUE), data = rdd_data)
  results_emp[i+2,1] <- round(lm_model$coefficients[2],2)
  results_emp[i+2,2] <- round(sqrt(diag(vcov(lm_model)))[2],2)
}

# calculate bandwidth using Imbens and Kalyanaraman Bandwidth Selector
bandwidth <- IKbandwidth(rdd_data$pretest, rdd_data$posttest, cutpoint = 0,
                         kernel = "triangular")
# # trim the sampled data based on bandwidth for local linear and 
# quadratic regressions
local_data <-  filter(rdd_data,0-bandwidth<rdd_data$pretest & 
                        rdd_data$pretest<0+bandwidth)

# for up to two order of local polynomial
for (i in 1:2) {
  # define linear regression model taking rdd in account to calculate 
  # estimates and standard errors of the treatment effect for local
  # polynomial
  lm_model <- lm(local_data$posttest ~ local_data$treat * 
                   poly(local_data$pretest, i,raw=TRUE), data = local_data)
  results_emp[i,1] <- round(lm_model$coefficients[2],2)
  results_emp[i,2] <- round(sqrt(diag(vcov(lm_model)))[2],2)
}

# print the results
View(results_emp)


# Argument 3 --------------------------------------------------------------

# read and extract the excel data, don't forget to change the PATH
PATH_emp = "../Life Expectancy Data.csv"
df_poly <- data.frame(read.csv(PATH_emp))

# store dependent and forcing variable and create a dataframe
y_var <- df_poly$Life.expectancy
x_var <- df_poly$Income.composition.of.resources
new_data <- data.frame(y_var,x_var)

# removing rows for which the values of forcing variable is zero
clean_data <- new_data[!new_data$x_var == 0, ]
clean_data$x_var[is.na(clean_data$x_var)] <- mean(clean_data$x_var,na.rm = TRUE)

# plot the variables to see the relation between them
plot(clean_data$x_var, clean_data$y_var)

# this simulation can take approximately 10 mins (or more depending on the 
# system specifications)

# using function 5 get the results
values <- rej_rate_emp(20000,200)

# storing the rejection rates and standard errors for global and local
# polynomials
values$rejectrate
round(values$median_se,3)

# end ---------------------------------------------------------------------
