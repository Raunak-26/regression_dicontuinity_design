#install.packages("dplyr")
library(dplyr)
library(boot)

set.seed(100)
X = rnorm(100000,0,5)
D = ifelse(X>=0,1,0)
Y = 2 + 2*D + 3*I(X-2) + rnorm(10000,0,1)
X_normalised = X-0
data = data.frame(Y,D,X_normalised)
maxdegree = 10
poly.mse =c()
for (d in 1:maxdegree){
  poly.fit = glm(Y~(D*poly(X_normalised,d,raw = TRUE)),data=data)
  mse = cv.glm(data,poly.fit,K=5)$delta[1]
  poly.mse = c(poly.mse,mse)
}
plot(poly.mse,type = 'b')
x=which.min(poly.mse)
x
points(x,poly.mse[x],pch=20,cex=2,col='red')

lm1 = lm(Y~D*poly(X_normalised,6))
summary(lm1)
plot(X,lm1$fitted.values)
lm1$coefficients[2]
lm3 = lm(Y~D*X_normalised+D*I(X_normalised^2)+D*I(X_normalised^3))
summary(lm3)
l = data.frame()
lm2$coefficients[2]
lm_result = matrix(NA,6,1)
for (i in 1:6){
  l=lm(Y~(D*I(X_normalised^i)))
  lm_result[i,]=l$coefficients[2]
}
summary(l)
list(lm1$coefficients[2],lm2$coefficients[2])
lm2$effects
D*X_normalised
mat = matrix(NA)
for (i in 1:6){
  matu = matrix(D*I(X_normalised^i))
  mat[,i]=matu
}

# Weights -----------------------------------------------------------------

#Treatment
k = 2
X_pos = sort(subset(X,X>0))
result = function(X_pos,k){

  
  
  SUMVec = matrix(rep(1,length(X_pos)),length(X_pos),1)
  e = matrix(c(1,rep(0,k)),k+1,length(X_pos))
  #t(e)

  lm=matrix(NA,k,length(X_pos))
  for (i in 1:k){
    lm[i,]=X_pos^i
  }
  Xvec = rbind(rep(1,length(X_pos)),lm)



  #data1 = matrix(NA,k+1,k+1)
   #for (i in 1:(k+1)){
    #data1[i,]=X_pos[1]^(i-1)
    #for(j in 1:(k)){
     #   data1[,j+1]= X_pos[1]*data1[,j]
      #}
  #}
  #View(data1)

  #data2 = matrix(NA,k+1,k+1)
  #for (i in 1:(k+1)){
  #  data2[i,]=X_pos[2]^(i-1)
  #  for(j in 1:(k)){
  #    data2[,j+1]= X_pos[2]*data1[,j]
  # }
  #}


  array1 = array(NA,dim=c(k+1,k+1,length(X_pos)))
  for(m in 1:length(X_pos)){
    for (i in 1:(k+1)){
      array1[,,m][i,]=X_pos[m]^(i-1)
      for(j in 1:(k)){
        array1[,,m][,j+1]= X_pos[m]*array1[,,m][,j]
      }
    }
  }
  
  #data1[6]
  #lol = array(data=c(data1,data2),dim=c(k+1,k+1,2))
  #lol[4,4,1]
  #array1[,,516]
  #array1[,,1]
  #lol[,,1]

  sum_all = apply(array1,1:2,sum)
  #sum_all
  solve(sum_all)
  weight = t(SUMVec)%*%t(e)%*%solve(sum_all)%*%Xvec
  return(weight)
}
mu1 = result(X_pos,k=1)
mu2 = result(X_pos,k=2)
mu3 = result(X_pos,k=3)
mu4 = result(X_pos,k=4)
mu5 = result(X_pos,k=5)
mu6 = result(X_pos,k=6)
mean(mu6)      
View(mu4)
data1 = data.frame(mu1,mu2,mu3,mu4,mu5,mu6,X_pos)
library(ggplot2)
ggplot(data1,aes(x=X_pos))+geom_line(aes(y=mu1))+
  geom_line(aes(y=mu2),color="blue")+geom_line(aes(y=mu3),color="red")+
  geom_line(aes(y=mu4),col="yellow")+geom_line(aes(y=mu5),col="magenta")+
  geom_line(aes(y=mu6),col="pink")+xlab("Forcing Variable")+ylab("Weights")+
  theme_classic()




plot(X_pos,mu1,type="l",ylim = c(-110,110),ylab = "Weights",xlab = "Forcing Variable")
lines(X_pos,mu2,type = "l",col="green")
lines(X_pos,mu3,col="blue")
lines(X_pos,mu4,col="red")
lines(X_pos,mu5,col="magenta")
lines(X_pos,mu6,col="purple")

hist(X_pos)
which.max(X_pos)
X_pos[19326]
