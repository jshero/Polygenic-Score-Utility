
##Converting percentage of variance explained by PGS into correlations by raking the square-root of the R2 value
sqrt(0.02)#0.141
sqrt(0.1)#0.316
sqrt(0.20)#0.447
sqrt(0.50)#0.707

#Loading required packages
library(MASS)
library(mvtnorm)

########Research Question 1######
#Creating a function to make matrices for use in subsequent analysis. r is correlation of PGS to end of year achievement
#s is correlation between progress monitoring and end of year achievement.
#since no distinction can be made between the correlations between PGS and end of year scores and between PGS and progress monitoring scores, r is used for both
matt <- function(r,s) {

  
  u=matrix(c(1,r,s,
             r,1,r,
             s,r,1),3)
  return(u)
}

#placing these values here again for reference
sqrt(0.02)#0.141
sqrt(0.1)#0.316
sqrt(0.20)#0.447
sqrt(0.50)#0.707

#creating the 12 possible matrices for the combinations of PGS and progress monitoring effectiveness levels
l1=matt(0.141,0.40)
m1=matt(0.141,0.55)
h1=matt(0.141,0.70)

l2=matt(0.316,0.40)
m2=matt(0.316,0.55)
h2=matt(0.316,0.70)

l3=matt(0.447,0.40)
m3=matt(0.447,0.55)
h3=matt(0.447,0.70)

l4=matt(0.707,0.40)
m4=matt(0.707,0.55)
h4=matt(0.707,0.70)

#Creating a vector of means. Setting all to zero standardizes and allows us to use correlations as we would covariances for data simulation.
mu=c(0,0,0)

#2% models
#Simulating dataset of n=1000000, with the means set to 0, and priorly specified l1 correlation matrix
l1d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = l1, empirical = T))
#running OLS regression of End of year scores predicted by PGS and progress monitoring 
l11=lm(V1~V2+V3,l1d)
summary(l11)
##running OLS regression of End of year scores predicted by just progress monitoring. Using the difference to find unique contribution of PGS to variance explained 
l12=lm(V1~V3,l1d)
summary(l12)

#The prior steps were used for all following analyses but with the correct correlation matrices specified
m1d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = m1, empirical = T))

m11=lm(V1~V2+V3,m1d)
summary(m11)

m12=lm(V1~V3,m1d)
summary(m12)

h1d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = h1, empirical = T))

h11=lm(V1~V2+V3,h1d)
summary(h11)

h12=lm(V1~V3,h1d)
summary(h12)

#10% models
l2d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = l2, empirical = T))

l21=lm(V1~V2+V3,l2d)
summary(l21)

l22=lm(V1~V3,l2d)
summary(l22)

m2d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = m2, empirical = T))

m21=lm(V1~V2+V3,m2d)
summary(m21)


m22=lm(V1~V3,m2d)
summary(m22)

h2d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = h2, empirical = T))

h21=lm(V1~V2+V3,h2d)
summary(h21)

h22=lm(V1~V3,h2d)
summary(h22)

#20% models
l3d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = l3, empirical = T))

l31=lm(V1~V2+V3,l3d)
summary(l31)

l32=lm(V1~V3,l3d)
summary(l32)

m3d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = m3, empirical = T))

m31=lm(V1~V2+V3,m3d)
summary(m31)

m32=lm(V1~V3,m3d)
summary(m32)

h3d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = h3, empirical = T))

h31=lm(V1~V2+V3,h3d)
summary(h31)

h32=lm(V1~V3,h3d)
summary(h32)

#50% models
l4d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = l4, empirical = T))

l41=lm(V1~V2+V3,l4d)
summary(l41)

l42=lm(V1~V3,l4d)
summary(l42)

m4d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = m4, empirical = T))

m41=lm(V1~V2+V3,m4d)
summary(m41)

m42=lm(V1~V3,m4d)
summary(m42)

h4d <- as.data.frame(mvrnorm(n = 1000000, mu = mu, Sigma = h4, empirical = T))

h41=lm(V1~V2+V3,h4d)
summary(h41)

h42=lm(V1~V3,h4d)
summary(h42)

########Research Question 2######
###Making a function that produces all of the information used in the paper
#e and f are desired correlations. Again since no distinction can be made between the correlations between PGS and end of year scores and between PGS and progress monitoring scores, e is used for both
#f then represents the correlation between progress monitoring scores and end of year achievement scores
#could be expanded to allow for a third correlation between variables by simply inputting another letter in the command, and changing out the two desired 'e's in the matrix with said letter
#c represents the cutoff for the outcome variable as a percentile
#d represents the cutoff for the predictor variables as a percentile
#if desired, could have 3 distinct cutoff variables, would simply call for adding another letter in the command line and creating another 'qnorm' line of code
per <- function(e,f,c,d) {
  mu1=c(0,0,0)####it's standardized so means all set to 0
  sigma <- matrix(c(1,e,f,
                    e,1,e,
                    f,e,1),
                  ncol = 3)###this is changing the input correlations (e and f) into a sigma matrix
  a=qnorm(c)###this is converting percentile (cutpoint) into a z-score
  b=qnorm(d)###Same as above
  ppp=pmvnorm(upper = c(a, b,b),#upper bounds set to the desired cutoff or +10 SD's when called for 
              lower=c(-10,-10,-10),#lower bounds set to the desire cutoff or -10S's when called for
              mean = mu1,
              sigma = sigma,
              algorithm = Miwa())#####Has the outcome (learning disability risk), both correctly predicted
  ppn=pmvnorm(upper = c(a, b,10),
              lower=c(-10,-10,b),
              mean = mu1,
              sigma = sigma,
              algorithm = Miwa())#####Has the outcome (learning disability risk), only the first variable correctly predicted
  pnp=pmvnorm(upper = c(a, 10,b),
              lower=c(-10,b,-10),
              mean = mu1,
              sigma = sigma,
              algorithm = Miwa())#####Has the outcome (learning disability risk), only the second variable correctly predicted
  pnn=pmvnorm(upper = c(a, 10,10),
              lower=c(-10,b,b),
              mean = mu1,
              sigma = sigma,
              algorithm = Miwa())#####Has the outcome (learning disability risk), neither variable correctly predicted it
  npp=pmvnorm(upper = c(10, b,b ),
              lower=c(a,-10,-10),
              mean = mu1,
              sigma = sigma,
              algorithm = Miwa())#####Does NOT have the outcome (learning disability risk), both incorrectly predicted that it would
  npn=pmvnorm(upper = c(10, b,10 ),
              lower=c(a,-10,b),
              mean = mu1,
              sigma = sigma,
              algorithm = Miwa())#####Does NOT have the outcome (learning disability risk), only the second variable correctly predicted
  nnp=pmvnorm(upper = c(10, 10,b ),
              lower=c(a,b,-10),
              mean = mu1,
              sigma = sigma,
              algorithm = Miwa())#####Does NOT have the outcome (learning disability risk), only the first variable correctly predicted
  nnn=pmvnorm(upper = c(10, 10,10 ),
              lower=c(a,b,b),
              mean = mu1,
              sigma = sigma,
              algorithm = Miwa())#####Does NOT have the outcome (learning disability risk), both variables correctly predicted
  ppv1=(ppp+pnp)/(ppp+pnp+npp+nnp)#Positive predictive value of variable 2 alone
  ppv2=(ppp+pnp+ppn)/(ppp+pnp+ppn+npp+npn+nnp)#positive predictive value when meeting either cutoff
  ppv3=(ppp)/(ppp+npp)#positive predictive value when meeting both cutoffs
  npv1=(nnn+npn)/(nnn+npn+pnn+ppn)#Negative predictive value of variable 2 alone
  npv2=(nnn)/(nnn+pnn)#Negative predictive value when meeting either cutoff
  npv3=(nnn)/(nnn+ppn+pnp+pnn)#Negative predictive value when meeting both cutoffs
  
  u=matrix(c('Shared True Positive',ppp,'V1 TP/ V2 FN',ppn,'V1 FN/ V2 TP',pnp,'Shared False Negative',pnn,'Shared False Positive',npp,
             'V1 FP/ V2 TN',npn,'V1 TN/ V2 FP',nnp,'Shared True Negative',nnn,'Positive Predictive Value of V2 Alone',ppv1,
             'Positive Predictive Value meeting Either Cutoff', ppv2,'Positive Predictive Value meeting both cutoffs',ppv3,
             'Negative Predictive Value of V2 Alone',npv1,'Negative Predictive Value meeting either Cutoff',npv2,'Negative Predictive Value meeting both cutoffs',npv3),2)
  return(u)
}

#All 8 different results (Shared true positive, only V1 true positive, etc...) were computed above by distinguishing... 
#which part of the probability density function we desired using different upper and lower bounds for each variable. 
#For example, if we wanted only PGS to accurately predict learning disability risk then we would look for... 
#all cases that had PGS between the bounds of -10SD's and the cutoff, end of year achievement scores between -10SD's and its cutoff, 
#and progress monitoring between the cutoff and +10SD's (i.e higher than the cutoff). 
#Here, progress monitoring scores were higher than the cutoff and did not accurately predict the learning disability risk that was present, but PGS was below the cutoff and accurately predicted this.


#Using this function, all 12 combinations of PGS and progress monitoring effectiveness were tested at the 4 different cutoffs