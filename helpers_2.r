
# June 13, 2013, same as helpers but uses cov in computation of the variance.


# These packages should be downloaded from R-CRAN 
library('MASS')
library('CompQuadForm')


# Regular analysis with rare variants
# Get p-values from CAST and C-alpha
#
##########################################
# This function should be in package #####
##########################################
#
# Input:
# Input values matrix of genotype calls
# Y - phenotype value, 1-cases and 0-controls
# X - matrix of genotypes 
# Matrix X - has next dimensions n by J, J - number of variants, n - sample size
# perm - number of permutations
# Resampling: Permutation
#
#
# Output:
# two p-values for CAST and C-alpha

get_pval = function(Y,X,perm){

S = calcS(X,Y)
Sigma = calcVS(X,Y)
w = rep(1,length(X[1,]))
SLobs = testL(S,Sigma,w)
SQobs = testC(S,Sigma)
Q = NULL
L = NULL

for (i in 1:perm){
k = sample(length(Y))
Y = Y[k]
S = calcS(X,Y)
Sigma = calcVS(X,Y)
w = rep(1,length(X[1,]))
L = c(L,testL(S,Sigma,w))
Q = c(Q,testC(S,Sigma))
}
pl = sum(L<=SLobs)/perm
pQ = sum(Q<=SQobs)/perm
return (c(pl,pQ))
}

# RVS analysis with rare variants
# Get p-values from RVS with CAST and C-alpha
#
##########################################
# This function should be in package #####
##########################################
#
# Input:
# Input values matrix of genotype calls
# Y - phenotype value, 1-cases and 0-controls
# X - matrix of conditional expectations 
# Matrix X - has next dimensions n by J, J - number of variants, n - sample size
# perm - number of permutations
# Resampling: bootstrap
#
#
# Output:
# two p-values for CAST and C-alpha

get_pval_b = function(Y,X,perm){
X1 = as.matrix(X[Y==1,])
X2 = as.matrix(X[Y==0,])
S = calcS(X,Y)
Sigma = calcVS_two(X1,X2,Y)
w = rep(1,length(X[1,]))
SLobs = testL(S,Sigma,w)
SQobs = testC(S,Sigma)
Q = NULL
L = NULL

for (i in 1:perm){
Xs = get_b_sample(X1,X2)
Xa = as.matrix(Xs[Y==1,])
Xb = as.matrix(Xs[Y==0,])
X = rbind(Xa,Xb)
S = calcS(X,Y)
Sigma = calcVS_two(Xa,Xb,Y)
w = rep(1,length(X[1,]))
L = c(L,testL(S,Sigma,w))
Q = c(Q,testC(S,Sigma))
}
pl = sum(L<=SLobs)/perm
pQ = sum(Q<=SQobs)/perm
return (c(pl,pQ))
}



# RVS analysis with rare variants
# Robust Variance Estimate
# Get p-values from RVS with CAST and C-alpha
#
##########################################
# This function should be in package #####
##########################################
#
#
# Input values matrix of expected values of genotypes given sequence data 
# Input Y - phenotype value, 1-cases and 0-controls
# X - matrix of conditional expected values
# Matrix X - has next dimensions n by J,- number of variants, n - sample size
# P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
# perm - number of permutations
# Resampling: Bootstrap
# Variance Estimate: Robust
# Output:
# two p-values for CAST and C-alpha

get_pval_b1 = function(Y,X,P,perm){
X1 = as.matrix(X[Y==1,])
X2 = as.matrix(X[Y==0,])
S = calcS(X,Y)
Sigma = calcVS_twoP(X1,X2,P,Y)
w = rep(1,length(X[1,]))
SLobs = testL(S,Sigma,w)
SQobs = testC(S,Sigma)
Q = NULL
L = NULL

for (i in 1:perm){
Xs = get_b_sample(X1,X2)
Xa = as.matrix(Xs[Y==1,])
Xb = as.matrix(Xs[Y==0,])
X = rbind(Xa,Xb)
S = calcS(X,Y)
Sigma = calcVS_two(Xa,Xb,Y)
w = rep(1,length(X[1,]))
L = c(L,testL(S,Sigma,w))
Q = c(Q,testC(S,Sigma))
}
pl = sum(L<=SLobs)/perm
pQ = sum(Q<=SQobs)/perm
return (c(pl,pQ))
}




# Helper functions

# Get score vector  S
#  

calcS = function(X,Y){ # NEED to check!!!
L = length(X[1,])
S = NULL
for (i in 1:L){
Yn = Y[!is.na(X[,i])]
Xn = X[!is.na(X[,i]),i]
xca = Xn[Yn==1]
xco = Xn[Yn==0]
my = mean(Yn)
s = sum(xca,na.rm=T)*(1-my) - my*sum(xco,na.rm=T)
S = c(S,s)
}
S = t(S)
return (S)
}
##
# Get Variance of vector S
# Regular variance 
calcVS = function(X,Y){

X1 = X[Y==1,]
X2 = X[Y==0,]
l1 = length((X1[,1]))
l2 = length(X2[,1])
J = length(X1[1,])

a =colSums(is.na(X1),na.rm=T)
b =colSums(is.na(X2),na.rm=T)

ncase = rep(l1,J) - a
ncont = rep(l2,J) - b
nn = ncase+ncont
L= length(Y)
Xm = cov(X,use="pairwise.complete.obs")
vs  = sqrt(diag(ncase*ncont/nn^2))
diag_S  = vs%*%Xm%*%vs*L
return (diag_S)
}

#centralize matrix X
mmatrix = function(X){
l = length(X[1,])
mX = NULL
for (i in 1:l){
mm = mean(X[,i],na.rm=T) 
mv = X[,i]-mm
mX = cbind(mX,mv)
}
return (mX)
}


# Linear test
testL = function(S,Sigma,w){
T = (w)%*%t(S)
Sg = t(w)%*%Sigma%*%w
linear = as.numeric(T)/sqrt(as.numeric(Sg))
p = 2*(1-pnorm(abs(linear)))
return (p)
}
#
#
# C-Alpha
testC = function(S,Sigma){
quad  =  (S)%*%t(S)
lambda = Re(eigen(Sigma)$values)
qp = davies(quad,lambda)
return (abs(qp$Qq))
}
#
#
# Hotelling
testH = function(S,Sigma){
dd  =  try(ginv(Sigma), silent = TRUE)
if ( class(dd) == 'try-error'){
cat('Inverse_error','\n')
Sigma = diag(diag(Sigma))
}
X =ginv(Sigma)
quad = S%*%X%*%t(S)
rank = qr(X)$rank
p = 1-pchisq(quad,rank)
return (p)
}


##
# Get Variance of vector S
# Robust variance estimate  
calcVS_two = function(X1,X2,Y){
l1 = length((X1[,1]))
l2 = length(X2[,1])
J = length(X1[1,])

a =colSums(is.na(X1),na.rm=T)
b =colSums(is.na(X2),na.rm=T)

ncase = rep(l1,J) - a
ncont = rep(l2,J) - b
nn = ncase+ncont

Yhat = mean(Y)
L= length(Y)
p = l2/l1
q = l1/l2
Xm1 = cov(X1,use="pairwise.complete.obs")*(l1-1)
Xm2 = cov(X2,use="pairwise.complete.obs")*(l2-1)

vs  = sqrt(diag(ncase*ncont/nn^2))
diag_S  = vs%*%(p*Xm1+ q*Xm2)%*%vs
return (diag_S)
}

#
# Bootstrap resampling
#
get_b_sample = function(X1,X2){
case = length(X1[,1])
cont = length(X2[,1])
X1 = mmatrix(X1)
X2 =  mmatrix(X2)
ca = sample(1:case,case,replace=TRUE)
co = sample(1:cont,cont,replace=TRUE)
Xca = as.matrix(X1[ca,])
Xco = as.matrix(X2[co,])
X = rbind(Xca,Xco)
return (X)
}

##
# Get Variance of vector S
# Robust variance estimate 
# P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
calcVS_twoP = function(X1,X2,P,Y){
L = length(P[,1])
V = NULL
for (i in 1:L){
V = c(V,vp(P[i,]))
}
V = diag(sqrt(V))
l1 = length(X1[,1])
l2 = length(X2[,1])
Yhat = mean(Y)
L= length(Y)
p = l2/l1
q = l1/l2
Sgcase = t(V)%*%cor(X1)%*%V*l1
Xm2 = mmatrix(X2)
vs  = var(Y)
diag_S  = vs*(p*Sgcase + q*t(Xm2)%*%Xm2)
return (diag_S)
}



