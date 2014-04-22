# This implementation of Lisa's project.
# This code is designed to mimic short and high read depth data
# We will simulate hard calls from


options(warn=0)

# EM algorithm to estimate genotype frequencies.
source('em.r')


#
# Generate sequence data under null hypothesis
#
##########################################
# This function should be in package #####
##########################################

## Input: 
#     L - number of variants (integer) >0
#     ncase - number of cases (integer) >0
#	  ncont - number of controls (integer) >0
#     mdcase - average read depth in cases (double) >0
#	  sddcase - standard deviation of read depth in cases  (double) >0
#     mdcont - average read depth in cases (double) >0
#	  sddcont - standard deviation of read depth in cases (double) >0  
#     em - average error rate, probability that sequence call is wrong (double) >0
#	  esd - standart deviation of error rate (double) >0
#     mmaf - MAF for L variants, HWE is assumed >0

####################################
## Output: Object with three fields
####################################
# $MM - conditional expected value E(Gij|Dij)
# $P -estimated genotype frequencies P(G=0), P(G=1), P(G=2)
# $G - true genotype from which sequence data generated
####
####
####

gen_matrix  = function(L,ncase,ncont,mdcase,sddcase,mdcont,sddcont,em,esd,mmaf){
MAF = NULL
A1 = 'C'
A2 = 'T'
MM = NULL
G = NULL
P = NULL
for (i in 1:L){
if (i %% 100 == 1) {cat(i,'\n')}
v = NULL
erv = NULL
rdv = NULL
gv = NULL
maf = mmaf[i]
MAF = c(MAF,maf)
for (j in 1:ncase){
rd = round(sddcase*rnorm(1) + mdcase)
if (rd <= 0){rd=1}
error = esd*rnorm(rd) + em
k = rbinom(1,2,maf)
gv = c(gv,k)
gen_vect = c('CC','CT','TT')
if (k==2){genotype=c('C','C')}
if (k==1){genotype=c('C','T')}
if (k==0){genotype=c('T','T')}
a = seq_call(genotype, error, gen_vect)
v = c(v,a)
erv = c(erv,error)
rdv = c(rdv,rd)
}
for (j in 1:ncont){
rd = round(sddcont*rnorm(1) + mdcont)
if (rd <= 0){rd=1}
error = esd*rnorm(rd) + em
k = rbinom(1,2,maf)
gv = c(gv,k)
gen_vect = c('CC','CT','TT')
if (k==2){genotype=c('C','C')}
if (k==1){genotype=c('C','T')}
if (k==0){genotype=c('T','T')}
a = seq_call(genotype, error, gen_vect)
v = c(v,a)
erv = c(erv,error)
rdv = c(rdv,rd)
}

G = cbind(G,gv)
M  = get_Mr(erv,v,rdv)
Mm = get_M(erv,v)
p = EM_calc(M)
#p = c(0.9^2,2*0.9*0.1,0.1^2)
P = rbind(P,p)
if ((sum(p< (-0.00001))>0) || (sum(p>1.00001)<0)){
cat('Problem!!!')
return (NULL)
}
EG = get_EG(Mm,p,rdv)
MM = cbind(MM,EG)
}
return (list(MM=MM,P=P,G=G))
}

#
# Generate sequence data under alternative hypothesis
#
##########################################
# This function should be in package #####
##########################################

## Input: 
#     L - number of variants (integer) >0
#     ncase - number of cases (integer) >0
#	  ncont - number of controls (integer) >0
#     mdcase - average read depth in cases (double) >0
#	  sddcase - standard deviation of read depth in cases  (double) >0
#     mdcont - average read depth in cases (double) >0
#	  sddcont - standard deviation of read depth in cases (double) >0  
#     em - average error rate, probability that sequence call is wrong (double) >0
#	  esd - standart deviation of error rate (double) >0
#     mmafCa - MAF for L variants in cases, HWE is assumed >0
#     mmafC0 - MAF for L variants in controls, HWE is assumed >0

####################################
## Output: Object with three fields
####################################
# $MM - conditional expected value E(Gij|Dij)
# $P -estimated genotype frequencies P(G=0), P(G=1), P(G=2)
# $G - true genotype from which sequence data generated
####
####
####


gen_matrix_c  = function(L,ncase,ncont,mdcase,sddcase,mdcont,sddcont,em,esd,mmafCa,mmafCo,R){
MAF = NULL
A1 = 'C'
A2 = 'T'
MM = NULL
MM1 = NULL
G = NULL
P = NULL
for (i in 1:L){
if (i %% 100 == 1) {cat(i,'\n')}
v = NULL
erv = NULL
rdv = NULL
gv = NULL
maf = mmafCa[i]
MAF = c(MAF,maf)
for (j in 1:ncase){
rd = round(sddcase*rnorm(1) + mdcase)
if (rd <= 0){rd=1}
error = esd*rnorm(rd) + em
k = rbinom(1,2,maf)
gv = c(gv,k)
gen_vect = c('CC','CT','TT')
if (k==2){genotype=c('C','C')}
if (k==1){genotype=c('C','T')}
if (k==0){genotype=c('T','T')}
a = seq_call(genotype, error, gen_vect, gen_freq,R)
v = c(v,a)
erv = c(erv,error)
rdv = c(rdv,rd)
}
maf = mmafCo[i]
for (j in 1:ncont){
rd = round(sddcont*rnorm(1) + mdcont)
if (rd <= 0){rd=1}
error = esd*rnorm(rd) + em
k = rbinom(1,2,maf)
gv = c(gv,k)
gen_vect = c('CC','CT','TT')
if (k==2){genotype=c('C','C')}
if (k==1){genotype=c('C','T')}
if (k==0){genotype=c('T','T')}
a = seq_call(genotype, error, gen_vect, gen_freq,R)
v = c(v,a)
erv = c(erv,error)
rdv = c(rdv,rd)
}

G = cbind(G,gv)
M  = get_Mr(erv,v,rdv)
Mm = get_M(erv,v)
p = EM_calc(M)
#p = c(0.9^2,2*0.9*0.1,0.1^2)
P = rbind(P,p)
if ((sum(p< (-0.00001))>0) || (sum(p>1.00001)<0)){
cat('Problem!!!')
return (NULL)
}
EG = get_EG(Mm,p,rdv)
EG2 = get_EG2(Mm,p,rdv)
MM = cbind(MM,EG)
MM1 = cbind(MM1,EG2)
}
return (list(MM=MM,P=P,G=G,t=MM1))
}



# Regular Score/Trend test
# Evaluation: by asymptotic distribution
#
##########################################
# This function should be in package #####
##########################################
#
# Input:
#
# Input values matrices of genotypes for cases and controls
# For cases: matrix M1 has next dimensions ncase by J, J - number of variants and ncase - number of cases (double)
# For controls: matrix M2 has next dimensions ncont by J, J - number of variants and ncont - number of controls (double)
# M1 - matrix with dimensions: number of cases by J (double)
# M2 - matrix with dimensions: number of controls by J (double)
#  
# Output: 
#
### vector of p-values for J variants (double)

test_EG = function(M1,M2){

L=length(M1[1,])
cc = rep(0,L)
for (i in 1:L){
X = c(M1[,i],M2[,i])
p = length(M1[,i])/length(X)
q = 1 - p
Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
Y = Y[!is.na(X)]
X = X[!is.na(X)]
a = (glm(Y~1,family='binomial'))
b = glm(Y~X,family='binomial')
cc[i] = 1-pchisq(anova(a,b,test='Rao')$Rao[2],1)
}
return (cc)
}

# RVS Score/Trend test
# Evaluation: by asymptotic distribution
#
##########################################
# This function should be in package #####
##########################################
#
# Input values matrices of expected values of genotypes given sequence data for  cases and controls
# Uses the Robust Variance Estimates described in manuscript
# Evaluation: by asymptotic distribution
#
# Input:
#
# For cases: matrix M1 (double) has next dimensions ncase by J, J - number of variants and ncase - number of cases
# For controls: matrix M2 (double) has next dimensions ncont by J, J - number of variants and ncont - number of controls
# M1 - matrix with dimensions: number of cases by J (double)
# M2 - matrix with dimensions: number of controls by J (double)
# P - matrix of estimated genotype frequency for J variants P(G=0), P(G=1) and P(G=2)
# P - J by 3 matrix (double) 
#
# Output:
#
### vector of p-values for J variants
### Evaluation by asymptotic distribution
test_EGC = function(M1,M2,P){
L=length(M1[1,])
cc = rep(0,L)
for (i in 1:L){
X = c(M1[,i],M2[,i])
mx = mean(X)
p = length(M1[,i])/length(X)
q = 1 - p
Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
Y = Y[!is.na(X)]
X = X[!is.na(X)]
a = (glm(Y~1,family='binomial'))
b = glm(Y~X,family='binomial')
p = length(X[Y==1])/length(X)
q = 1 - p 
v = q*vp(P[i,]) + p*var(X[Y==0],na.rm=TRUE)
x = anova(a,b,test='Rao')$Rao[2]*var(X,na.rm=TRUE)/v
cc[i] = 1-pchisq(x,1)
}
return (cc)
}

# RVS Score/Trend test
#
##########################################
# This function should be in package #####
##########################################
#
# Input values matrix of expected values of genotypes given sequence data
# Uses the Regular Variance Estimates for cases and controls
# Evaluation: by asymptotic distribution
#
# Input:
#
# For cases: matrix M1 (double) has next dimensions ncase by J, J - number of variants and ncase - number of cases
# For controls: matrix M2 (double) has next dimensions ncont by J, J - number of variants and ncont - number of controls
# M1 - number of cases by J (double)
# M2 - number of controls by J (double)
# P - matrix of estimated genotype frequency for J variants P(G=0), P(G=1) and P(G=2)
# P - J by 3 matrix (double) 
#
# Output:
#
### vector of p-values for J variants
### 
test_EGCV = function(M1,M2,P){
L=length(M1[1,])
cc = rep(0,L)
for (i in 1:L){
X = c(M1[,i],M2[,i])
Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
Y = Y[!is.na(X)]
X = X[!is.na(X)]
a = (glm(Y~1,family='binomial'))
b = glm(Y~X,family='binomial')
p = length(X[Y==1])/length(X)
q = 1 - p 
v = q*var(X[Y==1]) + p*var(X[Y==0])
x = anova(a,b,test='Rao')$Rao[2]*var(X)/v
cc[i] = 1-pchisq(x,1)
}
return (cc)
}


#
##########################################
# This function should be in package #####
##########################################
#
# RVS Score/Trend test
# Input values matrix of expected values of genotypes given sequence data
# Uses the Robust Variance Estimates for cases and controls
# Evaluation: bootstrap permutation
#
# Input:
#
# For cases: matrix M1 (double) has next dimensions ncase by J, J - number of variants and ncase - number of cases
# For controls: matrix M2 (double) has next dimensions ncont by J, J - number of variants and ncont - number of controls
# M1 - number of cases by J (double)
# M2 - number of controls by J (double)
# P - matrix of estimated genotype frequency for J variants P(G=0), P(G=1) and P(G=2)
# P - J by 3 matrix (double) 
#
# Output:
#
### vector of p-values for J variants


test_EGB = function(M1,M2,P,perm){

L=length(M1[1,])
cc = rep(0,L)
for (i in 1:L){
if (i %% 100 == 1) {cat(i,'\n')}
X = c(M1[,i],M2[,i]) 
Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
Y = Y[!is.na(X)]
X = X[!is.na(X)]
case = sum(Y==1)
cont = sum(Y==0)
p = length(X[Y==1])/length(X)
q = 1 - p 
vcase = vp(P[i,])
#vcase = var(M1[,i])
vcont = var(X[Y==0])
Tobs = (mean(X[Y==1]) - mean(X[Y==0]))/sqrt(vcase/case+vcont/cont) 
X1 = X[Y==1] - mean(X[Y==1])
X2 = X[Y==0] - mean(X[Y==0])
C = NULL
n = perm
for (j in 1:n){
Xca = sample(X1[],case,replace=TRUE)
Xco = sample(X2[],cont,replace=TRUE)
vcase = var(Xca)
vcont = var(Xco)
C =c(C,(mean(Xca) - mean(Xco))/sqrt(vcase/case+vcont/cont))
}
cc[i] = sum(abs(C)>=abs(Tobs))/n
}
return (cc)
}

#
##########################################
# This function should be in package #####
##########################################
#
# Regular Score/Trend test
# Evaluation: by permutation distribution
#
# Input:
#
# For cases: matrix M1 (double) has next dimensions ncase by J, J - number of variants and ncase - number of cases
# For controls: matrix M2 (double) has next dimensions ncont by J, J - number of variants and ncont - number of controls
# M1 - number of cases by J (double)
# M2 - number of controls by J (double)
# P - matrix of estimated genotype frequency for J variants P(G=0), P(G=1) and P(G=2)
# P - J by 3 matrix (double) 
# perm - number of permutations (integer)
#
# Output:
#
### vector of p-values for J variants

test_EGPv = function(M1,M2,P,perm){
case = length(M1[,1])
cont = length(M2[,1])
L=length(M1[1,])
cc = rep(0,L)
for (i in 1:L){
if (i %% 100 == 1) {cat(i,'\n')}
X = c(M1[,i],M2[,i]) 
Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
Y = Y[!is.na(X)]
X = X[!is.na(X)]
case = sum(Y==1)
cont = sum(Y==0)
test = calct(X[1:case],X[(case+1):(case+cont)])
C = NULL
n = perm
for (j in 1:n){
k = sample(length(X))
X = X[k]
a = calct(X[1:case],X[(case+1):(case+cont)])
C = c(C,a)
}
cc[i] = sum(C>=test)/n
}
return (cc)
}



#
# RVS Score/Trend test
# Evaluation: bootstrap distribution
#
##########################################
# This function should be in package #####
##########################################
#
# Input:
#
# Input values matrix of expected values of genotypes given sequence data
# Uses the Regular Variance Estimates for cases and controls
# For cases: matrix M1 has next dimensions ncase by J, J - number of variants and ncase - number of cases
# For controls: matrix M2 has next dimensions ncont by J, J - number of variants and ncont - number of controls
# M1 - matrix with dimensions: number of cases by J (double)
# M2 - matrix with dimensions: number of controls by J (double)
# P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2)
# P - matrix of J by 3 (double)
# perm - number of bootstrap permutations (integer)
#
# Output: 
#
### vector of p-values for J variants


test_EGBnew = function(M1,M2,P,perm){
case = length(M1[,1])
cont = length(M2[,1])
L=length(M1[1,])
cc = rep(0,L)
for (i in 1:L){
if (i %% 100 == 1) {cat(i,'\n')}
X = c(M1[,i],M2[,i]) 
Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
Y = Y[!is.na(X)]
X = X[!is.na(X)]
p = length(X[Y==1])/length(X)
q = 1 - p 
#vcase = vp(P[i,])
vcase = var(X[Y==1])
vcont = var(X[Y==0])
Tobs = (mean(X[Y==1]) - mean(X[Y==0]))/sqrt(vcase/case+vcont/cont) 
X1 = X[Y==1] - mean(X[Y==1])
X2 = X[Y==0] - mean(X[Y==0])
C = NULL
n = perm
for (j in 1:n){
Xca = sample(X1[],case,replace=TRUE)
Xco = sample(X2[],cont,replace=TRUE)
vcase = var(Xca)
vcont = var(Xco)
C =c(C,(mean(Xca) - mean(Xco))/sqrt(vcase/case+vcont/cont))
}
cc[i] = sum(abs(C)>=abs(Tobs))/n
}
return (cc)
}



###
# This is helper function
###
# 1) vector of base call errors.
# 2) genotype, vector of two letters.

## Input
# error - vector of errors (double) of length of debth
# debth - integer
# genotype - vector of size 2 (e.g. 'A','C')

## Output
# vector of row sequence reads ('A','G','T','C') for single variant (length debth) 
get_row_reads = function(error,genotype,debth){
vect_row_reads = NULL
etal = c('A','G','T','C')
for (i in 1:debth){
s = rbinom(1,1,0.5)+1
g = genotype[s]
ng = etal[etal!=g]
e3 = error[i]/3
e23 = e3*2
e = e3*3
r = runif(1)
if (r>e) {vect_row_reads = c(vect_row_reads,g)}
if (r<e3){vect_row_reads = c(vect_row_reads,ng[1])}
if (r>e3 & r<e23){vect_row_reads = c(vect_row_reads,ng[2])}
if (r>e23 & r<e){vect_row_reads = c(vect_row_reads,ng[3])}
}
return (vect_row_reads)
}


##
# This is help function
##
## 
# Calculate genotype likelihood for single sequece read

##Input
# A1, A2 - two reference alleles (e.g. A1-'A',A2='C')
# error - sequencing error (double)
# vect_row_read - single read (e.g. 'A')

##Output
# p - likelihood of A1 and A2 for single read

get_L = function(vect_row_read,A1,A2,error){
if (vect_row_read == A1){
p1 = 1-error
}else{
p1 = error/3
}

if (vect_row_read == A2){
p2 = 1-error
}else{
p2 = error/3
}
p = 1/2*p1 + 1/2*p2

return (p)
}


# This is helper function. It produces sequence call from given genotype
# Require vector of errors = Phen scores
# gen_vect = vector of prior genotypes
# Lastly vector of genotype frequencies.
seq_call = function(genotype, error, gen_vect){
debth = length(error)
vect_row_reads = get_row_reads(error,genotype,debth)
return (vect_row_reads)
}


#get conditional expected value of genotype
# This is helper function 
#
##Input:
#
# M - genotype likelihoods AA, Aa, aa, matrix sum(rdv) by 3 (double) 
# p - genotype frequencies AA, Aa, aa (double)
# rdv - read depth (vector of integers) for all samples

##Output:
# EG - conditional expectation ()

get_EG = function(M,p,rdv){
L = length(rdv)
S = 0
EG = NULL
g = c(0,1,2)
for (i in 1:L){
m = NULL
for  (j in 1:3){
L = 1
for (kk in 1:rdv[i]){

L = L*M[S + kk,j]
}
m = c(m,L*p[j])
}
S = S + rdv[i]
pm = sum(m/sum(m)*g)
EG = c(EG,pm)
}
return (EG)
}

#get conditional variance value of genotype
get_EG2 = function(M,p,rdv){
L = length(rdv)
S = 0
EG = NULL
g = c(0,1,2)
for (i in 1:L){
m = NULL
for  (j in 1:3){
L = 1
for (kk in 1:rdv[i]){

L = L*M[S + kk,j]
}
m = c(m,L*p[j])
}
S = S + rdv[i]
pm = sum(m/sum(m)*g^2) - sum(m/sum(m)*g)^2
EG = c(EG,pm)
}
return (EG)
}

# Get likelihood P(Di|G)
# This is helper function 

##Input:
#  Error - error rates
#  vect_row_reads - vector of reads
#

## Output:
# M - matrix of P(read|AA), P(read|Aa) and P(read|aa)
# M - matrix length(Error) by 3

get_M = function(Error,vect_row_reads){
L = length(Error)
M = NULL
LG  =c('TT','CT','CC')

for (i in 1:L){
m = NULL
for  (j in 1:3){
G = LG[j]
A1 = substring(G,1,1)
A2 = substring(G,2,2)
m = c(m,get_L(vect_row_reads[i],A1,A2,Error[i]))
}
M = rbind(M,m)
}
return (M)
}

# Get likelihood P(Di|G)
# This is helper function 

##Input:
#  Error - error rates
#  vect_row_reads - vector of reads
#

## Output:
# M - matrix of P(reads|AA), P(reads|Aa) and P(reads|aa)
# M - matrix length(rdv) by 3
#

get_Mr = function(Error,vect_row_reads,rdv){
L = length(rdv)
M = NULL
LG  =c('TT','CT','CC')
S = 0
for (i in 1:L){
m = NULL
for  (j in 1:3){
LL = 1
for (kk in 1:rdv[i]){
G = LG[j]
A1 = substring(G,1,1)
A2 = substring(G,2,2)
LL = LL*get_L(vect_row_reads[S+kk],A1,A2,Error[S+kk])
}
m =c(m,LL)
}
S =S+rdv[i]
M = rbind(M,m)
}
return (M)
}


# Variance of genotypes
vp = function(P){
Sq = 4*P[3] + P[2]
Sm = 2*P[3] + P[2]
S = Sq - Sm^2
return (S)
}

# Calculate score test
calct = function(M1,M2){
X = c(M1,M2)
Y = c(rep(1,length(M1)),rep(0,length(M2)))
p = length(M1)/length(X)
q = 1 - p
S = q*sum(M1)-p*sum(M2)
vs = p*q*(length(X))*var(X)
x = S^2/vs
return (x)
}
