# This part is new part of the analysis 
# We can use reads instead of hard calls
# For this to work we need to estimate genotype frequencies.

# let's use EM algorithm

# Data consists of matrix n by 3
# p = P(G=0)
# q = P(G=1)

EM_calc = function(M){
p_0 = 0.15
q_0 = 0.15

q_n = 1
p_n = 0
d_n = 0
k = 0
while ((p_n - p_0)^2 + (q_n - q_0)^2>0.000001){
d_0 = 1-p_0 - q_0
v = c(p_0,q_0,d_0)
p_D = M%*%(v)
E_p = M[,1]*p_0/p_D
E_q = M[,2]*q_0/p_D
p_n= p_0
q_n = q_0
d_n = 1-q_0-p_0
p_0 = sum(E_p)/length(E_p)
q_0 = sum(E_q)/length(E_q)
k = k+1
if (k==1000){
cat('hi','\n')
return (c(p_0,q_0,1-p_0-q_0))
}
}

return (c(p_0,q_0,1-p_0-q_0))

}