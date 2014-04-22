

# I think this is the hardest fucntion to explain
#
##########################################
# This function should be in package #####
##########################################

#
# This functions gets genotype likelihoods from vcf file
# VCF file should be in specific form



# These are helper functions for the analysis based on
# vcf files. 
# Step 1: read vcf files:
# 	Separate vcf file functions based on data
# 	Output: ref allele and alternative allele.
#	Output: variant locations.
#   Output: Indicator for PASS
#	Output: likelihood: vector (L(D|0),L(D|1),L(D|2)).
# VCF should be in the next form GT:AD:DP:GQ:PL 
# For exampe: 0/0:2,0:4:6:0,6,42
# PL = Phred-scaled likelihoods for genotypes as defined in the VCF specification in the format of 1000 Genome project
# Phred-scaled likelihoods: log10(L(D|0)),log10(L(D|1)),log10(L(D|2)).

## Input:
# Matrix M - vcf file without header 
# 2nd column - location
# 4th column - reference allele
# 5th column - alternative allele
# 7th column - Filter False/True
## Output:
# Object with three matrices:
# L0 - matrix of P(D|0), rows are individuals and columns are SNPS
# L1 - matrix of P(D|1), rows are individuals and columns are SNPS
# L2 - matrix of P(D|2), rows are individuals and columns are SNPS


get_L_vcf = function(M){
RA = NULL
AA = NULL
VL = NULL
VI = NULL
LD0 = NULL
LD1 = NULL
LD2 = NULL
Lv = length(M[,1])
Lind = length(M[1,])

for (i in 1:Lv){
RA = c(RA,as.vector(M[i,4]))
AA = c(AA,as.vector(M[i,5]))
VL = c(VL,M[i,2])
VI = c(VI,as.vector(M[i,7]))
l0 = NULL
l1 = NULL
l2 = NULL
for (j in 10:Lind){
LL = get_l(M[i,j])
l0 = c(l0,LL[1])
l1 = c(l1,LL[2])
l2 = c(l2,LL[3])
}
LD0 = cbind(LD0,l0)
LD1 = cbind(LD1,l1)
LD2 = cbind(LD2,l2)
}
return (list(ref_a = RA, alt_a = AA, var_loc=VL, var_indic = VI, L_matrix=list(L0 = LD0,L1 = LD1, L2 = LD2)))
}

# Helper function for 1000G data
# gives a vector of likelihood
# L0,L1,L2

get_l = function(M){
a = unlist(strsplit(as.vector(M),':'))
l = as.numeric(unlist(strsplit(a[3],',')))
l0 = l[1]
l1 = l[2]
l2 = l[3]
Pl0 = 10^l0
Pl1 = 10^l1
Pl2 = 10^l2
return (c(Pl0,Pl1,Pl2))
}
