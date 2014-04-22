#
# Read vcf file
#

a = 'C:/chr11_113low_56high/1g115low_1g56exomehigh_filtered.hg19.chr11.vcf'

#
# Read vcf helper functions
#

source('likelihood_vcf.r')
filen = a
filecon = file(filen, open='r')
#
# Skip header of vcf file.
# n = may be changed until reach header that contains list of samples
tt2 = readLines(filecon, n=128)
#
# S contains list of 169 samples
# One should change accordingly
#
S = unlist(strsplit(tt2[128],'\t'))[10:178]


# Genotype likelihoods
A0M = NULL
A1M = NULL
A2M = NULL
# Cordinates
Cord = NULL


s = TRUE
l = 0
while (s){
l=l+1
F = try(read.table(filecon, nrows = 10000, sep='\t'), silent = TRUE)
if ( class(F) == 'try-error'){break}
if ( length(F) == 0){break}

# Contains all information.
AA = get_L_vcf(F)

}




