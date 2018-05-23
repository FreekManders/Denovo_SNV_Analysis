#Remove white spaces from ensembl regbuild.

#zcat Ensemblregbuild.bed.gz | tr ' ' '_' | bgzip > Ensemblregbuild2.bed.gz
#tabix -p bed Ensemblregbuild.bed.gz
