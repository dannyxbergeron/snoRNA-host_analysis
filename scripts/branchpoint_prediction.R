# For the installation of branchpointer, for some reason it doesn't work
# with snakemake. Accordingly, it must be installed mannually by activating
# the conda environnement, running R, and run the 2 following commands:
#
# install.packages("BiocManager")
# BiocManager::install("branchpointer")

library(branchpointer)

exons <- gtfToExons(snakemake@input[["mini_gtf"]])

queryIntron <- readQueryFile(snakemake@input[["intron_file"]],
                             queryType = "region",
                             exons = exons)

branchpointPredictionsIntron <- predictBranchpoints(queryIntron,
                                                    queryType = "region",
                                                    genome=snakemake@input[["genome"]],
                                                    bedtoolsLocation=snakemake@params[["bedtools"]])
write.csv(branchpointPredictionsIntron, snakemake@output[["bp_distance"]], row.names = FALSE)
