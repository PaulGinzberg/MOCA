# Ensure that the directory containing "mocafunctions.R" is your working directory

source("mocafunctions.R")

whitelistinfo <- read.table("S1_Table.txt",header = TRUE)

# Run the MOCA pipeline on the TCGA breast cancer dataset.
out1 <- moca.pipeline("./example-data/brca-tcga.dat",gene.whitelist="./example-data/genewhitelist.txt",save.results="./example-output/brca-example",
              test.by.prefix=TRUE,max.merges=1e2,max.size=4)

# Create a vector of input files to use.
examplefiles <- list.files(path="./example-data/",pattern="tcga.dat",full.names=TRUE)
# This is equivalent to examplefiles <- c("./example-data/brca-tcga.dat","./example-data/prad-tcga.dat")

#Run the MOCA pipeline on the TCGA breast cancer dataset AND the TCGA prostate cancer dataset
out2 <- moca.pipeline(examplefiles,gene.whitelist="./example-data/genewhitelist.txt",save.results="./example-output/brcaANDprad-example",
              test.by.prefix=TRUE,merge.duplicates=TRUE,max.merges=1e2,max.size=10,max.p=0.01)










