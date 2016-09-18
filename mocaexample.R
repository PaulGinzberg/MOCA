
source("mocafunctions.R")

# A table with the gene symbol/name, EntrezID and source for all genes in the whitelist
whitelistinfo <- read.table("S1_Table.txt",header = TRUE)
# The file location of the whitelist (which uses the same naming convention as the data)
gene.whitelist <- "./example-data/genewhitelist.txt"

output1 <- moca.pipeline("./example-data/prad-tcga-somatic-mutations.dat",gene.whitelist=gene.whitelist,save.results="./example-output/brca-example",
                         test.by.prefix=TRUE,max.merges=1e2,max.size=4)
# (This should find 0 significant gene sets)

# Extract the data files from example-data.zip into the directory /MOCA/example-data/ before running the code below

# The following examples may take minutes to run.
# Run the MOCA pipeline on the TCGA breast cancer dataset.
output2 <- moca.pipeline("./example-data/brca-tcga-somatic-mutations.dat",gene.whitelist=gene.whitelist,save.results="./example-output/brca-example",
                         test.by.prefix=TRUE,max.merges=1e2,max.size=4)

# Run the MOCA pipeline on the TCGA breast cancer dataset AND the TCGA prostate cancer dataset, combining p-values.
examplefiles <- c("./example-data/brca-tcga-somatic-mutations.dat","./example-data/prad-tcga-somatic-mutations.dat")
output3 <- moca.pipeline(examplefiles,gene.whitelist=gene.whitelist,save.results="./example-output/brcaANDprad-example",
                         test.by.prefix=TRUE,merge.duplicates=TRUE,max.merges=1e2,max.size=10,max.p=0.01)
# Warning messages are due to the fact that some candidate sets get an automatic p-value of 1 because
# all alterations in the gene set are from the same gene. Ignore them safely.

# # To run on all datasets combined, and reproduce our results (may take hours to run), use instead
# examplefiles <- list.files(path="./example-data/",pattern="mutations.dat",full.names=TRUE)
# output4 <- moca.pipeline(moca.data.files,gene.whitelist=gene.whitelist,save.results="./example-output/pancancer-test2",
#                         test.by.prefix=TRUE,merge.duplicates=TRUE,max.merges=5e3,max.size=10,max.p=0.05,randomised=TRUE)
