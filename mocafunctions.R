

moca.pipeline <- function(mat, groups=NULL, gene.whitelist=NULL, gene.screening.function=potentially.significant,
                          max.size=5, max.merges=1e3, max.dist=1, max.p=0.05, alpha.bonferroni=max.p, max.interesting=10,
                          consider=c("subclusters","original","maximalonly"),
                          merge.types=FALSE,merge.duplicates=FALSE,test.by.prefix=FALSE,
                          randomised=FALSE,
                          save.results=NULL,na.action=NA.as.FALSE) {
  # This function runs the whole moca.pipeline starting from a binary matrix of mutations.
  # 
  # INPUT:
  # mat: Either 1) A binary matrix with genes/alterations in rows and samples in columns
  #             2) The filepath of a file containing such a matrix as a plain text table
  #             3) A character vector of filepaths from which to load multiple such matrices.
  #      Note that some features assume that the row names of mat are of the form prefix_suffix,
  #      Where prefix is the gene name / gene ID, and suffix is the alteration type (e.g. SNV/AMP/DEL)
  # groups: A factor with the same length as ncol(mat) which groups samples by cancer type. 
  #         Statistical analysis is perfomed on each cancer type separately and then combined using Stouffer's method.
  #         groups will be generated automatically if multiple mutation matrices are loaded.
  # gene.whitelist: Either 1) A vector of the rownames in mat which should be used (length>>1)
  #                        2) The filepath of a file containing such a vector in plain test (length==1)
  # gene.screening.function: A function to select genes based on their mutation counts (input=mat, output=mat with some rows removed)
  #                          By default only genes with a number of mutations greater than log2(number of genes kept) are kept.
  #                          Other options would be for example function(x){x[rowSums(x)>=10,,drop=FALSE]} or function(x){x[rowSums(x)>=0.1*ncol(x),,drop=FALSE]}.
  #                          For no screening, set gene.screening.function=identity.
  # max.size: The maximum allowed size for a gene set
  # max.merges: The maximum number of candidate gene sets to be generated greedily
  # max.dist: Maximum allowed pvalue for merges (irrelevant if >=1 as in default setting).
  #           A sensible lower choice for this value might be 1/(number of genes used), since no multiple hypothesis testing
  #           correction is performed at this stage.
  # max.p: The threshold for significance to report a gene set. (Note that a weighted Bonferroni correction is performed)
  # alpha.bonferroni: (An upper bound on) the threshold for significance to include an additional alteration in a gene set.
  #                   equal to max.p by default.
  # max.interesting: The maximum number of interesting examples to report.
  # consider: If "original", then all candidate gene sets are considered.
  #           If "maximalonly", only maximal candidate gene sets are considered. (not recommended)
  #           If "subclusters" (default) then subsets of the original candidate gene sets are also considered. This might be slow if max.size is large.
  # merge.types: If TRUE, suffixes in the row names of mat are ignored and rows with matching prefixes are merged. (can be used to ignore the distrinction between SNV, AMP and DEL alterations)
  # merge.duplicates: If TRUE, rows of mat that have the same mutation pattern are merged. The suffix ORO is used in naming any such merged rows.
  # test.by.prefix: If TRUE, when computing p-values the amount of cooccurrence between genes having the same prefix is ignored.
  #                 Although the default is FALSE, this should be set to TRUE if both amplifications and deletions of a gene are considered,
  #                 and the prefix_suffix naming convension is followed.
  # randomised: Whether to use randomised p-values to mitigate the conservative nature of discrete p-values. This does not apply to the generation of candidate gene sets.
  # save.results: (Optional) A name to use when saving intermediate computational results are saved to *.Rdata files
  # na.action: A function to handle missing values. By default NA.as.FALSE: missing values in mat will be treated as FALSE.
  #            Set this to na.omit to omit genes with missing values, and to na.pass to handle missing value in each later calculation
  # 
  # OUTPUT:
  # 
  # 
  elapsed <- proc.time()
  # Read and pre-process the data
  consider <- match.arg(consider)
  call <- match.call()
  cat("\nReading and pre-processing inputs...")
  if(!is.null(gene.whitelist)&&length(gene.whitelist)<2) {
    gene.whitelist <- as.character(read.table(gene.whitelist)[[1]])
  }
  if(is.character(mat)) {
    matdir <- mat
    matlist <- list()
    for(k in 1:length(matdir)) {
      mat <- read.table(matdir[k],header=TRUE,row.names=1)
      if(!is.null(gene.whitelist)) {
        mat <- mat[rownames(mat) %in% gene.whitelist,,drop=FALSE]
      }
      matlist <- c(matlist,list(mat))
    }
    if(length(matdir)>1) {
      matlist <- matlist2mat(matlist)
      mat <- matlist$mat
      if(!is.null(groups)) {warning("Argument 'groups' overwritten because multiple matrices were loaded.")}
      groups <- matlist$groups
    }
    rm(matlist)
  }
  mat <- as.matrix(mat)
  if(!is.logical(mat)) {
    if(any(range(mat,na.rm=TRUE)!=c(0,1))) {
      warning("mat converted to logical matrix from ",mode(mat),". Formatting of input data may not be supported.")
    }
    mat <- mat==1
  }
  mat <- na.action(mat)
  if(merge.types) { # Merge e.g. SNP AMP and DEL versions of each gene
    rownames(mat) <- suffixrm(rownames(mat))
    mat <- duplicatemerge(mat)
  }
  if(merge.duplicates) {
    # Merge genes that have exactly the same pattern of mutations
    matd<-duplicaterm(mat)
    #matorig <- mat
    mat <- matd$mat
    matd$mat <- NULL
  } else {
    matd <- list(dlistseed=list(),dlist=list())
  }
  mat <- gene.screening.function(mat)
  if(any(rowMeans(mat,na.rm=TRUE)%in%c(0,1))) { # Remove any genes which are never/always mutated.
    warning(paste(rownames(mat)[which(rowMeans(mat,na.rm=TRUE) %in% c(0,1))],collapse=", ")," is mutated in all/no samples and has been removed.")
    mat <- mat[!(rowMeans(mat,na.rm=TRUE)%in%c(0,1)),]
  }
  min.samples <- floor(log2(nrow(mat)))+1
  cat("DONE\n")
  cat("There are ",nrow(mat)," genes and ",ncol(mat),"samples after pre-processing\n")
  if(nrow(mat)<2||ncol(mat)<2) {stop("Dataset mat is too small (only ",nrow(mat)," by ",ncol(mat)," )")}
  
  # Generate candidate clusters
  clustout <- cooccurrence.clustering(mat,groups=groups,max.merges=max.merges,max.size=max.size,max.dist=max.dist)
  
  if(consider=="subclusters") {
    clusters <- maximal.clusters(clustout$merge,onlymaximal=TRUE)$clusters
    clusters <- all.subclusters(clusters)
  } else if(consider=="original") {
    clusters <- maximal.clusters(clustout$merge,onlymaximal=FALSE)$clusters
  } else if(consider=="maximalonly") {
    clusters <- maximal.clusters(clustout$merge,onlymaximal=TRUE)$clusters
  }
  cat(length(clusters)," candidate gene sets were generated.\n")
  
  cat("Computing p-values...")
  lengths<-sapply(clusters,length)
  if(test.by.prefix) {
    coverageall <- sapply(clusters,function(x)meta.coverage.test(duplicatemerge(mat[x,],names=suffixrm(rownames(mat)[x])),groups=groups,randomised=randomised))
  } else {
    coverageall<-sapply(clusters,function(x)meta.coverage.test(mat[x,],groups=groups,randomised=randomised))
  }
  # Apply weighted Bonferroni Correction
  coverageallc<-coverageall*sapply(lengths,subset.bonferroni,n=nrow(mat),order="size",kmax=max.size,alpha=alpha.bonferroni)
  cat("DONE\n")
  
  # Select Significant clusters
  significantsets<-representative.clusters(clusters,mat,n=Inf,max.p=max.p,pvals=coverageallc)
  # Select interesting clusters
  interestingsets <- significantsets
  maxoverlap <- max.size+1
  while(length(interestingsets$clusters)>max.interesting && maxoverlap>0) {
    maxoverlap <- maxoverlap-1
    interestingsets <- representative.clusters(clusters,mat,n=Inf,max.p=max.p,pvals=coverageallc,maxoverlap=maxoverlap)
  }
  if(length(interestingsets$clusters)==0) {
    interestingsets <- representative.clusters(clusters,mat,n=1,max.p=Inf,pvals=coverageallc,maxoverlap=maxoverlap)
  } else if(length(interestingsets$clusters)>max.interesting) {
    interestingsets <- representative.clusters(clusters,mat,n=max.interesting,max.p=Inf,pvals=coverageallc,maxoverlap=maxoverlap)
  }
  interestingsets$clusters <- relist(rownames(mat)[unlist(interestingsets$clusters)],skeleton = interestingsets$clusters)
  significantsets$clusters <- relist(rownames(mat)[unlist(significantsets$clusters)],skeleton = significantsets$clusters)
  connectedcomponents <- get.connected.components(significantsets$clusters)
  elapsed <- proc.time()-elapsed
  cat(length(significantsets$clusters)," significant gene sets were found, and combined they form a graph with ",length(connectedcomponents)," connected components.\n")
  output <- list(call=call,elapsed=elapsed,significantsets=significantsets,connectedcomponents=connectedcomponents,interestingsets=interestingsets,mat=mat,matd=matd,groups=groups,clustout=clustout,clusters=clusters,raw.p.values=coverageall,corrected.p.values=coverageallc,min.samples=min.samples,maxoverlap=maxoverlap)
  if(!is.null(save.results)) {
    cat("Writing output files ... ")
    #save(call,clustout,elapsed,clusters,significantsets,connectedcomponents,mat,matd,groups,coverageall,coverageallc,maxoverlap,interestingsets,min.samples, file=paste0(save.results,".Rdata"))
    save(output, file=paste0(save.results,".Rdata"))
    description <- data.frame(p.value=coverageall[interestingsets$selection],p.value.corrected=interestingsets$pvals,gene.set=sapply(interestingsets$clusters,paste,collapse=", "))
    write.table(description,paste0(save.results,"-interesting.txt"),row.names=FALSE)
    description <- data.frame(p.value=coverageall[significantsets$selection],p.value.corrected=significantsets$pvals,gene.set=sapply(significantsets$clusters,paste,collapse=", "))
    write.table(description,paste0(save.results,"-significant.txt"),row.names=FALSE)
    description <- data.frame(gene.set=sapply(connectedcomponents,paste,collapse=", "))
    write.table(description,paste0(save.results,"-connectedcomponents.txt"),row.names=FALSE)
    cat("DONE\n")
  }
  #return(list(call=call,elapsed=elapsed,significantsets=significantsets,connectedcomponents=connectedcomponents,interestingsets=interestingsets,mat=mat,matd=matd,groups=groups,clustout=clustout,clusters=clusters,raw.p.values=coverageall,corrected.p.values=coverageallc,min.samples=min.samples,maxoverlap=maxoverlap))
  return(output)
}

get.connected.components <- function(clusters) {
  # Gives the connected components of the graph obtained by adding an edge between two genes
  # whenever they share a cluster.
  # INPUT:
  # clusters: A list of gene sets.
  # OUTPUT:
  # components: A list of gene sets, where any gene sets with genes in common are merged.
  components<-list()
  while(length(clusters)>=1){
    current <- clusters[[1]]
    matches <- sapply(clusters,function(x)any(x %in% current))
    if(sum(matches)>1) {
      clusters[[1]] <- Reduce(union,clusters[matches])
      matches[1] <- FALSE
      clusters[matches] <- NULL
    } else {
      components <- c(components,clusters[1])
      clusters[[1]] <- NULL
    }
  }
  return(components)
}

NA.as.FALSE <- function(x){
  # Convert NA entries to FALSE
  x[is.na(x)]<-FALSE
  return(x)
}
potentially.significant <- function(mat,OtherGenesCoverage=0.5,min.genes=0,round.up=TRUE) {
  # Function to filter out genes which do not have a reasonable chance of being
  # significantly anti-co-occurring after Bonferroni correction.
  # 
  # INPUT:
  # mat: A binary matrix of mutations, where each row corresponds to a gene
  # OtherGenesCoverage: (default 50%) The hypothetical (maximum) coverage of the other gene(s) against which
  #                     each gene might be tested. See DETAILS.
  # min.genes: A minimum number of genes to keep in the output (overrides everything else)
  # round.up: It is possible that including some but not all of the genes with a certain minimum mutation
  #           count makes significance possible. If round.up=TRUE (default) none of them are included.
  #           If round.up=FALSE then all of those genes are included (even though they will not be significant).
  #
  # OUTPUT:
  # mat: same as imput mat, but without the rows which do not have a reasonable chance of being significant.
  #
  # DETAILS:
  # Only genes with mutation count greater or equal to min.samples will be kept.
  # To determin min.samples we use a heuristic based on the binomial approximation,
  # and the idea that it should be possible for the remaining genes to be significantly 
  # anti-co-occurring with respect to a gene that is mutated in a proportion
  # OtherGenesCoverage of the samples (half of the samples by default).
  # A Bonferroni correction of 1/(number of genes kept) is applied.
  samplespergene <- rowSums(mat,na.rm=TRUE)
  #samplespergene <- sort(samplespergene2,decreasing=TRUE)
  #plot(samplespergene,log="y"); lines(1:length(samplespergene),log2(1:length(samplespergene))); #log(1:length(samplespergene),1/(1-samplespergene/ncol(mat)))
  min.samples <- min(samplespergene[samplespergene>log(rank(-samplespergene,ties.method=ifelse(round.up,"max","first"))-1,1/(1-OtherGenesCoverage))])
  # Remove rare mutations
  if(sum(samplespergene>=min.samples)>=min.genes) {
    mat <- mat[samplespergene>=min.samples,,drop=FALSE]
  } else {
    mat <- mat[order(samplespergene,decreasing=TRUE)[1:min.genes],,drop=FALSE]
    warning("min.samples was not enforced due to too few genes remaining. Top ",min.genes," genes used instead.")
  }
  return(mat)
}


duplicaterm<- function(mat) {
  # Remove duplicate rows from mat, but keep track of their names.
  # The purpose of duplicaterm is to only keep one example 
  # when multiple genes have the same mutation pattern.
  # Genes that are only different in samples that have missing values are
  # considered as having the same mutation pattern.
  # The OUTPUT is a list containing
  # dlistseed: a list of the rownames in the new mat matrix which had duplicates in the old mat matrix.
  # Note that the 
  # dlist: a list of vectors of the rownames in the old mat matrix which had the same mutation pattern.
  # The genes in dlist[[i]] have been replaced by the single gene dlistseed[[i]].
  # mat: a new possibly smaller version of the input matrix mat, without duplicate rows.
  # Note that the name used for rows which were duplicated will end with _oro instead of _AMP/_DEL/_SNV
  # oro stands for "or other".
  
  matred <- mat[,!apply(mat,2,function(x)any(is.na(x)))]
  checkforcopies <- function(y) apply(matred,1,function(x)all(x==y,na.rm=TRUE))
  ncopies <- function(y) sum(checkforcopies(y))
  duplicates <- which(duplicated(matred) | duplicated(matred, fromLast = TRUE))
  dlist<-list()
  dlistseed<-list()
  if(length(duplicates)==0) {return(list(dlistseed=dlistseed,dlist=dlist,mat=mat))}
  while(length(duplicates)>0) {
    #print(length(duplicates))
    sel<-duplicates[1]
    aka<-setdiff(which(checkforcopies(mat[sel,])),sel)
    dlist<-c(dlist,list(rownames(mat)[aka]))
    dlistseed<-c(dlistseed,list(rownames(mat)[sel]))
    duplicates<-setdiff(duplicates,c(sel,aka))
  }
  mat<-mat[!(rownames(mat) %in% unlist(dlist)),]
  mn<-match(unlist(dlistseed),rownames(mat))
  namesplus<-strsplit(rownames(mat),split="_")
  namesplus[mn]<-lapply(namesplus[mn],function(x) return(c(x[1],"oro")))
  rownames(mat)<-sapply(namesplus,paste,collapse="_")
  dlist<-mapply(c,dlistseed,dlist)
  names(dlist) <- rownames(mat)[mn]
  dlistseed<-as.list(rownames(mat)[mn])
  return(list(dlistseed=dlistseed,dlist=dlist,mat=mat))
}

duplicatemerge <- function(mat,names=rownames(mat)) {
  # Merges rows of the matrix mat that have the same names
  # Can be used to merge SNP AMP and DEL versions of each gene,
  # so that the gene is considered mutated if any of those three types of alteration is present.
  if(anyDuplicated(names)!=0) {
    for(name in unique(names)) {
      mat[match(name,names),] <- apply(mat[names==name,,drop=FALSE],2,any)#,na.rm=TRUE)
    }
    mat <- mat[match(unique(names),names),,drop=FALSE]
  }
  return(mat)
}


matlist2mat <- function(matlist) {
  # convert a list of matrices/data.frames matlist to
  # a single matrix mat and factor groups indicating which columns
  # came from which matrix in matlist
  groupcounts <- sapply(matlist,ncol)
  groupnames <- names(matlist)
  if(is.null(groupnames)) {groupnames <- paste0("group",1:length(matlist))}
  groups <- factor(rep(groupnames,times=groupcounts))
  
  matlist <- lapply(matlist,function(x)cbind(data.frame(NAMEformerging=rownames(x)),x))
  mat <- matlist[[1]]
  if(length(matlist)>1) {
    for(k in 2:length(matlist)) {
      mat <- merge(mat,matlist[[k]],by=c("NAMEformerging"),all=TRUE)
    }
  }
  rownames(mat) <- mat$NAMEformerging
  mat$NAMEformerging <- NULL
  mat <- as.matrix(mat)
  return(list(mat=mat,groups=groups))
}

cooccurrence.clustering<-function(matlist, groups=NULL, max.merges=1e3, max.dist=1, max.size=Inf) {
  # Simplified / stripped down version.
  # Generates candidate sets of mutually exclusive mutations through a greedy algorithm similar to hierarchical clustering.
  # The output is a series of merges with heights (similar to hierarchical clustering) 
  # and a matrix containing the mutation pattern for each set considered.
  # INPUT
  # matlist: Either a binary matrix with genes in rows and samples in columns, or a list of such binary matrices having the same number of rows
  # groups: a factor indicating the source of the columns of matlist, which if provided is used to split the single matrix matlist into a list of matrices
  # Only anti-cooccurrence within each group will contribute to significance, not between groups.
  # max.merges: the maximum number of clusters to generate
  # max.dist: the maximum pvalue or height at which to merge two clusters. This should be set as a significance level.
  # max.size: the maximum allowed number of genes in each cluster.
  # OUTPUT
  # A list containing
  # merge: an n by 2 matrix, where n is the number of clusters generated. Row i of merge describes the merging of clusters at step i of the clustering. If an element j in the row is negative, then observation -j was merged at this stage. If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm. Thus negative entries in merge indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons.
  # height: a set of n pvalues corresponding to the n merges
  # extended.mat: an extended version of matlist with one additional row for each cluster
  
  if(class(matlist)=="matrix") {
    if(is.null(groups)) {
      matlist <- list(onlyone=matlist)
    } else {
      matlist <- by(t(matlist),INDICES=groups,FUN=t,simplify=FALSE)
    }
  }
  Nmats <- length(matlist)
  Nrows <- nrow(matlist[[1]])
  norig <- Nrows
  
  analysing.fun <- function(ij) {
    # Function which returns the p-value for a Fisher exact test
    # between the ij[1] and ij[2] rows. If Nmats>1 then the multiple tests are combined
    
    # Note that if any of the individual pvalues is 1 then the combined pvalue is 1.
    # This means that including information from
    # small samples has the potential of "breaking" the meta-analysis.
    # To avoid this, we use the mid-p-value instead of standard p-values.
    # Using the mid p-value also means lack of evidence against the null leads to lower
    # p-values than weak evidence in favour of the null.
    # Note however that this choice leads to type I error control being approximate.
    
    # Earlier solution: (implemented by commented-out code marked with "AVOIDONES")
    # To avoid this, a group is ignored entirely when
    # a rough estimate of the probability of getting a pvalue of 1 is larger than the 
    # group's normalised (squared) weight in the Stouffer method.
    # Also note that the weights are based on asymptotic theory,
    # so small samples will probably not be given an optimal weight and might create avoidable problems.
    
    #w2sum<-0 # Normalisation for weights
    #meta.p<-0
    at.least.one.computable <- FALSE
    #xytable <- matrix(0,nrow=2,ncol=2)
    
    # initialise quantities
    selectionscore <- numeric(0)
    weights <- numeric(0)
    pvec <- numeric(0)
    
    for(k in 1:Nmats) {
      #xytable <- table(matlist[[k]][ij[1],],matlist[[k]][ij[2],]) #table() is slow so I wrote a simplified version from scratch
      v1 <- matlist[[k]][ij[1],]
      v2 <- matlist[[k]][ij[2],]
      v12nas <- is.na(v1)|is.na(v2)
      v1[v12nas] <- NA
      v2[v12nas] <- NA
      v1s <- sum(v1,na.rm=TRUE)
      v2s <- sum(v2,na.rm=TRUE)
      v12s <- sum(v1&v2,na.rm=TRUE)
      #xytable[2,2] <- v12s
      #xytable[1,2] <- v2s - v12s
      #xytable[2,1] <- v1s - v12s
      #xytable[1,1] <- length(v1) - sum(v12nas) - sum(xytable)
      Ntot <- length(v1)-sum(v12nas)
      
      #if(length(dim(xytable))==2 && all(dim(xytable)==c(2,2)) && v1s>0 && v2s>0)
      if(v1s>0 && v2s>0) {
        at.least.one.computable <- TRUE
        if(Nmats==1) {
          return(phyper(v12s,v1s,Ntot-v1s,v2s))
        } else {
          ##w <- meta.weight(sum(xytable),c(xytable[2,1]+xytable[2,2],xytable[1,2]+xytable[2,2])) # weight
          w <- meta.weight(Ntot,c(v1s,v2s))
          weights <- c(weights,w)
          
          ##AVOIDONES
          #selectionscore <-c(selectionscore,w*(max(v1s,v2s)/Ntot)^-min(v1s,v2s))
          
        }
      } else {
        w <- 0
        #p <- 1
      }
      if(w>0) {
        #p <- fisher.test(xytable,alternative="less",conf.int=FALSE)$p.value # This was unnecessarily slow so replaced with direct p-value calculation
        #p <- phyper(xytable[2,2],xytable[2,1]+xytable[2,2],xytable[1,1]+xytable[1,2],xytable[1,2]+xytable[2,2])
        p <- phyper(v12s,v1s,Ntot-v1s,v2s)
        # Add mid p-value correction
        p <- 0.5*p + 0.5*phyper(v12s-1,v1s,Ntot-v1s,v2s)
        pvec <- c(pvec,p)
        #meta.p <- meta.p+qnorm(p,lower=FALSE)*w
        #w2sum <- w2sum + w^2
      }
    }
    if(!at.least.one.computable) {return(1)}
    
    ##AVOIDONES
    #ord <- order(selectionscore,decreasing=TRUE)
    #cumweight <- cumsum(weights[ord]^2)
    #cumweight[ord] <- cumweight
    #selection <- selectionscore>cumweight
    #return(pnorm(sum(qnorm(pvec[selection],lower=FALSE)*weights[selection]/sqrt(sum(weights[selection]^2))),lower=FALSE))
    
    return(pnorm(sum(qnorm(pvec,lower=FALSE)*weights/sqrt(sum(weights^2))),lower=FALSE))
    #return(pnorm(meta.p/sqrt(w2sum),lower=FALSE))
  }
  
  
  # Initialise quantities
  merge <- matrix(0,ncol=2,nrow=0)
  height <- vector("numeric",0)
  clusters <- as.list(1:norig)
  clustersizes <- rep(1,norig)
  isfull <- rep(FALSE,Nmats)
  
  # Do not consider genes which are always/never mutated.
  permanent.blacklist <- apply(sapply(matlist,rowMeans,na.rm=TRUE),1,function(x)all(x==0)|all(x==1))
  pair.sub <- all.pairs(Nrows)
  pair.sub <- pair.sub[,!((pair.sub[1,] %in% which(permanent.blacklist))|(pair.sub[2,] %in% which(permanent.blacklist))),drop=FALSE]
  # Compute initial distances
  cat("Computing pairwise pvalues ... ")
  pair.dist <- apply(pair.sub,2,analysing.fun)
  cat("DONE\n")
  pair.clustersizes <- rep(2L,length(pair.dist))
  
  #   if(size.downweighting){
  #     sizerank <- rank(rowSums(mat,na.rm=TRUE),ties.method="max")
  #     pair.dist <- pair.dist*pmax(sizerank[pair.sub[1,]],sizerank[pair.sub[2,]])
  #   }
  
  # Main Loop
  #clustersizescount <- rep(0L,max.size)
  #clustersizescount[1] <- norig
  counter <- 0
  temp.max.size <- 2
  cat("Generating clusters of size ",temp.max.size," ... ")
  while(counter<max.merges && temp.max.size<=max.size) {
    if(counter >= max.merges*(temp.max.size-1)/(max.size-1)) {
      temp.max.size <- temp.max.size+1
      cat("DONE\nGenerating clusters of (maximum) size ",temp.max.size," ... ")
    }
    select <- which.min(pair.dist+(pair.clustersizes>temp.max.size))
    new.height <- pair.dist[select]
    if(new.height >= max.dist) {
      temp.max.size <- temp.max.size+1
      cat("DONE\nGenerating clusters of (maximum) size ",temp.max.size," ... ")
      next
    }
    pair <- pair.sub[,select]
    # Make sure the selected pair will never be selected again
    pair.dist[select] <- Inf
    new.cluster <- c(clusters[[pair[1]]],clusters[[pair[2]]])
    # Avoid creating duplicate entries.
    if(any(sapply(clusters,setequal,y=new.cluster))) next
    # Merging of the new cluster with blacklisted clusters will not be attempted.
    #blacklist <- sapply(clusters,function(x,y)length(intersect(x,y)),y=new.cluster)>0
    blacklist <- sapply(clusters,function(x,y)any(y %in% x),y=new.cluster)
    blacklist <- blacklist|permanent.blacklist
    counter <- counter+1
    merge <- rbind(merge,pair.sub[,select])    
    clusters <- c(clusters,list(new.cluster))
    height <- c(height,new.height)
    clustersizes <- c(clustersizes,length(new.cluster))
    #clustersizescount[length(new.cluster)] <- clustersizescount[length(new.cluster)]+1
    # Extend inputs
    #mat<-rbind(mat,ifelse(is.na(mat[pair[1],])|is.na(mat[pair[2],]),NA,mat[pair[1],]|mat[pair[2],])) # old version
    Nrows <- Nrows+1
    for(k in 1:Nmats) {
      matlist[[k]] <- rbind(matlist[[k]],matlist[[k]][pair[1],]|matlist[[k]][pair[2],])
      # Note that | produces the right behaviour with missing data.
      isfull[k] <- all(matlist[[k]][Nrows,],na.rm=TRUE)
    }
    # List new combinations to try
    new.pair.sub <- rbind(Nrows,(1:(Nrows-1))[!blacklist])
    #if(max.size!=Inf) {
    # Avoid merges that would lead to excessively large clusters
    new.pair.clustersizes <- clustersizes[new.pair.sub[1,]]+clustersizes[new.pair.sub[2,]]
    new.pair.sub <- new.pair.sub[,new.pair.clustersizes<=max.size,drop=FALSE]
    new.pair.clustersizes <- new.pair.clustersizes[new.pair.clustersizes<=max.size]
    #}
    if(all(isfull)|clustersizes[Nrows]>=max.size) {
      # If all samples are already mutated for a cluster, it doesn't make sense to try further merging.
      permanent.blacklist<-c(permanent.blacklist,TRUE)
      next
    } else {
      permanent.blacklist<-c(permanent.blacklist,FALSE)
    }
    # Compute additional distances
    new.pair.dist <- apply(new.pair.sub,2,analysing.fun)
    pair.dist<-c(pair.dist,new.pair.dist)
    pair.sub<-cbind(pair.sub,new.pair.sub)
    pair.clustersizes <- c(pair.clustersizes,new.pair.clustersizes)
  }
  cat("DONE\n")
  merge[merge<=norig] <- -merge[merge<=norig]
  merge[merge>norig] <- merge[merge>norig]-norig
  return(list(merge=merge,height=height,extended.mat=matlist))
}

all.pairs <- function(nvar){
  # For nvar variables, output a matrix of all pairs of variables. 
  # Each column is of the form c(i,j) corresponding to the pair (i,j).
  ind2sub <- function(ind,nrow)  c(r = ((ind-1) %% nrow) + 1, floor((ind-1) / nrow) + 1) # Turn vector index into pair of (i,j) matrix indices
  type <- matrix(FALSE,nrow=nvar,ncol=nvar)
  pair.ind <- which(lower.tri(type))
  #rm(type)
  pair.sub <- apply(matrix(pair.ind),1,ind2sub,nrow=nvar) # each column corresponds to a pair of mutations
  dimnames(pair.sub) <- NULL
  #rownames(pair.sub) <- c("i","j");
  return(pair.sub)
}

pinterhyper <- function(q,m,N,ntries=1,tolerance=1e-10) {
  # Outputs the probability that the intersection of random sets of size m[1],m[2],etc 
  # in a universe of size N is less than or equal to q.
  # ntries is the number of genes for multiple hypothesis correction (work in progress)
  # tolerance is the error tolerance for warning about probabilities >1.
  
  if(length(m)==1) return(as.numeric(q>=m))
  if(length(m)==2) return(phyper(q,m[2],N-m[2],m[1]))
  #if(length(m)==2) return(-expm1(phyper(q,m[2],N-m[2],m[1],lower.tail=FALSE,log.p=TRUE)*ntries))
  k<-m[1]
  pk<-1
  for(iter in 2:(length(m)-1)){
    # Safety checks have been commented out
    #if(any(is.na(pk))) {print("x\nx");print(which(is.na(pk)));print("iter");print(iter);print("m");print(m[iter]);print("pk");print(rbind(k,pk))}
    k<-k[pk>0&!is.na(pk)]
    pk<-pk[pk>0&!is.na(pk)]
    #if(sum(pk)<0.999|sum(pk)>1.001) {print("pksum");print(sum(pk))}
    #if(any(is.na(k))|length(k)==0) {print("kkk");print(k);}
    #if(is.na(m[iter])) {print("mmm");print(m);print(iter)}
    x<-max(min(k)-(N-m[iter]),0):min(m[iter],max(k))
    if(ntries==1) {
      p<-dhyper(rep(x,each=length(k)),m[iter],N-m[iter],rep(k,times=length(x)))
    } else {
      p<-dhyper.corrected(rep(x,each=length(k)),m[iter],N-m[iter],rep(k,times=length(x)),ntries=ntries)
    }
    pk<-colSums(matrix(p,ncol=length(x))*pk)
    k<-x
  }
  if(ntries==1){
    p<-phyper(q,m[length(m)],N-m[length(m)],k)
  } else {
    p<- -expm1(phyper(q,m[length(m)],N-m[length(m)],k,lower.tail=FALSE,log.p=TRUE)*ntries)
  }
  pk<-sum(p*pk)
  if(pk>1) {
    if(pk>1+tolerance) {warning("Produced a probability of ",pk," this will be rounded to 1")}
    return(1)
  }
  if(pk<0) {
    warning("Produced a probability of ",pk," this will be rounded to 0")
    return(0)
  }
  return(pk)
}

coverage.test <- function(mat,ntries=1,mid.pvalue=FALSE,randomised=FALSE) {
  # Statistical test for whether the coverage of a gene set is significantly large
  # given mutation counts for each gene and total sample size
  # INPUT
  # mat: a binary matrix where each row corresponds to a gene from the gene set and each column to a sample
  # ntries: (not fully implemented / obsolete) the multiple hypothesis correction (e.g. number of total genes).
  # mid.pvalue: If TRUE compute the mid p-value for the test instead of the p-value
  # randomised: If TRUE randomise the p-value so that it follows a 
  # continuous U(0,1) null distribution instead of the default conservative 
  # discrete distribution.
  # Setting randomised=TRUE overrides mid.pvalue.
  # OUTPUT
  # The pvalue of the test
  if(any(is.na(mat))) {mat<-mat[,!apply(mat,2,function(x)any(is.na(x))),drop=FALSE]}
  q<-sum(apply(!mat,2,all))
  m<-rowSums(!mat)
  names(m)<-NULL
  p.value <- pinterhyper(q,m,ncol(mat),ntries)
  if(randomised) {
    U <- runif(1)
  } else if(mid.pvalue) {
    U <- 0.5
  } else {
    # Standard p-value
    return(pinterhyper(q,m,ncol(mat),ntries))
  }
  # Weighted p-value (either mid or randomised)
  return(U*pinterhyper(q,m,ncol(mat),ntries)+(1-U)*pinterhyper(q-1,m,ncol(mat),ntries))
}

# Now let's generalise coverage.test to allow for a meta-analysis of multiple cancer types.
meta.weight <- function(N,n) {
  # Compute the unnormalised weight for Stouffer's meta-analysis given to a certain group or cancer type Z-score
  # N is the total number of samples (ncol(mat))
  # n is a vector containing the number of mutations for each gene in the set. (rowSums(mat))
  # (Both N and n are based only on data from a specific cancer type / group)
  # The weight is based on a heuristic and asymptotic approximation of the power
  
  m <- N-n
  # For each pair of genes compute the asymptotic standard deviation of the sample log-odds-ratio.
  if(length(n)==2) {
    #z*sqrt(n) is the standard deviation
    z <- sqrt(1/(n[1]*n[2])+1/(n[1]*m[2])+1/(m[1]*n[2])+1/(m[1]*m[2]))
    w <- 1/z
  } else {
    #z <- sqrt(1/outer(n,n)+1/outer(n,m)+1/outer(m,n)+1/outer(m,m))
    #w <- sum(1/z[upper.tri(z)])
    z <- 1/outer(n,n)+1/outer(n,m)+1/outer(m,n)+1/outer(m,m)
    w <- sqrt(sum(1/z[upper.tri(z)]))
  }
  # Note that 1/0=Inf so 1/(1/0+...)=0. Luckily this is what we want.
  # (NO LONGER APPLIES TO THIS VERSION: 1/w is (proportional to) the harmonic mean of standard deviations.)
  # It would be (proportional to) a good approximation for the standard deviation of our test if the overlaps between pairs were independent.
  w <- w/sqrt(N) #*sqrt(2/length(n))
  return(w)
}
meta.coverage.test <- function(mat,groups=NULL,ntries=1,avoidones=FALSE,randomised=FALSE) {
  # Statistical test for whether the coverage of a gene set is significantly large
  # given mutation counts for each gene and total sample size
  # INPUT
  # mat: a binary matrix where each row corresponds to a gene from the gene set and each column to a sample
  # groups: a factor for column grouping (e.g. type of cancer for each sample).
  # We condition on the mutation counts in each group so that between-group patterns are not counted.
  # pvalues from different groups are combined using Stouffer's method.
  # ntries: (not fully implemented) the multiple hypothesis correction (e.g. number of total genes).
  # avoidones: if TRUE, a heuristic will be used to avoid pvalues of 1 from small groups. (additional comments in code)
  # OUTPUT
  # The pvalue of the test
  # We condition on the mutation counts in each group so that between-group patterns are not counted.
  if(is.null(groups)) {return(coverage.test(mat,ntries=ntries))}
  if(any(is.na(mat))) { #Handle missing values
    uselesssamples <- apply(mat,2,function(x)any(is.na(x)))
    # Note: for future work, also group by by is.na on each row of mat, to use some of the partially-missing data better.
    # uselesssamples <- apply(mat,2,function(x)sum(!is.na(x)))<2
    # etc
    mat <- mat[,!uselesssamples,drop=FALSE]
    groups <- groups[!uselesssamples,drop=TRUE]
  }
  
  Nsamplist <- tapply(rep(1,length(groups)),groups,sum)
  nmutlist <- by(t(mat),INDICES=groups,FUN=colSums,na.rm=TRUE)
  #print(nmutlist)
  weights <- mapply(meta.weight,N=Nsamplist,n=nmutlist)
  # Compute individual pvalues
  if(any(weights>0)) {
    if(avoidones) {
      # Avoid using groups where the probability of getting a pvalue of 1
      # is so large that it "outweighs" the contribution from that group.
      # In practice we check that the normalised square weight is larger then the probability
      # of getting a pvalue of 1 under the null, which is computed using a binomial approximation.
      # Note that this is just a heuristic.
      
      selection <- weights^2*(sapply(nmutlist,max)/Nsamplist)^-sapply(nmutlist,function(x)sum(x)-max(x))
      ord <- order(selection,decreasing=TRUE)
      cumweight <- cumsum(weights[ord]^2)
      cumweight[ord] <- cumweight
      selection <- selection>cumweight
    } else {
      selection <- weights>0
    }
    # Note that all the fancy weighting happens before accessing any information about the amount of overlap.
    pvec <- by(t(mat),INDICES=groups,FUN=function(x)coverage.test(t(x),randomised=randomised))
    # Do weighted Stouffer scoring. We ignore weights of 0 since they produce NaNs.
    p <- pnorm(sum(qnorm(pvec[selection],lower=FALSE)*weights[selection]/sqrt(sum(weights[selection]^2))),lower=FALSE)
    return(p)
  } else {
    warning("No useable data.")
    return(1)
  }
}

representative.clusters <- function(clusters,mat=NULL,n=Inf,max.p.value=0.05,pvals=NULL,groups=NULL,kmax=10,maxoverlap=Inf) {
  # This function selects a sublist of clusters from the input clusters.
  # Choose only the clusters with p-values (pvals) less than max.p.value.
  # If input pvals is missing, pvalues will be computed from mat (and groups)
  # If more than n clusters are obtained limit the output to n clusters.
  # Clusters outputted will be sorted from most to least significant.
  # maxoverlap controls the number of genes in a new cluster which is allowed to overlap with a previous cluster
  # Setting maxoverlap=0 will output only disjoint clusters.
  # If a cluster already has a superset/subset in the output (with a smaller p-value), then it will not be added to the output.
  # OUTPUT
  # A list containing
  # clusters: the selected clusters (ordered from smallest to largest pvalue)
  # selection: the indices of the selected clusters in the input list clusters
  # pvals: the pvalues of the selected clusters (taken from input pvals)
  lengths<-sapply(clusters,length)
  if(is.null(pvals)) {
    coverageall<-sapply(clusters,function(x)meta.coverage.test(mat[x,],groups=groups))
    pvals<-coverageall*sapply(lengths,subset.bonferroni,n=nrow(mat),order="size",kmax=kmax)
  }
  if(length(clusters) != length(pvals)) {stop("pvals and clusters should have the same length")}
  if(any(is.na(pvals))) {
    warning("pvals ",paste(which(pvals),collapse=", "),"are NA.")
    pvals[is.na(pvals)] <- Inf
  }
  ind <- seq(along.with=pvals)
  ind <- ind[pvals<=max.p.value]
  clusters <- clusters[pvals<=max.p.value]
  lengths <- lengths[pvals<=max.p.value]
  pvals <- pvals[pvals<=max.p.value]
  
  if(length(clusters)<1) {
    selection <- integer(0)
    return(list(clusters=clusters[selection],selection=selection,pvals=pvals[selection])) 
  }
  ord <- order(pvals)
  clusters <- clusters[ord]
  pvals <- pvals[ord]
  lengths <- lengths[ord]
  ind <- ind[ord]
  
  clustersmat<-clusters2mat(clusters,genes=unique(unlist(clusters)))
  overlap<-t(clustersmat) %*% clustersmat
  
  # Find whether there any subset of a gene set was already included as a more significant gene set
  # Also remove a set if it is a subset of a previously included gene set
  overlap1 <- overlap
  overlap1[lower.tri(overlap,diag=TRUE)] <- 0
  overlap2 <- t(overlap1)
  overlap1 <- overlap1/lengths
  overlap1 <- overlap1 == 1
  overlap2 <- overlap2/lengths
  overlap2 <- overlap2 == 1
  selection <- sort(which(!(apply(overlap1,2,any)|apply(overlap2,1,any))))
  if(maxoverlap<max(lengths)) { # No point in checking for overlaps if this is false.
    k <- 2
    while(k<=length(selection)) {
      if(any(overlap[selection[k],selection[1:(k-1)]]>maxoverlap)) {
        selection <- selection[-k]
      } else {
        k <- k+1
      }
    }
  }
  if(length(selection)>n) {selection <- selection[1:n]}
  return(list(clusters=clusters[selection],selection=ind[selection],pvals=pvals[selection]))  
}


subset.bonferroni <- function(n,k,nmax=n,kmin=2,kmax=k,order=c("size","stratifiedsize","genesize","simple"),alpha=0.05) {
  # We will assume that order="size" (default)
  # computes the weighted Bonferroni correction for a gene set of size k
  # n|nmax is the number of genes to choose from
  # kmin, kmax are the smallest (resp. largest) number of genes allowed in a gene set
  # the weight is chosen such that the probability of a passenger mutation being added
  # to an alteration set is less than alpha
  #(explanation of the order types:
  # simple: unweighted bonferroni correction
  # genesize: order the hypotheses by the coverage of the least mutated gene
  # stratified size: perform a Bonferroni correction separately for each set size, and then for the allowed range of gene set sizes.)
  order <- match.arg(order)
  if(order=="size") {
    #\frac{\sum_{\ell=2}^{\kmax}{m \choose \ell}\prod_{k=2}^{\ell}\left(1-(1-\alpha)^\frac{1}{m-k+1}\right)}{\prod_{k=2}^{|M|}\left(1-(1-\alpha)^\frac{1}{m-k+1}\right)}
    return(sum(choose(nmax,kmin:kmax)*cumprod(-expm1(log1p(-alpha)/(nmax+1-kmin:kmax))))/prod(-expm1(log1p(-alpha)/(nmax+1-kmin:k))))
  }
#   if(order=="size") {
#     return(choose(nmax,k)*factorial(k)/sum(((2*alpha)^(kmin:kmax)/factorial(kmin:kmax)))*(2*alpha)^(-k))
#   }
  if(order=="stratifiedsize") {
    return(choose(nmax,k)*(kmax-kmin+1))
    #return(sum(choose(nmax,kmin:k)*sum(choose(nmax,kmin:kmax)/cumsum(choose(nmax,kmin:kmax)))))
  }
  if(order=="genesize") {
    return(sum(choose(n,kmin:k))*(log(sum(choose(nmax,kmin:kmax)))+1))
  }
  if(order=="simple") {
    return(sum(choose(nmax,kmin:kmax)))
  }
}

suffixrm <- function(x,part=1,split="_") {
  # By default, assume x is a vector or (recursive) list containing character strings of the form
  # "prefix_suffix"
  # and return an object of the same length/shape, where only
  # "prefix"
  # remains.
  # Setting part=2 would keep the suffix instead of the prefix.
  y <- unlist(x)
  y <- strsplit(y,split=split,fixed=TRUE)
  y <- sapply(y,function(x)x[part])
  return(relist(y,skeleton=x))
}

clusters2mat <- function(clusters,genes=NULL,genenames=NULL) {
  # Convert a list of clusters to a binary matrix (genes in rows, clusters in columns)
  if(is.null(genes)) genes<-unique(unlist(clusters))
  if(is.null(genenames)) genenames<-genes
  genesvsclusters<-sapply(clusters,function(x,genes) genes %in% x,genes=genes)
  rownames(genesvsclusters)<-genenames
  colnames(genesvsclusters)<-names(clusters)
  return(genesvsclusters)
}

maximal.clusters <- function(merge,onlymaximal=TRUE) {
  # Generate the clusters based on the merges described in merge (see hclust)
  # OUTPUT:
  # A list containing:
  # clusters: a list of clusters (only the maximal clusters are retained if onlymaximal=TRUE (default))
  # ismaximal: a logical vector indicating which of the merges generate maximal clusters
  names(merge)<-NULL
  clusters<-list()
  for(k in 1:nrow(merge)){
    tmp<-sort(merge[k,])
    if(all(tmp<0)) {
      new.cluster <- -tmp
    } else {
      if(tmp[1]<0) {
        new.cluster <- c(-tmp[1],clusters[[tmp[2]]])
      } else {
        new.cluster <- c(clusters[[tmp[1]]],clusters[[tmp[2]]])
      }      
    }
    clusters<-c(clusters,list(new.cluster))
  }
  ismaximal1 <- !((1:nrow(merge)) %in% merge[merge>0])
  if(!onlymaximal) allclusters<-clusters
  clusters <- clusters[ismaximal1]
  ismaximal2 <- !sapply(clusters,function(x,clusters) sum(sapply(clusters,function(cc,x)all(x %in% cc),x=x))>1,clusters=clusters)
  clusters <- clusters[ismaximal2]
  ismaximal1[ismaximal1]<-ismaximal2
  if(!onlymaximal) clusters<-allclusters
  return(list(clusters=clusters,ismaximal=ismaximal1))
}

subsets <- function(set,min.size=1,max.size=Inf) {
  # List all subsets of a set (vector)
  # subsets of sizes between min.size and max.size only are considered.
  max.size <- min(max.size,length(set))
  if(min.size>max.size) {return(list())}
  out <- lapply(min.size:max.size,function(x)combn(set,x,simplify=FALSE))
  out <- unlist(out,recursive=FALSE)
  return(out)
}

all.subclusters <- function(clusters,min.size=2,max.size=Inf) {
  # Given a list of sets clusters,
  # produce a list of all subsets of the sets.
  # The size of subsets must be between min.size and max.size
  out <- lapply(clusters,subsets,min.size=min.size,max.size=Inf)
  out <- unlist(out,recursive=FALSE)
  out <- out[!duplicated(out)]
  return(out)
}

min.possible.coverage.test <- function(nsample,k,correct=FALSE,nmax=k,alpha=0.05) {
  # Returns the minimum possible p-value from coverage.test when there are
  # nsample: the number of samples
  # k: the number of genes in the gene set
  n <- nsample
  mat <- matrix(FALSE,nrow=k,ncol=n)
  g <- rep(floor(n/k),k)+c(rep(1,n%%k),rep(0,k-(n%%k)))
  cumg <- c(0L,cumsum(g))
  for(i in 1:k)  {
    mat[i,cumg[i]:cumg[i+1]] <- TRUE
  }
  out <- coverage.test(mat)
  if(correct) {out <- out*subset.bonferroni(n=nmax,kmax=k,k=k,alpha=alpha)}
  return(out)
}
get.kmax.upperbound <- function(nsample,nmax,alpha=0.05,correct=TRUE) {
  for(k in 2:nmax) {
    if(min.possible.coverage.test(nsample=nsample,k=k,correct=correct,alpha=alpha,nmax=nmax)>alpha) {
      return(k-1)
    }  
  }
  return(nmax)
}
