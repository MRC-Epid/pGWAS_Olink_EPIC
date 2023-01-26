rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## options
options(stringsAsFactors = F)

## --> packages required <-- ##
require(data.table)
require(susieR)
require(doMC)

## --> import parameters <-- ##
#olink = name of the Olink protein
#id = Olink protein id
#chr.s = chromosome of the cis-region 
#pos.s = start position of the cis-region
#pos.e = end position of the cis-region

olink <- args[1]
id    <- args[2]
chr.s <- args[3]
pos.s <- as.numeric(args[4])
pos.e <- as.numeric(args[5])

cat("Run Fine-Mapping using SuSiE with", olink, id, chr.s, pos.s, pos.e, "\n")

#-----------------------------------------#
##--     load regional assoc stats     --##
#-----------------------------------------#

cat("------------------------------\n")
cat("Reading summary statistics in \n")

## read the relevant GWAS data for the protein
res.olink        <- paste0("zcat ~/intert_path/",
                           olink,"/", olink,"GWAS.txt.gz",
                           " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e,
                           " '{if(($2 == chr && $3 >= low && $3 <= upp) || NR==1) print $0}' -")
## import
res.olink        <- data.table::fread(cmd = res.olink, sep = "\t", header = T, data.table = F)

## Create a MAF variable
res.olink$MAF    <- ifelse(res.olink$Freq1 > .5, 1 - res.olink$Freq1, res.olink$Freq1)

## create a common MarkerName 
res.olink$snp.id <- apply(res.olink, 1, function(x){
  paste0("chr", as.numeric(x[2]), ":", as.numeric(x[3]), "_", paste(toupper(sort(x[4:5])), collapse = "_"))
})

cat("Found", nrow(res.olink), "entries \n")
cat("------------------------------\n")

#--------------------------------------------------------------#
##--        load the individual level phenotype data        --##
#--------------------------------------------------------------#

cat("------------------------------\n")
cat("Reading phenotype variables in\n")

## import trait data
pheno      <- fread("EPIC.Olink.3K.txt", sep="\t", header=T, select = c("omicsid", id), data.table = F)

cat("Found", nrow(pheno), "samples \n")
cat("------------------------------\n")

#-----------------------------------------------------------------------------#
##--       import the individual level genotype data and dosage info       --##
#-----------------------------------------------------------------------------#

cat("------------------------------\n")
cat("Reading snp data in\n")

## write file to obtain SNP dosages
write.table(res.olink$MarkerName, paste("tmpdir/tmp", olink, chr.s, pos.s, pos.e, "lst", sep="."), row.names = F, col.names = F, quote = F)

## import function to do so
source("scripts/import_snp_dosages.R")
## import
ld        <- get.ld.info(olink, chr.s, pos.s, pos.e)
## ease downstream coding
snp.dat   <- ld[[1]]
snp.info  <- ld[[2]]
## delete and clean
rm(ld); gc(); gc()

## delete file no longer needed
system(paste("rm tmpdir/tmp", olink, chr.s, pos.s, pos.e, "dosage", sep="."))

## add rsid to Olink file
res.olink <- merge(res.olink, snp.info[, c("snp.id", "rsid")])

## SNPs of interest
snps      <- snp.info$id[snp.info$snp.id %in% res.olink$snp.id]

cat("Found", length(snps), "SNPs \n")
cat("------------------------------\n")

#-----------------------------------------#
##--       combine the data sets       --##
#-----------------------------------------#

## combine phenotype with SNPs
pheno <- merge(pheno, snp.dat, by.x="omicsid", by.y="ID_1")

## free some space
rm(snp.dat); gc()

#-----------------------------------------#
##--              run SuSiE            --##
#-----------------------------------------#

cat("------------------------------\n")
cat("Running SuSiE with L 2 - 10\n")

## set seed for reproducibility
set.seed(43)
registerDoMC(6)

## run across multiple cores

## --> run with an increasing number of potential causal variants <-- ##

res.susie <- mclapply(2:10, function(l){

  cat("Using L =", l, "\n")

  ## run with given L
  fitted   <- susie(as.matrix(pheno[,snps]), pheno[,id], L = l)

  ## get credible sets
  cred.set <- susie_get_cs(fitted, X = as.matrix(pheno[,snps]))

  ## get the elbow fit value
  ev       <- susie_get_objective(fitted)

  #-----------------------------------------#
  ##--        prepare final results      --##
  #-----------------------------------------#

  ## create data sets for credible sets
  if(length(cred.set$cs) > 0){
    tmp <- lapply(1:length(cred.set$cs), function(x){
      ## get all identifier
      jj <- cred.set$cs[[x]]
      ## map back to snps
      ii <- snps[jj]
      ## get pips
      pp <- fitted$pip[jj]
      ## identify top signal
      ts <- which.max(pp)
      ## identify in the GWAS results
      ts <- snp.info$snp.id[which(snp.info$id == ii[ts])]
      ## store information
      return(data.frame(cset=x, snp.id=ts, c.length=length(jj), pip.top=max(pp, na.rm=T), cred.set$purity[x,,drop=F],
                        c.vars=paste(snp.info$snp.id[which(snp.info$id %in% ii)], collapse = ", ")))
    })
    tmp      <- do.call(rbind, tmp)
    ## add information from the GWAS
    tmp      <- merge(tmp, res.olink, all.x=T)

    ## get the data needed
    ts       <- res.olink[which.max(res.olink$Effect/res.olink$StdErr),"snp.id"]
    as       <- sapply(c(ts, tmp$snp.id), function(x) snp.info$id[which(snp.info$snp.id == x)])
    as       <- data.frame(snp.id=names(as[-1]), ld.lead.ma=cor(pheno[, as[1]], pheno[, as[-1]])^2)

    ## add to the output
    tmp      <- merge(tmp, as)

    ## add elbow value
    tmp$elbo <- ev
    tmp$L    <- l

    #-----------------------------------------#
    ##--          run joint model          --##
    #-----------------------------------------#

    ## define formula to run a joint model of the top variant (pip) from each credible set
    tmp.snps  <- unique(tmp$snp.id)
    tmp.snps  <- snp.info$id[which(snp.info$snp.id %in% tmp.snps)]
    ## generate the formula
    ff        <- paste0(id, " ~ ", paste(tmp.snps, collapse = " + "))
    ff        <- data.frame(summary(lm(ff, data=pheno))$coefficients)
    ## add id to merge with the output from SuSiE
    ff$id     <- rownames(ff)
    ## rename
    names(ff) <- c("Effect.joint", "StdErr.joint", "tval.joint", "Pvalue.joint", "id")
    ff        <- merge(ff, snp.info[, c("id", "snp.id", "alleleA", "alleleB")])

    ## add to the data (drops also intercept from the linear model)
    tmp       <- merge(tmp, ff, all.x=T, by="snp.id")
    return(tmp)

  }else{
    cat("no signals found\n")
    return(NA)
  }
}, mc.cores=6)

cat("Done \n")
cat("------------------------------\n")

#-----------------------------------------#
##--  delete NA entries from the list  --##
#-----------------------------------------#

## drop empty entries
res.susie <- Filter(function(x) !all(is.na(x)), res.susie)

## proceed only if any left
if(length(res.susie) > 0){

  #-----------------------------------------#
  ##--  go through the list and select   --##
  #-----------------------------------------#

  cat("------------------------------\n")
  cat("Select L for final run \n")

  ## go through and count number of csets
  c.num  <- unlist(lapply(res.susie, nrow))
  ## go through and count how often top varaints are in LD (R2>.1) with lead variant
  c.lead <- unlist(lapply(res.susie, function(x) sum(x$ld.lead.ma > .1)))

  ## choose lowest if c.lead > 1
  if(length(c.lead)  > 2){
    c.lead <- c.lead[which(c.lead < 3)]
    c.num  <- c.num[1:length(c.lead)]
  }

  ## choose appropriate
  if(sum(!is.na(c.num))>0 & length(c.num) > 1){
    ## go through and select either the largest value or if two are the same for increasing L
    ii <- c.num[1]
    jj <- c.num[2]
    if(ii == jj){
      tmp <- res.susie[[ii]]
    }else{
      k <- 2
      while(ii < jj & (k < length(c.num)+1)){
        k  <- k+1
        ii <- c.num[k-1]
        jj <- c.num[k]
      }
      # print(c.num[ii-1])
      tmp <- res.susie[[ii-1]]
    }
  }else{
    tmp    <- NA
  }

  cat("Found\n")
  print(tmp)
  cat("------------------------------\n")

  #-----------------------------------------#
  ##--  fit a last time to obtain PIPs   --##
  #-----------------------------------------#

  cat("------------------------------\n")
  cat("Running very last round for SuSiE\n")

  if(!is.null(nrow(tmp))){

    ## fit the model
    cat("Fit final model with L =", tmp$L[1], "\n")

    ## set seed for reproducibility
    set.seed(43)

    ## run with given L
    fitted       <- susie(as.matrix(pheno[,snps]), pheno[,id], L = tmp$L[1])

    ## add PIP to results (may drop some rare SNPs with no LD information)
    tmp          <- data.frame(id=snps, pip=fitted$pip)
    res.olink    <- merge(res.olink, snp.info[, c("snp.id", "id")])
    res.olink    <- merge(res.olink, tmp, all.x = T)

    ## add credible set information
    tmp          <- fitted$sets$cs
    ## use LD matrix to map names back, since ordering in the results file has changed
    tmp          <- lapply(tmp, function(x) snps[x])

    #------------------------------------------------#
    ##-- run joint model to ensure GWAS threshold --##
    #------------------------------------------------#

    ## identify top SNP for each credible set
    top.snp         <- lapply(tmp, function(x){
      ## ensure that exact mapping is preserved
      foo <- res.olink[which(res.olink$id %in% x), ]
      return(foo$id[which.max(foo$pip)])
    })

    ## run model
    m.joint         <- summary(lm(paste(id, "~", paste(unlist(top.snp), collapse = " + ")), pheno))$coefficients[-1, ,drop=F]
    ## get indicator whether threshold is met
    m.joint         <- ifelse(m.joint[,4] < 5e-8, T, F)
    ## subset credible sets
    tmp             <- tmp[m.joint]

    ## only proceed it al least one remaining
    if(length(tmp) > 0){
      ## very last round
      top.snp      <- lapply(tmp, function(x){
        ## ensure that exact mapping is preserved
        foo <- res.olink[which(res.olink$id %in% x), ]
        return(foo$id[which.max(foo$pip)])
      })

      ## add info to the results file
      res.olink$cs   <- NA

      ## loop through all credible sets
      for(j in 1:length(tmp)){

        ## set all variants to j if included in the current credible set
        res.olink$cs[which(res.olink$id %in% tmp[[j]])] <- j

        ## add LD with top variant in the set
        ii                            <- res.olink[which(res.olink$cs == j),]
        ii                            <- ii$id[which.max(ii$pip)]
        ## add LD
        res.olink[, paste0("R2.", j)] <- t(cor(pheno[,ii], pheno[, res.olink$id])^2)
      }

      ## store final results
      write.table(res.olink, paste("output/fine.mapped", olink, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)

      ## run final joint model
      m.joint        <- data.frame(summary(lm(paste(id, "~", paste(unlist(top.snp), collapse = " + ")), pheno))$coefficients[-1, ,drop=F])
      ## keep what is needed
      names(m.joint) <- c("Effect.joint", "StdErr.joint", "tval.joint", "Pvalue.joint")
      ## add id
      m.joint$id     <- rownames(m.joint)
      ## add credible set identifier
      m.joint$cs     <- 1:nrow(m.joint)
      ## add allele coding
      m.joint        <- merge(snp.info[, c("id", "alleleA", "alleleB")], m.joint)
      ## combine with marginal stats
      m.joint        <- merge(res.olink, m.joint)

      ## store results
      write.table(m.joint, paste("joint_models/fine.mapped", olink, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)

    }else{
      ## place empty file
      write.table(data.frame(snp.id=NA), paste("output/fine.mapped", olink, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)
    }

  }else{
    ## place empty file
    write.table(data.frame(snp.id=NA), paste("output/fine.mapped", olink, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)
  }

}else{

  ## place empty file
  write.table(data.frame(snp.id=NA), paste("output/fine.mapped", olink, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)

}


cat("Done!")