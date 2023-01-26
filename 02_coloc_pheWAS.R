rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

##options
options(stringsAsFactors = F)


## --> import parameters <-- ##

#olink = name of the Olink protein
#chr.s = chromosome of the cis-region 
#pos.s = start position of the cis-region
#pos.e = end position of the cis-region

olink <- args[1]
chr.s <- as.numeric(args[2])
pos.s <- as.numeric(args[3])
pos.e <- as.numeric(args[4])


#------------------------------------------------#
##--           import protein GWAS data       --##
#------------------------------------------------#

## read the relevant data
res.olink        <- paste0("zcat ~/insert_path/",
                           olink,"/", olink,"GWAS.txt.gz",
                           " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e,
                           " '{if(($2 == chr && $3 >= low && $3 <= upp) || NR==1) print $0}' -")
## import
res.olink        <- data.table::fread(cmd = res.olink, sep = "\t", header = T, data.table = F)


## create MAF variable
res.olink$MAF    <- ifelse(res.olink$Freq1 > .5, 1 - res.olink$Freq1, res.olink$Freq1)


## create new MarkerName 
res.olink$snp.id <- apply(res.olink, 1, function(x){
  paste0("chr", as.numeric(x[2]), ":", as.numeric(x[3]), "_", paste(toupper(sort(x[4:5])), collapse = "_"))
})

## define break criteria here (min p-value!)
p.min            <- min(as.numeric(res.olink$Pvalue), na.rm = T)

#-----------------------------------------#
##--   only do if suggestive evidence  --##
#-----------------------------------------#

if(p.min < 1e-6){

  #-------------------------------------------------------------#
  ##--       import the indivudal level genotype data        --##
  #-------------------------------------------------------------#

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

  #-----------------------------------------#
  ##--  define top SNP and get proxies   --##
  #-----------------------------------------#

  ## identify top SNP by Z-score
  top.snp    <- res.olink$snp.id[which.max(abs(res.olink$Effect/res.olink$StdErr))]

  ## get proxies (use snp.id)
  proxy.snps <- cor(snp.dat[, snp.info$id[which(snp.info$snp.id == top.snp)]], snp.dat[, snp.info$id])^2
  proxy.snps <- data.frame(id=colnames(proxy.snps), R2=t(proxy.snps))
  proxy.snps <- subset(proxy.snps, R2 >= .8)
  ## ease merging
  proxy.snps <- merge(proxy.snps, snp.info)
  ## ensure overlap with summary stats
  gc()

  #-----------------------------------------#
  ##--  perform PheWAS for all proxies   --##
  #-----------------------------------------#

  ## load package to for API to https://gwas.mrcieu.ac.uk/
  require(ieugwasr)

  ## query everything (omit FinnGen and BBJ as sample sizes are missing)
  tmp.phewas   <- phewas(variants=proxy.snps$rsid, pval=1e-6, batch=c("ebi-a", "ieu-a", "ieu-b", "met-a", "ubm-a", "ukb-b"))

  ## drop phenotypes with missing information on sample size
  tmp.phewas$n <- as.numeric(tmp.phewas$n)
  tmp.phewas   <- subset(tmp.phewas, !is.na(n))

  ## define break criteria here as well
  if(nrow(tmp.phewas)>0){

    ## to ease some downstream analysis
    require(data.table)

    ## take the strongest association for each phenotype
    tmp     <- as.data.table(tmp.phewas[order(tmp.phewas$trait, abs(tmp.phewas$beta/tmp.phewas$se), decreasing = T),])
    tmp[, ind:=1:.N, by="trait"]
    tmp     <- subset(tmp, ind == 1)
    tmp$ind <- NULL

    ## test whether any remaining phenotypes may exist
    if(nrow(tmp)>0){
      tmp2 <- subset(tmp.phewas, rsid %in% proxy.snps$rsid & !(trait %in% tmp$trait))
    }else{
      tmp2 <- subset(tmp.phewas, rsid %in% proxy.snps$rsid)
    }

    ## do iterative subsetting to keep unqie entries only
    if(nrow(tmp2)>1){
      ## take the strongest association for each phenotype
      tmp2     <- as.data.table(tmp2[order(tmp2$trait, abs(tmp2$beta/tmp2$se), decreasing = T),])
      tmp2[, ind:=1:.N, by="trait"]
      tmp2     <- subset(tmp2, ind == 1)
      tmp2$ind <- NULL
    }

    ## combine both
    tmp                 <- rbind(tmp, tmp2)
    ## add SNP information
    tmp                 <- merge(proxy.snps[, c("rsid", "id", "snp.id", "R2")], tmp, by="rsid", suffixes = c(".snp", ".ieu"))
    ## rename and add lead SNP
    tmp$snp.id.lead     <- top.snp
    tmp$rsid.lead       <- snp.info$rsid[which(snp.info$snp.id == top.snp)]

    ## add protein stats
    phewas.results      <- merge(res.olink, tmp, by=c("snp.id", "chr"), suffix=c(".pQTL", ".trait"))

    #-----------------------------------------#
    ##--    run coloc for all results      --##
    #-----------------------------------------#

    require(coloc)

    res.coloc <- lapply(1:nrow(phewas.results), function(x){

      ## obtain summary stats and trait info
      reg             <- paste0(chr.s, ":", pos.s,"-", pos.e)
      res.trait       <- associations(reg, phewas.results$id.ieu[x])
      ## just to make sure everything will work fine
      res.trait       <- subset(res.trait, !is.na(beta))
      ## make unique, since some data processing error
      res.trait       <- unique(res.trait)

      ## get meta information (careful some might miss those)
      tr.info         <- tryCatch(
        {
          gwasinfo(phewas.results$id.ieu[x])
        }, error=function(e){
          return(NA)
        })

      ## merge (this might ignore INDELS)
      res.all                   <- merge(res.olink, res.trait, by="rsid")
      ## remove non-biallelelic variants
      ii                        <- table(res.all$rsid)
      res.all                   <- subset(res.all, rsid %in% names(ii[ii==1]))
      ## account for possible INDELs
      res.all[, c("ea", "nea")] <- t(apply(res.all[, c("ea", "nea")], 1, function(x){

        if(nchar(x[1]) > 1 | nchar(x[2]) > 1){
          ## replace
          if(nchar(x[1]) > nchar(x[2])){
            return(c("I", "D"))
          }else{
            return(c("D", "I"))
          }
        }else{
          return(x)
        }

      }))


      ## align effect estimates
      res.all$beta.aligned <- ifelse(toupper(res.all$Allele1) == res.all$ea, res.all$beta, -res.all$beta)

      #-----------------------------------------#
      ##--               sanity check            --##
      #-----------------------------------------#

      ## top signal for Olink in overlap
      it    <- res.all$snp.id[which.max(abs(res.all$Effect/res.all$StdErr))]
      ## get the top SNP for the outcome in overlap only
      io    <- res.all$snp.id[which.max(abs(res.all$beta/res.all$se))]

      ## keep names
      isnps <- sapply(c(top.snp, it, io), function(x) snp.info$id[which(snp.info$snp.id == x)])

      ## LD for check
      it    <- cor(snp.dat[, isnps[1]], snp.dat[, isnps[2]])^2
      ## LD within
      is    <- cor(snp.dat[, isnps[2]], snp.dat[, isnps[3]])^2
      ## LD top
      io    <- cor(snp.dat[, isnps[1]], snp.dat[, isnps[3]])^2

      #-----------------------------------------#
      ##--                  run coloc            --##
      #-----------------------------------------#

      ## add to ease mapping of LD matrix
      res.all     <- merge(res.all, snp.info[, c("snp.id", "id")], by="snp.id", suffix=c(".trait", ".snp"))

      ## order by position
      res.all     <- res.all[order(res.all$pos),]

      ## prepare input
      D1          <- list(beta=res.all$Effect, varbeta=res.all$StdErr^2,
                          type="quant",
                          N=max(res.all$TotalSampleSize, na.rm=T),
                          sdY=1,
                          MAF=res.all$MAF,
                          snp=res.all$id.snp,
                          position=1:nrow(res.all))

      ## try out
      if(!("ncase" %in% names(tr.info))){
        print("use quantitative")
        D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2,
                            type="quant",
                            # N=tr.info$sample_size,
                            N=max(res.all$n, na.rm=T),
                            sdY=1,
                            MAF=res.all$MAF,
                            snp=res.all$id.snp,
                            position=1:nrow(res.all))
      }else{
        ## binary outcome
        D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2,
                            type="cc",
                            s=tr.info$ncase/(tr.info$ncontrol+tr.info$ncase),
                            N=tr.info$sample_size,
                            MAF=res.all$MAF,
                            snp=res.all$id.snp,
                            position=1:nrow(res.all))
      }

      ## do naive coloc as well
      naive.coloc                              <- coloc.signals(D1, D2, method="single", p12=5e-6)

      ## add checks to the data
      naive.coloc$summary$ld.preserved.overlap <- it
      naive.coloc$summary$ld.check.top.overlap <- io
      naive.coloc$summary$ld.check.top.general <- is

      ## add the trait id
      naive.coloc$summary$id.ieu               <- phewas.results$id.ieu[x]

      #-----------------------------------------#
      ##--              draw selected            --##
      #-----------------------------------------#

      if(naive.coloc$summary$PP.H4.abf > .7 | naive.coloc$summary$ld.check.top.overlap > .8){
        source("scripts/plot_coloc_results.R")
        png(paste0("graphics/", olink, ".", phewas.results$id.ieu[x], ".", chr.s, ".", pos.s, ".", pos.e, ".png"), width=16, height=8, units="cm", res=150)
        par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, cex.main=.6, font.main=2)
        ## more complex layout for gene assignment
        layout(matrix(c(1,1,1,2,3,4),3,2), heights = c(.43,.37,.2))
        plot.regional.coloc(res.all, naive.coloc$summary, snp.dat, phewas.results$id.ieu[x], olink, snp.info)
        dev.off()
      }


      ## write results to file
      return(naive.coloc$summary)

    })
    res.coloc <- do.call(rbind, res.coloc)

    #-----------------------------------------#
    ##--        combine and store the data     --##
    #-----------------------------------------#

    ## combine both
    phewas.results       <- merge(phewas.results, res.coloc[, c("id.ieu", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
                                                                "ld.preserved.overlap", "ld.check.top.overlap", "ld.check.top.general")])

    ## add phenotype (be precise in region tested to account for proteins encoded by multiple genes)
    phewas.results$pheno <- paste(olink, chr.s, pos.s, pos.e, sep=".")

    ## write to file
    write.table(phewas.results, paste("output/results.3k.phewas", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)


  }else{
    cat("found no evidence for", olink, "on chr", chr.s, "between", pos.s, "and", pos.e, "\n")
    ## write to file
    write.table(data.frame(pheno=paste(olink, chr.s, pos.s, pos.e, sep="."), id.ieu=NA),
                paste("output/results.3k.phewas", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
  }

}else{
  cat("found no evidence for", olink, "on chr", chr.s, "between", pos.s, "and", pos.e, "\n")
  ## write to file
  write.table(data.frame(pheno=paste(olink, chr.s, pos.s, pos.e, sep="."), id.ieu=NA),
              paste("output/results.3k.phewas", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
}
