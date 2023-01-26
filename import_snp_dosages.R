########################################
## function to obtain LD information
## by importing SNPs

get.ld.info <- function(olink, chr.s, pos.s, pos.e){
  
  ## 'olink' -- identifier for the protein
  ## 'chr.s' -- chromosome
  ## 'pos.s' -- start position
  ## 'pos.e' -- end position
  
  ## run as as a bash job
  system(paste("./scripts/get_SNP_dosages.sh", olink, chr.s, pos.s, pos.e))
  
  ## read the dosage file
  require(data.table)
  tmp           <- fread(paste("tmpdir/tmp", olink, chr.s, pos.s, pos.e, "dosage", sep="."), sep=" ", header=T, data.table=F)
  ## transpose
  rownames(tmp) <- tmp$rsid
  ## store allele information to realign effect estimates afterwards
  tmp.info      <- tmp[, 1:6]
  tmp           <- t(tmp[,-c(1:6)])
  ## retransform to data set (keep X in mind for variable names)
  tmp           <- data.frame(ID_1=rownames(tmp), tmp)
  
  ## create another column to info to map names from the SNP data set
  tmp.info$id         <- sapply(tmp.info$rsid, function(x) ifelse(substr(x, 1, 2) == "rs", x, paste0("X", gsub(":", ".", x))))
  ## edit some IDs (X-chromosome)
  tmp.info$id         <- gsub("XX", "X", tmp.info$id)
  tmp.info$id         <- gsub("XAffx-", "Affx.", tmp.info$id)
  tmp.info$id         <- gsub(",", ".", tmp.info$id)
  ## set to those included in the data set, just to be sure
  tmp.info            <- subset(tmp.info, id %in% names(tmp))
  ## create markername to ease merging with GWAS results
  tmp.info$snp.id     <- apply(tmp.info, 1, function(x){
    paste0("chr", as.numeric(x[1]), ":", as.numeric(x[4]), "_", paste(toupper(sort(x[5:6])), collapse = "_"))
  })
  
  
  return(list(tmp, tmp.info))
  
}
