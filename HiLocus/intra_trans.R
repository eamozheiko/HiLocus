## functions
chrom_pair <- function(path_to_hic, chri, binsize, input_format){
  chrj = chri
  ## read data
  if(input_format == ".hic"){
    # read hic format file
    bed <- strawr::straw("NONE", path_to_hic, as.character(chri), as.character(chri), "BP", binsize)
    names(bed) = c("posi", "posj", "counts")
  } else {
    # read pairs format file 
    command = paste("cat ", path_to_hic, " | sed \'s/chr//g\' | ", "awk", "\ '{ if($1 == chri && $3 == chrj){ print $0 } }\' ", "chri=", chri, " chrj=", chrj,   sep="")
    bed = read.table(pipe(command), stringsAsFactors = F, head=F)
    
    # rm chr cols
    bed = bed[, -c(1, 3)]
    if(length(bed) == 3){
      names(bed) = c("posi", "posj", "counts")
    } else {
      names(bed) = c("posi", "posj")
    }
  }
  
  ## rm de
  w = which(abs(bed$posi - bed$posj) >= 10^4)
  bed = bed[w, ]
  
  
  
  return(bed)
}

bed_to_sparse <- function(bed, binsize, binsize_frame, nrows, ncols){
  if(ncol(bed) == 2){
    x = rep(1, nrow(bed))
  } else {
    x = bed$counts
  }

  i = c(bed$posi, bed$posj)%/%binsize + 1
  j = c(bed$posj, bed$posi)%/%binsize_frame + 1
  x = c(x, x)
  
  sparse.hic = sparseMatrix(i = i, j = j, x = x, dims = c(nrows, ncols))
  
  return(sparse.hic)
}

normalization <- function(m){
  # by cols
  m = t(t(m)/colSums(m))
  
  # by rows
  m =  m/rowSums(m, na.rm=T)
  
  return(m)
}

filter_intra_trans <- function(vals_res, sd_cis, mean_cis, probability_threshold){
  N = nrow(vals_res)
  q = as.numeric(quantile(vals_res[, 6], 0:10/10))

  if(N < 100){
    coverage = vals_res[w, 6]
    coverage[which(coverage > 200)] = 200
    
    sdev = sd_cis[coverage]
    m = mean_cis[coverage]
    
    val = vals_res[w, 4]
    pr = pnorm(val, m, sdev, lower.tail = F)
    vals_res[w, 5] = pr
  } else {
    for (i in 1:10) {
      w = which(vals_res[, 6] >= q[i] & vals_res[, 6] <= q[i + 1])
      coverage = vals_res[w, 6]
      
      sdev = sd(vals_res[w, 4])
      m = median(vals_res[w, 4])
      
      val = vals_res[w, 4]
      pr = pnorm(val, m, sdev, lower.tail = F)
      vals_res[w, 5] = pr
    }
  }
  
  
  w = which(vals_res[, 6] > 50 & vals_res[, 7] > 50 & vals_res[, 5] < probability_threshold)
  
  res = data.frame(rbind(vals_res[w, ]))
  names(res) = c("CHR_FROM", "POS_FROM_START", "POS_FROM_END", "VAL", "PROBABILITY", "COVERAGE_CASE_CIS", "COVERAGE_CONTROL_CIS")
  
  return(res)
}




hilocus_intra_trans <- function(args){
  DIR = args[1]
  outdir_path = args[2]
  vp_case_path = args[3]
  vp_control_path = args[4]
  binsize = as.numeric(args[5])
  sample_name = args[6]
  probability_threshold = as.numeric(args[7])
  input_format_case = args[8]
  input_format_control = args[9]
  
  ## output file name
  filename = paste(outdir_path, "/intra.translocations.binsize", binsize, ".", sample_name, ".tsv", sep = "")
  
  ## interchromosomal translocation algorithm
  vals_res0 = matrix(0, nrow = 1, ncol = 7)
  
  for (chr in 1:22) {
    ## read data
    #print(chr)
    gc()
    ######################################################################################################################
    ## read hic intra
    case.bed = try(chrom_pair(vp_case_path, chr, binsize, input_format_case),  silent = TRUE)
    control.bed = try(chrom_pair(vp_control_path, chr, binsize, input_format_control),  silent = TRUE)
    if(inherits(case.bed, "try-error") || inherits(control.bed, "try-error") )
    {
      #error handling code, maybe just skip this iteration using
      next
    }
    m = c(case.bed$posi, case.bed$posj, control.bed$posi, control.bed$posj)
    nrows = max(m%/%binsize + 1)
    ncols = max(m%/%binsize_frame + 1)
    
    case.sparse = bed_to_sparse(case.bed, binsize, binsize_frame, nrows, ncols)
    control.sparse = bed_to_sparse(control.bed, binsize, binsize_frame, nrows, ncols)
    
    ## normalization
    case.sparse.norm = normalization(case.sparse)
    control.sparse.norm  = normalization(control.sparse)
    
    ## result
    vals = rowSums(abs(case.sparse.norm - control.sparse.norm), na.rm = TRUE)
    N = length(vals)
    
    vals_res = matrix(0, nrow = N, ncol = 7)
    vals_res[, 1] = chr
    vals_res[, 2] = seq(from = 0, to = (N - 1)*binsize, by = binsize )
    vals_res[, 3] = seq(from = binsize, to = N*binsize, by = binsize )
    vals_res[, 4] = vals
    vals_res[, 6] = rowSums(case.sparse)[1:N]
    vals_res[, 7] = rowSums(control.sparse)[1:N]
    
    w = which(vals_res[, 4] != Inf & vals_res[, 4] != 0 & !is.na(vals_res[, 4] & vals_res[, 4] != 2))
    vals_res0 = rbind(vals_res0, vals_res[w,])
  }
  vals_res = vals_res0[-1, ]
  
  ## filtration
  # read pre-calculated standard deviation and mean depending on coverage
  sd_cis = read.table(paste(DIR,"/data/sd_cis.tsv", sep = ""), stringsAsFactors = F, head=F)
  sd_cis = sd_cis[, 1]
  
  mean_cis = read.table(paste(DIR,"/data/mean_cis.tsv", sep = ""), stringsAsFactors = F, head=F)
  mean_cis = mean_cis[, 1]
  
  # filter
  intra_translocations = filter_intra_trans(vals_res, sd_cis, mean_cis, probability_threshold)
  
  ## save results
  write.table(intra_translocations, filename, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
}

library(Matrix)
library(strawr)

options(scipen=999)
args = commandArgs(trailingOnly=TRUE)

## default variables
#outdir_path = getwd()
#binsize = 10000
binsize_frame = 5*10^6
#probability_threshold = 10^(-6)



## hilocus_inter_trans algorithm
hilocus_intra_trans(args)

