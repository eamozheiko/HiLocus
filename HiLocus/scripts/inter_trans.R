## functions
chrom_pair <- function(path_to_hic, chri, chrj, binsize, input_format){
  ## read data
  if(input_format == ".hic"){
    # read hic format file
    bed <- strawr::straw("NONE", path_to_hic, as.character(chri), as.character(chrj), "BP", binsize)
    names(bed) = c("posi", "posj", "counts")
  } else {
    # read pairs format file 
    if(chri > chrj){
      command = paste("cat ", path_to_hic, " | sed \'s/chr//g\' | ", "awk", "\ '{ if($1 == chri && $3 == chrj){ print $0 } }\' ", "chri=", chrj, " chrj=", chri,   sep="")
      bed = read.table(pipe(command), stringsAsFactors = F, head=F)
    } else {
      command = paste("cat ", path_to_hic, " | sed \'s/chr//g\' | ", "awk", "\ '{ if($1 == chri && $3 == chrj){ print $0 } }\' ", "chri=", chri, " chrj=", chrj,   sep="")
      bed = read.table(pipe(command), stringsAsFactors = F, head=F)
    }
    
    # rm chr cols
    bed = bed[, -c(1, 3)]
    if(length(bed) == 3){
      names(bed) = c("posi", "posj", "counts")
    } else {
      names(bed) = c("posi", "posj")
    }
  }
  
  ## rm de
  if(chri == chrj){
    w = which(abs(bed$posi - bed$posj) >= 10^6)
    bed = bed[w, ]
  }
  
  ## switch chromosomes
  if(chri > chrj){
    names(bed)[1:2] = c("posj", "posi")
  }
  
  ## index positions
  bed$posi = bed$posi%/%binsize + 1
  bed$posj = bed$posj%/%binsize + 1
  
  return(bed)
}

bed_to_rowsums <- function(bed, nrows, ncols, chr_intra, chr_inter){
  if(ncol(bed) == 2){
    x = rep(1, nrow(bed))
  } else {
    x = bed$counts
  }
  if(chr_intra == chr_inter){
    i = c(bed$posi, bed$posj)
    j = c(bed$posj, bed$posi)
    x = c(x, x)
    nrows = max(nrows, ncols)
    ncols = nrows
  } else {
    i = bed$posi
    j = bed$posj
  }
  sparse.hic = sparseMatrix(i = i, j = j, x = x, dims = c(nrows, ncols))
  
  return(rowSums(sparse.hic))
}

normalization <- function(vals_res, cis_coverage_case, cis_coverage_control, trans_coverage_case, trans_coverage_control){
  ## normalization
  for (chr_intra in 1:22) {
    for (chr_inter in 1:22) {
      if(chr_intra != chr_inter){
        w = which(vals_res[, 1] == chr_intra & vals_res[, 2] == chr_inter)
        if(length(w) != 0){
          ## calculate cis/trans normalization coefficient
          cis_sum_case = sum(cis_coverage_case[-chr_intra])
          cis_sum_control = sum(cis_coverage_control[-chr_intra])
          trans_sum_case = sum(trans_coverage_case[-chr_inter, -chr_inter])
          trans_sum_control = sum(trans_coverage_control[-chr_inter, -chr_inter])
          cis_trans_norm_coefficient = (trans_sum_case/cis_sum_case)/(trans_sum_control/cis_sum_control)
          
          if(is.na(cis_trans_norm_coefficient) || cis_trans_norm_coefficient == Inf || cis_trans_norm_coefficient == 0){
            cis_trans_norm_coefficient = 1
          }
          
          ## calculate coverage normalization coefficient
          case_coverage_chr_inter = sum(trans_coverage_case[chr_inter, -chr_inter])
          case_coverage_chr_intra = sum(trans_coverage_case[chr_intra, -chr_inter])
          control_coverage_chr_inter = sum(trans_coverage_control[chr_inter, -chr_inter])
          control_coverage_chr_intra = sum(trans_coverage_control[chr_intra, -chr_inter])
          
          coverage_norm_coefficient = (case_coverage_chr_inter/case_coverage_chr_intra)/(control_coverage_chr_inter/control_coverage_chr_intra)
          if(is.na(coverage_norm_coefficient) || coverage_norm_coefficient == Inf || coverage_norm_coefficient == 0){
            coverage_norm_coefficient = 1
          }
          ## cis/trans normalization
          vals_res[w, 5] = vals_res[w, 5]/cis_trans_norm_coefficient
          
          ## coverage normalization
          vals_res[w, 5] = vals_res[w, 5]/coverage_norm_coefficient 
        }
      }
    }
  }
  return(vals_res)
}

only_max_vals <- function(vals_res){
  ## extract only max val for every bin
  ord = order(vals_res[, 1], vals_res[, 3], vals_res[, 5], decreasing = T)
  vals_res = vals_res[ord, ]
  
  vals_res_pre = rbind(vals_res[1, ])
  w0 = c()
  for (i in 2:nrow(vals_res)) {
    if(vals_res[i, 1]!=vals_res_pre[, 1] | vals_res[i, 3]!=vals_res_pre[, 3]){
      vals_res_pre[] = vals_res[i, ]
      w0 = c(w0, i)
    }
  }
  vals_res = vals_res[w0,]
  
  ord = order(vals_res[, 1], vals_res[, 2], vals_res[, 3])
  vals_res = vals_res[ord, ]
  
  return(vals_res)
}

filter_inter_trans <- function(vals_res, probability_threshold){
  N = nrow(vals_res)
  q = as.numeric(quantile(vals_res[, 8], 0:10/10))
  
  if(N < 100){    
    val = log(vals_res[w, 5])
    pr = pnorm(val, 0, 0.45, lower.tail = F)
    vals_res[w, 6] = pr
  } else {
    for (i in 1:10) {
      w = which(vals_res[, 8] >= q[i] & vals_res[, 8] <= q[i + 1])
      coverage = vals_res[w, 8]
      sdev = sd(log(vals_res[w, 5]))
      val = log(vals_res[w, 5])
      pr = pnorm(val, 0, sdev, lower.tail = F)
      vals_res[w, 6] = pr
    }
  }
  
  ## extract only max val for every bin
  vals_res = only_max_vals(vals_res)
  
  w = which(vals_res[, 7] > 5 & vals_res[, 8] > 5 & vals_res[, 9] > 5 & vals_res[, 10] > 5 & vals_res[, 6] < probability_threshold)
  
  res = data.frame(rbind(vals_res[w, ]))
  names(res) = c("CHR_FROM", "CHR_TO", "POS_FROM_START", "POS_FROM_END", "VAL", "PROBABILITY", "COVERAGE_CASE_CIS", "COVERAGE_CASE_TRANS", "COVERAGE_CONTROL_CIS", "COVERAGE_CONTROL_TRANS")
  
  return(res)
}

hilocus_inter_trans <- function(args){
  outdir_path = args[1]
  vp_case_path = args[2]
  vp_control_path = args[3]
  binsize = as.numeric(args[4])
  sample_name = args[5]
  probability_threshold = as.numeric(args[6])
  input_format_case = args[7]
  input_format_control = args[8]
  

  ## output file name
  filename = paste(outdir_path, "/inter.translocations.binsize", binsize, ".", sample_name, ".tsv", sep = "")
  
  ## interchromosomal translocation algorithm
  vals_res0 = matrix(0, nrow = 1, ncol = 10)
  cis_coverage_case = rep(0, 22)
  cis_coverage_control = rep(0, 22)
  trans_coverage_case = matrix(0, nrow = 22, ncol =  22)
  trans_coverage_control = matrix(0, nrow = 22, ncol =  22)
  
  for (chr_intra in 1:22) {
    ## read data
    #print(chr_intra)
    gc()
    ######################################################################################################################
    ## read hic intra
    case.intra.bed = try(chrom_pair(vp_case_path, chr_intra, chr_intra, binsize, input_format_case),  silent = TRUE)
    control.intra.bed = try(chrom_pair(vp_control_path, chr_intra, chr_intra, binsize, input_format_control),  silent = TRUE)
    if(inherits(case.intra.bed, "try-error") || inherits(control.intra.bed, "try-error") )
    {
      #error handling code, maybe just skip this iteration using
      next
    }
    nrows = max(c(case.intra.bed$posi, control.intra.bed$posi))
    ncols = max(c(case.intra.bed$posj , control.intra.bed$posj))
    
    case.intra.rowsums = bed_to_rowsums(case.intra.bed, nrows, ncols, chr_intra, chr_intra)
    control.intra.rowsums = bed_to_rowsums(control.intra.bed, nrows, ncols, chr_intra, chr_intra)
    rm(case.intra.bed)
    rm(control.intra.bed)
    
    ## cis sum
    cis_coverage_case[chr_intra] = cis_coverage_case[chr_intra] + sum(case.intra.rowsums)
    cis_coverage_control[chr_intra] = cis_coverage_control[chr_intra] + sum(control.intra.rowsums)
    ######################################################################################################################
    
    for (chr_inter in 1:22) {
      if(chr_inter!=chr_intra){
        #print(chr_inter)
        ######################################################################################################################
        ## read hic inter
        case.inter.bed = try(chrom_pair(vp_case_path, chr_intra, chr_inter, binsize, input_format_case),  silent = TRUE)
        control.inter.bed = try(chrom_pair(vp_control_path, chr_intra, chr_inter, binsize, input_format_control),  silent = TRUE)
        if(inherits(case.inter.bed, "try-error") || inherits(control.inter.bed, "try-error") )
        {
          #error handling code, maybe just skip this iteration using
          next
        }
        
        nrows = max(c(case.inter.bed$posi, control.inter.bed$posi))
        ncols = max(c(case.inter.bed$posj , control.inter.bed$posj))
        
        case.inter.rowsums = bed_to_rowsums(case.inter.bed, nrows, ncols, chr_intra, chr_inter)
        control.inter.rowsums = bed_to_rowsums(control.inter.bed, nrows, ncols, chr_intra, chr_inter)
        
        ## trans sum
        trans_coverage_case[chr_intra, chr_inter] = trans_coverage_case[chr_intra] + sum(case.inter.rowsums)
        trans_coverage_control[chr_intra, chr_inter] = trans_coverage_control[chr_intra] + sum(control.inter.rowsums)
        trans_coverage_case[chr_intra, chr_inter] = trans_coverage_case[chr_inter, chr_intra]
        trans_coverage_control[chr_intra, chr_inter] = trans_coverage_control[chr_inter, chr_intra]
        ######################################################################################################################
        ## result
        N = min(c(length(case.inter.rowsums), length(case.intra.rowsums), length(control.inter.rowsums), length(control.intra.rowsums)))
        vals = (case.inter.rowsums[1:N]/case.intra.rowsums[1:N])/(control.inter.rowsums[1:N]/control.intra.rowsums[1:N])
        
        vals_res = matrix(0, nrow = N, ncol = 10)
        vals_res[, 1] = chr_intra
        vals_res[, 2] = chr_inter
        vals_res[, 3] = seq(from = 0, to = (N - 1)*binsize, by = binsize )
        vals_res[, 4] = seq(from = binsize, to = N*binsize, by = binsize )
        vals_res[, 5] = vals
        vals_res[, 7] = case.intra.rowsums[1:N]
        vals_res[, 8] = case.inter.rowsums[1:N]
        vals_res[, 9] = control.intra.rowsums[1:N]
        vals_res[, 10] = control.inter.rowsums[1:N]
        
        w = which(vals_res[, 5]!=Inf & vals_res[,5]!=0 & !is.na(vals_res[,5]))
        vals_res0 = rbind(vals_res0, vals_res[w,])
      }
    }
  }
  if(nrow(vals_res0!=1)){
    vals_res0 = vals_res0[-1, ]
    
    ## normalization (coverage + cis/trans)
    vals_res = normalization(vals_res0, cis_coverage_case, cis_coverage_control, trans_coverage_case, trans_coverage_control)
    
    ## filtration
    inter_translocations = filter_inter_trans(vals_res, probability_threshold)
    
    ## save results
    write.table(inter_translocations, filename, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
    #write.table(getwd(), paste(filename, "_wd", sep = ""), row.names = F, col.names = F, append = F, quote = F, sep = "\t") 
  } else {
    message("No matching data, please try another binsize or input")
  }
}


library(Matrix)
library(strawr)


options(scipen=999)
args = commandArgs(trailingOnly=TRUE)

## default variables
#binsize = 10000
#binsize_frame = 5000000
#probability_threshold = 10^(-6)



## hilocus_inter_trans algorithm
hilocus_inter_trans(args)

