options(scipen=999)

thresholds = NA

# probability
filter_intra_trans <- function(s, mean_intra_trans, sd_intra_trans, thresholds, filename, thr_frame){
  if(!is.na(thresholds)){
    print("")
  } else {
    s[,14] = s[,11]%/%10
    s[which(s[,14]>100), 14] = 100
    s = s[which(s[,14]>0),]
    m = mean_intra_trans[s[,14],1]
    v = sd_intra_trans[s[,14],1]
    val = log(s[, 4])
    pr = pnorm(val, m, v, lower.tail = F)
    s[,15] = pr
    
    w = which(s[,15]<10^-5 & s[,11] > 25)
    s = s[w,]
    q = s[, 13]
    s[, 13:14] = s[, 14:15]
    s[, 15] = q
    
    w = which(s[, 15] > thr_frame)
    s = s[w, ]
    q = data.frame(s)
    names(q) = c("CHR_FROM", "POS_FROM_ST", "POS_FROM_EN", "VALS", "IS_TRANS_IN_CASE", "HALF_TRANS_DIST(5MB_BIN_COORDINATE)", "MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL_IN_5MB_FRAME",
                 "MIN_VALUE_OF_CASE_DIVIDED_BY_CONTROL_IN_5MB_FRAME", "COORDINATE_OF_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL_IN_5MB_FRAME",
                 "COORDINATE_OF_MIN_VALUE_OF_CASE_DIVIDED_BY_CONTROL_IN_5MB_FRAME", "COVERAGE_CONTROL_INTRA", "COVERAGE_CASE_INTRA", "COVERAGE_CONTROL_INTRA_ROUNDED", 
                 "PROBABILITY", "ABSOLUTE_CASE_VALUE_OF_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL_IN_5MB_FRAME")
    write.table(q, paste(filename, "_filtred", sep = ""), row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  }
}

get_matrix <- function(s, genome_size_chr, genome_size_chr1, binsize, binsize_frame){
  Z1 <- matrix(0 , nrow = genome_size_chr, ncol = genome_size_chr1)
  Z2 <- matrix(0 , nrow = genome_size_chr, ncol = genome_size_chr1)
  
  ##
  x = s[, 2]%/%binsize + 1
  y = s[, 4]%/%binsize_frame + 1
  n = nrow(Z1)
  v = x + n*(y - 1)
  for (i in 1:length(v)) {
    Z1[v[i]] = Z1[v[i]] + 1
  }
  
  ##
  x = s[, 4]%/%binsize + 1
  y = s[, 2]%/%binsize_frame + 1
  n = nrow(Z2)
  v = x + n*(y - 1)
  for (i in 1:length(v)) {
    Z2[v[i]] = Z2[v[i]] + 1
  }
  
  ##

  Z = Z1 + Z2
  
  Z[which(is.na(Z))] = 0
  
  return(Z)
}

vp_convert <- function(s, binsize){
  s = s[which(s[, 1]!="X" & s[, 1]!="Y" & s[, 3]!="X" & s[, 3]!="Y"), ]
  s = s[which((s[, 4] - s[, 2])>binsize | s[, 1]!=s[, 3]), ]
  s = s[, 1:4]
  s[, 1] = as.numeric(s[, 1])
  s[, 3] = as.numeric(s[, 3])
  s = s[which(s[,1] == s[,3]),]
  return(s)
}

get_genome_size_all_chr <- function(chr_s, binsize){
  genome_size = 0
  for (i in 1:22) {
    genome_size = genome_size + chr_s[i,2]%/%binsize + 1
  }
  return(genome_size)
}

get_chr_index_correction <- function(chr_s, binsize){
  ic = matrix(0, nrow = 24, ncol = 1)
  ic[1,1] = 0
  for (i in 2:24) {
    ic[i,1] = sum(chr_s[1:(i - 1),2]%/%binsize + 1)
  }
  return(ic)
}


options(scipen=999)
args = commandArgs(trailingOnly=TRUE)


DIR = args[1]
outdir_path = args[2]
vp_case_path = args[3]
vp_control_path = args[4]
binsize = args[5]
sample_name = args[6]
thr_frame = args[7]

binsize = as.numeric(binsize)
thr_frame = as.numeric(thr_frame)

chr_s = read.table(paste(DIR,"/data/chr_sizes_hg19", sep = ""), stringsAsFactors = F, head=F)
chr_s[,2] = as.numeric(chr_s[,2])

mean_intra_trans = read.table(paste(DIR,"/data/cis_m", sep = ""), stringsAsFactors = F, head=F, dec = ",")
mean_intra_trans[,1] = as.numeric(mean_intra_trans[,1])

sd_intra_trans = read.table(paste(DIR,"/data/cis_sd", sep = ""), stringsAsFactors = F, head=F, dec = ",")
sd_intra_trans[,1] = as.numeric(sd_intra_trans[,1])


binsize_frame = 5000000
genome_size = get_genome_size_all_chr(chr_s, binsize)
ic = get_chr_index_correction(chr_s, binsize)



vp_case = read.table(vp_case_path, stringsAsFactors = F, head=F)
vp_control = read.table(vp_control_path, stringsAsFactors = F, head=F)

vp_case = vp_convert(vp_case, binsize)
vp_control = vp_convert(vp_control, binsize)

vals_res = matrix(0, nrow = 1, ncol = 15)

filename = paste(outdir_path, "/small_intra_translocations_", sample_name, sep = "")

if (file.exists(filename)) {
  #Delete file if it exists
  file.remove(filename)
}
### for all chromosomes implement cis variations algoritm


for (chr in 1:22) {
  ## initialize varibles
  # genome sizes for different binsizes
  genome_size_chr = chr_s[chr, 2]%/%binsize + 1
  genome_size_chr1 = chr_s[chr, 2]%/%binsize_frame + 1
  
  # filter trans contacts
  vp_case_intra_chr = vp_case[which(vp_case[, 1]==chr & vp_case[, 3]==chr), ]
  vp_control_intra_chr = vp_control[which(vp_control[, 1]==chr & vp_control[, 3]==chr), ]
  
  # calculate contacts matrix
  a = get_matrix(vp_case_intra_chr, genome_size_chr, genome_size_chr1, binsize, binsize_frame)
  b = get_matrix(vp_control_intra_chr, genome_size_chr, genome_size_chr1, binsize, binsize_frame)
  
  
  
  ## cis variations algoritm
  
  # initialize varibles
  vals = rep(0, genome_size_chr)
  vals1 = rep(0, genome_size_chr)
  coords = rep(0, genome_size_chr)
  max_dev_val = rep(0, genome_size_chr)
  min_dev_val = rep(0, genome_size_chr)
  max_dev_coord = rep(0, genome_size_chr)
  min_dev_coord = rep(0, genome_size_chr)
  max_dev_abs_val = rep(0, genome_size_chr)
  a_norm = a/rowSums(a)
  b_norm = b/rowSums(b)
  
  # algoritm
  for (bin in 1:nrow(a)) {
    if(sum(a[bin,])!=0 && sum(b[bin,])!=0){
      a1 = a_norm[bin,]
      b1 = b_norm[bin,]
      
      case_left = 0
      case_right = sum(a1)
      control_left = 0
      control_right = sum(b1)
      
      c1 = a1/b1
      w = which(!is.na(c1) & c1!=Inf & c1!=0)
      c1 = c1[w]
      if(length(c1!=0)){
        max_dev_val[bin] = max(c1)
        min_dev_val[bin] = min(c1)
        
        n1 = 1:genome_size_chr1
        n2 = 1:genome_size_chr1
        n1 = n1[w]
        n2 = n2[w]
        
        max_dev_coord[bin] = n1[which.max(c1)]
        min_dev_coord[bin] = n2[which.min(c1)]
        max_dev_abs_val[bin] = a[bin, n1[which.max(c1)]]
      }

      
      for (bin1 in 1:ncol(a)) {
      
        case_left = case_left + a1[bin1]
        case_right = case_right - a1[bin1]
        control_left = control_left + b1[bin1]
        control_right = control_right - b1[bin1]
        l = (case_left - control_left)
        r = (case_right - control_right)
        val = abs(r - l)
        val1 = 0
        ql = (bin1 - 1)*(binsize_frame/binsize)
        #qr = bin1*(binsize_frame/binsize)
        if(r > l && bin >= ql){
          val1 = 1
        }
        if(l > r && bin < ql){
          val1 = 1
        }

        if(vals[bin] < val){
          vals[bin] = val
          vals1[bin] = val1
          coords[bin] = bin1
        }
        
        
      
        
      }
    }
  }
  vals_res_chr = matrix(0, nrow = genome_size_chr, ncol = 15)
  vals_res_chr[,1] = chr
  vals_res_chr[,2] = seq(from = 0, to = (genome_size_chr - 1)*binsize, by = binsize )
  vals_res_chr[,3] = seq(from = binsize, to = genome_size_chr*binsize, by = binsize )
  vals_res_chr[,4] = vals
  vals_res_chr[,5] = vals1
  vals_res_chr[,6] = coords
  vals_res_chr[,7] = max_dev_val
  vals_res_chr[,8] = min_dev_val
  vals_res_chr[,9] = max_dev_coord
  vals_res_chr[,10] = min_dev_coord
  vals_res_chr[,11] = rowSums(b)*(sum(a)/sum(b))
  vals_res_chr[,12] = rowSums(a)
  vals_res_chr[,13] = max_dev_abs_val
  vals_res = rbind(vals_res, vals_res_chr)
  
}
vals_res = vals_res[-1, ]

#vals_res = data.frame(vals_res)
#names(vals_res) = c("CHR_FROM", "POS_FROM_ST", "POS_FROM_EN", "VALS", "COVERAGE_CASE_INTRA", "COVERAGE_CASE_INTER",
#             "CHR_TO", "COVERAGE_CONTROL_INTRA", "COVERAGE_CONTROL_INTER", "ROUND_COVERAGE_CONTROL_INTER", "PROBABILITY", 
#             "BIN_NUMBER_CHR_TO_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL", "CHR_TO_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL", 
#             "BIN_NUMBER_CHR_TO_MIN_VALUE_OF_CASE_DIVIDED_BY_CONTROL", "CHR_TO_MIN_VALUE_OF_CASE_DIVIDED_BY_CONTROL", 
#             "ABSOLUTE_CASE_VALUE_OF_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL")

# threshold
w = which(vals_res[, 12] > 0)
if(length(w) > 0){
  if(length(w) == 1){
    write.table(rbind(vals_res[w, ]), filename, row.names = F, col.names = T, append = F, sep = "\t")
  } else {
    write.table(vals_res[w, ], filename, row.names = F, col.names = T, append = F, sep = "\t")
  }
}

vals_res = vals_res[w, ]


## extract high confidence translocations

filter_intra_trans(vals_res, mean_intra_trans, sd_intra_trans, thresholds, filename, thr_frame)


