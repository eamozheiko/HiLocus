
case_in_control = 0
thresholds = NA
thr_frame = 0
binsize = 10000

# probability
filter_inter_trans <- function(s, sd_inter_trans, thresholds, filename, thr_frame){
  if(!is.na(thresholds)){
    print("")
  } else {
    s[, 9] = as.numeric(s[, 9])
    s[, 4] = as.numeric(s[, 4])
    s[which(is.na(s[, 9])), 9] = 0
    
    s[,10] = s[,9]%/%1 + 1
    
    s[which(s[,10]>64), 10] = 64
    v = sd_inter_trans[s[,10],1]
    val = log(s[, 4])
    pr = pnorm(val, 0, v, lower.tail = F)
    s[,11] = pr
    
    w = which(s[,9]>10 & s[,11]<0.000001)
    s = s[w,]
    
    prev_chr = 0
    prev_pos = 0
    prev_val = 2
    
    ord = order(s[,1], s[, 2])
    s = s[ord, ]
    s[, 17] = 1
    for (i in 1:nrow(s)) {
      if(prev_chr == s[i, 1] && prev_pos == s[i, 2]){
        if(prev_val > s[i, 11]){
          s[i - 1, ] = s[i, ]
        } else {
          s[i, ] = s[i - 1, ]
        }
        s[i - 1, 17] = 0
      }
      prev_chr = s[i, 1]
      prev_pos = s[i, 2]
      prev_val = s[i, 11]
    }
    
    w = which(s[, 17] == 1 & s[, 16] > thr_frame)
    s = s[w, 1:16]
    q = data.frame(s)
    names(q) = c("CHR_FROM", "POS_FROM_ST", "POS_FROM_EN", "VALS", "COVERAGE_CASE_INTRA", "COVERAGE_CASE_INTER",
                 "CHR_TO", "COVERAGE_CONTROL_INTRA", "COVERAGE_CONTROL_INTER", "ROUND_COVERAGE_CONTROL_INTER", "PROBABILITY", 
                 "BIN_NUMBER_CHR_TO_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL", "CHR_TO_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL", 
                 "BIN_NUMBER_CHR_TO_MIN_VALUE_OF_CASE_DIVIDED_BY_CONTROL", "CHR_TO_MIN_VALUE_OF_CASE_DIVIDED_BY_CONTROL", 
                 "ABSOLUTE_CASE_VALUE_OF_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL")
    write.table(q, paste(filename, "_filtred", sep = ""), row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  }
}

get_matrix <- function(s, chr_s, chr1, chr2, binsize, binsize_frame){
  genome_size_chr = chr_s[chr1, 2]%/%binsize + 1
  gonome_size_chr_frame = chr_s[chr2, 2]%/%binsize_frame + 1
  
  Z <- matrix(0 , nrow = genome_size_chr, ncol = gonome_size_chr_frame)
  
  if(chr1!=chr2){
    x = s[, 2]%/%binsize + 1
    y = s[, 4]%/%binsize_frame + 1
  } else {
    x = c(s[, 2]%/%binsize + 1, s[, 4]%/%binsize + 1)
    y = c(s[, 4]%/%binsize_frame + 1, s[, 2]%/%binsize_frame + 1)
  }
  ##
  n = nrow(Z)
  v = x + n*(y - 1)
  
  for (i in 1:length(v)) {
    Z[v[i]] = Z[v[i]] + 1
  }
  
  return(Z)
}

vp_convert <- function(s, binsize){
  s = s[which(s[, 1]!="X" & s[, 1]!="Y" & s[, 3]!="X" & s[, 3]!="Y"), ]
  s = s[which((s[, 4] - s[, 2])>binsize | s[, 1]!=s[, 3]), ]
  s = s[, 1:4]
  s[, 1] = as.numeric(s[, 1])
  s[, 3] = as.numeric(s[, 3])
  #s = s[which(s[,1] == s[,3]),]
  return(s)
}

get_genome_size <- function(chr_s, binsize){
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

cis_trans_norm <- function(coverage_track_case, coverage_track_control, chr_s, binsize, ic){
  d1 = matrix(0, nrow = 22, ncol = 22)
  d2 = matrix(0, nrow = 22, ncol = 22)
  d3 = matrix(0, nrow = 22, ncol = 22)
  d4 = matrix(0, nrow = 22, ncol = 22)
  d5 = matrix(0, nrow = 22, ncol = 22)
  d6 = matrix(0, nrow = 22, ncol = 22)
  d7 = matrix(0, nrow = 22, ncol = 22)
  d8 = matrix(0, nrow = 22, ncol = 22)
  d9 = matrix(0, nrow = 22, ncol = 22)
  d10 = matrix(0, nrow = 22, ncol = 22)
  for (chr_intra in 1:22) {
    print(chr_intra)
    
    genome_size_chr = chr_s[chr_intra, 2]%/%binsize + 1
    
    case = coverage_track_case[1:genome_size_chr + ic[chr_intra, 1], ]
    control = coverage_track_control[1:genome_size_chr + ic[chr_intra, 1], ]
    
    #inter_case = case[, -c(chr_inter + 3)]
    #inter_case = rowSums(inter_case[,4:24])
    
    #inter_control = control[, -c(chr_intra + 3)]
    #inter_control = rowSums(inter_control[,4:24])
    #cis_trans_coef
    
    a1 = case
    a1 = a1[, -c(chr_intra)]
    a2 = case[, chr_intra]
    b1 = control
    b1 = b1[, -c(chr_intra)]
    b2 = control[, chr_intra]
    
    for (chr_inter in 1:22) {
      if(chr_inter!=chr_intra){
        
        d1[chr_intra, chr_inter] = sum(case[, chr_inter])
        d2[chr_intra, chr_inter] = sum(case[, chr_intra])
        d3[chr_intra, chr_inter] = sum(control[, chr_inter])
        d4[chr_intra, chr_inter] = sum(control[, chr_intra])
        d5[chr_intra, chr_inter] = sum(case)
        d6[chr_intra, chr_inter] = sum(control)
        d7[chr_intra, chr_inter] = sum(a1)
        d8[chr_intra, chr_inter] = sum(a2)
        d9[chr_intra, chr_inter] = sum(b1)
        d10[chr_intra, chr_inter] = sum(b2)
      }
    }
  }
  d = (d1/d2)/(d3/d4)
  d1 = (d7/d8)/(d9/d10)
  d1[which(is.na(d1))] = 0
  d[which(d==Inf | d==-Inf | is.na(d))] = 0
  md = colSums(d)/21
  md1 = rowSums(d1)/21
  md11 = md1/mean(md1)
  ctm = cbind(md11)%*%rbind(md)
  return(1/ctm)
}


DIR = "/home/evgeniy/clu/202204142105ws21/ExoC"
outdir_path = "/home/evgeniy/clu/202204142105ws21/fs_norm/s141/"
bic_path = "/home/evgeniy/clu/202204142105ws21/bic_out/BIC_out_cnv_s141_fs_norm"
vp_case_path = "/home/evgeniy/clu/202204142105ws21/allValidPairs/t"
vp_control_path = "/home/evgeniy/clu/202204142105ws21/allValidPairs/t1"
chr_sizes_path = "/mnt/scratch/ws/eamozheiko/202204142105ws21/FSNorm/chr_sizes_hg19"

DIR = "/home/evgeniy/clu/202204142105ws21/FSNorm"
outdir_path = "/home/evgeniy/clu/202204142105ws21/fs_norm/s97_10kb/"
bic_path = "/home/evgeniy/clu/202204142105ws21/bic_out/BIC_out_cnv_s97_fs_norm"
vp_case_path = "/home/evgeniy/clu/202204142105ws21/allValidPairs/allValidPairs97_10kb"
vp_control_path = "/home/evgeniy/clu/202204142105ws21/allValidPairs/allValidPairs97_10kb_control"
chr_sizes_path = "/mnt/scratch/ws/eamozheiko/202204142105ws21/FSNorm/chr_sizes_hg19"

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

sd_inter_trans = read.table(paste(DIR,"/data/trans_sd", sep = ""), stringsAsFactors = F, head=F, dec = ",")

sd_inter_trans[,1] = as.numeric(sd_inter_trans[,1])

#binsize = 10000
binsize_frame = 5000000

genome_size = get_genome_size(chr_s, binsize)
ic = get_chr_index_correction(chr_s, binsize)


vp_case = read.table(vp_case_path, stringsAsFactors = F, head=F)
vp_control = read.table(vp_control_path, stringsAsFactors = F, head=F)

vp_case = vp_convert(vp_case, binsize)
vp_control = vp_convert(vp_control, binsize)


### for all chromosomes implement cis variations algoritm
coverage_track_case = matrix(0, nrow = genome_size, ncol = 22)
coverage_track_control = matrix(0, nrow = genome_size, ncol = 22)
track_which_max = matrix(0, nrow = genome_size, ncol = 22)
track_max = matrix(0, nrow = genome_size, ncol = 22)
track_which_min = matrix(0, nrow = genome_size, ncol = 22)
track_min = matrix(0, nrow = genome_size, ncol = 22)
track_abs_val = matrix(0, nrow = genome_size, ncol = 22)
vals_res = matrix(0, nrow = 1, ncol = 17)

### calc row sums by all chr pairs
for (chr1 in 1:22) {
  print(chr1)
  ## initialize varibles
  # genome sizes for different binsizes
  genome_size_chr = chr_s[chr1, 2]%/%binsize + 1
  
  
  for (chr2 in 1:22) {
    w_case = which((vp_case[, 1]==chr1 & vp_case[, 3]==chr2) | (vp_case[, 1]==chr2 & vp_case[, 3]==chr1))
    w_control = which((vp_control[, 1]==chr1 & vp_control[, 3]==chr2) | (vp_control[, 1]==chr2 & vp_control[, 3]==chr1))
    if(chr2 < chr1){
      vp_case_chr = vp_case[w_case, c(2,4,1,2)]
      vp_control_chr = vp_control[w_control, c(2,4,1,2)]
    } else {
      vp_case_chr = vp_case[w_case, c(1:4)]
      vp_control_chr = vp_control[w_control, c(1:4)]
    }
    # calculate contacts matrix
    a = get_matrix(vp_case_chr, chr_s, chr1, chr2, binsize, binsize_frame)
    b = get_matrix(vp_control_chr, chr_s, chr1, chr2, binsize, binsize_frame)
    bx = b
    bx[which(bx==0)] = 1
    c = a/bx
    c[which(is.na(c) | c==Inf | c==-Inf) ] = 0
    for (r in 1:nrow(c)) {
      wm = which.max(c[r, ])
      m = max(c[r, ])
      abs_val_p = a[r, wm]
      
      track_which_max[r + ic[chr1, 1], chr2] = wm
      track_max[r + ic[chr1, 1], chr2] = m
      track_abs_val[r + ic[chr1, 1], chr2] = abs_val_p
      
      ##
      c1 = c[r,]
      c1[which(c1==0)] = Inf
      wm = which.min(c1)
      m = min(c1)
      
      track_which_min[r + ic[chr1, 1], chr2] = wm
      track_min[r + ic[chr1, 1], chr2] = m
      
    }
    coverage_track_case[1:genome_size_chr + ic[chr1, 1], chr2] = rowSums(a)
    coverage_track_control[1:genome_size_chr + ic[chr1, 1], chr2] = rowSums(b)
  }
}


### cis trans normalization coefficient
cis_trans_norm_coef = cis_trans_norm(coverage_track_case, coverage_track_control, chr_s, binsize, ic)


## del res file

filename = paste(outdir_path, "/small_inter_translocations_", sample_name, sep = "")

if (file.exists(filename)) {
  #Delete file if it exists
  file.remove(filename)
}


### interchromasomal translocation algoritm
for (chr_intra in 1:22) {
  ## read data
  genome_size_chr = chr_s[chr_intra, 2]%/%binsize + 1
  
  
  case = coverage_track_case[1:genome_size_chr + ic[chr_intra, 1], ]
  control = coverage_track_control[1:genome_size_chr + ic[chr_intra, 1], ]
  wmax = track_which_max[1:genome_size_chr + ic[chr_intra, 1], ]
  mmax = track_max[1:genome_size_chr + ic[chr_intra, 1], ]
  wmin = track_which_min[1:genome_size_chr + ic[chr_intra, 1], ]
  mmin = track_min[1:genome_size_chr + ic[chr_intra, 1], ]
  abs_val = track_abs_val[1:genome_size_chr + ic[chr_intra, 1], ]
  
  
  ## control correction
  if(case_in_control == 1){
    control = control - case
  }
  
  for (chr_inter in 1:22) {
    print(chr_inter)
    if(chr_inter!=chr_intra){
      case_intra = case[, chr_intra]
      case_inter = case[, chr_inter]
      control_intra = control[, chr_intra]
      control_inter = control[, chr_inter]
      wmax_chr = wmax[, chr_inter]
      mmax_chr = mmax[, chr_inter]
      wmin_chr = wmin[, chr_inter]
      mmin_chr = mmin[, chr_inter]
      aabs_val = abs_val[, chr_inter]
        
      vals = (case_inter/case_intra)/(control_inter/control_intra)
      
      
      
      
      # cis trans correction
      vals = vals*cis_trans_norm_coef[chr_intra, chr_inter]
      
      # NA and Inf removing
      vals[which((vals == Inf) | is.na(vals))] = 0
      
      

      vals_res_chr = matrix(0, nrow = genome_size_chr, ncol = 17)
      
      
      
      vals_res_chr[, 1] = chr_intra
      vals_res_chr[, 2] = seq(from = 0, to = (genome_size_chr - 1)*binsize, by = binsize )
      vals_res_chr[, 3] = seq(from = binsize, to = genome_size_chr*binsize, by = binsize )
      vals_res_chr[, 4] = vals
      vals_res_chr[, 5] = case[, chr_intra]
      vals_res_chr[, 6] = case[, chr_inter]
      vals_res_chr[, 7] = chr_inter
      vals_res_chr[, 8] = control[, chr_intra]*(rowSums(case)/rowSums(control))
      vals_res_chr[, 9] = control[, chr_inter]*(rowSums(case)/rowSums(control))
      a = case[, chr_inter]
      b = control[, chr_inter]*(rowSums(case)/rowSums(control))

      vals_res_chr[, 12] = wmax_chr
      vals_res_chr[, 13] = mmax_chr
      vals_res_chr[, 14] = wmin_chr
      vals_res_chr[, 15] = mmin_chr
      vals_res_chr[, 16] = aabs_val

      vals_res = rbind(vals_res, vals_res_chr)
      

    }
  }
}
vals_res = vals_res[-1, ]

#vals_res = data.frame(vals_res)
#names(vals_res) = c("CHR_FROM", "POS_FROM_ST", "POS_FROM_EN", "VALS", "COVERAGE_CASE_INTRA", "COVERAGE_CASE_INTER",
#             "CHR_TO", "COVERAGE_CONTROL_INTRA", "COVERAGE_CONTROL_INTER", "ROUND_COVERAGE_CONTROL_INTER", "PROBABILITY", 
#             "BIN_NUMBER_CHR_TO_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL", "CHR_TO_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL", 
#             "BIN_NUMBER_CHR_TO_MIN_VALUE_OF_CASE_DIVIDED_BY_CONTROL", "CHR_TO_MIN_VALUE_OF_CASE_DIVIDED_BY_CONTROL", 
#             "ABSOLUTE_CASE_VALUE_OF_MAX_VALUE_OF_CASE_DIVIDED_BY_CONTROL")

# threshold
w = which(vals_res[, 6] > 0)
if(length(w) > 0){
  if(length(w) == 1){
    write.table(rbind(vals_res[w, ]), filename, row.names = F, col.names = T, append = F, sep = "\t")
  } else {
    write.table(vals_res[w, ], filename, row.names = F, col.names = T, append = F, sep = "\t")
  }
}

vals_res = vals_res[w, ]


## extract high confidence translocations

filter_inter_trans(vals_res, sd_inter_trans, thresholds, filename, thr_frame)
  

#s2 = read.table("/home/evgeniy/clu/202204142105ws21/small_inter_translocations_sample1", stringsAsFactors = F, head=T, dec = ".")
#s = matrix(unlist(s), nrow = nrow(s))




      
      
      
