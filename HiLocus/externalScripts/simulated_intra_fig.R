## model_intra() - modeling of intrachromosomal translocations depending on the coverage
## and implementation of the HiLocus trans algorithm
model_intra <- function(j, coverage, qq, vp_case, vp_control, chr_s, binsize, binsize_frame){
  #print(i)
  ## extract probabilities for a given distance over which locus was translocated
  s1 = vp_case[which(vp_case[, 6] == qq[j]), ]
  chr1 = as.numeric(s1[1, 1])
  t = table(c(s1[,2], s1[,4]))
  pos1 = as.numeric(names(t[which.max(t)]))
  chr2 = as.numeric(s1[1, 3])
  pos2 = as.numeric(names(t[which.max(t)]))
  vp_control1 = vp_control[which(vp_control[, 6] == qq[j]), ]
  vp_case1 = vp_case[which(vp_case[, 6] == qq[j]), ]
  
  genome_size_chr = chr_s[chr1, 2]%/%binsize + 1
  genome_size_chr1 = chr_s[chr1, 2]%/%binsize_frame + 1

  ## calculate contacts vector of translocated locus
  a = getvector_model(vp_case1, genome_size_chr1, binsize, binsize_frame, pos1)
  b = getvector(vp_control1, genome_size_chr1, binsize, binsize_frame, pos1)
  case_coverage = sum(a)
  control_coverage = sum(b)
  
  ## depcov - simulated overall sequencing depth
  p0 = control_coverage/14873176
  depcov = round(coverage/p0)
  
  ## simulate contacts for a given locus
  a_binom = rep(0, length(a))
  for (j in 1:length(a_binom)) {
    p = a[j]/14873176
    a_binom[j] = rbinom(1, depcov, p)
  }
  
  ## HiLocus trans algorithm for cis tranclocations
  if(sum(a_binom)==2){
    vals = NA
  } else {
    
    
    # algoritm
    if(case_coverage!=0 && control_coverage!=0){
      vals = 0
      
      # normalization
      a1 = a_binom/sum(a_binom)
      b1 = b/control_coverage
      
      # algorithm
      case_left = 0
      case_right = sum(a1)
      control_left = 0
      control_right = sum(b1)
      
      
      for (bin1 in 1:genome_size_chr1) {
        
        case_left = case_left + a1[bin1]
        case_right = case_right - a1[bin1]
        control_left = control_left + b1[bin1]
        control_right = control_right - b1[bin1]
        l = (case_left - control_left)
        r = (case_right - control_right)
        val = abs(r - l)
        
        if(vals < val){
          vals = val
        }
      }
    } else {
      vals = NA
    }
  }
  
  return(vals)
}

getvector_model <- function(s, genome_size_chr1, binsize, binsize_frame, pos1){
  Z <- matrix(0 , nrow = 1, ncol = genome_size_chr1)
  x <- matrix(0 , nrow = 1, ncol = genome_size_chr1)
  
  x = c(s[,2], s[, 4])
  y = c(s[, 5], s[, 5])
  w = which(x!=pos1)
  x = x[w]
  y = y[w]
  ##
  v = x%/%binsize_frame + 1
  for (i in 1:length(v)) {
    Z[1, v[i]] = Z[v[i]] + y[i]
  }
  
  
  return(Z[1,])
}


getvector <- function(s, genome_size_chr1, binsize, binsize_frame, pos1){
  Z <- matrix(0 , nrow = 1, ncol = genome_size_chr1)
  x <- matrix(0 , nrow = 1, ncol = genome_size_chr1)
  
  x = c(s[,2], s[, 4])
  w = which(x>=pos1 & x<(pos1 + binsize))
  x = x[-w]
  ##
  v = x%/%binsize_frame + 1
  for (i in 1:length(v)) {
    Z[1, v[i]] = Z[v[i]] + 1
  }
  
  
  return(Z[1,])
}

## precalculate sd and mean per coverage
i = 1
sample_hilocus_trans_result = read.table(paste("/home/evgeniy/clu/202212241802ws26/fs_norm_old_real_intra/s", i, "_10kb/small_intra_translocations_v6", sep = ""), stringsAsFactors = F, head=F)
samples_set = c(2,3,5,7,8,9,10,14,41,42,43,44,45,46,47,48)

for (i in samples_set) {
  print(i)
  sx = read.table(paste("/home/evgeniy/clu/202212241802ws26/fs_norm_old_real_intra/s", i, "_10kb/small_intra_translocations_v6", sep = ""), stringsAsFactors = F, head=F)
  sample_hilocus_trans_result = rbind(sample_hilocus_trans_result, sx)
}


## reading simulated probabilities for translocations
vp_case = read.table("/home/evgeniy/clu/202212241802ws26/tmp_intra_1diag_final/vp_case", stringsAsFactors = F, head=F)
vp_control = read.table("/home/evgeniy/clu/202212241802ws26/tmp_intra_1diag_final/vp_control", stringsAsFactors = F, head=F)
chr_s = read.table("/home/evgeniy/clu/202212241802ws26/FSNorm/chr_sizes_hg19", stringsAsFactors = F, head=F)
binsize = 10000
binsize_frame = 5000000

## reading simulated translocations metadata
s = read.table("/home/evgeniy/clu/202212241802ws26/mutvar_7.txt", stringsAsFactors = F, head=F)

## sort simulated translocations distance
q1 = rep(0,length(s[,1]))
for (i in 1:length(s[,1])) {
  q1[i] = unlist(gregexpr("0-5M", s[i,1]))
}
q2 = rep(0,length(s[,1]))
for (i in 1:length(s[,1])) {
  q2[i] = unlist(gregexpr("5-10M", s[i,1]))
}
q3 = rep(0,length(s[,1]))
for (i in 1:length(s[,1])) {
  q3[i] = unlist(gregexpr("10-50", s[i,1]))
}
q4 = rep(0,length(s[,1]))
for (i in 1:length(s[,1])) {
  q4[i] = unlist(gregexpr("50-", s[i,1]))
}
q5 = rep(0,length(s[,1]))
for (i in 1:length(s[,1])) {
  q5[i] = unlist(gregexpr("100-", s[i,1]))
}
s[which(q1!=-1), 1] = "0-5Mb"
s[which(q2!=-1), 1] = "5-10Mb"
s[which(q3!=-1), 1] = "10-50Mb"
s[which(q4!=-1), 1] = "50-100Mb"
s[which(q5!=-1), 1] = "100Mb"


## filter simulated translocations
q = s[,4]%%10000
w = which(q==0 & (q1!=-1 | q2!=-1 | q3!=-1 | q4!=-1 | q5!=-1))
w = w[which(w<1014)]
s = s[w,]
sample_hilocus_trans_result = sample_hilocus_trans_result[which(sample_hilocus_trans_result[,4]!=0 & sample_hilocus_trans_result[,4]!=2),]

## calculate ROC depending on coverage
f0 = 1:20*(200/20)
ff = c("0-5Mb", "5-10Mb", "10-50Mb", "50-100Mb", "100Mb")
faction_of_true_translocations_all_distance = matrix(0, nrow = length(f0), ncol = length(ff))
nn = 0
for (kl in ff) {
  print(kl)
  nn = nn + 1
  qq = as.numeric(s[which(s[,1]==kl), 2])
  
  roc = rep(0, length(f0))
  faction_of_true_translocations = rep(0, length(f0))
  for (i in 2:length(f0)) {
    print(i)
    coverage = f0[i]
    
    w = which(sample_hilocus_trans_result[, 11]>=f0[i - 1] & sample_hilocus_trans_result[, 11]<f0[i])
    vals_real = sample_hilocus_trans_result[w, 4]
    if(length(vals_real>1000)){
      vals_simulated = rep(0, length(qq))
      for (j in 1:length(qq)) {
        if(length(which(vp_case[, 6] == qq[j])) != 0){
          vals_simulated[j] = model_intra(j, coverage, qq, vp_case, vp_control, chr_s, binsize, binsize_frame) 
        }
      }
      vals_simulated = vals_simulated[which(!is.na(vals_simulated))]
      ##############
      m = median(vals_real)
      sdev = sd(vals_real)
      vals_probabilities = pnorm(vals_simulated, m, sdev, lower.tail = F)

      faction_of_true_translocations[i] = length(which(vals_probabilities < 10^(-8)))/length(vals_probabilities)
      #vals_real = vals_real[order(vals_real, decreasing = T)]
      #faction_of_true_translocations[i] = length(which(vals_simulated>vals_real[10]))/length(vals_simulated)
    }
    
  }
  faction_of_true_translocations_all_distance[, nn] = faction_of_true_translocations
}

pal <- colorRampPalette(c("blue", "red"))
png("faction_of_true_translocations_intra.png",width = 8, height = 5, units="in", res=300)
matplot(f0, faction_of_true_translocations_all_distance,type="l", col = pal(20)[1:5*2], lty=c(1,1),  ylim = c(0,1), lwd = 2, xlab = "Coverage", ylab = "Fraction of called true translocations")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 
legend("topleft", legend = rev(ff), col = rev(pal(20)[1:5*2]), pch = 19, bty = "n")
dev.off()


###############################################################################################################################
###############################################################################################################################
## precalculate sd and mean per coverage
coverage = 1:201
sd_cis = 1:200
mean_cis = 1:200
mean_cis1 = 1:200
for (i in 2:201) {
  w = which(sample_hilocus_trans_result[, 12] >= coverage[i - 1] & sample_hilocus_trans_result[, 12] < coverage[i])
  sd_cis[i - 1] = sd(sample_hilocus_trans_result[w, 4])
  mean_cis[i - 1] = median(sample_hilocus_trans_result[w, 4])
}
write.table(sd_cis, "sd_cis", row.names =F, col.names = F, append = F)
write.table(mean_cis, "mean_cis", row.names =F, col.names = F, append = F)