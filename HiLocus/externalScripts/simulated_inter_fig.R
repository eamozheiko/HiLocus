## model_inter() - modeling of interchromosomal translocations depending on the coverage
## and implementation of the HiLocus trans algorithm without normalization
library(Matrix)
library(ggplo2)

model_inter <- function(coverage, case_probabilities, control_probabilities, chr_s, binsize){
  ## coef - coefficient used at the first stage of intrachromosomal translocation modeling
  coef = 14873176
  
  ## evaluate case probabilities
  case_cis_probability = case_probabilities[, 1]/coef
  case_trans_probability = case_probabilities[, 2]/coef
  control_cis_probability = control_probabilities[, 1]
  control_trans_probability = control_probabilities[, 2]
  
  ## depcov - simulated overall sequencing depth
  control_coverage = control_probabilities[,2]
  p0 = control_coverage/coef
  depcov = round(coverage/p0)
  
  case_simulated_cis_contacts_count = rep(0, length(case_cis_probability))
  case_simulated_trans_contacts_count = rep(0, length(case_trans_probability))
  
  ## simulate contacts count
  for (i in 1:length(case_cis_probability)) {
    p1 = case_cis_probability[i]
    p2 = case_trans_probability[i]
    case_simulated_cis_contacts_count[i] = rbinom(1, depcov[i], p1)
    case_simulated_trans_contacts_count[i] = rbinom(1, depcov[i], p2)
  }
  
  ## HiLocus trans algorithm without normalization
  control_val = control_trans_probability/control_cis_probability
  case_val = case_simulated_trans_contacts_count/case_simulated_cis_contacts_count
  val = case_val/control_val
  
  return(val)
}


## reading simulated probabilities for translocations
vp_case = read.table("/mnt/scratch/ws/eamozheiko/202212241802ws26/tmp_inter_1diag_final/vp_case", stringsAsFactors = F, head=F)
vp_control = read.table("/mnt/scratch/ws/eamozheiko/202212241802ws26/tmp_inter_1diag_final/vp_control", stringsAsFactors = F, head=F)
chr_s = read.table("/mnt/scratch/ws/eamozheiko/202212241802ws26/FSNorm/chr_sizes_hg19", stringsAsFactors = F, head=F)
binsize = 10000

## get number of simulated translocations
u = unique(vp_case[,6])

## calculate cis- and trans-contact probabilities for translocated loci
case_probabilities = matrix(0, nrow = length(u), ncol = 5)
control_probabilities = matrix(0, nrow = length(u), ncol = 3)
n = 0
for (i in u) {
  n = n + 1
  case = vp_case[which(vp_case[,6] == i), ]
  case_probabilities[n, 1] = sum(case[which(case[,1] == case[,3]), 5])
  case_probabilities[n, 2] = sum(case[which(case[,1] != case[,3]), 5])
  
  control = vp_control[which(vp_control[,6] == i), ]
  control_probabilities[n, 1] = sum(control[which(control[,1] == control[,3]), 5])
  control_probabilities[n, 2] = sum(control[which(control[,1] != control[,3]), 5])
}



## calculate ROC depending on coverage
coverages = 1:20*(50/20)
samples_set = c(1,2,3,5,8,9,10,14,41,42,43,44,45,46,47,48)
faction_of_true_translocations_all = matrix(0, nrow = length(coverages), ncol = length(samples_set))

n = 0
for (sample in samples_set) {
  n = n + 1
  print(sample)
  sample_hilocus_trans_result = read.table(paste("/mnt/scratch/ws/eamozheiko/202212241802ws26/fs_norm_old_real_inter/s", sample, "_10kb/small_inter_translocations_v8", sep = ""), stringsAsFactors = F, head=F)
  faction_of_true_translocations_sample = rep(0, length(coverages))

  for (i in 2:length(coverages)) {
    #print(i)
    coverage = coverages[i]
    
    w = which(sample_hilocus_trans_result[, 6]>=coverages[i - 1] & sample_hilocus_trans_result[, 6]<coverages[i])
    vals_real = sample_hilocus_trans_result[w, 4]
    if(length(vals_real)>100){
      vals_simulated = model_inter(coverage, case_probabilities, control_probabilities, chr_s, binsize)
      vals_simulated = vals_simulated[which(!is.na(vals_simulated) & vals_simulated!=Inf & vals_simulated!=-Inf & vals_simulated!=0)]
      
      vals_real = vals_real[order(vals_real, decreasing = T)]
      faction_of_true_translocations_sample[i] = length(which(vals_simulated>vals_real[10]))/length(vals_simulated)
    }
  }
  faction_of_true_translocations_all[, n] = faction_of_true_translocations_sample
}

#rocpp = read.table("/mnt/scratch/ws/eamozheiko/202212241802ws26/rocpp_inter6", stringsAsFactors = F, head=F)

pal <- colorRampPalette(c("blue", "red"))
w = which(faction_of_true_translocations_all==0)
faction_of_true_translocations_all[w] = NA
png("faction_of_true_translocations_inter.png",width = 8, height = 5, units="in", res=300)
matplot(coverages, faction_of_true_translocations_all,type="l", col = pal(length(coverages)), lty=c(1,1),  ylim = c(0,1), lwd = 2, xlab = "Coverage", ylab = "Fraction of called true translocations")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 
dev.off()

###############################################################################################################################
###############################################################################################################################
## precalculate sd per coverage
sample = 1
s = read.table(paste("/mnt/scratch/ws/eamozheiko/202212241802ws26/fs_norm_old_real_inter/s", sample, "_10kb/inter_translocations_", sample, sep = ""), stringsAsFactors = F, head=F)

samples_set = c(2,3,5,7,8,9,10,14,41,42,43,44,45,46,47,48)

for (sample in samples_set) {
  print(sample)
  sx = read.table(paste("/mnt/scratch/ws/eamozheiko/202212241802ws26/fs_norm_old_real_inter/s", sample, "_10kb/inter_translocations_", sample, sep = ""), stringsAsFactors = F, head=F)
  s = rbind(s, sx)
}
coverage = 1:51
sd_trans = 1:50
prev = 0
for (i in 2:51) {
  w = which(s[, 5] >= coverage[i - 1] & s[, 8] < coverage[i])
  sd_trans[i - 1] = sd(log(s[w, 5]))
}
# correct
for (i in 1:49) {
  if(sd_trans[i + 1] > sd_trans[i]){
    sd_trans[i + 1] = sd_trans[i]
  }
}
write.table(sd_trans, "sd_trans", row.names =F, col.names = F, append = F)






