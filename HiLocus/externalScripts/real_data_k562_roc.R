#rm(list = ls())
#.rs.restartR()
#install.packages("Matrix")
#install.packages("ggplo2")
library(strawr)
library(Matrix)
library(ggplo2)

convert_hic <- function(hic, binsize){
  names(hic) = c("x", "y", "counts")
  hic = hic[which(abs(hic$x - hic$y) >= 10^6), ]
  hic$x = hic$x%/%binsize + 1
  hic$y = hic$y%/%binsize + 1
  hict = hic
  hict$x = hic$y
  hict$y = hic$x
  
  hic_full = rbind(hic, hict)
  return(hic_full)
}


binsize = 10^6
path_to_hic_case = "/media/evgeniy/OS/fastq/GSE63525_K562_combined_30.hic"
path_to_hic_control = "/media/evgeniy/OS/fastq/GSE63525_GM12878_insitu_primary_30.hic"


vals_res0 = matrix(0, nrow = 1, ncol = 10)
cis_coverage_case = rep(0, 22)
cis_coverage_control = rep(0, 22)
trans_coverage_case = rep(0, 22)
trans_coverage_control = rep(0, 22)


### interchromasomal translocation algoritm
for (chr_intra in 1:22) {
  ## read data
  print(chr_intra)
  ######################################################################################################################
  ## read hic intra
  case.intra.bed <- strawr::straw("NONE", path_to_hic_case, as.character(chr_intra), as.character(chr_intra), "BP", binsize)
  control.intra.bed <- strawr::straw("NONE", path_to_hic_control, as.character(chr_intra), as.character(chr_intra), "BP", binsize)
  
  case.intra.bed = convert_hic(case.intra.bed, binsize)
  control.intra.bed  = convert_hic(control.intra.bed, binsize)
  
  nrows = max(c(case.intra.bed$x, control.intra.bed$x))
  ncols = max(c(case.intra.bed$y, control.intra.bed$y))
  
  case.intra.sparse = sparseMatrix(i = case.intra.bed$x, j = case.intra.bed$y, x = case.intra.bed$counts, dims = c(nrows, ncols))
  control.intra.sparse = sparseMatrix(i = control.intra.bed$x, j = control.intra.bed$y, x = control.intra.bed$counts, dims = c(nrows, ncols))
  case_intra = rowSums(case.intra.sparse)
  control_intra = rowSums(control.intra.sparse)
  
  
  ## cis sum
  cis_coverage_case[chr_intra] = cis_coverage_case[chr_intra] + sum(case_intra)
  cis_coverage_control[chr_intra] = cis_coverage_control[chr_intra] + sum(control_intra)
  ######################################################################################################################
  
  for (chr_inter in 1:22) {
    print(chr_inter)
    if(chr_inter!=chr_intra){
      ######################################################################################################################
      ## read hic inter
      case.inter.bed <- strawr::straw("NONE", path_to_hic_case, as.character(chr_intra), as.character(chr_inter), "BP", binsize)
      control.inter.bed <- strawr::straw("NONE", path_to_hic_control, as.character(chr_intra), as.character(chr_inter), "BP", binsize)
      
      case.inter.bed = convert_hic(case.inter.bed, binsize)
      control.inter.bed  = convert_hic(control.inter.bed, binsize)
      
      nrows = max(c(case.inter.bed$x, control.inter.bed$x))
      ncols = max(c(case.inter.bed$y, control.inter.bed$y))
      
      case.inter.sparse = sparseMatrix(i = case.inter.bed$x, j = case.inter.bed$y, x = case.inter.bed$counts, dims = c(nrows, ncols))
      control.inter.sparse = sparseMatrix(i = control.inter.bed$x, j = control.inter.bed$y, x = control.inter.bed$counts, dims = c(nrows, ncols))
      case_inter = rowSums(case.inter.sparse)
      control_inter = rowSums(control.inter.sparse)
      
      ## trans sum
      trans_coverage_case[chr_intra] = trans_coverage_case[chr_intra] + sum(case_inter)
      trans_coverage_control[chr_intra] = trans_coverage_control[chr_intra] + sum(control_inter)
      ######################################################################################################################
      ## result
      N = min(c(length(case_inter), length(case_intra), length(control_inter), length(control_intra)))
      vals = (case_inter[1:N]/case_intra[1:N])/(control_inter[1:N]/control_intra[1:N])
      
      vals_res = matrix(0, nrow = N, ncol = 10)
      vals_res[, 1] = chr_intra
      vals_res[, 2] = chr_inter
      vals_res[, 3] = seq(from = 0, to = (N - 1)*binsize, by = binsize )
      vals_res[, 4] = seq(from = binsize, to = N*binsize, by = binsize )
      vals_res[, 5] = vals
      vals_res[, 6] = 0
      vals_res[, 7] = case_intra[1:N]
      vals_res[, 8] = control_intra[1:N]
      vals_res[, 9] = case_inter[1:N]
      vals_res[, 10] = control_inter[1:N]

      w = which(vals_res[, 5]!=Inf & vals_res[,5]!=0 & !is.na(vals_res[,5]))
      vals_res0 = rbind(vals_res0, vals_res[w,])
    }
  }
}

## rm outliers
w0 = c()
for (i in 7:10) {
  q_up = as.numeric(quantile(vals_res0[which(vals_res0[, i]!=0), i], 0.99))
  q_down = as.numeric(quantile(vals_res0[which(vals_res0[, i]!=0), i], 0.01))
  w = which(vals_res0[, i] > q_down & vals_res0[, i] < q_up)
  w0 = c(w0, w)
}
w = unique(w0)
vals_res0 = vals_res0[w, ]


## normalization (coverage + cis/trans)
cis_sum_case = sum(cis_coverage_case)
cis_sum_control = sum(cis_coverage_control)
trans_sum_case = sum(trans_coverage_case)
trans_sum_control = sum(trans_coverage_control)
cis_trans_norm_coefficient = (trans_sum_case/cis_sum_case)/(trans_sum_control/cis_sum_control)

for (chr_intra in 1:22) {
  for (chr_inter in 1:22) {
    if(chr_intra != chr_inter){
      w = which(vals_res0[, 1] == chr_intra & vals_res0[, 2] == chr_inter)
      ## cis/trans normalization
      vals_res0[w, 4] = vals_res0[w, 4]/cis_trans_norm_coefficient
      
      ## coverage normalization
      coverage_norm_coefficient = (trans_coverage_case[chr_inter]/trans_coverage_case[chr_intra])/(trans_coverage_control[chr_inter]/trans_coverage_control[chr_intra])
      vals_res0[w, 4] = vals_res0[w, 4]/coverage_norm_coefficient 
    }
  }
}


## extract only max val for every bin
vals = vals_res0
ord = order(vals[, 1], vals[, 3], vals[, 5], decreasing = T)
vals = vals[ord, ]

valspre[] = vals[1, ]
w0 = c()
for (i in 2:nrow(vals)) {
  if(vals[i, 1]!=valspre[, 1] | vals[i, 3]!=valspre[, 3]){
    valspre[] = vals[i, ]
    w0 = c(w0, i)
  }
}

## res = vals_res0 and sort
res = vals[w0,]
res = res[order(res[, 1], res[, 2], res[, 3]), ]
res = res[-1,]


## eval roc/auc
q = rep(0, nrow(res))
for (i in 1:nrow(res)) {
  if(res[i, 1] < res[i, 2]){
    q[i] = paste(as.character(res[i, 1]), as.character(res[i, 2]), sep = "_")
  } else {
    q[i] = paste(as.character(res[i, 2]), as.character(res[i, 1]), sep = "_")
  }
}

res = list(res[,1], res[,2], res[,3], q, res[,5], res[,6], res[,7], res[,8], res[,9], res[,10])
names(res) <- c("chri", "chrj", "pos", "chr_pair", "val", "val1", "case_intra_coverage", "control_intra_coverage", "case_inter_coverage", "control_inter_coverage")

## tp_set - validated intra-chromosomal translocations (unique chromosome pairs) in K562 https://www.nature.com/articles/s41588-018-0195-8#Sec36
tp_set = c("1_6","1_18","1_20","3_10","3_18","5_6","9_22","12_21","9_13","13_22","6_16","9_17","10_17","6_18","2_22")

f = 20*0:1000/100
tpr = rep(0, length(f))
fpr = rep(0, length(f))
Nfp = (22*22 - 22)/2
for (i in f) {
  w = which(res$val > f[i])
  called = unique(res$chr_pair[w])
  TP = length(called[called %in% tp_set])
  FP = length(called) - TP
  tpr[i] = TP/length(tp_set)
  fpr[i] = FP/Nfp
}

tpr = c(1, tpr)
fpr = c(1, fpr)
auc = 0
for (i in 2:length(fpr)) {
  auc = auc + tpr[i]*abs(fpr[i] - fpr[i - 1])
}
png("k562_roc.png",width = 4, height = 4, units="in", res=300)
matplot(fpr, tpr,type="l", col = "red", lty=c(1,1),  ylim = c(0,1), lwd = 5, xlab = "False positive rate", ylab = "True positive rate")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 
dev.off()

