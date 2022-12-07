#!/bin/bash


#####################################################################################################
### PATHS TO INPUT FILES
#####################################################################################################
#PATH_TO_MODELS=
#PATH_TO_MODELS_LIST=
#PATH_TO_CONTROL_VALID_PAIRS=
#OUTDIR_FOLDER_PATH=

#####################################################################################################
## preparations
#####################################################################################################

mkdir ${OUTDIR_FOLDER_PATH}
cd ${OUTDIR_FOLDER_PATH}


NROW=$(cat ${PATH_TO_MODELS_LIST} | wc -l)
FROM=1
TO=$NROW

rm ${OUTDIR_FOLDER_PATH}/vp_case
rm ${OUTDIR_FOLDER_PATH}/vp_control


for ROW in $(seq $FROM $TO)
do 
#echo ${ROW}

echo ${ROW} > row

cat ${PATH_TO_MODELS_LIST} | awk 'BEGIN{getline row < "row"}{
if(NR==row){
print $0
}
}' > sample_inf

cat sample_inf | awk '{
print $2
}' > sample

cat sample_inf | awk '{
print $3
}' > chr1

cat sample_inf | awk '{
print $7
}' > chr2

## get only positions that are divisible by the binsize (10000 kb)
cat sample_inf | awk '{
if(($4 % 10000)!=5000){
  print $4
} else {
  print "-1"
}
}' > pos

SAMPLE=$(cat sample)
#CHR1=$(cat chr1)
#CHR2=$(cat chr2)
#POS=$(cat pos)

#####################################################################################################
## extract all contacts from translocated loci
#####################################################################################################

if [ $POS -ne  -1 ]
then
zcat ${PATH_TO_MODELS}/mut_${SAMPLE}.summ.pre.short.gz | sed 's/chr//g' | awk '
BEGIN{getline sample < "sample"; getline chr1 < "chr1"; getline chr2 < "chr2"; getline pos < "pos";}
{
if($1==chr1 && $3==chr2 && $2>=pos && $2<(pos + 10000)){
print $0 " " sample
}
if($1==chr2 && $3==chr1 && $4>=pos && $4<(pos + 10000)){
print $0 " " sample
}
if($1==chr1 && $3==chr1 && (($4>=pos && $4<(pos + 10000)) || ($2>=pos && $2<(pos + 10000)))){
print $0 " " sample
}
}' >> ${OUTDIR_FOLDER_PATH}/vp_case

cat ${PATH_TO_CONTROL_VALID_PAIRS} | awk '
BEGIN{getline sample < "sample"; getline chr1 < "chr1"; getline chr2 < "chr2"; getline pos < "pos";}
{
if($1==chr1 && $3==chr2 && $2>=pos && $2<(pos + 10000)){
print $0 " " sample
}
if($1==chr2 && $3==chr1 && $4>=pos && $4<(pos + 10000)){
print $0 " " sample
}
if($1==chr1 && $3==chr1 && (($4>=pos && $4<(pos + 10000)) || ($2>=pos && $2<(pos + 10000)))){
print $0 " " sample
}
}' >> ${OUTDIR_FOLDER_PATH}/vp_control
fi

done

rm sample_inf
rm sample
rm chr1
rm chr2
rm pos

exit
