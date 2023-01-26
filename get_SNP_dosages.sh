#!/bin/bash


## assign directories used for the analysis
BGENDIR=~/insert_path/
PHENOFILE=~/insert_path/file.sample
BGENIX=~/insert_path/Programs/bgen_tools/bgenix

## define coordinates 
olink=${1}
chr=${2}
LOWPOS=${3}
UPPOS=${4}

echo "Protein ${olink} : Chromosome ${chr} : Locus start ${LOWPOS} : Locus end ${UPPOS}"

#------------------------------#
## --> file for LD-matrix <-- ##
#------------------------------#

## Create temporary bgen file
$BGENIX -g ${BGENDIR}/chr${chr}.bgen \
-incl-rsids tmpdir/tmp.${olink}.${chr}.${LOWPOS}.${UPPOS}.lst > tmpdir/${olink}.${chr}.${LOWPOS}.${UPPOS}.bgen

## convert to dosage file
~/insert_path/Programs/qctool_v2.0.2/qctool_v2.0.2 \
-g tmpdir/${olink}.${chr}.${LOWPOS}.${UPPOS}.bgen \
-s ${PHENOFILE} \
-og - \
-ofiletype dosage > tmpdir/tmp.${olink}.${chr}.${LOWPOS}.${UPPOS}.dosage

rm tmpdir/${olink}.${chr}.${LOWPOS}.${UPPOS}.bgen

