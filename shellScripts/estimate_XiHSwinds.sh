#!/bin/bash
#$ -N XiHSwind
#$ -wd /ndata/home/genomeindia/gi-advance/gi3/dipto/aNalysis/selectSign/selscan/resultS/iHS/XiHSwinds/
#$ -pe smp 36
#$ -l h_vmem=500G
#$ -o estXiHSwind_SASfull.out
#$ -e estXiHSwind_SASfull.error

module load R-3.6.3 

Rscript /ndata/home/genomeindia/gi-advance/gi3/dipto/aNalysis/selectSign/selSweeps/iHSres_postProc/iHS_plots.R
