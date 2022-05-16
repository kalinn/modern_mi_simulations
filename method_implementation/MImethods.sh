module load R
Rscript MImethods.R $LSB_JOBINDEX MCAR
Rscript MImethods.R $LSB_JOBINDEX MAR
Rscript MImethods.R $LSB_JOBINDEX MNAR1
Rscript MImethods.R $LSB_JOBINDEX MNAR2
