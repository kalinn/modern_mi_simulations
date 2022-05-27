module load R
Rscript MImethods.R $LSB_JOBINDEX MCAR TRUE
Rscript MImethods.R $LSB_JOBINDEX MAR TRUE
Rscript MImethods.R $LSB_JOBINDEX MNAR1 TRUE
Rscript MImethods.R $LSB_JOBINDEX MNAR2 TRUE
