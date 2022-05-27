module load python/3.8
python aeLowDim.py $LSB_JOBINDEX MCAR 
python aeLowDim.py $LSB_JOBINDEX MAR 
python aeLowDim.py $LSB_JOBINDEX MNAR1 
python aeLowDim.py $LSB_JOBINDEX MNAR2 
