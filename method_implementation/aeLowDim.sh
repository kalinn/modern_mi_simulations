module load python/3.8
python aeLowDim.py $LSB_JOBINDEX MCAR 10 
python aeLowDim.py $LSB_JOBINDEX MAR 10 
python aeLowDim.py $LSB_JOBINDEX MNAR1 10
python aeLowDim.py $LSB_JOBINDEX MNAR2 10
python aeLowDim.py $LSB_JOBINDEX MCAR 30
python aeLowDim.py $LSB_JOBINDEX MAR 30
python aeLowDim.py $LSB_JOBINDEX MNAR1 30
python aeLowDim.py $LSB_JOBINDEX MNAR2 30
python aeLowDim.py $LSB_JOBINDEX MCAR 50
python aeLowDim.py $LSB_JOBINDEX MAR 50
python aeLowDim.py $LSB_JOBINDEX MNAR1 50
python aeLowDim.py $LSB_JOBINDEX MNAR2 50
