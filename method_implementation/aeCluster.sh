module load python/3.8
python aeCluster.py $LSB_JOBINDEX MCAR TRUE
python aeCluster.py $LSB_JOBINDEX MAR TRUE
python aeCluster.py $LSB_JOBINDEX MNAR1 TRUE
python aeCluster.py $LSB_JOBINDEX MNAR2 TRUE
