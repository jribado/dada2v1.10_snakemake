source activate /home/jribado/miniconda3/envs/dada2v1.6/
snakemake -p \
--jobs 100 \
--configfile config.yaml \
--cluster-config clusterconfig.yaml \
--restart-times 0 \
-s Snakefile \
--profile scg
source deactivate
