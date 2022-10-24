# Simulations

Generate simulated sequence pairs:

```shell
module load cmake gcc/10.2.0
RUN_GEN_DATASET=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/generate_dataset

mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/

for error in 0.001 0.01 0.05 0.10 0.20 0.30; do
  for len in `seq 250 10 20000`; do
    echo $error $len $n
    seq 1 1 | while read n; do
      PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/seq_l${len}_e${error}_n${n}.seq

      $RUN_GEN_DATASET -l $len -e $error -n 1 > $PATH_SEQ
    done
  done
done
```

Align sequences:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/

module load cmake gcc/10.2.0
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_alignments
mkdir -p $DIR_OUTPUT

sbatch -c 128 --exclusive --job-name biwfa-sim    --wrap 'echo biwfa-sim;    ls /gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.biwfa.log;    \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode ultralow --check correct --output ${PREFIX}.biwfa.out    2>> ${PREFIX}.$i.biwfa.log; done; done'
sbatch -c 128 --exclusive --job-name wfa-high-sim --wrap 'echo wfa-high-sim; ls /gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.wfa-high.log; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high     --check correct --output ${PREFIX}.wfa-high.out 2>> ${PREFIX}.$i.wfa-high.log; done; done'
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/
mkdir -p seq_statistics/

rm seq_statistics/statistics.simulations.tsv
rm seq_statistics/scores.simulations.tsv

for error in 0.001 0.01 0.05 0.10 0.20 0.30; do
  for len in `seq 250 10 20000`; do
    echo $error $len $n
    seq 1 1 | while read n; do
      PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/seq_l${len}_e${error}_n${n}.seq
      NAME=$(basename $PATH_SEQ .seq);
     
      cat seq_alignments/$NAME/*.log | python3 ../log2info.ns.py | awk -v OFS='\t' -v error=$error -v len=$len '{print(error,len,$0)}' >> seq_statistics/statistics.simulations.tsv
      
      grep '^-' seq_alignments/$NAME/*.out | cut -f 1 | cut -f 3 -d '/' | while read ROW; do
        MODE=$( echo $ROW | cut -f 3 -d '.' )
        SCORE=$( echo $ROW | cut -f 2 -d ':' )

        echo $error $len $NAME $MODE $SCORE | tr ' ' '\t' >> seq_statistics/scores.simulations.tsv
      done
     done
  done
done
```

pigz seq_statistics/statistics.simulations.tsv
pigz seq_statistics/scores.simulations.tsv
