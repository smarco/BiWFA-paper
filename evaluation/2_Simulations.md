# Simulations

Generate simulated sequence pairs:

```shell
module load cmake gcc/10.2.0
RUN_GEN_DATASET=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/generate_dataset

mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/

for error in 0.001 0.01 0.05 0.10 0.20 0.40; do
  for len in `seq 100 100 1000; seq 2000 1000 10000; seq 20000 10000 1000000`; do
    echo $error $len $n
    seq 1 1 | while read n; do
      PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/seq_l${len}_e${error}_n${n}.seq

      $RUN_GEN_DATASET -l $len -e $error -n 1 > $PATH_SEQ
    done
  done
done
```

Prepare batches (for better parallelization):
```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/

# Divide in chunks to avoid too long jobs
ls /gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/*.seq | split -l 10 - chunk_
```


Align sequences:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_alignments
mkdir -p $DIR_OUTPUT
SET=sim

module load cmake gcc/10.2.0
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark

ls chunk_* | while read PATH_CHUNK; do
  echo $PATH_CHUNK

  sbatch -c 128 --exclusive --job-name biwfa-$SET      --wrap 'echo biwfa-'$SET';     cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.biwfa.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode biwfa  --check correct --output ${PREFIX}.biwfa.out          --repeat 100000 2>> ${PREFIX}.biwfa.log; done'
  sbatch -c 128 --exclusive --job-name wfa-high-$SET   --wrap 'echo wfa-high-'$SET';  cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-high.log;       \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high   --check correct --output ${PREFIX}.wfa-high.out       --repeat 100000 2>> ${PREFIX}.wfa-high.log; done'
done
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/
mkdir -p seq_statistics/

rm seq_statistics/statistics.simulations.tsv
rm seq_statistics/scores.simulations.tsv

for error in 0.001 0.01 0.05 0.10 0.20 0.40; do
  for len in `seq 100 100 1000; seq 2000 1000 10000; seq 20000 10000 1000000`; do
    seq 1 1 | while read n; do
      echo $error $len $n

      cat seq_alignments/seq_l${len}_e${error}_n${n}.*.log | python3 ../log2info.py | awk -v OFS='\t' -v error=$error -v len=$len '{print(error,len,$0)}' >> seq_statistics/statistics.simulations.tsv
      
      grep '^-' seq_alignments/seq_l${len}_e${error}_n${n}.*.out | uniq | cut -f 1 | cut -f 2 -d '/' | while read ROW; do
        NAME=$( echo $ROW | cut -f 1,2 -d '.' )
        MODE=$( echo $ROW | cut -f 3 -d '.' )
        SCORE=$( echo $ROW | cut -f 2 -d ':' )

        echo $error $len $NAME $MODE $SCORE | tr ' ' '\t' >> seq_statistics/scores.simulations.tsv
      done
    done
  done
done

pigz seq_statistics/statistics.simulations.tsv -f
pigz seq_statistics/scores.simulations.tsv -f
```
