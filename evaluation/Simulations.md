# Simulations

Generate simulated sequence pairs:

```shell
module load cmake gcc/10.2.0
RUN_GEN_DATASET=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/generate_dataset

mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/

for error in 0.01 0.05 0.10 0.20 0.30; do
  for len in 50000 100000 200000 500000 1000000 2000000; do
    seq 1 10 | while read n; do
      echo $error $len $n
      PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/seq_l${len}_e${error}_n${n}.seq

      $RUN_GEN_DATASET -l $len -e $error -n 1 > $PATH_SEQ
    done
  done
done
```

Align sequences:

```shell
module load cmake gcc/10.2.0
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_alignments
mkdir -p $DIR_OUTPUT

for error in 0.01 0.05 0.10 0.20 0.30; do
  for len in 50000 100000 200000 500000 1000000 2000000; do    
    seq 1 10 | while read n; do
      echo $error $len $n
      PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/seq_l${len}_e${error}_n${n}.seq
      NAME=$(basename $PATH_SEQ .seq);
      PREFIX=$DIR_OUTPUT/$NAME;
      echo $NAME;
  
      sbatch -c 128 --exclusive --job-name biwfa     --wrap 'echo biwfa;     echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-bidirectional    --check correct --output '${PREFIX}'.biwfa.out          2> '${PREFIX}'.biwfa.log'
      sbatch -c 128 --exclusive --job-name wfa-med   --wrap 'echo wfa-med;   echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med  --check correct --output '${PREFIX}'.wfa-med.out        2> '${PREFIX}'.wfa-med.log'
      sbatch -c 128 --exclusive --job-name wfa-low   --wrap 'echo wfa-low;   echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low  --check correct --output '${PREFIX}'.wfa-low.out        2> '${PREFIX}'.wfa-low.log'
      sbatch -c 128 --exclusive --job-name wfa-high  --wrap 'echo wfa-high;  echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high --check correct --output '${PREFIX}'.wfa-high.out       2> '${PREFIX}'.wfa-high.log'
        
      sbatch -c 128 --exclusive --job-name wfalm     --wrap 'echo wfalm;     echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm                                 --check correct --output '${PREFIX}'.wfalm.out          2> '${PREFIX}'.wfalm.log'
      sbatch -c 128 --exclusive --job-name wfalm-low --wrap 'echo wfalm-low; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm-lowmem                          --check correct --output '${PREFIX}'.wfalm-lowmem.out   2> '${PREFIX}'.wfalm-lowmem.log'
      sbatch -c 128 --exclusive --job-name wfalm-rec --wrap 'echo wfalm-rec; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm-rec                             --check correct --output '${PREFIX}'.wfalm-rec.out      2> '${PREFIX}'.wfalm-rec.log'
      
      sbatch -c 128 --exclusive --job-name edlib     --wrap 'echo edlib;     echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a edlib                                 --check correct --output '${PREFIX}'.edlib.out          2> '${PREFIX}'.edlib.log'
      sbatch -c 128 --exclusive --job-name bitpal    --wrap 'echo bitpal;    echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a bitpal-scored                         --check correct --output '${PREFIX}'.bitpal-scored.out  2> '${PREFIX}'.bitpal-scored.log'
      sbatch -c 128 --exclusive --job-name ksw2      --wrap 'echo wfalm;     echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                        --check correct --output '${PREFIX}'.ksw2-extz2-sse.out 2> '${PREFIX}'.ksw2-extz2-sse.log'
    done
  done
done
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations/
mkdir -p seq_statistics/

rm seq_statistics/statistics.simulations.tsv
rm seq_statistics/scores.simulations.tsv

for error in 0.01 0.05 0.10 0.20 0.30; do
  for len in 50000 100000 200000 500000 1000000 2000000; do 
    echo $error $len
    seq 1 10 | while read n; do
      PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/simulations/seq_pairs/seq_l${len}_e${error}_n${n}.seq
      NAME=$(basename $PATH_SEQ .seq);
     
      cat seq_alignments/${NAME}.*.log | python3 ../log2info.py | awk -v OFS='\t' -v error=$error -v len=$len '{print(error,len,$0)}' >> seq_statistics/statistics.simulations.tsv
      
      grep '^-' seq_alignments/${NAME}.*.out | cut -f 1 | cut -f 2 -d '/' | while read ROW; do
        MODE=$( echo $ROW | cut -f 3 -d '.' )
        SCORE=$( echo $ROW | cut -f 2 -d ':' )

        echo $error $len $NAME $MODE $SCORE | tr ' ' '\t' >> seq_statistics/scores.simulations.tsv
      done
    done
  done
done
```
