# Simulations

Generate simulated sequence pairs:

```shell
module load cmake gcc/10.2.0

mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/simulations_table/
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations_table/

mkdir -p datasets

RUN_GEN_DATASET=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/generate_dataset
for error in 0.40 0.20 0.10 0.05 0.01 0.001; do
  $RUN_GEN_DATASET -o datasets/sim.l100.n1M.$error.seq  -n 1000000 -l 100     -e $error
  $RUN_GEN_DATASET -o datasets/sim.l1K.n100K.$error.seq -n 100000  -l 1000    -e $error
  $RUN_GEN_DATASET -o datasets/sim.l10K.n1K.$error.seq  -n 1000    -l 10000   -e $error
  $RUN_GEN_DATASET -o datasets/sim.l100K.n10.$error.seq -n 10      -l 100000  -e $error
  $RUN_GEN_DATASET -o datasets/sim.l1M.n1.$error.seq    -n 1       -l 1000000 -e $error
  $RUN_GEN_DATASET -o datasets/sim.l2M.n1.$error.seq    -n 1       -l 2000000 -e $error
done
```

Align sequences:

```shell
module load cmake gcc/10.2.0

DIR_RESULTS=/gpfs/projects/bsc18/bsc18995/biwfa/simulations_table/results
mkdir -p $DIR_RESULTS

for input in "sim.l100.n1M" "sim.l1K.n100K" "sim.l10K.n1K" "sim.l100K.n10" "sim.l1M.n1" "sim.l2M.n1"; do
  for error in 0.40 0.20 0.10 0.05 0.01 0.001
    PREFIX=${DIR_RESULTS}/$input.$error
    PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/simulations_table/datasets/$input.$error.seq
    
    RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark
    sbatch -c 128 --exclusive --job-name Bitpal          --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a bitpal-scored                                              -i '${PATH_SEQ}' -o '$PREFIX'.bitpal.alg      &> '$PREFIX'.bitpal.log'
    sbatch -c 128 --exclusive --job-name edlib           --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a edlib                                                      -i '${PATH_SEQ}' -o '$PREFIX'.edlib.alg       &> '$PREFIX'.edlib.log'
    sbatch -c 128 --exclusive --job-name ksw2-extz2-sse  --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a ksw2-extz2-sse                                             -i '${PATH_SEQ}' -o '$PREFIX'.ksw2.alg        &> '$PREFIX'.ksw2.log'

    sbatch -c 128 --exclusive --job-name wfalm           --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a wfalm                                                      -i '${PATH_SEQ}' -o '$PREFIX'.wfalm.alg       &> '$PREFIX'.wfalm.log'
    sbatch -c 128 --exclusive --job-name wfalm-lowmem    --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a wfalm-lowmem                                               -i '${PATH_SEQ}' -o '$PREFIX'.wfalmlow.alg    &> '$PREFIX'.wfalmlow.log'
    sbatch -c 128 --exclusive --job-name wfalm-rec       --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a wfalm-rec                                                  -i '${PATH_SEQ}' -o '$PREFIX'.wfalmrec.alg    &> '$PREFIX'.wfalmrec.log'
    
    RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark
    sbatch -c 128 --exclusive --job-name wfa-high        --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode=high                      -i '${PATH_SEQ}' -o '$PREFIX'.wfahigh.alg     &> '$PREFIX'.wfahigh.log'
    sbatch -c 128 --exclusive --job-name wfa-med         --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode=med                       -i '${PATH_SEQ}' -o '$PREFIX'.wfamed.alg      &> '$PREFIX'.wfamed.log'
    sbatch -c 128 --exclusive --job-name wfa-low         --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode=low                       -i '${PATH_SEQ}' -o '$PREFIX'.wfalow.alg      &> '$PREFIX'.wfalow.log'

    sbatch -c 128 --exclusive --job-name biwfa           --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode=ultralow                  -i '${PATH_SEQ}' -o '$PREFIX'.biwfa.alg       &> '$PREFIX'.biwfa.log'
    sbatch -c 128 --exclusive --job-name biwfa           --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode=ultralow --wfa-score-only -i '${PATH_SEQ}' -o '$PREFIX'.biwfa.score.alg &> '$PREFIX'.biwfa.score.log'
  done
done
```

Collect statistics:

```shell
for input in "sim.l100.n1M" "sim.l1K.n100K" "sim.l10K.n1K" "sim.l100K.n10" "sim.l1M.n1" "sim.l2M.n1"; do
  for error in 0.40 0.20 0.10 0.05 0.01 0.001; do
    #for algo in bitpal edlib ksw2 wfalm wfalmlow wfalmrec wfahigh wfamed wfalow biwfa biwfa.score; do
    for algo in wfahigh biwfa biwfa.score; do
      FILE=results/${input}.$error.$algo.log
      TIME=$(cat $FILE | grep "Time.Alignment" | awk '{print $3}')
      EXIT=$(cat $FILE | grep "Exit status:" | awk '{print $3}')
      MEM=$(cat $FILE | grep "Maximum resident set size" | awk '{print $6}')
    
      echo "$input $error $algo $EXIT $TIME $MEM"
    done
  done
done | column -t
```
