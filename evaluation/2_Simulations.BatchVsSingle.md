# Simulations

Generate simulated sequence pairs:

```shell
module load cmake gcc/10.2.0

mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/simulations_batch_vs_single/
cd /gpfs/projects/bsc18/bsc18995/biwfa/simulations_batch_vs_single/

mkdir -p datasets

RUN_GEN_DATASET=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/generate_dataset
for error in 0.01; do
  $RUN_GEN_DATASET -o datasets/sim.l500K.n100.$error.seq -n 100 -l 500000 -e $error
done


mkdir -p datasets_sep/

for input in sim.l500K.n100; do
  for error in 0.01; do
    echo $input $error

    PATH_SEQ=datasets/$input.$error.seq

    N=$(echo "$(wc -l $PATH_SEQ | cut -f 1 -d ' ')" | bc)

    seq 1 2 $N | while read i; do
      j=$(echo "$i + 1" | bc)
      sed -n "${i},${j}p" $PATH_SEQ > datasets_sep/$input.$error.${i}_${j}.seq
    done
  done
done
```

Align sequences:

```shell
module load cmake gcc/10.2.0

DIR_RESULTS=/gpfs/projects/bsc18/bsc18995/biwfa/simulations_batch_vs_single/results
mkdir -p $DIR_RESULTS

for input in sim.l500K.n100; do
  for error in 0.01; do
    PREFIX=${DIR_RESULTS}/$input.$error
    PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/simulations_batch_vs_single/datasets/$input.$error.seq
    
    RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark
    sbatch -c 128 --exclusive --job-name wfa-high        --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode high                      -i '${PATH_SEQ}' -o '$PREFIX'.wfahigh.alg     &> '$PREFIX'.wfahigh.log'
    sbatch -c 128 --exclusive --job-name biwfa           --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode ultralow                  -i '${PATH_SEQ}' -o '$PREFIX'.biwfa.alg       &> '$PREFIX'.biwfa.log'
  done
done


DIR_RESULTS=/gpfs/projects/bsc18/bsc18995/biwfa/simulations_batch_vs_single/results_sep
mkdir -p $DIR_RESULTS

for input in sim.l500K.n100; do
  for error in 0.01; do
    ls /gpfs/projects/bsc18/bsc18995/biwfa/simulations_batch_vs_single/datasets_sep/$input.$error.*.seq | while read PATH_SEQ; do
      NAME=$(basename $PATH_SEQ .seq);
      PREFIX=${DIR_RESULTS}/$NAME
    
      RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark
      sbatch -c 128 --exclusive --job-name wfa-high        --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode high                      -i '${PATH_SEQ}' -o '$PREFIX'.wfahigh.alg     &> '$PREFIX'.wfahigh.log'
      sbatch -c 128 --exclusive --job-name biwfa           --time 48 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa --wfa-memory-mode ultralow                  -i '${PATH_SEQ}' -o '$PREFIX'.biwfa.alg       &> '$PREFIX'.biwfa.log'
    done
  done
done
```

Collect statistics:

```shell
for input in sim.l500K.n100; do
  for error in 0.01; do
    cat results/*.log | python3 ../log2info.ns.py > batch.statistics.tsv
  done
done

for input in sim.l500K.n100; do
  for error in 0.01; do
    cat results_sep/*.log | python3 ../log2info.ns.py > single.statistics.tsv
  done
done
```
