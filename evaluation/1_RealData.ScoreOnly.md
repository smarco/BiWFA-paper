# Real data

## ONT MinION

We use the sequence pairs prepared in the `1_RealData.md` workflow.

Align sequences:

```shell
module load cmake gcc/10.2.0

# Other reads (<= 500 kbps)
cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul_500kbps/

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/ont_ul_500kbps/seq_alignments_score_only
mkdir -p $DIR_OUTPUT

# NOTE: keep the single ', to avoid bash variables interpolation before sending the jobs

RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark
mkdir -p $DIR_OUTPUT/biwfa
sbatch -c 128 --exclusive --job-name biwfa-ont_ul    --wrap 'echo biwfa-ont_ul;    ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul_500kbps/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/biwfa/$NAME;          TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode ultralow --wfa-score-only --check correct --output ${PREFIX}.biwfa.out         2> ${PREFIX}.biwfa.log; done'
mkdir -p $DIR_OUTPUT/wfa-high
sbatch -c 128 --exclusive --job-name wfa-high-ont_ul --wrap 'echo wfa-high-ont_ul; ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul_500kbps/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/wfa-high/$NAME;       TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high     --wfa-score-only --check correct --output ${PREFIX}.wfa-high.out      2> ${PREFIX}.wfa-high.log; done'

RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark
mkdir -p $DIR_OUTPUT/edlib
sbatch -c 128 --exclusive --job-name edlib-ont_ul    --wrap 'echo edlib-ont_ul;    ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul_500kbps/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/edlib/$NAME;         TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                                       --check correct --output ${PREFIX}.edlib.out         2> ${PREFIX}.edlib.log; done'
mkdir -p $DIR_OUTPUT/bitpal-scored
sbatch -c 128 --exclusive --job-name bitpal-ont_ul   --wrap 'echo bitpal-ont_ul;   ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul_500kbps/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/bitpal-scored/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                                               --check correct --output ${PREFIX}.bitpal-scored.out 2> ${PREFIX}.bitpal-scored.log; done'
```

Collect statistics:

```shell
# Other reads (<= 500 kbps)
cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul_500kbps/
mkdir -p seq_statistics/

rm seq_statistics/statistics.score_only.tsv
rm seq_statistics/lengths.score_only.tsv
rm seq_statistics/scores.score_only.tsv

for mode in biwfa wfa-high edlib bitpal-scored; do
  echo $mode
  cat seq_alignments_score_only/$mode/*.log | python3 ../log2info.py | awk -v OFS='\t' -v SET='ONT_UL_OTHER_SCORE_ONLY' '{print(SET,$0)}' >> seq_statistics/statistics.score_only.tsv
  grep '^-' seq_alignments_score_only/$mode/*.out | cut -f 1 | cut -f 2 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET="ONT_UL_OTHER_SCORE_ONLY" '{print(SET,$0)}' >> seq_statistics/scores.score_only.tsv
done

ls seq_pairs/*.seq | while read PATH_SEQ; do
  NAME=$(basename $PATH_SEQ .seq)
  PREFIX=seq_alignments_score_only/$SET/$NAME
    
  cat $PATH_SEQ | tr '\n' ' ' | less -S | awk -v OFS='\t' -v SET="ONT_UL_OTHER_SCORE_ONLY" -v NAME=$NAME '{print(SET, NAME, length($1), length($2))}' >> seq_statistics/lengths.score_only.tsv
done
```

## Statistics

Merge all statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa

grep wfa ont_ul_500kbps/seq_statistics/statistics.tsv > statistics.alignment_vs_score_only.tsv
grep wfa ont_ul_500kbps/seq_statistics/statistics.score_only.tsv >> statistics.alignment_vs_score_only.tsv

cp ont_ul_500kbps/seq_statistics/lengths.tsv lengths.alignment_vs_score_only.tsv
```

Plot the statistics with the `scripts/plotsscore_only.R` script.
