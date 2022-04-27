# Real data

## ONT PromethION

Prepare the data:

```shell
# Create the folder on the cluster
mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens

# Move the data on the cluster
scp -r mat_*_regions bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens
scp -r pat_*_regions bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens
scp -r ont_*_fastas bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens
scp -r chr13_hsat1_*_fastas bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens
```

In the `hsapiens` dataset:
- `mat_hprc_regions` and `mat_chm_regions` folders contain maternal HG002 contigs from HPRC and a sequence from the 
CHM13 v1.1 assembly that they align to;
- `pat_hprc_regions` and `pat_chm_regions` folders contain maternal HG002 contigs from HPRC and a sequence from the 
CHM13 v1.1 assembly that they align to;
- `ont_mapping_ref_fastas` and `ont_read_fastas` folders contain promethION reads from HG002 and the sequences they map 
to in CHM13 v1.1;
- `chr13_hsat1_read_fastas` and `chr13_hsat1_ref_mapping_fastas` folders contain minION reads from CHM13 and the sequences 
they map to in the v1.1 assembly, subsetted to reads contained in the HSAT1 array on chromosome 13 (highly repetitive region).

The estimated identity is 95% for ONT data and 99.5% for the contigs.


Prepare sequence pairs:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens

mkdir -p seq_pairs/mat_regions/
ls mat_hprc_regions/mat_hprc_regions/*.fa | while read PATH_QUERY_FA; do
  NAME=$(basename $PATH_QUERY_FA .fa)
  echo $NAME
  
  PATH_TARGET_FA=mat_chm_regions/mat_chm_regions/$NAME.fa
  PATH_SEQ=seq_pairs/mat_regions/$NAME.seq

  grep '^>' -v $PATH_QUERY_FA | tr -d '\n' | sed 's/^/>/' | sed 's/$/\n/' > $PATH_SEQ
  grep '^>' -v $PATH_TARGET_FA | tr -d '\n' | sed 's/^/</' >> $PATH_SEQ
done

mkdir -p seq_pairs/pat_regions/
ls pat_hprc_regions/pat_hprc_regions/*.fa | while read PATH_QUERY_FA; do
  NAME=$(basename $PATH_QUERY_FA .fa)
  echo $NAME
  
  PATH_TARGET_FA=pat_chm_regions/pat_chm_regions/$NAME.fa
  PATH_SEQ=seq_pairs/pat_regions/$NAME.seq

  grep '^>' -v $PATH_QUERY_FA | tr -d '\n' | sed 's/^/>/' | sed 's/$/\n/' > $PATH_SEQ
  grep '^>' -v $PATH_TARGET_FA | tr -d '\n' | sed 's/^/</' >> $PATH_SEQ
done

mkdir -p seq_pairs/ont_regions/
ls ont_read_fastas/ont_read_fastas/*.fa | while read PATH_QUERY_FA; do
  NAME=$(basename $PATH_QUERY_FA .fa)
  echo $NAME
  
  PATH_TARGET_FA=ont_mapping_ref_fastas/ont_mapping_ref_fastas/$NAME.fa
  PATH_SEQ=seq_pairs/ont_regions/$NAME.seq

  grep '^>' -v $PATH_QUERY_FA | tr -d '\n' | sed 's/^/>/' | sed 's/$/\n/' > $PATH_SEQ
  grep '^>' -v $PATH_TARGET_FA | tr -d '\n' | sed 's/^/</' >> $PATH_SEQ
done

mkdir -p seq_pairs/hsat1_regions/
ls chr13_hsat1_read_fastas/chr13_hsat1_read_fastas/*.fa | while read PATH_QUERY_FA; do
  NAME=$(basename $PATH_QUERY_FA .fa)
  echo $NAME
  
  PATH_TARGET_FA=chr13_hsat1_ref_mapping_fastas/chr13_hsat1_ref_mapping_fastas/$NAME.fa
  PATH_SEQ=seq_pairs/hsat1_regions/$NAME.seq

  grep '^>' -v $PATH_QUERY_FA | tr -d '\n' | sed 's/^/>/' | sed 's/$/\n/' > $PATH_SEQ
  grep '^>' -v $PATH_TARGET_FA | tr -d '\n' | sed 's/^/</' >> $PATH_SEQ
done
```

Align sequences:

```shell
module load cmake gcc/10.2.0
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark

ls seq_pairs/ | while read SET; do
  echo $SET

  DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_alignments/$SET
  mkdir -p $DIR_OUTPUT

  sbatch -c 128 --exclusive --job-name biwfa-$SET      --wrap 'echo biwfa-'$SET';     ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-bidirectional    --check correct --output ${PREFIX}.biwfa.out          2> ${PREFIX}.biwfa.log; done'
  sbatch -c 128 --exclusive --job-name wfa-med-$SET    --wrap 'echo wfa-med-'$SET';   ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med  --check correct --output ${PREFIX}.wfa-med.out        2> ${PREFIX}.wfa-med.log; done'
  sbatch -c 128 --exclusive --job-name wfa-low-$SET    --wrap 'echo wfa-low--'$SET';  ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low  --check correct --output ${PREFIX}.wfa-low.out        2> ${PREFIX}.wfa-low.log; done'
  sbatch -c 128 --exclusive --job-name wfa-high-$SET   --wrap 'echo wfa-high-'$SET';  ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high --check correct --output ${PREFIX}.wfa-high.out       2> ${PREFIX}.wfa-high.log; done'
  
  sbatch -c 128 --exclusive --job-name wfalm-$SET      --wrap 'echo wfalm-'$SET';     ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                 --check correct --output ${PREFIX}.wfalm.out          2> ${PREFIX}.wfalm.log; done'
  sbatch -c 128 --exclusive --job-name wfalm-low-$SET  --wrap 'echo wfalm-low-'$SET'; ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                          --check correct --output ${PREFIX}.wfalm-lowmem.out   2> ${PREFIX}.wfalm-lowmem.log; done'
  sbatch -c 128 --exclusive --job-name wfalm-rec-$SET  --wrap 'echo wfalm-rec-'$SET'; ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-rec                             --check correct --output ${PREFIX}.wfalm-rec.out      2> ${PREFIX}.wfalm-rec.log; done'

  sbatch -c 128 --exclusive --job-name edlib-$SET      --wrap 'echo edlib-'$SET';     ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                 --check correct --output ${PREFIX}.edlib.out          2> ${PREFIX}.edlib.log; done'
  sbatch -c 128 --exclusive --job-name bitpal-$SET     --wrap 'echo bitpal-'$SET';    ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                         --check correct --output ${PREFIX}.bitpal-scored.out  2> ${PREFIX}.bitpal-scored.log; done'
  sbatch -c 128 --exclusive --job-name ksw2-$SET       --wrap 'echo ksw2-'$SET';      ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                        --check correct --output ${PREFIX}.ksw2-extz2-sse.out 2> ${PREFIX}.ksw2-extz2-sse.log; done'
done
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens
mkdir -p seq_statistics/

rm seq_statistics/statistics.tsv
rm seq_statistics/scores.tsv
rm seq_statistics/lengths.tsv

ls seq_pairs/ | while read SET; do
  echo $SET
 
  cat seq_alignments/$SET/*.log | python3 ../log2info.py | awk -v OFS='\t' -v SET=$SET '{print(SET,$0)}' >> seq_statistics/statistics.tsv
  grep '^-' seq_alignments/$SET/*.out | cut -f 1 | cut -f 3 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0)}' >> seq_statistics/scores.tsv
  
  ls seq_pairs/$SET/*.seq | while read PATH_SEQ; do
    NAME=$(basename $PATH_SEQ .seq)
    PREFIX=seq_alignments/$SET/$NAME
    
    cat $PATH_SEQ | tr '\n' ' ' | less -S | awk -v OFS='\t' -v SET=$SET -v NAME=$NAME '{print(SET, NAME, length($1), length($2))}' >> seq_statistics/lengths.tsv

#    grep 'Time.Alignment' $PREFIX.*.log | grep call -v | cut -f 3 -d '/' | cut -f 1,2,4,5 -d '.' | tr -s ' ' | sed 's/ (s)//' | sed 's/.Alignment//' | sed 's/\./ /' | tr ' ' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"time_s")}' >> seq_statistics/statistics.tsv
#    grep 'Maximum resident' $PREFIX.*.log  | cut -f 3 -d '/' | sed 's/.log://' | sed 's/Maximum resident set size (kbytes): //g' | tr '.' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"memory_kb")}' >> seq_statistics/statistics.tsv
#    grep '^-' $PREFIX.*.out | cut -f 1 | cut -f 3 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"memory_kb")}' >> seq_statistics/scores.tsv
  done
done
```


## ONT MinION > 500 kbps

Prepare sequence pairs:

```shell
mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs

# Long reads
c=1
cat /gpfs/projects/bsc18/bsc18571/wfa2/datasets/real/Nanopore.UL.100K.seq | awk 'length($1) >= 500000' | while read -r line1; do
  read -r line2

  echo $line1 > /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/seq$c.seq
  echo $line2 >> /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/seq$c.seq
        
  c=$((c + 1))
done
```

You can find the sequence pairs in the `~/BiWFA-paper/evaluation/data/ONT_MinION_UL.500kbps.zip` file.

Align sequences:

```shell
module load cmake gcc/10.2.0
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark

cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_alignments
mkdir -p $DIR_OUTPUT

ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do
  NAME=$(basename $PATH_SEQ .seq);
  PREFIX=$DIR_OUTPUT/$NAME;
  echo $NAME;
  
  sbatch -c 128 --exclusive --job-name biwfa-ont_ul     --wrap 'echo biwfa-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-bidirectional    --check correct --output '${PREFIX}'.biwfa.out          2> '${PREFIX}'.biwfa.log'
  sbatch -c 128 --exclusive --job-name wfa-med-ont_ul   --wrap 'echo wfa-med-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med  --check correct --output '${PREFIX}'.wfa-med.out            2> '${PREFIX}'.wfa-med.log'
  sbatch -c 128 --exclusive --job-name wfa-low-ont_ul   --wrap 'echo wfa-low-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low  --check correct --output '${PREFIX}'.wfa-low.out            2> '${PREFIX}'.wfa-low.log'
  sbatch -c 128 --exclusive --job-name wfa-high-ont_ul  --wrap 'echo wfa-high-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high --check correct --output '${PREFIX}'.wfa-high.out           2> '${PREFIX}'.wfa-high.log'
    
  sbatch -c 128 --exclusive --job-name wfalm-ont_ul     --wrap 'echo wfalm-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm                                 --check correct --output '${PREFIX}'.wfalm.out          2> '${PREFIX}'.wfalm.log'
  sbatch -c 128 --exclusive --job-name wfalm-low-ont_ul --wrap 'echo wfalm-low-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm-lowmem                          --check correct --output '${PREFIX}'.wfalm-lowmem.out   2> '${PREFIX}'.wfalm-lowmem.log'
  sbatch -c 128 --exclusive --job-name wfalm-rec-ont_ul --wrap 'echo wfalm-rec-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm-rec                             --check correct --output '${PREFIX}'.wfalm-rec.out      2> '${PREFIX}'.wfalm-rec.log'
  
  sbatch -c 128 --exclusive --job-name edlib-ont_ul     --wrap 'echo edlib-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a edlib                                 --check correct --output '${PREFIX}'.edlib.out          2> '${PREFIX}'.edlib.log'
  sbatch -c 128 --exclusive --job-name bitpal-ont_ul    --wrap 'echo bitpal-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a bitpal-scored                         --check correct --output '${PREFIX}'.bitpal-scored.out  2> '${PREFIX}'.bitpal-scored.log'
  sbatch -c 128 --exclusive --job-name ksw2-ont_ul      --wrap 'echo wfalm-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                        --check correct --output '${PREFIX}'.ksw2-extz2-sse.out 2> '${PREFIX}'.ksw2-extz2-sse.log'
done
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/
mkdir -p seq_statistics/

rm seq_statistics/lengths.tsv

cat seq_alignments/*.log | python3 ../log2info.py | awk -v OFS='\t' -v SET='ONT_UL' '{print(SET,$0)}' > seq_statistics/statistics.tsv
grep '^-' seq_alignments/*.out | cut -f 1 | cut -f 2 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET="ONT_UL" '{print(SET,$0)}' > seq_statistics/scores.tsv

ls seq_pairs/*.seq | while read PATH_SEQ; do
  NAME=$(basename $PATH_SEQ .seq)
  PREFIX=seq_alignments/$SET/$NAME
    
  cat $PATH_SEQ | tr '\n' ' ' | less -S | awk -v OFS='\t' -v SET="ONT_UL" -v NAME=$NAME '{print(SET, NAME, length($1), length($2))}' >> seq_statistics/lengths.tsv
done
```


## PacBio HiFi

Prepare sequence pairs:

```shell
mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs

# First 1000 pairs
c=1
head -n 2000 /gpfs/projects/bsc18/bsc18571/wfa2/datasets/real/PacBio.HF.1M.seq | while read -r line1; do
    read -r line2
    
    echo $line1 > /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/seq$c.seq
    echo $line2 >> /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/seq$c.seq
    
    c=$((c + 1))
done
```

Align sequences:

```shell
module load cmake gcc/10.2.0
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark

cd /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_alignments
mkdir -p $DIR_OUTPUT

sbatch -c 48 --job-name biwfa-pacbio      --wrap 'echo biwfa-ont_ul;     ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-bidirectional    --check correct --output ${PREFIX}.biwfa.out          2> ${PREFIX}.biwfa.log; done'
sbatch -c 48 --job-name wfa-med-pacbio    --wrap 'echo wfa-med-ont_ul;   ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med  --check correct --output ${PREFIX}.wfa-med.out        2> ${PREFIX}.wfa-med.log; done'
sbatch -c 48 --job-name wfa-low-pacbio    --wrap 'echo wfa-low-ont_ul;   ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low  --check correct --output ${PREFIX}.wfa-low.out        2> ${PREFIX}.wfa-low.log; done'
sbatch -c 48 --job-name wfa-high-pacbio   --wrap 'echo wfa-high-ont_ul;  ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high --check correct --output ${PREFIX}.wfa-high.out       2> ${PREFIX}.wfa-high.log; done'

sbatch -c 48 --job-name wfalm-pacbio      --wrap 'echo wfalm-ont_ul;     ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                 --check correct --output ${PREFIX}.wfalm.out          2> ${PREFIX}.wfalm.log; done'
sbatch -c 48 --job-name wfalm-low-pacbio  --wrap 'echo wfalm-low-ont_ul; ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                          --check correct --output ${PREFIX}.wfalm-lowmem.out   2> ${PREFIX}.wfalm-lowmem.log; done'

sbatch -c 48 --job-name edlib-pacbio      --wrap 'echo edlib-ont_ul;     ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                 --check correct --output ${PREFIX}.edlib.out          2> ${PREFIX}.edlib.log; done'
sbatch -c 48 --job-name bitpal-pacbio     --wrap 'echo bitpal-ont_ul;    ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                         --check correct --output ${PREFIX}.bitpal-scored.out  2> ${PREFIX}.bitpal-scored.log; done'
sbatch -c 48 --job-name ksw2-pacbio       --wrap 'echo ksw2-ont_ul;      ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                        --check correct --output ${PREFIX}.ksw2-extz2-sse.out 2> ${PREFIX}.ksw2-extz2-sse.log; done'
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/
mkdir -p seq_statistics/
   
cat seq_alignments/*.log | python3 ../log2info.py | awk -v OFS='\t' -v SET='ONT_UL' '{print(SET,$0)}' > seq_statistics/statistics.tsv
grep '^-' seq_alignments/*.out | cut -f 1 | cut -f 2 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET="ONT_UL" '{print(SET,$0)}' > seq_statistics/scores.tsv

ls seq_pairs/*.seq | while read PATH_SEQ; do
  NAME=$(basename $PATH_SEQ .seq)
  PREFIX=seq_alignments/$SET/$NAME
    
  cat $PATH_SEQ | tr '\n' ' ' | less -S | awk -v OFS='\t' -v SET="PACBIO" -v NAME=$NAME '{print(SET, NAME, length($1), length($2))}' >> seq_statistics/lengths.tsv
done
```

## Statistics

Merge all statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa

cat hsapiens/seq_statistics/statistics.tsv > statistics_all.tsv
cat ont_ul/seq_statistics/statistics.tsv >> statistics_all.tsv
cat hsapiens/seq_statistics/scores.tsv > scores_all.tsv
cat ont_ul/seq_statistics/scores.tsv >> scores_all.tsv
cat hsapiens/seq_statistics/lengths.tsv > lengths_all.tsv
cat ont_ul/seq_statistics/lengths.tsv >> lengths_all.tsv
```

Plot the statistics with the `scripts/plots.R` script.