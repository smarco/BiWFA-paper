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
ls datasets/mat_hprc_regions/mat_hprc_regions/*.fa | while read PATH_QUERY_FA; do
  NAME=$(basename $PATH_QUERY_FA .fa)
  echo $NAME
  
  PATH_TARGET_FA=datasets/mat_chm_regions/mat_chm_regions/$NAME.fa
  PATH_SEQ=seq_pairs/mat_regions/$NAME.seq

  grep '^>' -v $PATH_QUERY_FA | tr -d '\n' | sed 's/^/>/' | sed 's/$/\n/' > $PATH_SEQ
  grep '^>' -v $PATH_TARGET_FA | tr -d '\n' | sed 's/^/</' >> $PATH_SEQ
  
  #A_LEN=$(head $PATH_SEQ -n 1 | awk '{ print length-1 }')
  #B_LEN=$(tail $PATH_SEQ -n 1 | awk '{ print length-1 }')
  #if [ "$A_LEN" -le 10000 ] && [ "$B_LEN" -le 10000 ]; then
  #  cp $PATH_SEQ seq_pairs/mat_regions_10kbps/$NAME.seq
  #fi
done

mkdir -p seq_pairs/pat_regions/
ls datasets/pat_hprc_regions/pat_hprc_regions/*.fa | while read PATH_QUERY_FA; do
  NAME=$(basename $PATH_QUERY_FA .fa)
  echo $NAME
  
  PATH_TARGET_FA=datasets/pat_chm_regions/pat_chm_regions/$NAME.fa
  PATH_SEQ=seq_pairs/pat_regions/$NAME.seq

  grep '^>' -v $PATH_QUERY_FA | tr -d '\n' | sed 's/^/>/' | sed 's/$/\n/' > $PATH_SEQ
  grep '^>' -v $PATH_TARGET_FA | tr -d '\n' | sed 's/^/</' >> $PATH_SEQ
  
  #A_LEN=$(head $PATH_SEQ -n 1 | awk '{ print length-1 }')
  #B_LEN=$(tail $PATH_SEQ -n 1 | awk '{ print length-1 }')
  #if [ "$A_LEN" -le 10000 ] && [ "$B_LEN" -le 10000 ]; then
  #  cp $PATH_SEQ seq_pairs/pat_regions_10kbps/$NAME.seq
  #fi
done

mkdir -p seq_pairs/ont_regions/
ls datasets/ont_read_fastas/ont_read_fastas/*.fa | while read PATH_QUERY_FA; do
  NAME=$(basename $PATH_QUERY_FA .fa)
  echo $NAME
  
  PATH_TARGET_FA=datasets/ont_mapping_ref_fastas/ont_mapping_ref_fastas/$NAME.fa
  PATH_SEQ=seq_pairs/ont_regions/$NAME.seq

  grep '^>' -v $PATH_QUERY_FA | tr -d '\n' | sed 's/^/>/' | sed 's/$/\n/' > $PATH_SEQ
  grep '^>' -v $PATH_TARGET_FA | tr -d '\n' | sed 's/^/</' >> $PATH_SEQ
  
  #A_LEN=$(head $PATH_SEQ -n 1 | awk '{ print length-1 }')
  #B_LEN=$(tail $PATH_SEQ -n 1 | awk '{ print length-1 }')
  #if [ "$A_LEN" -le 10000 ] && [ "$B_LEN" -le 10000 ]; then
  #  cp $PATH_SEQ seq_pairs/ont_regions_10kbps/$NAME.seq
  #fi
done

mkdir -p seq_pairs/hsat1_regions/
ls datasets/chr13_hsat1_read_fastas/chr13_hsat1_read_fastas/*.fa | while read PATH_QUERY_FA; do
  NAME=$(basename $PATH_QUERY_FA .fa)
  echo $NAME
  
  PATH_TARGET_FA=datasets/chr13_hsat1_ref_mapping_fastas/chr13_hsat1_ref_mapping_fastas/$NAME.fa
  PATH_SEQ=seq_pairs/hsat1_regions/$NAME.seq

  grep '^>' -v $PATH_QUERY_FA | tr -d '\n' | sed 's/^/>/' | sed 's/$/\n/' > $PATH_SEQ
  grep '^>' -v $PATH_TARGET_FA | tr -d '\n' | sed 's/^/</' >> $PATH_SEQ
  
  #A_LEN=$(head $PATH_SEQ -n 1 | awk '{ print length-1 }')
  #B_LEN=$(tail $PATH_SEQ -n 1 | awk '{ print length-1 }')
  #if [ "$A_LEN" -le 10000 ] && [ "$B_LEN" -le 10000 ]; then
  #  cp $PATH_SEQ seq_pairs/hsat1_regions_10kbps/$NAME.seq
  #fi
done
```

Prepare batches (for better parallelization):

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens

# Divide in chunks to avoid too long jobs
ls seq_pairs/ | grep ont | while read SET; do
  echo $SET
  ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/$SET/*.seq | split -l 200 - ${SET}_chunk_
done
```

Align sequences:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens

module load cmake gcc/10.2.0
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark

ls seq_pairs/ | grep ont | while read SET; do
  echo $SET
  
  DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_alignments/$SET
  mkdir -p $DIR_OUTPUT

  ls ${SET}_chunk_* | while read PATH_CHUNK; do
    echo $PATH_CHUNK

    # Write just one .out file, not one for each replicate
    # NOTE: keep the single ', to avoid bash variables interpolation before sending the jobs

    sbatch -c 128 --exclusive --job-name biwfa-$SET      --wrap 'echo biwfa-'$SET';     cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.biwfa.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode biwfa  --check correct --output ${PREFIX}.biwfa.out          --repeat 100000 2>> ${PREFIX}.biwfa.log; done'
    sbatch -c 128 --exclusive --job-name wfa-med-$SET    --wrap 'echo wfa-med-'$SET';   cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-med.log;        \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med    --check correct --output ${PREFIX}.wfa-med.out        --repeat 100000 2>> ${PREFIX}.wfa-med.log; done'
    sbatch -c 128 --exclusive --job-name wfa-low-$SET    --wrap 'echo wfa-low-'$SET';   cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-low.log;        \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low    --check correct --output ${PREFIX}.wfa-low.out        --repeat 100000 2>> ${PREFIX}.wfa-low.log; done'
    sbatch -c 128 --exclusive --job-name wfa-high-$SET   --wrap 'echo wfa-high-'$SET';  cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-high.log;       \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high   --check correct --output ${PREFIX}.wfa-high.out       --repeat 100000 2>> ${PREFIX}.wfa-high.log; done'

    sbatch -c 128 --exclusive --job-name wfalm-$SET      --wrap 'echo wfalm-'$SET';     cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                   --check correct --output ${PREFIX}.wfalm.out          --repeat 100000 2>> ${PREFIX}.wfalm.log; done'
    sbatch -c 128 --exclusive --job-name wfalm-low-$SET  --wrap 'echo wfalm-low-'$SET'; cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-lowmem.log;   \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                            --check correct --output ${PREFIX}.wfalm-lowmem.out   --repeat 100000 2>> ${PREFIX}.wfalm-lowmem.log; done'
    sbatch -c 128 --exclusive --job-name wfalm-rec-$SET  --wrap 'echo wfalm-rec-'$SET'; cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-rec.log;      \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-rec                               --check correct --output ${PREFIX}.wfalm-rec.out      --repeat 100000 2>> ${PREFIX}.wfalm-rec.log; done'

    sbatch -c 128 --exclusive --job-name edlib-$SET      --wrap 'echo edlib-'$SET';     cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.edlib.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                   --check correct --output ${PREFIX}.edlib.out          --repeat 100000 2>> ${PREFIX}.edlib.log; done'
    sbatch -c 128 --exclusive --job-name bitpal-$SET     --wrap 'echo bitpal-'$SET';    cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.bitpal-scored.log;  \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                           --check correct --output ${PREFIX}.bitpal-scored.out  --repeat 100000 2>> ${PREFIX}.bitpal-scored.log; done'
    sbatch -c 128 --exclusive --job-name ksw2-$SET       --wrap 'echo ksw2-'$SET';      cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.ksw2-extz2-sse.log; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                          --check correct --output ${PREFIX}.ksw2-extz2-sse.out --repeat 100000 2>> ${PREFIX}.ksw2-extz2-sse.log; done'
  done
done

# UNUSED CODE
## Short reads (<= 10 kbps) for ont_regions
## Small runtimes require multiple runs to take their averages (for the time) and maximum (for the memory)
#SET=ont_regions_10kbps
#DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_alignments/$SET
#mkdir -p $DIR_OUTPUT
## Write just one .out file, not one for each replicate
#RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark
#sbatch -c 128 --exclusive --job-name biwfa-$SET      --wrap 'echo biwfa-'$SET';     ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.biwfa.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode biwfa    --check correct --output ${PREFIX}.biwfa.out          2>> ${PREFIX}.$i.biwfa.log; done; done'
#sbatch -c 128 --exclusive --job-name wfa-med-$SET    --wrap 'echo wfa-med-'$SET';   ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.wfa-med.log;        \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med      --check correct --output ${PREFIX}.wfa-med.out        2>> ${PREFIX}.$i.wfa-med.log; done; done'
#sbatch -c 128 --exclusive --job-name wfa-low-$SET    --wrap 'echo wfa-low-'$SET';   ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.wfa-low.log;        \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low      --check correct --output ${PREFIX}.wfa-low.out        2>> ${PREFIX}.$i.wfa-low.log; done; done'
#sbatch -c 128 --exclusive --job-name wfa-high-$SET   --wrap 'echo wfa-high-'$SET';  ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.wfa-high.log;       \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high     --check correct --output ${PREFIX}.wfa-high.out       2>> ${PREFIX}.$i.wfa-high.log; done; done'
#RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark
#sbatch -c 128 --exclusive --job-name wfalm-$SET      --wrap 'echo wfalm-'$SET';     ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.wfalm.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                     --check correct --output ${PREFIX}.wfalm.out          2>> ${PREFIX}.$i.wfalm.log; done; done'
#sbatch -c 128 --exclusive --job-name wfalm-low-$SET  --wrap 'echo wfalm-low-'$SET'; ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.wfalm-lowmem.log;   \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                              --check correct --output ${PREFIX}.wfalm-lowmem.out   2>> ${PREFIX}.$i.wfalm-lowmem.log; done; done'
#sbatch -c 128 --exclusive --job-name wfalm-rec-$SET  --wrap 'echo wfalm-rec-'$SET'; ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.wfalm-rec.log;      \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-rec                                 --check correct --output ${PREFIX}.wfalm-rec.out      2>> ${PREFIX}.$i.wfalm-rec.log; done; done'
#sbatch -c 128 --exclusive --job-name edlib-$SET      --wrap 'echo edlib-'$SET';     ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.edlib.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                     --check correct --output ${PREFIX}.edlib.out          2>> ${PREFIX}.$i.edlib.log; done; done'
#sbatch -c 128 --exclusive --job-name bitpal-$SET     --wrap 'echo bitpal-'$SET';    ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.bitpal-scored.log;  \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                             --check correct --output ${PREFIX}.bitpal-scored.out  2>> ${PREFIX}.$i.bitpal-scored.log; done; done'
#sbatch -c 128 --exclusive --job-name ksw2-$SET       --wrap 'echo ksw2-'$SET';      ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); seq 1 100 | while read i; do PREFIX='$DIR_OUTPUT'/$NAME/$NAME; mkdir -p '$DIR_OUTPUT'/$NAME/; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo "Replicate $i" > ${PREFIX}.$i.ksw2-extz2-sse.log; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                            --check correct --output ${PREFIX}.ksw2-extz2-sse.out 2>> ${PREFIX}.$i.ksw2-extz2-sse.log; done; done'
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens
mkdir -p seq_statistics/

rm seq_statistics/statistics.tsv
rm seq_statistics/scores.tsv
rm seq_statistics/lengths.tsv

ls seq_alignments/ | while read SET; do
  echo $SET
 
  cat seq_alignments/$SET/*.log | python3 ../log2info.py | awk -v OFS='\t' -v SET=$SET '{print(SET,$0)}' >> seq_statistics/statistics.tsv
  grep '^-' seq_alignments/$SET/*.out | cut -f 1 | cut -f 3 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0)}' >> seq_statistics/scores.tsv
  
  ls seq_pairs/$SET/*.seq | while read PATH_SEQ; do
    NAME=$(basename $PATH_SEQ .seq)
    PREFIX=seq_alignments/$SET/$NAME
    
    cat $PATH_SEQ | tr '\n' ' ' | less -S | awk -v OFS='\t' -v SET=$SET -v NAME=$NAME '{print(SET, NAME, length($1)-1, length($2)-1)}' >> seq_statistics/lengths.tsv

#    grep 'Time.Alignment' $PREFIX.*.log | grep call -v | cut -f 3 -d '/' | cut -f 1,2,4,5 -d '.' | tr -s ' ' | sed 's/ (s)//' | sed 's/.Alignment//' | sed 's/\./ /' | tr ' ' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"time_s")}' >> seq_statistics/statistics.tsv
#    grep 'Maximum resident' $PREFIX.*.log  | cut -f 3 -d '/' | sed 's/.log://' | sed 's/Maximum resident set size (kbytes): //g' | tr '.' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"memory_kb")}' >> seq_statistics/statistics.tsv
#    grep '^-' $PREFIX.*.out | cut -f 1 | cut -f 3 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"memory_kb")}' >> seq_statistics/scores.tsv
  done
done
```


## ONT MinION

Prepare the sequence pairs (from [Bowden et. al, 2019](https://doi.org/10.1038/s41467-019-09637-5)):

```shell
FILTER_SEQ_PY=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/evaluation/scripts/filter_seq.py

cd /gpfs/projects/bsc18/bsc18995/biwfa
mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs

# Long reads (> 500 kbps)
c=1
cat /gpfs/projects/bsc18/bsc18571/wfa2/datasets/real/Nanopore.UL.100K.seq | python3 $FILTER_SEQ_PY 500001 9999999999 | \
  while read -r line1; do
  read -r line2

  echo $line1 > /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/seq$c.seq
  echo $line2 >> /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/seq$c.seq
        
  c=$((c + 1))
done

# Other reads (<= 500 kbps)
c=$(ls ont_ul/seq_pairs/ | wc -l)
c=$((c + 1))
cat /gpfs/projects/bsc18/bsc18571/wfa2/datasets/real/Nanopore.UL.100K.seq | python3 $FILTER_SEQ_PY 0 500000 | head -n 3904 | \
  while read -r line1; do
  read -r line2

  echo $line1 > /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/seq$c.seq
  echo $line2 >> /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/seq$c.seq
        
  c=$((c + 1))
done
```

You can find the sequence pairs in the `~/BiWFA-paper/evaluation/data/ONT_MinION_UL.zip` file. 

Prepare batches (for better parallelization):
```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/

# Divide in chunks to avoid too long jobs
ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | split -l 200 - chunk_
```

Align sequences:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_alignments
mkdir -p $DIR_OUTPUT
SET=ont_ul

module load cmake gcc/10.2.0
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark

ls chunk_* | while read PATH_CHUNK; do
  echo $PATH_CHUNK

  # Write just one .out file, not one for each replicate
  # NOTE: keep the single ', to avoid bash variables interpolation before sending the jobs

  sbatch -c 128 --exclusive --job-name biwfa-$SET      --wrap 'echo biwfa-'$SET';     cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.biwfa.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode biwfa  --check correct --output ${PREFIX}.biwfa.out          --repeat 100000 2>> ${PREFIX}.biwfa.log; done'
  sbatch -c 128 --exclusive --job-name wfa-med-$SET    --wrap 'echo wfa-med-'$SET';   cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-med.log;        \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med    --check correct --output ${PREFIX}.wfa-med.out        --repeat 100000 2>> ${PREFIX}.wfa-med.log; done'
  sbatch -c 128 --exclusive --job-name wfa-low-$SET    --wrap 'echo wfa-low-'$SET';   cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-low.log;        \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low    --check correct --output ${PREFIX}.wfa-low.out        --repeat 100000 2>> ${PREFIX}.wfa-low.log; done'
  sbatch -c 128 --exclusive --job-name wfa-high-$SET   --wrap 'echo wfa-high-'$SET';  cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-high.log;       \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high   --check correct --output ${PREFIX}.wfa-high.out       --repeat 100000 2>> ${PREFIX}.wfa-high.log; done'

  sbatch -c 128 --exclusive --job-name wfalm-$SET      --wrap 'echo wfalm-'$SET';     cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                   --check correct --output ${PREFIX}.wfalm.out          --repeat 100000 2>> ${PREFIX}.wfalm.log; done'
  sbatch -c 128 --exclusive --job-name wfalm-low-$SET  --wrap 'echo wfalm-low-'$SET'; cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-lowmem.log;   \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                            --check correct --output ${PREFIX}.wfalm-lowmem.out   --repeat 100000 2>> ${PREFIX}.wfalm-lowmem.log; done'
  sbatch -c 128 --exclusive --job-name wfalm-rec-$SET  --wrap 'echo wfalm-rec-'$SET'; cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-rec.log;      \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-rec                               --check correct --output ${PREFIX}.wfalm-rec.out      --repeat 100000 2>> ${PREFIX}.wfalm-rec.log; done'

  sbatch -c 128 --exclusive --job-name edlib-$SET      --wrap 'echo edlib-'$SET';     cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.edlib.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                   --check correct --output ${PREFIX}.edlib.out          --repeat 100000 2>> ${PREFIX}.edlib.log; done'
  sbatch -c 128 --exclusive --job-name bitpal-$SET     --wrap 'echo bitpal-'$SET';    cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.bitpal-scored.log;  \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                           --check correct --output ${PREFIX}.bitpal-scored.out  --repeat 100000 2>> ${PREFIX}.bitpal-scored.log; done'
  sbatch -c 128 --exclusive --job-name ksw2-$SET       --wrap 'echo ksw2-'$SET';      cat '$PATH_CHUNK' | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.ksw2-extz2-sse.log; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                          --check correct --output ${PREFIX}.ksw2-extz2-sse.out --repeat 100000 2>> ${PREFIX}.ksw2-extz2-sse.log; done'
done

# UNUSED CODE
#RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark_wfa
#sbatch -c 128 --exclusive --job-name biwfa-ont_ul     --wrap 'echo biwfa-ont_ul;     ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.biwfa.log;          \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode biwfa --check correct --output '${PREFIX}'.biwfa.out          2>> '${PREFIX}'.biwfa.log; done'
#sbatch -c 128 --exclusive --job-name wfa-med-ont_ul   --wrap 'echo wfa-med-ont_ul;   ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfa-med.log;        \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med   --check correct --output '${PREFIX}'.wfa-med.out        2>> '${PREFIX}'.wfa-med.log; done'
#sbatch -c 128 --exclusive --job-name wfa-low-ont_ul   --wrap 'echo wfa-low-ont_ul;   ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfa-low.log;        \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low   --check correct --output '${PREFIX}'.wfa-low.out        2>> '${PREFIX}'.wfa-low.log; done'
#sbatch -c 128 --exclusive --job-name wfa-high-ont_ul  --wrap 'echo wfa-high-ont_ul;  ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfa-high.log;       \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high  --check correct --output '${PREFIX}'.wfa-high.out       2>> '${PREFIX}'.wfa-high.log; done'
#
#RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark
#sbatch -c 128 --exclusive --job-name wfalm-ont_ul     --wrap 'echo wfalm-ont_ul;     ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfalm.log;          \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm                                  --check correct --output '${PREFIX}'.wfalm.out          2>> '${PREFIX}'.wfalm.log; done'
#sbatch -c 128 --exclusive --job-name wfalm-low-ont_ul --wrap 'echo wfalm-low-ont_ul; ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-lowmem.log;   \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm-lowmem                           --check correct --output '${PREFIX}'.wfalm-lowmem.out   2>> '${PREFIX}'.wfalm-lowmem.log; done'
#sbatch -c 128 --exclusive --job-name wfalm-rec-ont_ul --wrap 'echo wfalm-rec-ont_ul; ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-rec.log;      \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm-rec                              --check correct --output '${PREFIX}'.wfalm-rec.out      2>> '${PREFIX}'.wfalm-rec.log; done'
#
#sbatch -c 128 --exclusive --job-name edlib-ont_ul     --wrap 'echo edlib-ont_ul;     ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.edlib.log;          \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a edlib                                  --check correct --output '${PREFIX}'.edlib.out          2>> '${PREFIX}'.edlib.log; done'
#sbatch -c 128 --exclusive --job-name bitpal-ont_ul    --wrap 'echo bitpal-ont_ul;    ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.bitpal-scored.log;  \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a bitpal-scored                          --check correct --output '${PREFIX}'.bitpal-scored.out  2>> '${PREFIX}'.bitpal-scored.log; done'
#sbatch -c 128 --exclusive --job-name ksw2-ont_ul      --wrap 'echo ksw2-ont_ul;      ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.ksw2-extz2-sse.log; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                         --check correct --output '${PREFIX}'.ksw2-extz2-sse.out 2>> '${PREFIX}'.ksw2-extz2-sse.log; done'

# UNUSED CODE
#ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do
#  NAME=$(basename $PATH_SEQ .seq);
#  PREFIX=$DIR_OUTPUT/$NAME;
#  echo $NAME;
#
#  sbatch -c 128 --exclusive --job-name biwfa-ont_ul     --wrap 'echo biwfa-ont_ul;     echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.biwfa.log;          \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode biwfa  --check correct --output '${PREFIX}'.biwfa.out          2>> '${PREFIX}'.biwfa.log'
#  sbatch -c 128 --exclusive --job-name wfa-med-ont_ul   --wrap 'echo wfa-med-ont_ul;   echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfa-med.log;        \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med    --check correct --output '${PREFIX}'.wfa-med.out        2>> '${PREFIX}'.wfa-med.log'
#  sbatch -c 128 --exclusive --job-name wfa-low-ont_ul   --wrap 'echo wfa-low-ont_ul;   echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfa-low.log;        \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low    --check correct --output '${PREFIX}'.wfa-low.out        2>> '${PREFIX}'.wfa-low.log'
#  sbatch -c 128 --exclusive --job-name wfa-high-ont_ul  --wrap 'echo wfa-high-ont_ul;  echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfa-high.log;       \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high   --check correct --output '${PREFIX}'.wfa-high.out       2>> '${PREFIX}'.wfa-high.log'
#
#  sbatch -c 128 --exclusive --job-name wfalm-ont_ul     --wrap 'echo wfalm-ont_ul;     echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfalm.log;          \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm                                   --check correct --output '${PREFIX}'.wfalm.out          2>> '${PREFIX}'.wfalm.log'
#  sbatch -c 128 --exclusive --job-name wfalm-low-ont_ul --wrap 'echo wfalm-low-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-lowmem.log;   \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm-lowmem                            --check correct --output '${PREFIX}'.wfalm-lowmem.out   2>> '${PREFIX}'.wfalm-lowmem.log'
#  sbatch -c 128 --exclusive --job-name wfalm-rec-ont_ul --wrap 'echo wfalm-rec-ont_ul; echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-rec.log;      \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a wfalm-rec                               --check correct --output '${PREFIX}'.wfalm-rec.out      2>> '${PREFIX}'.wfalm-rec.log'
#  
#  sbatch -c 128 --exclusive --job-name edlib-ont_ul     --wrap 'echo edlib-ont_ul;     echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.edlib.log;          \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a edlib                                   --check correct --output '${PREFIX}'.edlib.out          2>> '${PREFIX}'.edlib.log'
#  sbatch -c 128 --exclusive --job-name bitpal-ont_ul    --wrap 'echo bitpal-ont_ul;    echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.bitpal-scored.log;  \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a bitpal-scored                           --check correct --output '${PREFIX}'.bitpal-scored.out  2>> '${PREFIX}'.bitpal-scored.log'
#  sbatch -c 128 --exclusive --job-name ksw2-ont_ul      --wrap 'echo wfalm-ont_ul;     echo '$NAME'; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo -n "InfoAlignment: " > ${PREFIX}.ksw2-extz2-sse.log; \time -v '${RUN_BENCHMARK}' -i '${PATH_SEQ}' --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                          --check correct --output '${PREFIX}'.ksw2-extz2-sse.out 2>> '${PREFIX}'.ksw2-extz2-sse.log'
#done
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/
mkdir -p seq_statistics/

rm seq_statistics/lengths.tsv

cat seq_alignments/*.log | python3 ../log2info.py | awk -v OFS='\t' -v SET='ONT' '{print(SET,$0)}' > seq_statistics/statistics.tsv
grep '^-' seq_alignments/*.out | cut -f 1 | cut -f 2 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET="ONT" '{print(SET,$0)}' > seq_statistics/scores.tsv

ls seq_pairs/*.seq | while read PATH_SEQ; do
  NAME=$(basename $PATH_SEQ .seq)
  PREFIX=seq_alignments/$SET/$NAME
    
  cat $PATH_SEQ | tr '\n' ' ' | less -S | awk -v OFS='\t' -v SET="ONT" -v NAME=$NAME '{print(SET, NAME, length($1)-1, length($2)-1)}' >> seq_statistics/lengths.tsv
done
```


## PacBio HiFi

Prepare sequence pairs:

```shell
mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs

# First 2000 pairs
c=1
head -n 4000 /gpfs/projects/bsc18/bsc18571/wfa2/datasets/real/PacBio.HF.1M.seq | while read -r line1; do
    read -r line2
    
    echo $line1 > /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/seq$c.seq
    echo $line2 >> /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/seq$c.seq
    
    c=$((c + 1))
done
```

Align sequences:

```shell
module load cmake gcc/10.2.0

cd /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_alignments
mkdir -p $DIR_OUTPUT

RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark_wfa
sbatch -c 48 --job-name biwfa-pacbio      --wrap 'echo biwfa-pacbio;     ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.biwfa.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode biwfa --check correct --output ${PREFIX}.biwfa.out          2>> ${PREFIX}.biwfa.log; done'
sbatch -c 48 --job-name wfa-med-pacbio    --wrap 'echo wfa-med-pacbio;   ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-med.log;        \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med   --check correct --output ${PREFIX}.wfa-med.out        2>> ${PREFIX}.wfa-med.log; done'
sbatch -c 48 --job-name wfa-low-pacbio    --wrap 'echo wfa-low-pacbio;   ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-low.log;        \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low   --check correct --output ${PREFIX}.wfa-low.out        2>> ${PREFIX}.wfa-low.log; done'
sbatch -c 48 --job-name wfa-high-pacbio   --wrap 'echo wfa-high-pacbio;  ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfa-high.log;       \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high  --check correct --output ${PREFIX}.wfa-high.out       2>> ${PREFIX}.wfa-high.log; done'

RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark
sbatch -c 48 --job-name wfalm-pacbio      --wrap 'echo wfalm-pacbio;     ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                  --check correct --output ${PREFIX}.wfalm.out          2>> ${PREFIX}.wfalm.log; done'
sbatch -c 48 --job-name wfalm-low-pacbio  --wrap 'echo wfalm-low-pacbio; ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-lowmem.log;   \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                           --check correct --output ${PREFIX}.wfalm-lowmem.out   2>> ${PREFIX}.wfalm-lowmem.log; done'
sbatch -c 48 --job-name wfalm-rec-pacbio  --wrap 'echo wfalm-rec-pacbio; ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.wfalm-rec.log;      \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-rec                              --check correct --output ${PREFIX}.wfalm-rec.out      2>> ${PREFIX}.wfalm-rec.log; done'

sbatch -c 48 --job-name edlib-pacbio      --wrap 'echo edlib-pacbio;     ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.edlib.log;          \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                  --check correct --output ${PREFIX}.edlib.out          2>> ${PREFIX}.edlib.log; done'
sbatch -c 48 --job-name bitpal-pacbio     --wrap 'echo bitpal-pacbio;    ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.bitpal-scored.log;  \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                          --check correct --output ${PREFIX}.bitpal-scored.out  2>> ${PREFIX}.bitpal-scored.log; done'
sbatch -c 48 --job-name ksw2-pacbio       --wrap 'echo ksw2-pacbio;      ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; echo $NAME; echo -n "InfoAlignment: " > ${PREFIX}.ksw2-extz2-sse.log; \time -v '${RUN_BENCHMARK}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                         --check correct --output ${PREFIX}.ksw2-extz2-sse.out 2>> ${PREFIX}.ksw2-extz2-sse.log; done'
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/
mkdir -p seq_statistics/
   
cat seq_alignments/*.log | python3 ../log2info.py | awk -v OFS='\t' -v SET="HIFI" '{print(SET,$0)}' > seq_statistics/statistics.tsv
grep '^-' seq_alignments/*.out | cut -f 1 | cut -f 2 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET="HIFI" '{print(SET,$0)}' > seq_statistics/scores.tsv

ls seq_pairs/*.seq | while read PATH_SEQ; do
  NAME=$(basename $PATH_SEQ .seq)
  PREFIX=seq_alignments/$SET/$NAME
    
  cat $PATH_SEQ | tr '\n' ' ' | less -S | awk -v OFS='\t' -v SET="HIFI" -v NAME=$NAME '{print(SET, NAME, length($1)-1, length($2)-1)}' >> seq_statistics/lengths.tsv
done
```


## Statistics

Merge all statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa

cat hsapiens/seq_statistics/statistics.tsv > statistics_all.tsv
cat hsapiens/seq_statistics/statistics.10kbps.tsv >> statistics_all.tsv
cat ont_ul/seq_statistics/statistics.tsv >> statistics_all.tsv
cat ont_ul_10kbps/seq_statistics/statistics.10kbps.tsv >> statistics_all.tsv
cat ont_ul_50kbps/seq_statistics/statistics.tsv >> statistics_all.tsv

cat hsapiens/seq_statistics/scores.tsv > scores_all.tsv
cat hsapiens/seq_statistics/scores.10kbps.tsv >> scores_all.tsv
cat ont_ul/seq_statistics/scores.tsv >> scores_all.tsv
cat ont_ul_10kbps/seq_statistics/scores.10kbps.tsv >> scores_all.tsv
cat ont_ul_50kbps/seq_statistics/scores.tsv >> scores_all.tsv

cat hsapiens/seq_statistics/lengths.tsv > lengths_all.tsv
cat hsapiens/seq_statistics/lengths.10kbps.tsv >> lengths_all.tsv
cat ont_ul/seq_statistics/lengths.tsv >> lengths_all.tsv
cat ont_ul_10kbps/seq_statistics/lengths.10kbps.tsv >> lengths_all.tsv
cat ont_ul_50kbps/seq_statistics/lengths.tsv >> lengths_all.tsv

pigz statistics_all.tsv -11 -f
pigz scores_all.tsv -f
pigz lengths_all.tsv -f
```

Plot the statistics with the `scripts/plots.R` script.