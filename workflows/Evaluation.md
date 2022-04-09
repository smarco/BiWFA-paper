# Evaluation

Prepare the tools:

```shell
git clone --recursive https://github.com/smarco/BiWFA-paper.git
git clone --recursive https://github.com/smarco/WFA2-lib.git

# Move the repositories on the cluster
scp -r BiWFA-paper bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/
scp -r WFA2-lib bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/

# Build
module load cmake gcc/10.2.0

cd /gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper
make clean all
cd ..

cd WFA2-lib
git checkout benchmark
# Remove block-aligner
make clean all
```

## Dataset 1

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
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark
RUN_BENCHMARK_EXT=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark

ls seq_pairs/ | while read SET; do
  echo $SET

  DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_alignments/$SET
  mkdir -p $DIR_OUTPUT

  sbatch -c 48 --job-name biwfa-$SET     --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-bidirectional    --check correct --output ${PREFIX}.biwfa.out          2> ${PREFIX}.biwfa.log; done'
  sbatch -c 48 --job-name wfa-med-$SET   --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med  --check correct --output ${PREFIX}.wfa-mem.out        2> ${PREFIX}.wfa-mem.log; done'
  sbatch -c 48 --job-name wfa-low-$SET   --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low  --check correct --output ${PREFIX}.wfa-low.out        2> ${PREFIX}.wfa-low.log; done'
  sbatch -c 48 --job-name wfa-high-$SET  --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high --check correct --output ${PREFIX}.wfa-high.out       2> ${PREFIX}.wfa-high.log; done'
  
  sbatch -c 48 --job-name wfalm-$SET     --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                                 --output ${PREFIX}.wfalm.out          2> ${PREFIX}.wfalm.log; done'
  sbatch -c 48 --job-name wfalm-low-$SET --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                                          --output ${PREFIX}.wfalm-lowmem.out   2> ${PREFIX}.wfalm-lowmem.log; done'

  sbatch -c 48 --job-name edlib-$SET      --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                                 --output ${PREFIX}.edlib.out          2> ${PREFIX}.edlib.log; done'
  sbatch -c 48 --job-name bitpal-$SET     --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                                         --output ${PREFIX}.bitpal-scored.out  2> ${PREFIX}.bitpal-scored.log; done'
  sbatch -c 48 --job-name ksw2-$SET       --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_pairs/'$SET'/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                                        --output ${PREFIX}.ksw2-extz2-sse.out 2> ${PREFIX}.ksw2-extz2-sse.log; done'

#  ls seq_pairs/$SET/*.seq | while read PATH_SEQ; do
#  cat ONTaa | while read PATH_SEQ; do
#    NAME=$(basename $PATH_SEQ .seq)
#    PREFIX=$DIR_OUTPUT/$NAME
#
#    FULL_PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/$PATH_SEQ
#    
#    #sbatch -c 1  --job-name $SET --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${FULL_PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-bidirectional    --check correct --output '${PREFIX}'.bid.out  2> '${PREFIX}'.bid.log'
#    #sbatch -c 1  --job-name $SET --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${FULL_PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med  --check correct --output '${PREFIX}'.med.out  2> '${PREFIX}'.med.log'
#    #sbatch -c 1  --job-name $SET --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${FULL_PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low  --check correct --output '${PREFIX}'.low.out  2> '${PREFIX}'.low.log'
#    # -c 48 because the memory consumption can be high
#    #sbatch -c 48 --job-name $SET --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -i '${FULL_PATH_SEQ}' --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high --check correct --output '$PREFIX'.high.out 2> '$PREFIX'.high.log'
#  done
done
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/hsapiens
mkdir -p seq_statistics/

rm seq_statistics/statistics.tsv
rm seq_statistics/scores.tsv

ls seq_pairs/ | while read SET; do
  echo $SET
      
  ls seq_pairs/$SET/*.seq | while read PATH_SEQ; do
    NAME=$(basename $PATH_SEQ .seq)
    PREFIX=seq_alignments/$SET/$NAME

    grep 'Time.Alignment' $PREFIX.*.log | grep call -v | cut -f 3 -d '/' | cut -f 1,2,4,5 -d '.' | tr -s ' ' | sed 's/ (s)//' | sed 's/.Alignment//' | sed 's/\./ /' | tr ' ' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"time_s")}' >> seq_statistics/statistics.tsv
    grep 'Maximum resident' $PREFIX.*.log  | cut -f 3 -d '/' | sed 's/.log://' | sed 's/Maximum resident set size (kbytes): //g' | tr '.' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"memory_kb")}' >> seq_statistics/statistics.tsv

    grep '^-' $PREFIX.*.out | cut -f 1 | cut -f 3 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' -v SET=$SET '{print(SET,$0,"memory_kb")}' >> seq_statistics/scores.tsv
  done
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
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark
RUN_BENCHMARK_EXT=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark

cd /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_alignments
mkdir -p $DIR_OUTPUT

sbatch -c 48 --job-name biwfa-pacbio      --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-bidirectional    --check correct --output ${PREFIX}.biwfa.out          2> ${PREFIX}.biwfa.log; done'
sbatch -c 48 --job-name wfa-med-pacbio    --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med  --check correct --output ${PREFIX}.wfa-mem.out        2> ${PREFIX}.wfa-mem.log; done'
sbatch -c 48 --job-name wfa-low-pacbio    --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low  --check correct --output ${PREFIX}.wfa-low.out        2> ${PREFIX}.wfa-low.log; done'
sbatch -c 48 --job-name wfa-high-pacbio   --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high --check correct --output ${PREFIX}.wfa-high.out       2> ${PREFIX}.wfa-high.log; done'

sbatch -c 48 --job-name wfalm-pacbio      --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                                 --output ${PREFIX}.wfalm.out          2> ${PREFIX}.wfalm.log; done'
sbatch -c 48 --job-name wfalm-low-pacbio  --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                                          --output ${PREFIX}.wfalm-lowmem.out   2> ${PREFIX}.wfalm-lowmem.log; done'

sbatch -c 48 --job-name edlib-pacbio      --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                                 --output ${PREFIX}.edlib.out          2> ${PREFIX}.edlib.log; done'
sbatch -c 48 --job-name bitpal-pacbio     --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                                         --output ${PREFIX}.bitpal-scored.out  2> ${PREFIX}.bitpal-scored.log; done'
sbatch -c 48 --job-name ksw2-pacbio       --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                                        --output ${PREFIX}.ksw2-extz2-sse.out 2> ${PREFIX}.ksw2-extz2-sse.log; done'
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/pacbio/
mkdir -p seq_statistics/

rm seq_statistics/statistics.tsv
rm seq_statistics/scores.tsv
   
ls seq_pairs/*.seq | while read PATH_SEQ; do
  NAME=$(basename $PATH_SEQ .seq)
  PREFIX=seq_alignments/$NAME

  grep 'Time.Alignment' $PREFIX.*.log | grep call -v | cut -f 2 -d '/' | cut -f 1,2,4,5 -d '.' | tr -s ' ' | sed 's/ (s)//' | sed 's/.Alignment//' | sed 's/\./ /' | tr ' ' '\t' | awk -v OFS='\t' '{print("pacbio",$0,"time_s")}' >> seq_statistics/statistics.tsv
  
  grep 'Maximum resident' $PREFIX.*.log  | cut -f 2 -d '/' | sed 's/.log://' | sed 's/Maximum resident set size (kbytes): //g' | tr '.' '\t' | awk -v OFS='\t' '{print("pacbio",$0,"memory_kb")}' >> seq_statistics/statistics.tsv

  grep '^-' $PREFIX.*.out | cut -f 1 | cut -f 2 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' '{print("pacbio",$0,"memory_kb")}' >> seq_statistics/scores.tsv
done
```

## ONT UL

Prepare sequence pairs:

```shell
mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs

# First 1000 pairs
c=1
head -n 2000 /gpfs/projects/bsc18/bsc18571/wfa2/datasets/real/Nanopore.UL.100K.seq | while read -r line1; do
    read -r line2
    
    echo $line1 > /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/seq$c.seq
    echo $line2 >> /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/seq$c.seq
    
    c=$((c + 1))
done
```

Align sequences:

```shell
RUN_BENCHMARK=/gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper/bin/align_benchmark
RUN_BENCHMARK_EXT=/gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib/bin/align_benchmark

cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/

DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_alignments
mkdir -p $DIR_OUTPUT

sbatch -c 48 --job-name biwfa-ont_ul      --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-bidirectional    --check correct --output ${PREFIX}.biwfa.out          2> ${PREFIX}.biwfa.log; done'
sbatch -c 48 --job-name wfa-med-ont_ul    --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode med  --check correct --output ${PREFIX}.wfa-mem.out        2> ${PREFIX}.wfa-mem.log; done'
sbatch -c 48 --job-name wfa-low-ont_ul    --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode low  --check correct --output ${PREFIX}.wfa-low.out        2> ${PREFIX}.wfa-low.log; done'
sbatch -c 48 --job-name wfa-high-ont_ul   --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}'     -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a gap-affine-wfa --wfa-memory-mode high --check correct --output ${PREFIX}.wfa-high.out       2> ${PREFIX}.wfa-high.log; done'

sbatch -c 48 --job-name wfalm-ont_ul      --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm                                                 --output ${PREFIX}.wfalm.out          2> ${PREFIX}.wfalm.log; done'
sbatch -c 48 --job-name wfalm-low-ont_ul  --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a wfalm-lowmem                                          --output ${PREFIX}.wfalm-lowmem.out   2> ${PREFIX}.wfalm-lowmem.log; done'

sbatch -c 48 --job-name edlib-ont_ul      --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a edlib                                                 --output ${PREFIX}.edlib.out          2> ${PREFIX}.edlib.log; done'
sbatch -c 48 --job-name bitpal-ont_ul     --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a bitpal-scored                                         --output ${PREFIX}.bitpal-scored.out  2> ${PREFIX}.bitpal-scored.log; done'
sbatch -c 48 --job-name ksw2-ont_ul       --wrap 'ls /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/seq_pairs/*.seq | while read PATH_SEQ; do NAME=$(basename $PATH_SEQ .seq); PREFIX='$DIR_OUTPUT'/$NAME; TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK_EXT}' -i ${PATH_SEQ} --affine-penalties 0,4,6,2 -a ksw2-extz2-sse                                        --output ${PREFIX}.ksw2-extz2-sse.out 2> ${PREFIX}.ksw2-extz2-sse.log; done'
```

Collect statistics:

```shell
cd /gpfs/projects/bsc18/bsc18995/biwfa/ont_ul/
mkdir -p seq_statistics/

rm seq_statistics/statistics.tsv
rm seq_statistics/scores.tsv
   
ls seq_pairs/*.seq | while read PATH_SEQ; do
  NAME=$(basename $PATH_SEQ .seq)
  PREFIX=seq_alignments/$NAME

  grep 'Time.Alignment' $PREFIX.*.log | grep call -v | cut -f 2 -d '/' | cut -f 1,2,4,5 -d '.' | tr -s ' ' | sed 's/ (s)//' | sed 's/.Alignment//' | sed 's/\./ /' | tr ' ' '\t' | awk -v OFS='\t' '{print("ont_ul",$0,"time_s")}' >> seq_statistics/statistics.tsv
  
  grep 'Maximum resident' $PREFIX.*.log  | cut -f 2 -d '/' | sed 's/.log://' | sed 's/Maximum resident set size (kbytes): //g' | tr '.' '\t' | awk -v OFS='\t' '{print("ont_ul",$0,"memory_kb")}' >> seq_statistics/statistics.tsv

  grep '^-' $PREFIX.*.out | cut -f 1 | cut -f 2 -d '/' | sed 's/out://' | tr '.' '\t' | awk -v OFS='\t' '{print("ont_ul",$0,"memory_kb")}' >> seq_statistics/scores.tsv
done
```
