# Homo sapiens

Prepare the tool:

```shell
git clone --recursive https://github.com/smarco/BiWFA-paper.git

# Move the repository on the cluster
scp -r BiWFA-paper bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/

# Build
module load cmake gcc/10.2.0

cd BiWFA-paper
make clean all
```

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

In the dataset:
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

ls seq_pairs/ | while read SET; do
  echo $SET
  
  DIR_OUTPUT=/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/seq_alignments/$SET
  mkdir -p $DIR_OUTPUT
    
  ls seq_pairs/$SET/*.seq | while read PATH_SEQ; do
    NAME=$(basename $PATH_SEQ .seq)
    PREFIX=$DIR_OUTPUT/$NAME
       
    FULL_PATH_SEQ=/gpfs/projects/bsc18/bsc18995/biwfa/hsapiens/$PATH_SEQ
    
    sbatch -c 1 --job-name $SET --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa -i '${FULL_PATH_SEQ}' --affine-penalties 0,4,6,2 --wfa-bidirectional    --check correct --output '${PREFIX}'.bid.out  2> '${PREFIX}'.bid.log'
    sbatch -c 1 --job-name $SET --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa -i '${FULL_PATH_SEQ}' --affine-penalties 0,4,6,2 --wfa-memory-mode med  --check correct --output '${PREFIX}'.med.out  2> '${PREFIX}'.med.log'
    sbatch -c 1 --job-name $SET --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa -i '${FULL_PATH_SEQ}' --affine-penalties 0,4,6,2 --wfa-memory-mode low  --check correct --output '${PREFIX}'.low.out  2> '${PREFIX}'.low.log'
  
    # -c 48 because the memory consumption can be high
    sbatch -c 48 --job-name $SET --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; \time -v '${RUN_BENCHMARK}' -a gap-affine-wfa -i '${FULL_PATH_SEQ}' --affine-penalties 0,4,6,2 --wfa-memory-mode high --check correct --output '$PREFIX'.high.out 2> '$PREFIX'.high.log'
  done
  break
done
```
