# Preparation

Prepare the tools:

```shell
git clone --recursive https://github.com/smarco/BiWFA-paper.git

# Move the repositories on the cluster
scp -r BiWFA-paper bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/

# Build
module load cmake gcc/10.2.0

cd /gpfs/projects/bsc18/bsc18995/biwfa/BiWFA-paper
git checkout benchmark
make clean all


git clone --recursive https://github.com/smarco/WFA2-lib.git

# Move the repositories on the cluster
scp -r WFA2-lib bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/

# Build
module load cmake gcc/10.2.0

cd /gpfs/projects/bsc18/bsc18995/biwfa/WFA2-lib
git checkout 4600d5d42ea0ee1cd7c346fc3861777f4492e961
sed -i 's/BUILD_EXAMPLES=1/BUILD_EXAMPLES=0/g' Makefile
sed -i 's/BUILD_WFA_PARALLEL=1/BUILD_WFA_PARALLEL=0/g' Makefile
make clean all
```
