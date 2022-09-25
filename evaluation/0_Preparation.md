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
```
