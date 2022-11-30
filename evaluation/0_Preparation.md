# Preparation

Prepare the tools:

```shell
git clone --recursive https://github.com/smarco/BiWFA-paper
tar -zcvf BiWFA-paper.tar.gz BiWFA-paper

# Compress and move the repository on the cluster
scp -r BiWFA-paper.tar.gz bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/biwfa/

# Decompress and checkout the benchmark branch
mkdir -p /gpfs/projects/bsc18/bsc18995/biwfa/
cd /gpfs/projects/bsc18/bsc18995/biwfa/

tar -xf BiWFA-paper.tar.gz
rm BiWFA-paper.tar.gz
cd BiWFA-paper
git checkout benchmark

# Build
module load cmake gcc/10.2.0
make clean all
```
