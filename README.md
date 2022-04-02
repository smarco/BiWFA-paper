# Bidirectional WFA

### 1 What is this?

A WFA2-lib development branch implementing the bidirectional WFA alignment which allows `O(ns+s^2)` alignment time using `O(s)` memory.

Usage: As WFA2-lib but adding the option `--wfa-bidirectional`

```
./bin/align_benchmark -a gap-affine-wfa -i sequences.seq -o out.alg --affine-penalties 0,3,6,1 --wfa-bidirectional
```