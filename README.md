
# PaRCA
## Pathogen Research in Clinical Applications

---

## **The pipeline**

### Stage 1: Quality control and error correction

### Stage 2: Assembly

### Stage 3: Kraken and Kaiju

### Stage 4: Parse hits

### Stage 5: Blast processing

### Stage 6: Blast sliced database

### Stage 7: Blast remaining reads

### Stage 8: Format results

## Usage

```
snakemake \
    -rp \
    -s main.smk \
    --use-singularity \
    --use-conda \
    --cluster-config cluster.yaml \
    --drmaa "qsub -S /bin/bash -pe mpi {cluster.cores} -q {cluster.queue} -N {rule} -l excl=1 -S /bin/bash" \
    --jobs 999 \
    --latency-wait 60
```

### To-Do
- Add qsub_profile
