
# PaRCA
## Pathogen detection for Research and Clinical Applications

---
## Prerequisites:
* Databases for Kraken and Kaiju are currently manually downloaded 
* Pollux and Fiona are not available from conda and has to be manually downloaded

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
    --cluster-config config/cluster.yaml \
    --profile qsub_profile
```

### To-Do
