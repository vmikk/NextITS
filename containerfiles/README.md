# NextITS container images

Reproducible computational environments are essential for scientific workflows.  
NextITS provides pre-built container images to ensure consistent software versions and dependencies across different computing platforms, eliminating "it works on my machine" issues and enabling reproducible bioinformatics analyses.  

In general, NextITS will pull the required images automatically (e.g., when providing the `-profile singularity` or `-profile docker` flag to the Nextflow command).  

However, if you want to pull or build the container images manually, you can do so using the following instructions.

## Pull pre-built images

### Docker Hub

```bash
# Pull specific version
docker pull vmikk/nextits:1.1.0
```

### Singularity library

```bash
# Pull specific version
singularity pull library://vmiks/nextits/nextits:1-1-0
```

