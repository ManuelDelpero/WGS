# Enhanced Scalable WGS Pipeline in Snakemake

This pipeline is designed for efficient processing and analysis of Whole Genome Sequencing (WGS) data. Utilizing renowned tools like `GATK4`, `hisat2`, `samtools`, `fastqc`, and `htseq-count`, this pipeline streamlines your WGS data workflow. Configure the pipeline easily through the `config.yaml` file, which holds all essential parameters.

## Getting Started: Prerequisites

Before diving in, ensure that Mamba is installed on your system. If not, visit [Mamba Installation Guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) for easy installation instructions.

## Setting Up the Pipeline

1. **Obtain the Pipeline**: Clone or download the pipeline code from our repository.
2. **Configuration**: Verify the presence of `config.yaml` in the main directory, alongside the Snakefile.

## Creating Your Conda Environment

1. **Environment Creation**: Formulate a new environment (e.g., named "WGS" or another name of your choosing) and install dependencies from the `environment.yaml` file:
   
   ```bash
   mamba env create -f envs/WGS.yaml
   ```

   This command establishes a new environment named "WGS" and installs all necessary dependencies.
   
2. **Customization**: Modify the config file with your specific parameters.

## Activating the Conda Environment and Running the Pipeline

Execute the following commands to activate your environment and start the pipeline:

```bash
conda activate WGS
snakemake --cores <num_cores>
```

## Optional: Tertiary Analysis Setup for WGS Pipeline

For those interested in tertiary analysis using VEP and ClinVar, set up the VEP cache directory and download the ClinVar VCF:

```bash
mkdir WGS/VEP_cache
cd WGS/VEP_cache
curl -O https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
tar xzf homo_sapiens_vep_110_GRCh38.tar.gz
cd ..
cd clinvar
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

## Docker Integration

Prefer Docker? Utilize our Dockerfile to build the image and run the pipeline within a Docker container. The built image includes all dependencies and the GRCh38 reference genome. The configfile can be downaloded from the repository.
If using Docker the only paramenters to modify in the config file are the samples, computing_threads, and tertiary.

Build and run the pipeline using Docker:

```bash
docker build -t wgs .

docker run \
  -v /path/to/config.yaml:/app/config.yaml \
  -v /path/to/raw_data/:/data \
  -v /path/to/results/:/results \
  wgs \
  conda run -n WGS snakemake --cores 20
```
