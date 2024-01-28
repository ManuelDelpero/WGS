# Use the Miniforge3 image as a base, which comes with conda
FROM condaforge/miniforge3:latest

# Set the working directory in the container to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Create the conda environment from the YAML file
RUN mamba env create -f /app/envs/WGS.yaml

# Volumes for input data, configuration, and output results
VOLUME ["/data", "/results"]

# Create a reference directory
RUN mkdir -p /app/reference

# Download and decompress the reference genome
RUN wget -q -O /app/reference/homo_sapiens.GRCh38.dna.primary_assembly.fa.gz https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    && gunzip /app/reference/homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Set an environment variable for the default number of cores
ENV SNAKEMAKE_CORES=1

# Ensure the container is executable
CMD ["/bin/bash"]

docker run \
  -v /Data/WGS_test/data/config.yaml:/app/config.yaml \
  -v /Data/WGS_test/data/:/data \
  -v /Data/WGS_test/data/results/:/results \
  -e SNAKEMAKE_CORES=20 \
  wgs \
  conda run -n WGS snakemake --cores ${SNAKEMAKE_CORES}


  ### REQUIRED ARGUMENTS BY USER
rawdir:
  '/data'
output_dir:
  '/results'

samples: [
  'G5',
  'G6'
]

ref:
  '/app/reference/homo_sapiens.GRCh38.dna.primary_assembly.fa'

computing_threads:
  15

Tertiary:
  1
