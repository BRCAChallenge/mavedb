## Overview

[MaveDB](https://www.mavedb.org/) is a biological database for Multiplex Assays of Variant Effect (MAVE) datasets, such as those generated by Deep Mutational Scanning (DMS) or Massively Parallel Reporter Assays (MPRA).  The software in this repo generates content for the **Variant Effect Maps** track hub, which renders MaveDB assay results mapped against the Hg38 reference genome.  This software relies on the mappings of the raw MaveDB scores from the original reference sequences to the canonical human reference sequences, building upon mappings generated by [Arbesfeld et al](https://www.biorxiv.org/content/10.1101/2023.06.20.545702v1).

## Prerequisites

This software relies on the [cool-seq-tool](https://pypi.org/project/cool-seq-tool/) python package to map protein to genomic coordinates.  The cool-seq-tool package in turn relies on the [Universal Transcript Archive (UTA)](https://github.com/biocommons/uta) database, the [SeqRepo](https://github.com/biocommons/biocommons.seqrepo) sequence mapping resource, and other resources.  See the installation instructions for [cool-seq-tool](https://pypi.org/project/cool-seq-tool/) for more information.

This software also requires the UCSC Genome Browser executables `bedSort` and `bedToBigBed`.  These should be installed in your path prior to execution.

## Installation 

1. Install the [cool-seq-tool](https://pypi.org/project/cool-seq-tool/) package, following the package installation instructions.
   
2. Install DynamoDBLocal_lib

## Execution

1. Launch Postgresql
   
2. In the *DynamoDBLocal_lib* directory, execute the command `java -Djava.library.path=./DynamoDBLocal_lib -jar DynamoDBLocal.jar -sharedDb`

3. In a separate shell, launch the cool-seq-tool API service as follows:
   
   1. Navigate to the cool_seq_tool source directory (which contains the file `api.py`)
  
   2. Execute the command `uvicorn cool_seq_tool.api:app --reload`
   
5. In a third shell:

    1. Activate the cool-seq-tool virtual environment
 
    2. Navigate to the cool_seq_tool source code directory.  
      
    4. Execute the command `pipenv shell`
      
    5. Set the following environment variables
       
        * `export UTA_VERSION=uta_20210129.pgd.gz`
   
        * `export UTA_PASSWORD=uta`
      
        * `export SEQREPO_ROOT_DIR=/usr/local/share/seqrepo/latest`
      
        * `export TRANSCRIPT_MAPPINGS_PATH=${PWD}/cool_seq_tool/data/transcript_mapping.tsv`
      
        * `export AWS_ACCESS_KEY_ID=fakeMyKeyId`
      
        * `export AWS_SECRET_ACCESS_KEY=fakeSecretAccessKey`
      
  6. Navigate to the `src` directory under `mavedb`. 
    
  8. Execute with a command as follows:

     ```src/mavedb_to_trackhub.py \
            -i input/00000098-a-1.json \
            -n "Variant Effect Maps" \
            -t trackDb.txt \
            -b output/bed \
            -l output/lm
            -d 1
     ```
     
     where
     
        `--input (-i)` specifies the name of one or more input json files (wildcards accepted)
     
        `--track_name (-n)` specifies the name of the track (default: "Variant Effect Maps")
     
        `--trackDb (-t)` specifies the pathname of the output trackDb file
     
        `--bed_dir (-b)` specifies a subdirectory for the output bed files
     
        `--location_matrix_dir (-l)` specifies a subdirectory for the output location matrix files
     
        `--debug (-d)` turns on debugging information


To Be Continued

src/bigHeat.py  output/bed/00000098-a-1.bed output/lm/00000098-a-1.lm ../hg38.chrom.sizes output/bb 
Processing output/lm/00000098-a-1.lm to output directory output/bb

