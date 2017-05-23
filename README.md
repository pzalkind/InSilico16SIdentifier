# InSilico16SIdentifier
Open Source In Silico 16S Identification of prokaryotic organism. Takes FASTA as input and return the BLAST hits closest to your sequence. Based on BLAST + RDP database


## Requirements

* docker :)

## How to Install

docker pull pzalkind/ssu_identifier

## How to use
### Identification
You only need a FASTA file :) Mount the folder containing the FASTA file when running the container, so it gets access to the file.

Cmd:
docker run --rm -v <path_to_your_folder_with_fasta_files>:/data ssu_identifier <file>.fasta

Ideally move to the folder containing the fasta file and type (Unix systems):
docker run --rm -v $(pwd):/data ssu_identifier <file>.fasta
