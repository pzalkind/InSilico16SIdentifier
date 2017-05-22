# InSilico16SIdentifier
Open Source In Silico 16S Identification of prokaryotic organism. Takes FASTA as input and return the BLAST hits closest to your sequence. Based on BLAST + RDP database


# Identification script Documentation



## Requirements

* ncbi-blast+
* python3 (not necessary to replace your default python version)
* python3-pandas (sudo apt-get install)
* GenostarMicrob-x.y.z.tar.gz (> 1.0.7)
* 16S BLAST bank

## Software installation

* Copy the complete folder wherever you need it
  * `cd <path_to>/Identification_script`
  * `tar zxvf GenostarMicrob-x.y.z.tar.gz`
  * `export PATH="${PATH}:<path_to>/ncbi-blast-2.3.0+/bin:<path_to>/jre/bin"`

## Prepare input data

* Download and prepare RDP data

Identification database is a subset of the bacterial 16S RDP database, current release 11.4 (https://rdp.cme.msu.edu/download/releaseREADME.txt)
Original database can be found at https://rdp.cme.msu.edu/misc/resources.jsp, section Release Alignment, bacterial 16S, unaligned (fasta format)

`wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz`
`gunzip current_Bacteria_unaligned.fa.gz`

A few changes have been made to the original base, 
1. to keep only entries with names compliant to standard taxonomy (*Genus species*) :
 - Headers have been formatted to allow direct extraction of taxonomy from blast output [removed phylogenetic part of the headers, reformatted name and id to: >ID|Genus_species_(subsp._subspecies_str._strain)]
 - Entries with following pattern have been excluded from the database (case insensitive) :
  - 'unidentified'
  - 'uncultured'
  - 'uncultivated'
 - Entries with name not consistent with nomenclatural taxonomy were also removed, meaning :
  - names with first word not starting with an uppercase
  - names being less than 2 word long
  - names with 2nd word (species) being either 'sp.' or 'bacterium'

2. to keep only Type strains (representative genome) :
 - Headers have been formatted to allow direct extraction of taxonomy from blast output [removed phylogenetic part of the headers, reformatted name and id to: >ID|Genus_species_(subsp._subspecies_str._strain)]
 - Only entries with following pattern have been included in the refined database (case insensitive) :
  - '(T)' (marker of type strain genomes)
 - Entries with name not consistent with nomenclatural taxonomy were also removed, meaning :
  - names with first word not starting with an uppercase
  - names being less than 2 word long
  - names with 2nd word (species) being either 'sp.' or 'bacterium'
  
All necessary steps have been compiled in a script:

`python3 refineRDPdatabase.py current_Bacteria_unaligned.fa`

* Create the 16S BLAST bank:
  * `makeblastdb -in 16S_RDP_named.fasta -out 16S_RDP_named -title 16S_RDP_named -dbtype nucl`
  * `makeblastdb -in 16S_RDP_type_strains.fasta -out 16S_RDP_type_strains -title 16S_RDP_type_strains -dbtype nucl`
  * `export BLASTDB=$BLASTDB:$(pwd)`


* Download and prepare microB data

The microB rRNA16S bank contains all the 16S rRNA sequences from all organisms from the microB Reference database, the bank is provided with MicroB release.

* Export files path:
  * `export BLASTDB=$BLASTDB:$(pwd)` (or =$BLASTDB:$(pwd))

## Configure and run

* Complete the configuration file (identification.conf)
  !! WARNING : due to the use of an automatic parser, parameters should now be filled without single quotes
  * sample_list_path = #<path_to>/<SampleList>.csv
  * microB_params = #52.31.123.234:6000:microB_nestle:microb-user:microb-user:5435:listeria
  * contigs_super_folder = #<path_to>/data/
  * metadata_file = #<path_to>/Metadata.csv
  * bank_type = #RDP or microB (or RDP,microB)
  * bank_name = #16S_RDP_named or microB_bank_rRNA16S (or 16S_RDP_named,microB_bank_rRNA16S)
  * other parameters have already been filled : * min_identity = 99.5
                                                * min_align_len = 1400
                                                * pattern = filtered_contigs
                                                * subdirectory = assembling
                                                * num_threads = 8


* Script should be ready to run !
  * `python3 findOrganismAndReferences.py identification.conf range [typeStrainsList] [doNotCompareWithMetadataList]`
  where range is compliant to any of the following format :
    *  x..y : identify from line x to line y
    * x,y,z : identify lines x, y and z (any number of lines)
  typeStrainsList is an optional argument, used to force identification when known in advance.
  If provided it has to be a list of organism, formatted as a table with a header, with 1 row by organism with following column : 
  1=id;2=Genus;3=Species;4=subspecies(optional); e.g. NCC51;Lactobacillus;helveticus;;
  doNotCompareWithMetadataList is an optional argument, used to avoid comparison with Metadata (keep 16S results even if it is not consistent with metadata)
  If provided it has to be a list of organism, with one ID (only the number, no prefix) by row, and nothing else.


* Output data

All outputs are saved in the script main folder
 - identification hits table : 'identification_results_table_<bank_name>.csv' (one for each blast bank used)
 - completed sample list : 'completed_SampleList.csv'
 - log file : 'log_identification.txt'

(Warnings are all printed in the log file)



## Specs


### 16S identification procedure:

 1. check the format of the sample list, add a species_name column if not already present
  (content = full name with genus, species, and subsp. if needed)
 2. For each blast banks specified in the configuration file :
      For each organism_id in the sample list :
       * if strain_name is not already known
         * try to get contigs ('filtered_contigs.fasta'), if it exists (in subfolders '<id>/assembling'):
           * BLASTN contigs v database of 16S rRNA (RDP)
           * Filter output, need % of identity >= 99.5 (hits between 97.0 and 99.5 % are stored in the identification hits table)
                                 alignment length >= 1400 (hits between 500 and 1400 bp are stored in the identification hits table)
           * sort the different hits for each organism by descending % of identity, and by alignment length if there's a tie
           * Warning if several equally scoring hits or no hits at all (if no hit leave row empty and change status to 'SKIP_NO_IDENT'), else complete sample list with identified genus and species
       * save all identification hits
       * Compare identification results with already given organism's names
           * if there's a conflict between identification results and metadata :
              * leave row empty (status = 'SKIP_CONFLICT_IDENT')
           * a table referencing type strains can be supplied as an optional argument to the script, if present:
              * when there's a conflict :
               search for data in the list of type strains, and if corresponding entry is found, update species_name with data from the list of type strains
           * a list of ids not to compare can be supplied as an optional argument to the script, if present:
              * for each organism in the list, avoid comparison with metadata (keep 16S results, even if it is inconsistent with metadata)
 

### Reference search procedure :

 4. For each organism_id in the sample list, if organism is identified (genus + species) :
   * if taxonomic node for species, genus and family are not already known :
     * Find corresponding taxonomic node for species, genus and family
   * Check the presence of reference organisms in microB
     * if reference organisms exist for species taxon :
       * ref taxon = species taxon
     * else if reference organisms exist for genus taxon :
       * ref_taxon = genus_taxon + Warning if multiple reference organisms at genus level
      * else :
        * Warning & leave row empty (status = 'SKIP_NO_IDENT')
 5. save the completed sample list