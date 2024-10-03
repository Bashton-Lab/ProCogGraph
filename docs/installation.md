# Installation

## Docker

ProCogGraph is both a pipeline for analysis of structures and a database of cognate ligand-domain mappings. To get started, the easiest method, described below, is to run ProCogGraph in a Docker container - for installation instructions for the database on bare metal, and for running the Nextflow pipeline see the [installation](docs/installation.md) guide. As part of the installation process, the latest flat files are downloaded from Zenodo. These files are currently 175.9 MB, and so the download may take some time depending on your internet connection. The total size of the database once built is approximately 4GB - ensure you have sufficient disk space available before beginning the install.

NOTE: Currently, the NeoDash Docker image does not contain a build for arm based Mac devices. There is an [open issue](https://github.com/neo4j-labs/neodash/issues/754) in NeoDash related to this, and until it is fixed by the developers, ProCogGraph cannot be setup via Docker on arm-based Mac devices. Therefore execution of Docker steps is limited to x86 Mac devices. ProCogGraph can still be installed directly on arm-based Mac devices by following the steps in the [installation](docs/installation.md) guide and using the web-hosted (by Neo4j) NeoDash web app.

1. Download and install Docker from the [Docker website](https://www.docker.com/get-started)

2. Clone the ProCogGraph repository:

    ``` bash
    git clone bashton-lab/ProCogGraph
    cd ProCogGraph
    ```

3. Run the setup script to download the latest flat files and create the necessary directories and Docker compose files if running on Linux/OSX:

    ``` bash
    ./setup_docker_linux.sh
    ```

    or for Windows (in Powershell with administrative access)

    ``` powershell
    Set-ExecutionPolicy Unrestricted
    ./setup_docker_windows.ps1
    Set-ExecutionPolicy Restricted
    ```

    This script creates the necessary directories for setting up the database, downloads the latest flat files from Zenodo and produces two yaml files, one to build the database (run first time only) and one to run the database (run each time you want to start the database).

4. Run the build command:

    ``` bash
    docker compose -f compose-build.yml up
    ```

5. Run the database:

    ``` bash
    docker compose -f compose-run.yml up
    ```

    After running the Docker Compose script, three containers are started, one for the Neo4j database, one for the NeoDash dashboard and an Nginx server which serves the iframe visualisations available within the dashboard. The database can be accessed by navigating to `http://localhost:7474` in a web browser to access the neo4j browser tool or connecting to ProCogDash via [localhost:5005](http://localhost:5005/). The compose-run.yml file can be modified to specify memory allocation for the Neo4j database, which can be adjusted as necessary for your system. Currently, these are not set by the install script, and so will operate with the memory configured in docker. To adjust these parameters add the following lines to the environment section of the compose_run.yaml file:

    ``` yaml
      - NEO4J_server_memory_heap_initial__size=3600m
      - NEO4J_server_memory_heap_max__size=3600m
      - NEO4J_server_memory_pagecache_size=2g
      - NEO4J_server_jvm_additional=-XX:+ExitOnOutOfMemoryError
    ```

6. Access the dashboard. The ProCogDash dashboard is built using NeoDash, a Neo4j plugin. The dashboard can be accessed by connecting to a running instance of the database in Docker at [localhost:5005](localhost:5005). The dashboard requires a username and password, which are set to `neo4j` and `procoggraph` by default.

## Neo4j

Installation instructions for running the database on bare metal, rather than Docker, are described below.

1. Download the latest database flat files from Zenodo [here](https://zenodo.org/records/13165852) and clone the ProCogGraph repository:

    ``` bash
    git clone Bashton-Lab/ProCogGraph
    curl https://zenodo.org/records/13165852/files/procoggraph_flat_files_v1-0.zip?download=1 -o /PATH/TO/DATABASE_FLAT_FILES/procoggraph_flat_files_v1-0.zip
    unzip /PATH/TO/DATABASE_FLAT_FILES/procoggraph_flat_files_v1-0.zip
    ```

2. Download and install Neo4j community edition from the [Neo4j website](https://neo4j.com/download/). The database was built using Neo4j version 5.

3. Copy the build script from the repository to the Neo4j database directory (e.g. neo4j-5.4.0) and the database flat files to the import directory:

    ``` bash
    cp -r /PATH/TO/PROCOGGRAPH_REPOSITORY/nextflow/bin/import_neo4j_data.sh /PATH/TO/NEO4J_DATABASE/
    cp -r /PATH/TO/DATABASE_FLAT_FILES/* /PATH/TO/NEO4J_DATABASE/import/
    ```

4. Run the build script:

    ``` bash
        cd /PATH/TO/NEO4J_DATABASE/
        ./import_neo4j_data.sh
    ```

5. Start the Neo4j database:

    ``` bash
    bin/neo4j start
    ```

6. Access the database by navigating to `http://localhost:7474` in a web browser and update the default password (set to user `neo4j` and password `neo4j` by default). 

7. Access ProCogDash via [NeoDash](http://neodash.graphapp.io/). The dashboard can be loaded into Neodash by expanding the menu option in the bottom left of the screen, clicking the `+` icon and importing the dashboard from a JSON file. Upload the file from the repository at `procogdash/dashboard.json`.

## ProCogGraph Pipeline

The ProCogGraph pipeline is built using Nextflow for workflow management. To run the pipeline, follow these steps:

1. The pipeline utilises data from a number of different sources to build the ProCogGraph database. To begin, prepare a data files directory with the following:

    | File | Description | Download |
    | ---- | ---- | ---- |
    | pdb_chain_enzyme.tsv.gz | Protein chain EC ID annotation from SIFTS for PDe structures. | [SIFTS](https://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/pdb_chain_enzyme.csv) |
    | assemblies_data.csv.gz | Assembly data for PDBe structures from PDBe-KB | [PDBe-KB](https://ftp.ebi.ac.uk/pub/databases/pdbe-kb/complexes/assemblies_data.csv) |
    | enzclass.txt | Enzyme classification hierarchy | [ExPASy](https://ftp.expasy.org/databases/enzyme/) |
    | enzyme.dat | ENZYME database records | [ExPASy](https://ftp.expasy.org/databases/enzyme/) |
    | cath-names.txt | CATH domain names | [CATH](http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/) |
    | cath-domain-description-file.txt | CATH domain descriptions | [CATH](http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/) |
    | dir.des.scop.1_75.txt | SCOP domain descriptions | [SCOP](https://scop.berkeley.edu/downloads/) |
    | dir.cla.scop.1_75.txt | SCOP domain classifications | [SCOP](https://scop.berkeley.edu/downloads/) |
    | clan_membership.txt.gz | Pfam clan membership | [InterPro](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/) |
    | clan.txt.gz | Pfam clan descriptions | [InterPro](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/) |
    | interpro.xml.gz | InterPro domain annotations | [InterPro](https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/) |
    | rhea-reaction-smiles.tsv | RHEA reaction smiles strings | [RHEA](https://www.rhea-db.org/help/download) |
    | rhea2ec.tsv | RHEA to EC number mappings | [RHEA](https://www.rhea-db.org/help/download) |
    | rhea-directions.tsv | RHEA reaction directions | [RHEA](https://www.rhea-db.org/help/download) |
    | chebi_names.tsv.gz | ChEBI names | [ChEBI](https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/) |
    | relation.tsv | ChEBI relations | [ChEBI](https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/) |
    | ChEBI_Results.tsv | ChEBI records with database cross references to KEGG GLYCAN and KEGG COMPOUND, where a structure exists for the record, generated with advanced search function. | [ChEBI](https://www.ebi.ac.uk/chebi/advancedSearchForward.do) |
    | scop2-cla-latest.txt | SCOP2 domain classifications | [EBI](https://www.ebi.ac.uk/pdbe/scop/download) |
    | scop2-des-latest.txt | SCOP2 domain descriptions | [EBI](https://www.ebi.ac.uk/pdbe/scop/download) |
    | ccd.cif | Chemical Component Dictionary Structures | [CCD](https://www.wwpdb.org/data/ccd) |
    | pubchem_substance_id_mapping.txt | PubChem substance ID mappings from PubChem search for KEGG data source. | [PubChem](https://www.ncbi.nlm.nih.gov/pcsubstance?term=%22KEGG%22%5BSourceName%5D%20AND%20hasnohold%5Bfilt%5D) |

2. Clone this repository and install dependencies:

    ``` bash
    git clone m-crown/ProCogGraph
    cd ProCogGraph
    conda env create -f nextflow/envs/environment.yml
    ```

3. Preprocess RHEA reaction files:

    ``` bash
    cd /PATH/TO/DATA_DIR/
    python3 preprocess_rhea.py --rhea_ec_mapping rhea2ec.tsv --rhea_reaction_directions rhea-directions.tsv --rd_dir rd/ --outdir . --chebi_names chebi_names.tsv.gz
    ```

4. Produce final manifest file of structures to be processed:

    ``` bash
    python3 download_mmcif.py --sifts_file /PATH/TO/DATA_DIR/pdb_chain_enzyme.tsv.gz --assemblies_file /PATH/TO/DATA_DIR/assemblies_data.csv.gz --chunk_size 100 --output_dir /PATH/TO/STRUCTURES_DIR
    ```

5. Run the nextflow pipeline:

    To configure the nextflow pipeline, the nextflow.config file within the repository should be modified. A SLURM cluster profile, specific for the development of the pipeline within the Bashton Group at Northumbria, is included called 'crick'. The standard profile is designed for running the pipeline on a local machine, and is by default configured with a large amount of memory and CPU resources. This should be adjusted before running.

    Four additional parameters must be set specific to the user's environment:

    - params.data_dir - the path to the data directory created including data files described above.

    - params.cache_in - the path to the cache directory for the pipeline, if pipeline has been run previously.

    - params.output_dir - the desired output directory.

    - params.manifest - the path to the manifest file created in step 3.

    ``` bash
    cd /PATH/TO/PROCOGGRAPH_REPOSITORY/nextflow
    nextflow run main.nf -resume -profile standard
    ```
