# Installation on Bare Metal

ProCogGraph is both a pipeline for analysis of structures and a database of cognate ligand-domain mappings. Installation instructions for running the database on bare metal, rather than Docker, are described below.

1. Download the latest database flat files from Zenodo [here](https://ZENODOIDHERE) and clone the ProCogGraph repository:

    ``` bash
    git clone m-crown/ProCogGraph
    wget https://zenodo.org/record/IDHERE -O /PATH/TO/DATABASE_FLAT_FILES
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



    NOTE: Running on Docker for Windows requires specific formatting of file paths for volume mappings. An example valid command is shown below:

    ``` bash
    docker run --volume=C:/Users/Matt/ProCogGraph/neo4j_docker/data:/data `
    --volume=C:/Users/Matt/ProCogGraph/neo4j_docker/logs:/logs `
    --volume=C:/Users/Matt/ProCogGraph/neo4j_docker/conf:/conf `
    --volume=C:/Users/Matt/ProCogGraph/neo4j_docker/plugins:/plugins `
    --volume=C:/Users/Matt/ProCogGraph/neo4j_docker/import:/var/lib/neo4j/import `
    --volume=C:/Users/Matt/ProCogGraph/nextflow/bin/import_neo4j_data.sh:/import_neo4j_data.sh `
    neo4j:latest /import_neo4j_data.sh
    ```

### ProCogGraph Pipeline

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

#### From Source

If running ProCogGraph from source, the database can be started by running the following command:

``` bash
cd /PATH/TO/NEO4J_DATABASE/
bin/neo4j start
```

The database can then be accessed by navigating to `http://localhost:7474` in a web browser or connecting to ProCogDash via [NeoDash](http://neodash.graphapp.io/).