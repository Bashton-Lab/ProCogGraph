# ProCogGraph

![ProCogGraph Logo](images/PROCOGGRAPH%20full%20logo%20v1%20small.png)

A graph based cognate ligand-domain interaction database for exploring and mining domain-cognate ligand interactions.

## Table of Contents

- [Features](#features)
  - [Database Schema](#database-schema)
- [Installation](#installation)
  - [Docker Image](#docker-image)
  - [Building directly](#building-directly)
  - [ProCogDash Dashboard](#procogdash-dashboard)
  - [ProCogGraph Pipeline](#procoggraph-pipeline)
- [Usage](#usage)
  - [Starting the Database](#starting-the-database)
    - [Docker](#docker)
    - [From Source](#from-source)
    - [Dashboard](#dashboard)
      - [PDB Search Mode](#pdb-search-mode)
      - [Cognate/Bound Entity Search Mode](#cognatebound-entity-search-mode)
      - [EC Search Mode](#ec-search-mode)
    - [Custom Queries](#custom-queries)
- [Cognate Ligands](#cognate-ligands)
  - [Sources](#sources)
  - [Similarity](#similarity)
- [Domains](#domains)
  - [Domain Annotation Sources](#domain-annotation-sources)
  - [Interaction Modes](#interaction-modes)
- [License](#license)
- [Citations](#citations)

## Features

ProCogGraph is a graph database which maps domains - cognate ligand interactions in enzyme structures in the PDBe. The database builds upon principles originally described by Bashton et al. in the PROCOGNATE database, and expands upon this database by expanding the domain databases used, including a wider range of cognate ligands, and updating the domain-ligand interaction mode and ligand similarity scoring methods.

To learn more, check out the ProCogGraph preprint on bioRxiv [here](LINKHEREWHENSUBMITTED).

### Database Schema

Figure X below shows the schema of the ProCogGraph database, which is built using Neo4j. The database is built around the following key nodes:

![ProCogGraph Schema](images/ProCogGraphSchemaLatest_schema.png)

- Entry: A PDB structure, which contains one or more protein chains and bound entities.

- Domain: An annotated sequence of a protein chain, which is classified into a domain database.

- Bound Entity: A small molecule or oligosaccharide which interacts with a domain.

- Cognate Ligand: Represents a ligand whcih is part of an enzyme reaction, and is mapped to one or more EC numbers.

A dashboard has been built to access the data stored within the database. To learn more see [the docs](#dashboard)

## Installation

ProCogGraph is both a pipeline for analysis of structures and a database of cognate ligand-domain mappings. Installation instructions for the database as a Docker image and from source, and the pipeline are provided below.

### ProCogGraph Database

1. For both the Docker image and building from source, the first step is to download the latest database flat files from Zenodo [here](https://ZENODOIDHERE) and clone the ProCogGraph repository:

    ``` bash
    mkdir 
    wget https://zenodo.org/record/IDHERE -O /PATH/TO/DATABASE_FLAT_FILES
    git clone m-crown/ProCogGraph
    ```

#### Docker Image

To setup the ProCogGraph database as a Docker image, follow these steps:

1. Download and install Docker from the [Docker website](https://www.docker.com/get-started) and pull the latest Neo4j image:

    ``` bash
    docker pull neo4j:latest
    ```

2. Prepare a database directory for persistant storage of data from the database:

    ``` bash
    mkdir /PATH/TO/NEO4J_DOCKER_DIR
    mkdir /PATH/TO/NEO4J_DOCKER_DIR/data
    mkdir /PATH/TO/NEO4J_DOCKER_DIR/logs
    mkdir /PATH/TO/NEO4J_DOCKER_DIR/conf
    mkdir /PATH/TO/NEO4J_DOCKER_DIR/plugins
    mkdir /PATH/TO/NEO4J_DOCKER_DIR/import
    ```

3. Run the import command:

    ``` bash
    docker run --volume=/PATH/TO/NEO4J_DOCKER_DIR/data:/data \                        
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/data/logs:/logs \
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/data/conf:/conf \
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/data/plugins:/plugins \
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/data/import:/var/lib/neo4j/import \
    --volume=/PATH/TO/PROCOGGRAPH_REPOSITORY/nextflow/bin/import_neo4j_data.sh:/import_neo4j_data.sh \       
    neo4j:latest /import_neo4j_data.sh
    ```

#### Building directly

To build the ProCogGraph database from source, follow these steps:

1. Download and install Neo4j community edition from the [Neo4j website](https://neo4j.com/download/). The database was built using Neo4j version 5.

2. Copy the build script from the repository to the Neo4j database directory (e.g. neo4j-5.4.0) and the database flat files to the import directory:

    ``` bash
    cp -r /PATH/TO/PROCOGGRAPH_REPOSITORY/nextflow/bin/import_neo4j_data.sh /PATH/TO/NEO4J_DATABASE/
    cp -r /PATH/TO/DATABASE_FLAT_FILES/* /PATH/TO/NEO4J_DATABASE/import/
    ```

3. Run the build script:

    ``` bash
        cd /PATH/TO/NEO4J_DATABASE/
        ./import_neo4j_data.sh
    ```

#### ProCogDash Dashboard

The ProCogDash dashboard is built using NeoDash, a Neo4j plugin, and is stored as a node within the database itself. The dashboard can be accessed by connecting to a running instance of the database from the [NeoDash](http://neodash.graphapp.io/) website or by setting up and running a local instance of NeoDash with Docker ([NeoDash Docker Guide](https://neo4j.com/labs/neodash/2.1/developer-guide/build-and-run/)).

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

2. Clone this repository:

    ``` bash
    git clone m-crown/ProCogGraph
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

## Usage

### Starting the Database

#### Docker

After setting up the Docker image as described above, the ProCogGraph database can be started by running the following command:

``` bash
docker run --name procoggraph_docker -p7474:7474 -p7687:7687 --env NEO4J_PLUGINS='["apoc"]' \
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/data:/data \
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/logs:/logs \
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/conf:/conf \
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/plugins:/plugins \
    --volume=/PATH/TO/NEO4J_DOCKER_DIR/import:/var/lib/neo4j/import \
    neo4j:latest
```

The database can then be accessed by navigating to `http://localhost:7474` in a web browser or connecting to ProCogDash via [NeoDash](http://neodash.graphapp.io/).

Depending on your system memory, the Neo4j database memory settings may need to be adjusted. The following command can be run to get suggestions for your specific system (e.g. here with 8gb):

``` bash
neo4j_docker % docker run -p7474:7474 -p7687:7687 \
--volume=/PATH/TO/NEO4J_DOCKER_DIR/data:/data \
--volume=/PATH/TO/NEO4J_DOCKER_DIR/logs:/logs \
--volume=/PATH/TO/NEO4J_DOCKER_DIR/conf:/conf \
--volume=/PATH/TO/NEO4J_DOCKER_DIR/plugins:/plugins \
--env NEO4J_PLUGINS='["apoc"]' \
neo4j:latest bin/neo4j-admin server memory-recommendation --memory=8g --docker
```

This will output a number of memory settings which can be added to the `docker run` command to adjust the memory settings for the Neo4j database. The overall command will look something like this:

``` bash
neo4j_docker % docker run \
-p7474:7474 -p7687:7687 \
--volume=/PATH/TO/NEO4J_DOCKER_DIR/data:/data \
--volume=/PATH/TO/NEO4J_DOCKER_DIR/logs:/logs \
--volume=/PATH/TO/NEO4J_DOCKER_DIR/conf:/conf \
--volume=/PATH/TO/NEO4J_DOCKER_DIR/plugins:/plugins \
--env NEO4J_PLUGINS='["apoc"]' \
--env NEO4J_server_memory_heap_initial__size='3600m' \
--env NEO4J_server_memory_heap_max__size='3600m' \
--env NEO4J_server_memory_pagecache_size='2g' \
--env NEO4J_server_jvm_additional='-XX:+ExitOnOutOfMemoryError' \
neo4j:latest
```

#### From Source

If running ProCogGraph from source, the database can be started by running the following command:

``` bash
cd /PATH/TO/NEO4J_DATABASE/
bin/neo4j start
```

The database can then be accessed by navigating to `http://localhost:7474` in a web browser or connecting to ProCogDash via [NeoDash](http://neodash.graphapp.io/).

### Dashboard

The dashboard contains five key visualisation modes: PDB, Cognate Ligand, PDB Ligand, Domain and EC. The homepage provides summary statistics for the number of structures and ligands represented in the current version of the graph, as well as the number of cognate ligand matches for the currently specified cutoff.

From the search page, global and visualisation specific settings can be specified:

| Parameter | Description | Options |
| ---- | ---- | ---- |
| Cutoff | The minimum similarity score for a cognate ligand match to be considered. | 0.0 - 1.0 |
| Domain Database | The domain database to be searched against. | CATH, Gene3D, Pfam, SCOP, SUPERFAMILY, SCOP2-SF , SCOP2-FA |
| Cognate Ligand Filter | The type of cognate ligand matches to be considered. | All , Best , Any (see below) |

When filtering cognate ligand matches in the database, users can select All, Best or Any. These parameters are described below:

- All: all bound entities for a PDB are displayed in the interaction table, regardless of whether they have a mapping to a cognate ligand.
- Any: all cognate ligand matches for a bound entity which are above the scoring threshold are presented.
- Best: only cognate ligands with the highest score for a given bound entity are shown.

It should be noted that even when set to “Best”, a bound entity may have matches to multiple cognate ligands with the same maximum score. ProCogGraph is designed to serve as an information source, and so does not make an effort to select a particular best match as the “Best” best match, instead leaving this up to the user.

#### PDB Search Mode

To search for a structure, the PDB search box is used, and a PDB ID can be matched from any partial searches via a dropdown list. A clickable link is then presented next to the PDB search box, which takes you to the PDB visualisation mode (see image above). The PDB exploration page results are described in the table and image below:

![PDB Search Results](images/procoggraph_dashboard_pdb_view_formatted.png)

| Section | Ref | Description |
| ---- | ---- | ---- |
| Summary Report | A | A summary of the PDB structure, including the number of chains, ligands, and domains present. |
| Domain Interaction Table | B | A table of interactions between domains and bound entities in the structure. |
| PDB Ligand Table | C | A table of bound entities present in the structure. |
| Domain Interaction Visualisation | D | An embedded iframe visualisation of the interacting residues between bound entities and domain residues in the structure,  using PDBe-Molstar. Residues from the currently focussed domain are highlighted in blue, and other domain residues highlighted in purple |
| PARITY Score Visualisation | E | PARITY score between cognate ligands and bound entities can also be viewed, showing both the MCS match and the atom matches, making up the PARITY score - visualised with RDKit JS. |

The iframe visualisations are loaded from GitHub pages (https://m-crown.github.io/ProCogGraph/) if using the online version of the graph or from an Nginx server running in the distributed Docker image if running locally.

For each domain listed in the PDB structure page, breakout links are accessible to a domain summary page (domains can also be searched for directly from the search page using the domain search box). This page includes a summary report of the number of ligands the domain is known to interact with, together with links to the external domain annotation. The report also summarises interactions for a domain at a “group” level which varies depending on the domain database being examined, for example, in the CATH/Gene3D/SCOP/SUPERFAMILY, the group level is Superfamily, and for Pfam, it is the family level. Summaries are presented on the following:
Group interactions table: lists all cognate ligands a group level are known to interact with in the database, together with the number of domains that interact with the ligand.
Domain Contexts: This query lists the contexts in which a domain interacts with a ligand i.e, the other domains involved in the interaction and their interaction modes.
Domain Cognate Ligand breakdowns: table lists the cognate ligands the specific domain searched for interacts with, and the percentage of the overall group the domain belongs to which also interact with the ligand. This is useful for identify if the cognate ligand a domain binds to reflects typical superfamily activity or is an outlier.

#### Cognate/Bound Entity Search Mode

Cognate and PDB ligands can be viewed in detail through breakout from the PDB structure page view. Both pages contain a similar set of results tables including an RDKit.js visualisation of the ligand structure, a summary report detailing the cognate ligand database cross-references or the number of times a PDB ligand has been observed in the database, and the domain interactions that are observed for a ligand. Additionally, PDB ligands contain link tables to cognate ligand mappings. As with domains, PDB and cognate ligands can also be searched for directly from the search page, either by hetcode or name for PDBligands, or database ID (format DB:ID) or name for cognate ligand.

#### EC Search Mode

When searching for an EC number, results are aggregated and links presented to relevent structures, cognate and PDB ligands, together with a summary of domains known to interact wit ligands In this reaction. In addition the reactions associated with the EC number within the RHEA database are visualised using RDKit-js, allowing for dynamic generation based on the reaction smiles strings associated with the EC nodes in the graph.
Every result presented in ProCogDash is generated using a Cypher query, which are contained within the neodash dashboard, and which is stored as a node within the database, allowing it to be versioned and distributed alongside the database itself. In addition to this, all queries are also made available within the repository as a single YAML file, where each query contains additional comments describing the underlying process.

### Custom Queries

Custom queries can be executed using the Cypher query language in the Neo4j browser (`http://localhost:7474` when using a local instance of the database). For example, to EXAMPLE QUERY HERE, the following query can be executed:

``` cypher
MATCH EXAMPLE QUERY HERE
RETURN EXAMPLE
```

## Cognate Ligands

### Sources

Cognate ligands in ProCogGraph are aggregated from the following sources:

- KEGG
- ChEBI
- RHEA
- PubChem
- GlyTouCan

SMILES representations are obtained for each ligand, and each cognate ligand is mapped to one or more EC IDs. Cognate ligands are processed using the RDKit library in Python, with structures neutralised and re-canonicalised to reduce the number of duplicate structures. A total of XXX cognate ligands are currently represented in the database.

### Similarity

ProCogGraph defines cognate ligand similarity using the PARITY method, which uses a permissive maximum common substructure and measures the number of matching atoms between ligands. The score ranges from 0-1 with 1 representing identical ligands and 0 representing no similarity.

EXAMPLE HERE?

A threshold value for defining a cognate ligand match is set at 0.40, based on the mean 95th percentile score for 5 sets of 2000 randomly paired ligands. Bound entity - cognate ligand pairs with a score below the threshold are not included in the database.

## Domains

### Domain Annotation Sources

ProCogGraph uses domain annotations from SIFTS to describe domains. The following domain databases are included:

- CATH

- Gene3D

- Pfam

- SCOP

- SUPERFAMILY

- SCOP2 (Split into SCOP2-Superfamily and SCOP2-Family)

Interactions between domains and ligands are considered indepdently for each domain database included. The dashboard allows users to select a specific domain database to search against.

### Interaction Modes

In ProCogGraph, domain ownership is represented in three different ways based on domain contact percentage (as determined using PDBe-Arpeggio), depending on the number of interaction domains (single, dual and multi-domain interactions):

| # Domains | Contact % | Type | Description |
| ---- | ---- | ---- | ---- |
| 1 | 100% | Exclusive | A single domain contacts the ligand. |
| 2+ | 90+% | Dominant | Two or more domains interact, and this domain dominates in terms of contact % |
| 2+ | <10% | Minor | Two or more domains interact, and this domain only plays a minor role in the interaction interface |
| 2+ | 10-90% Contacts | Major | Two + domains interact, and this domain is the only domain with more than 10% contacts (i.e. other domains are all minor) |
| 2+ | 50-90% Contacts | Major Partner | Two or more domains interact, and there is more than 1 non-minor domains. This domain contributes more than half of the total contacts measured for all domains. |
| 2+ | 10-50% Contacts | Partner | Two or more domains interact, and there is more than 1 non-minor domains. This domain has a non-minor contribution, but provides less than 50% of the total contacts measured for all domains. |

## License

ProCogGraph is licensed under the MIT License.
