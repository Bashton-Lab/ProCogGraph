# ProCogGraph

![ProCogGraph Logo](images/PROCOGGRAPH%20full%20logo%20v1%20small_white.png)

A graph based cognate ligand-domain interaction database for exploring and mining domain-cognate ligand interactions.

## Table of Contents

- [Features](#features)
- [Cloud Instance](#Cloud-Instance)
- [Quick Start](#quick-start)
- [Database Schema](#database-schema)
- [Dashboard](#dashboard)
- [Custom Queries](#custom-queries)
- [Tutorial](#tutorial)
- [Database Information](#database-information)
  - [Cognate Ligands](#cognate-ligands)
  - [Domains](#domains)
- [License](#license)

## Features

ProCogGraph is a graph database which maps domains - cognate ligand interactions in enzyme structures in the PDBe. The database builds upon principles originally described by Bashton et al. in the PROCOGNATE database, and expands upon this database by expanding the domain databases used, including a wider range of cognate ligands, and updating the domain-ligand interaction mode and ligand similarity scoring methods.

To learn more, check out the ProCogGraph preprint on bioRxiv [here](https://www.biorxiv.org/content/10.1101/2024.08.08.607191v1).

## Cloud Instance

ProCogGraph is currently hosted at [http://procoggraphprod.uksouth.cloudapp.azure.com:5005/](http://procoggraphprod.uksouth.cloudapp.azure.com:5005/) (use password: `procoggraph`) although we can't guarantee the longevity of this instance long term. Hence we created the demo mode and end user based installation options below.

## Demo Mode

To facilitate understanding of the types of data available in ProCogGraph, a demo database is available to access online. This database can be accessed from the [NeoDash](http://neodash.graphapp.io/) web app. You will need to download the `dashboard_remote.json` file from the `procogdash` directory, or copy its contents, in order to use this mode. To access the  database, navigate to the [NeoDash](http://neodash.graphapp.io/) site and use the following parameters on the splash screen:

- Protocol: `neo4j+s`
- Hostname: `neo4j+s://4135bfe1.databases.neo4j.io`
- Port: `7687`
- Username: `neo4j`
- Password: `hSIIDVT4T8KDFMHYuDSoxKovBSiKk5xlSSZhv8n4TfQ`

Once connected, you can import the dashboard by clicking the expand menu icon in the bottom left of the screen, and selecting the `+` icon in the top right of this menu. Here, you have the option to import a dashboard, either by pasting the code from the `dashboard_remote.json` file or uploading it as a file. This will load the dashboard into the NeoDash web app, allowing you to explore the database.

A video walkthrough of this process, and the features of the dashboard, is available [here](https://www.youtube.com/watch?v=qObmQNvKPnw).

Additionally, should the free hosted version of the demo database be unavailable, the dump file for this database is available in this repo (`neo4j_demo.dump`), and can be imported into a local or cloud instance of Neo4j following the instructions available at the [Neo4j website](https://neo4j.com).

## Quick Start

ProCogGraph is both a pipeline for analysis of structures and a database of cognate ligand-domain mappings. To get started, the easiest method, described below, is to run ProCogGraph in a Docker container - for installation instructions for the database on bare metal, and for running the Nextflow pipeline see the [installation](docs/installation.md) guide. As part of the installation process, the latest flat files are downloaded from Zenodo. These files are currently 175.9 MB, and so the download may take some time depending on your internet connection. The total size of the database once built is approximately 4GB - ensure you have sufficient disk space available before beginning the install.

NOTE: Currently, the NeoDash Docker image does not contain a build for arm based Mac devices. There is an [open issue](https://github.com/neo4j-labs/neodash/issues/754) in NeoDash related to this, and until it is fixed by the developers, ProCogGraph cannot be setup via Docker on arm-based Mac devices. Therefore execution of Docker steps is limited to x86 Mac devices. ProCogGraph can still be installed directly on arm-based Mac devices by following the steps in the [installation](docs/installation.md) guide and using the web-hosted (by Neo4j) NeoDash web app.

1. Download and install Docker from the [Docker website](https://www.docker.com/get-started)

2. Clone the ProCogGraph repository:

    ``` bash
    git clone Bashton-Lab/ProCogGraph
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

    This script creates the necessary directories for setting up the database, downloads the latest flat files from Zenodo and produces two yaml files on Linux/MACOS, and two powershell files on windows, one to build the database (run first time only) and one to run the database (run each time you want to start the database).

4. Run the build command:
    Linux/MACOS:

    ``` bash
    docker compose -f compose-build.yml up
    ```

    Windows:

    ``` powershell
    ./run_build.ps1
    ```

5. Run the database:
    Linux/MACOS:

    ``` bash
    docker compose -f compose-run.yml up
    ```

    Windows:

    ``` powershell
    Set-ExecutionPolicy Unrestricted
    ./run_services_.ps1
    ```

    After running the Docker Compose script, three containers are started, one for the Neo4j database, one for the NeoDash dashboard and an Nginx server which serves the iframe visualisations available within the dashboard. The database can be accessed by navigating to `http://localhost:7474` in a web browser to access the neo4j browser tool or connecting to ProCogDash via [localhost:5005](http://localhost:5005/). On linux, the compose-run.yml file can be modified (or the `run_services.ps1` on Windows) to specify memory allocation for the Neo4j database, which can be adjusted as necessary for your system. Currently, these are not set by the install script, and so will operate with the memory configured in docker. To adjust these parameters add the following lines to the environment section of the compose_run.yaml file (or add as environment parameters in the `run_services.ps1` file on Windows):

    ``` yaml
      - NEO4J_server_memory_heap_initial__size=3600m
      - NEO4J_server_memory_heap_max__size=3600m
      - NEO4J_server_memory_pagecache_size=2g
      - NEO4J_server_jvm_additional=-XX:+ExitOnOutOfMemoryError
    ```

6. Access the dashboard. The ProCogDash dashboard is built using NeoDash, a Neo4j plugin. The dashboard can be accessed by connecting to a running instance of the database in Docker at [localhost:5005](localhost:5005). The dashboard requires a username and password, which are set to `neo4j` and `procoggraph` by default.

7. To stop the database, run the following command:

    Linux/MACOS:

    ``` bash
    docker compose -f compose-run.yml down
    ```

    Windows:

    ``` powershell
    Set-ExecutionPolicy Unrestricted
    ./stop_services.ps1
    Set-ExecutionPolicy Restricted
    ```

## Database Schema

The image below shows the schema of the ProCogGraph database, which is built using Neo4j. The database is built around the following key nodes:

![ProCogGraph Schema](images/ProCogGraphSchema_schema.png)

- Entry: A PDB structure, which contains one or more protein chains and bound entities.

- Domain: An annotated sequence of a protein chain, which is classified into a domain database.

- Bound Entity: A small molecule or oligosaccharide which interacts with a domain.

- Cognate Ligand: Represents a ligand whcih is part of an enzyme reaction, and is mapped to one or more EC numbers.

## Dashboard

*NOTE: A visual bug is currently affecting the Dashboard. When navigating through the cognate ligand page, reports from this page will remain visible when navigating to other pages. This is a known issue and is being addressed with the Neodash authors see [this issue](https://github.com/neo4j-labs/neodash/issues/936).*

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

### PDB Search Mode

To search for a structure, the PDB search box is used, and a PDB ID can be matched from any partial searches via a dropdown list. A clickable link is then presented next to the PDB search box, which takes you to the PDB visualisation mode (see image above). The PDB exploration page results are described in the table and image below:

![PDB Search Results](images/procoggraph_dashboard_pdb_view_formatted.png)

| Section | Ref | Description |
| ---- | ---- | ---- |
| Summary Report | A | A summary of the PDB structure, including the number of chains, ligands, and domains present. |
| PDB Ligand Table | B | A table of bound entities present in the structure and their cognate ligand mappings (if applicable). |
| Domain Interaction Table | C | A table of domain-bound entity interactions, the interaction mode and the number of residues involved. |
| Domain Interaction Visualisation | D | An embedded iframe visualisation of the interacting residues between bound entities and domain residues in the structure,  using PDBe-Molstar. Residues from the currently focussed domain are highlighted in blue, and other domain residues highlighted in purple |
| PARITY Score Visualisation | D | PARITY score between cognate ligands and bound entities can also be viewed, showing both the MCS match and the atom matches, making up the PARITY score - visualised with RDKit JS. |

For each domain listed in the PDB structure page, breakout links are accessible to a domain summary page (domains can also be searched for directly from the search page using the domain search box). This page includes links to the external domain annotation. The report also summarises interactions for a domain at a “group” level which varies depending on the domain database being examined, for example, in the CATH/Gene3D/SCOP/SUPERFAMILY, the group level is Superfamily, and for Pfam, it is the family level. 

The domain combinatorial interactions table lists the contexts in which a domain interacts with a ligand i.e, the other domains involved in the interaction and their interaction modes.

### Cognate/Bound Entity Search Mode

Cognate and PDB ligands can be viewed in detail through breakout from the PDB structure page view. Both pages contain a similar set of results tables including an RDKit.js visualisation of the ligand structure, a summary report detailing the cognate ligand database cross-references or the number of times a PDB ligand has been observed in the database, and the domain interactions that are observed for a ligand. Additionally, PDB ligands contain link tables to cognate ligand mappings. As with domains, PDB and cognate ligands can also be searched for directly from the search page, either by hetcode or name for PDBligands, or database ID (format DB:ID) or name for cognate ligand.

### EC Search Mode

When searching for an EC number, results are aggregated and links presented to relevent structures, cognate and PDB ligands, together with a summary of domains known to interact wit ligands In this reaction. In addition the reactions associated with the EC number within the RHEA database are visualised using RDKit-js, allowing for dynamic generation based on the reaction smiles strings associated with the EC nodes in the graph.
Every result presented in ProCogDash is generated using a Cypher query, which are contained within the neodash dashboard, and which is stored as a node within the database, allowing it to be versioned and distributed alongside the database itself. In addition to this, all queries are also made available within the repository as a single YAML file, where each query contains additional comments describing the underlying process.

## Custom Queries

Custom queries can be executed using the Cypher query language in the Neo4j browser (`http://localhost:7474` when using a local instance of the database). For example, to match all domains belonging to protein chains with the EC ID 3.2.1.1 and an exclusive interaction mode, the following query can be executed:

``` cypher
MATCH (pc)<-[:IS_IN_PROTEIN_CHAIN]-(d:domain)-[int:INTERACTS_WITH_LIGAND]->(be:boundEntity) WHERE "3.2.1.1" in pc.ecList and int.interactionMode = "exclusive"
RETURN *
```
## Tutorial

A [tutorial](docs/tutorial.md) with a walkthrough of some of these features is available in the [docs](docs/) directory.

## Database Information

### Cognate Ligands

#### Sources

Cognate ligands in ProCogGraph are aggregated from the following sources:

- [KEGG](https://kegg.jp/)
- [ChEBI](https://www.ebi.ac.uk/chebi/)
- [RHEA](https://www.rhea-db.org/)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [GlyTouCan](https://glytoucan.org/)

SMILES representations are obtained for each ligand, and each cognate ligand is mapped to one or more EC IDs. Cognate ligands are processed using the RDKit library in Python, with structures neutralised and re-canonicalised to reduce the number of duplicate structures. A total of XXX cognate ligands are currently represented in the database.

#### Similarity

ProCogGraph defines cognate ligand similarity using the PARITY method, which uses a permissive maximum common substructure and measures the number of matching atoms between ligands. The score ranges from 0-1 with 1 representing identical ligands and 0 representing no similarity.

A threshold value for defining a cognate ligand match is set at 0.40, based on the mean 95th percentile score for 5 sets of 2000 randomly paired ligands. Bound entity-cognate ligand pairs with a score below the threshold are not included in the database.

### Domains

#### Domain Annotation Sources

ProCogGraph uses domain annotations from SIFTS to describe domains. The following domain databases are included:

- [CATH](https://www.cathdb.info/)

- [Gene3D](http://gene3d.biochem.ucl.ac.uk/)

- [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/#table)

- [SCOP](https://scop.berkeley.edu/ver=1.75)

- [SUPERFAMILY](https://supfam.org/SUPERFAMILY/)

- [SCOP2](https://www.ebi.ac.uk/pdbe/scop) (Split into SCOP2-Superfamily and SCOP2-Family)

Interactions between domains and ligands are considered independently for each domain database included. The dashboard allows users to select a specific domain database to search against.

#### Interaction Modes

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
