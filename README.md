# ProCogGraph

![ProCogGraph Logo](images/PROCOGGRAPH%20full%20logo%20v1%20small.png)

A graph based cognate ligand-domain interaction database for exploring and mining domain-cognate ligand interactions.

## Table of Contents

- [Features](#features)
    - [Database Schema](#database-schema)
- [Installation](#installation)
  - [Docker Image](#docker-image)
  - [Building from Source](#building-from-source)
  - [ProCogDash Dashboard](#procogdash-dashboard)
  - [ProCogGraph Pipeline](#procoggraph-pipeline)
- [Usage](#usage)
  - [Starting the Database](#starting-the-database)
    - [Docker](#docker)
    - [From Source](#from-source)
    - [ProCogDash Dashboard](#procogdash-dashboard-1)
    - [Custom Queries](#custom-queries)
- [Cognate Ligands](#cognate-ligands)
  - [Aggregation](#aggregation)
  - [Similarity](#similarity)
- [Domains](#domains)
  - [Aggregation](#aggregation-1)
  - [Interaction Modes](#interaction-modes)
- [License](#license)
- [Citations](#citations)

## Features

KEY FEATURES HERE

### Database Schema

## Installation

ProCogGraph is both a pipeline for analysis of structures and a database of cognate ligand-domain mappings. Installation instructions for the database as a Docker image and from source, and the pipeline are provided below.

### ProCogGraph Database

#### Docker Image

To setup the ProCogGraph database as a Docker image, follow these steps:

1. Download the latest database dump from Zenodo [here](https://ZENODOIDHERE).

2. Download and install Docker from the [Docker website](https://www.docker.com/get-started).

3. 

#### Building from Source

To build the ProCogGraph database from source, follow these steps:

1. Download the latest database flat files from Zenodo [here](https://ZENODOIDHERE).

2. Download and install Neo4j community edition from the [Neo4j website](https://neo4j.com/download/). The database was built using Neo4j 5.4.0.

3. Clone this repository:

    ``` bash
    git clone m-crown/ProCogGraph
    ```

4. Copy the build script from the repository to the database flat file directory (e.g. ProCogGraphData).

    ``` bash
    cp -r /PATH/TO/ProCogGraph/nextflow/bin/import_neo4j_data.sh /PATH/TO/ProCogGraphData/
    ```

5. Run the build script:

    ``` bash
        cd ProCogGraphData
        ./import_neo4j_data.sh
    ```

### ProCogDash Dashboard

The ProCogDash dashboard is built using NeoDash, a Neo4j plugin, and is stored as a node within the database itself. The dashboard can be accessed by connecting to a running instance of the database from the [NeoDash](http://neodash.graphapp.io/) website or by setting up and running a local instance of NeoDash with Docker ([NeoDash Docker Guide](https://neo4j.com/labs/neodash/2.1/developer-guide/build-and-run/)).

### ProCogGraph Pipeline

The ProCogGraph pipeline is built using Nextflow for workflow management. To run the pipeline, follow these steps:

1. Clone this repository:

    ``` bash
    git clone m-crown/ProCogGraph
    ```

2. 

## Usage

### Starting the Database

#### Docker

After setting up the Docker image as described above, the ProCogGraph database can be started by running the following command:

``` bash
docker run -d -p 7474:7474 -p 7687:7687 --name procoggraph procoggraph:latest
```

The database can then be accessed by navigating to `http://localhost:7474` in a web browser.

#### From Source

If running ProCogGraph from source, the database can be started by running the following command:

``` bash
cd /path/to/neo4j
neo4j start
```

The database can then be accessed by navigating to `http://localhost:7474` in a web browser.

### ProCogDash Dashboard



### Custom Queries

Custom queries can be executed using the Cypher query language in the Neo4j browser (`http://localhost:7474` when an instance of the database is running). For example, to EXAMPLE QUERY HERE, the following query can be executed:

``` cypher
MATCH EXAMPLE QUERY HERE
RETURN EXAMPLE
```

## Cognate Ligands 

### Aggregation

COGNATE LIGANDS SOURCING HERE

### Similarity

ProCogGraph defines cognate ligand similarity using the PARITY method, which uses a permissive maximum common substructure and measures the number of matching atoms between ligands. The score ranges from 0-1 with 1 representing identical ligands and 0 representing no similarity.

EXAMPLE HERE?

A threshold value for defining a cognate ligand match is set at 0.40, based on the mean 95th percentile score for 5 sets of 2000 randomly paired ligands. Bound entity - cognate ligand pairs with a score below the threshold are not included in the database.

## Domains 

### Aggregation

DOMAIN SOURCES HERE

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

## Citations

CITATIONS HERE