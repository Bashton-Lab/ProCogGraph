#!/usr/bin/env bash

# Make the directories required for running and persisting data in the neo4j database.
mkdir -p neo4j_docker/data
mkdir -p neo4j_docker/logs
mkdir -p neo4j_docker/conf
mkdir -p neo4j_docker/plugins
mkdir -p neo4j_docker/import

cd neo4j_docker/import

# Download the zenodo database flat zip.
# ZENODO_ID="IDHERE"
# ZENODO_URL="https://zenodo.org/record/${ZENODO_ID}"
# wget "${ZENODO_URL}"

# Unzip the database flat files into the import directory, to be used by the import script.
unzip procoggraph_flat_files_v1-0.zip
cd ../../

# Specify fully qualified paths for building the compose yaml file (required for running docker volumes).
NEO4J_DOCKER_DIR=$(pwd)/neo4j_docker
PROCOGGRAPH_REPOSITORY=$(pwd)

# Create the docker-compose.yaml file
cat <<EOF > compose_build.yaml
services:
  neo4j_build:
    image: neo4j:latest
    container_name: neo4j
    ports:
      - "7474:7474"
      - "7687:7687"
    volumes:
      - ${NEO4J_DOCKER_DIR}/data:/data
      - ${NEO4J_DOCKER_DIR}/logs:/logs
      - ${NEO4J_DOCKER_DIR}/conf:/conf
      - ${NEO4J_DOCKER_DIR}/plugins:/plugins
      - ${NEO4J_DOCKER_DIR}/import:/var/lib/neo4j/import
      - ${PROCOGGRAPH_REPOSITORY}/nextflow/bin/import_neo4j_data.sh:/import_neo4j_data.sh
    environment:
      - NEO4J_PLUGINS=["apoc"]
      - NEO4J_server_memory_heap_initial__size=3600m
      - NEO4J_server_memory_heap_max__size=3600m
      - NEO4J_server_memory_pagecache_size=2g
      - NEO4J_server_jvm_additional=-XX:+ExitOnOutOfMemoryError
    entrypoint: ["/bin/bash", "/import_neo4j_data.sh"]
EOF

cat <<EOF > compose_run.yaml

services:
  neo4j_run:
    image: neo4j:latest
    container_name: neo4j
    ports:
      - "7474:7474"
      - "7687:7687"
    volumes:
      - /Users/matthewcrown/GitHub/ProCogGraph/neo4j_docker/data:/data
      - /Users/matthewcrown/GitHub/ProCogGraph/neo4j_docker/logs:/logs
      - /Users/matthewcrown/GitHub/ProCogGraph/neo4j_docker/conf:/conf
      - /Users/matthewcrown/GitHub/ProCogGraph/neo4j_docker/plugins:/plugins
      - /Users/matthewcrown/GitHub/ProCogGraph/neo4j_docker/import:/var/lib/neo4j/import
      - /Users/matthewcrown/GitHub/ProCogGraph/nextflow/bin/import_neo4j_data.sh:/import_neo4j_data.sh
    environment:
      - NEO4J_PLUGINS=["apoc"]
      - NEO4J_server_memory_heap_initial__size=4600m
      - NEO4J_server_memory_heap_max__size=4600m
      - NEO4J_server_memory_pagecache_size=3g
      - NEO4J_server_jvm_additional=-XX:+ExitOnOutOfMemoryError

  nginx:
    image: nginx:latest
    container_name: nginx
    ports:
      - "8080:80"
    volumes:
      - /Users/matthewcrown/GitHub/ProCogGraph/procogdash/web/nginx.conf:/etc/nginx/conf.d/default.conf:ro
      - /Users/matthewcrown/GitHub/ProCogGraph/procogdash/web:/usr/share/nginx/html:ro

  neodash:
    image: neo4jlabs/neodash
    container_name: neodash
    ports:
      - "5005:5005"
    environment:
      - ssoEnabled=false
      - standalone=true
      - standaloneProtocol=bolt
      - standaloneHost=localhost
      - standalonePort=7687
      - standaloneDatabase=neo4j
      - standaloneDashboardName=http://localhost:8080/dashboard.json
    stdin_open: true
    tty: true
EOF

echo "Successfully downloaded database files and generated docker compose files. To get started, run 'docker compose -f compose_build.yaml up', then 'docker compose -f compose_run.yaml up'."