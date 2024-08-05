# Create directories for Neo4j data persistence
$directories = @("neo4j_docker/data", "neo4j_docker/logs", "neo4j_docker/conf", "neo4j_docker/plugins", "neo4j_docker/import")
foreach ($dir in $directories) {
    New-Item -ItemType Directory -Force -Path $dir
}

cd "neo4j_docker/import"

# Download the zenodo database flat zip (uncomment and update ZENODO_ID if needed)
# $ZENODO_ID = "IDHERE"
# $ZENODO_URL = "https://zenodo.org/record/$ZENODO_ID"
# Invoke-WebRequest -Uri $ZENODO_URL -OutFile "zenodo_database.zip"

# Unzip the database flat files into the import directory

Expand-Archive -Path "../procoggraph_flat_files_v1-0.zip" -DestinationPath .
cd "../.."

# Get the fully qualified paths for building the compose yaml file
$NEO4J_DOCKER_DIR = (Get-Location).Path + "\neo4j_docker"
$PROCOGGRAPH_REPOSITORY = Get-Location

# Create the docker-compose build YAML file
$composeBuildContent = @"
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
"@
$composeBuildContent | Out-File -FilePath "compose_build.yaml" -Encoding utf8

# Create the docker-compose run YAML file
$composeRunContent = @"
services:
  neo4j_run:
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
      - ${PROCOGGRAPH_REPOSITORY}/procogdash/web/nginx.conf:/etc/nginx/conf.d/default.conf:ro
      - ${PROCOGGRAPH_REPOSITORY}/procogdash/web:/usr/share/nginx/html:ro

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
"@
$composeRunContent | Out-File -FilePath "compose_run.yaml" -Encoding utf8

Write-Output "Successfully downloaded database files and generated docker compose files. To get started, run 'docker compose -f compose_build.yaml up', then 'docker compose -f compose_run.yaml up'."
