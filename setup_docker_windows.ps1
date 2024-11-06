# Create directories for Neo4j data persistence
$directories = @("neo4j_docker/data", "neo4j_docker/logs", "neo4j_docker/conf", "neo4j_docker/plugins", "neo4j_docker/import")
foreach ($dir in $directories) {
    $null = New-Item -ItemType Directory -Force -Path $dir
}

cd "neo4j_docker/import"

# Download the zenodo database flat zip
$ZENODO_URL = "https://zenodo.org/records/13929716/files/procoggraph_flat_files_v1-0-2.zip?download=1"
Invoke-WebRequest -Uri $ZENODO_URL -OutFile "procoggraph_flat_files_v1-0-2.zip"

# Unzip the database flat files into the import directory
Expand-Archive -Force -Path "procoggraph_flat_files_v1-0-2.zip" -DestinationPath .
cd "../.."

# Get the fully qualified paths for building the docker run commands
$NEO4J_DOCKER_DIR = (Get-Location).Path + "\neo4j_docker"
$PROCOGGRAPH_REPOSITORY = Get-Location

# Create the PowerShell script for the build command
$buildScriptContent = @"
docker run --name neo4j_build -p 7474:7474 -p 7687:7687 -v ${NEO4J_DOCKER_DIR}/data:/data -v ${NEO4J_DOCKER_DIR}/logs:/logs -v ${NEO4J_DOCKER_DIR}/conf:/conf -v ${NEO4J_DOCKER_DIR}/plugins:/plugins -v ${NEO4J_DOCKER_DIR}/import:/var/lib/neo4j/import -v ${PROCOGGRAPH_REPOSITORY}/nextflow/bin/import_neo4j_data.sh:/import_neo4j_data.sh -e NEO4J_AUTH=neo4j/procoggraph neo4j:latest /import_neo4j_data.sh
"@
$buildScriptContent | Out-File -FilePath "run_build.ps1" -Encoding utf8

# Create the PowerShell script for the run commands
$runScriptContent = @"
docker run -d --name neo4j_run -p 7474:7474 -p 7687:7687 -v ${NEO4J_DOCKER_DIR}/data:/data -v ${NEO4J_DOCKER_DIR}/logs:/logs -v ${NEO4J_DOCKER_DIR}/conf:/conf -v ${NEO4J_DOCKER_DIR}/plugins:/plugins -v ${NEO4J_DOCKER_DIR}/import:/var/lib/neo4j/import -v ${PROCOGGRAPH_REPOSITORY}/nextflow/bin/import_neo4j_data.sh:/import_neo4j_data.sh -e NEO4J_PLUGINS=[""apoc""] -e NEO4J_AUTH=neo4j/procoggraph neo4j:latest

docker run -d --name nginx -p 8080:80 -v ${PROCOGGRAPH_REPOSITORY}/procogdash/nginx.conf:/etc/nginx/conf.d/default.conf:ro -v ${PROCOGGRAPH_REPOSITORY}/procogdash:/usr/share/nginx/html:ro nginx:latest

docker run -d --name neodash -p 5005:5005 -e ssoEnabled=false -e standalone=true -e standaloneProtocol=bolt -e standaloneHost=localhost -e standalonePort=7687 -e standaloneDatabase=neo4j -e standaloneDashboardName=ProCogGraph -e standaloneDashboardDatabase=neo4j -e standaloneDashboardURL=http://localhost:8080/dashboard.json neo4jlabs/neodash
"@
$runScriptContent | Out-File -FilePath "run_services.ps1" -Encoding utf8

Write-Output "`nSuccessfully downloaded database files and generated PowerShell scripts. To get started, run 'run_build.ps1', then 'run_services.ps1'`n."
