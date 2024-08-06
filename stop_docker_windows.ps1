# Stop and remove the running containers

$containers = @("neo4j_run", "nginx", "neodash")

foreach ($container in $containers) {
    try {
        Write-Output "Stopping container: $container"
        docker stop $container
        Write-Output "Removing container: $container"
        docker rm $container
    } catch {
        Write-Output "Error stopping or removing container: $container. It might not be running."
    }
}

Write-Output "`nAll specified containers have been stopped and removed.`n"