#/usr/bin/env python
from ftplib import FTP
import pandas as pd
import argparse
import requests
from gemmi import cif
import tarfile
from io import BytesIO
import io
import time
import os
import math
from pathlib import Path

def process_pdb_chunk(pdb_ids):
    cwd = os.getcwd()

    prexisting_updated = []
    updated_mmcif_ids = []
    print(pdb_ids)
    for pdb in pdb_ids:
        if os.path.exists(f"{pdb}_updated.cif") or os.path.exists(f"{pdb}_updated.cif.gz"): #if the file is already downloaded, we don't need to download it again
            prexisting_updated.append(pdb)
    print(f"Prexisting updated structures: {len(prexisting_updated)}")
    pdb_ids = list(set(pdb_ids) - set(prexisting_updated))
    print(f"Missing updated structures: {len(pdb_ids)}")
    if len(pdb_ids) > 0:
        print("Downloading updated mmcif files")
        query_list = [f"id={pdb_id}" for pdb_id in pdb_ids]

        query_string = "&".join(query_list)
        updated_mmcif_url = f"https://www.ebi.ac.uk/pdbe/download/api/pdb/entry/updated?{query_string}"

        response = requests.get(updated_mmcif_url)
        while response.status_code == 202:
            time.sleep(1)
            response = requests.get(updated_mmcif_url)
        if response.status_code == 200:
            data = response.json()
            download_url = data['url']  #URl to the archive file
            file_response = requests.get(download_url)
            while file_response.status_code == 202:
                time.sleep(10) #wait 10 seconds before trying again - to give time for download to be ready (10 seconds to prevent spamming the server with requests)
                file_response = requests.get(download_url)
            if file_response.status_code == 200:
                with tarfile.open(fileobj=BytesIO(file_response.content), mode='r:gz') as tar:
                    tar.extractall()
                print("Extraction completed successfully.")
            else:
                print("Failed to download the tar archive.")
        else:
            print("Failed to get data from the API.")
            
        with open('contains.txt', 'r') as file:
            updated_mmcif_ids = [line.strip() for line in file]
            if len(updated_mmcif_ids) < len(pdb_ids):
                print(f"Only {len(updated_mmcif_ids)} out of {len(pdb_ids)} updated mmcif files were downloaded")
            else:
                print(f"All {len(updated_mmcif_ids)} updated mmcif files were downloaded")
        os.remove('contains.txt')

    updated_mmcif_ids = updated_mmcif_ids + prexisting_updated

    manifest = pd.DataFrame([{"pdb_id": pdb, "updated_mmcif": f"{cwd}/{pdb}_updated.cif"} for pdb in updated_mmcif_ids])

    #remove the contains.txt file after reading 
    
    #before running the modelserver api call, check if the protonated structure already exists
    #with large (200 record) chunks we experience 504 errors, so to reduce the chance of time outs, we sub chunk this modelserver query into chunks of 50 records at a time
    prexisting_protonated = []
    for pdb_id in updated_mmcif_ids:
        if os.path.exists(f"{pdb_id}_bio-h.cif") or os.path.exists(f"{pdb_id}_bio-h.cif.gz"): #if the file is already downloaded, we don't need to download it again
            prexisting_protonated.append(pdb_id)
    print(f"Prexisting protonated structures: {len(prexisting_protonated)}")
    updated_mmcif_ids = list(set(updated_mmcif_ids) - set(prexisting_protonated))
    print(f"Missing protonated structures: {len(updated_mmcif_ids)}")
    protonated_ids = make_modelserver_query(updated_mmcif_ids, chunk_size = 50)
    if len(protonated_ids) < len(updated_mmcif_ids):
        print(f"Only {len(protonated_ids)} out of {len(updated_mmcif_ids)} protonated structures were downloaded")
    protonated_ids.extend(prexisting_protonated)

    protonated_manifest = pd.DataFrame([{"pdb_id": pdb, "protonated_assembly": f"{cwd}/{pdb}_bio-h.cif"} for pdb in protonated_ids])

    manifest = manifest.merge(protonated_manifest, on="pdb_id", how="left")

    #get the sifts xml files.
    host = 'ftp.ebi.ac.uk'
    remote_dir = '/pub/databases/msd/sifts/xml'
    local_dir = os.getcwd()
    xml_manifest = download_sifts_xml(host, remote_dir, pdb_ids, local_dir)

    manifest = manifest.merge(xml_manifest, on="pdb_id", how="left")

    return manifest


def make_modelserver_query(pdb_ids, chunk_size = 50):
    protonated_id_list = []

    for i in range(0, len(pdb_ids), chunk_size):
        protonated_query_ids = pdb_ids[i:i+chunk_size]
        print(f"Protonated Chunk {i}") 
        protonated_assembly_query_list = [{"entryId": pdb, "query": "full", "data_source":"pdb-h"} for pdb in protonated_query_ids]
        protonated_assembly_query_json = {"queries": protonated_assembly_query_list}

        protonated_assembly_response = requests.post('https://www.ebi.ac.uk/pdbe/model-server/v1/query-many',
                            json=protonated_assembly_query_json,
                            headers={'accept': 'text/plain', 'Content-Type': 'application/json'})

        if protonated_assembly_response.status_code == 200:
            file_object = protonated_assembly_response.text
            cif_file = cif.read_string(file_object)
            failed_ids = []
            for block in cif_file:
                if block.find_pair("_model_server_error.message") is not None:
                    print(f"Error in {block.name}:\n{block.find_pair('_model_server_error.message')[1]}")
                    failed_ids.append(block.name) #we need to do something with these ids
                else:
                    pdb_id = block.find_pair("_entry.id")[1].lower()
                    protonated_id_list.append(pdb_id)
                    block.write_file(f"{pdb_id}_bio-h.cif")
            print("Response processed")
        elif protonated_assembly_response.status_code == 503:
            #if the request fails, try again at normal chunk size
            protonated_ids = make_modelserver_query(protonated_query_ids, chunk_size = chunk_size)
            protonated_id_list.extend(protonated_ids)
        elif protonated_assembly_response.status_code in [502, 504]:
            #process the chunk in smaller groups
            print(f"API Query failed error {protonated_assembly_response.status_code}, trying new chunk size {math.ceil(chunk_size/2)}")
            sub_ids = list(set(protonated_query_ids) - set(protonated_id_list))
            protonated_ids = make_modelserver_query(list(sub_ids), chunk_size = math.ceil(chunk_size/2))
            protonated_id_list.extend(protonated_ids)
        else:
            # Print an error message if the request failed
            print(f"Request failed with status code {protonated_assembly_response.status_code} on ids: {protonated_query_ids}")

    return protonated_id_list

def download_sifts_xml(host, remote_dir, pdb_ids, local_dir):
    ftp = FTP(host)
    ftp.login()

    # Change to the remote directory
    ftp.cwd(remote_dir)
    pdb_file_paths = []
    pdb_file_ids = []
    # Iterate over pdb ids and download each file
    for pdb in pdb_ids:
        file_name = f"{pdb}.xml.gz"  # Modify the file extension as needed
        local_file_path = os.path.join(local_dir, file_name)
        if os.path.exists(local_file_path):
            pdb_file_paths.append(local_file_path)
            pdb_file_ids.append(pdb)
            continue
        try:
            with open(local_file_path, 'wb') as local_file:
                ftp.retrbinary(f"RETR {file_name}", local_file.write)
            pdb_file_paths.append(local_file_path)
            pdb_file_ids.append(pdb)
        except:
            print(f"Error downloading {file_name}")
            pdb_file_paths.append(None)
            pdb_file_ids.append(pdb)

    ftp.quit()

    xml_manifest = pd.DataFrame({"pdb_id": pdb_file_ids, "sifts_xml": pdb_file_paths})

    return xml_manifest

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sifts_file', type=str, help='')
    parser.add_argument('--assemblies_file', type=str, help='')
    parser.add_argument('--chunk_size', type=int, help='number of pdb ids to process at a time')
    parser.add_argument('--output_dir', default = ".", type=str, help='output directory')
    args = parser.parse_args()
    #specify an output dir in the args, and switch to this output dir before starting.
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    os.chdir(args.output_dir)

    sifts_file = pd.read_csv(args.sifts_file, sep = "\t", skiprows = 1)
    assemblies_file = pd.read_csv(args.assemblies_file, usecols = ["ASSEMBLIES","PREFERED_ASSEMBLIES"]) #from https://ftp.ebi.ac.uk/pub/databases/pdbe-kb/complexes/, see https://www.biorxiv.org/content/10.1101/2023.05.15.540692v1.full
    assemblies_file[["PDB", "ASSEMBLY_ID"]] = assemblies_file["ASSEMBLIES"].str.split("_", expand = True)
    assemblies_file_prefered = assemblies_file.loc[assemblies_file.PREFERED_ASSEMBLIES == 1] #bool true
    sifts_file_ec = sifts_file.loc[~sifts_file.EC_NUMBER.isin(["?"]), ["PDB"]].drop_duplicates()
    sifts_assemblies = sifts_file_ec.merge(assemblies_file_prefered, how = "left", on = "PDB", indicator = True)
    assert(len(sifts_assemblies.loc[sifts_assemblies._merge != "both"]) == 0)
    sifts_assemblies.drop(columns = "_merge", inplace = True)

    query_list = sifts_assemblies.PDB.tolist()
    manifests = []
    print(f"Processing {len(query_list)} PDB entries in {args.chunk_size} chunks")
    for i in range(0, len(query_list), args.chunk_size):
        print(f"Processing structures {i}-{i+args.chunk_size} of {len(query_list)}")
        chunk = query_list[i:i+args.chunk_size]
        max_retries = 3
        retry = 1
        while retry <= max_retries:
            try:
                manifest = process_pdb_chunk(chunk)
                manifests.append(manifest)
                retry = max_retries + 1
            except Exception as e:
                print(f"Error processing chunk {i}-{i+args.chunk_size}: {e}")
                retry += 1
    
    combined_manifest = pd.concat(manifests)
    final_manifest = sifts_assemblies.merge(combined_manifest, left_on="PDB", right_on="pdb_id", how="left")
    final_manifest.drop(columns = ["pdb_id"], inplace = True)
    final_manifest.to_csv("final_manifest.csv", index=False)

    
if __name__ == "__main__":
    main()