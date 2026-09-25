#!/gpfs/data/yarmarkovichlab/Frank/pan_cancer/antigen_portal/spectrum_env/bin/python3.8

import sys
import os
import numpy as np
import pandas as pd
import requests
import io
import subprocess

'''
RAW ZIP MZML FASTA
SPECTRA PROTEOME
'''


API_KEY = "d839b8b34ecc41a19624bc38c05a64f3e1224e7ad25113a6bb8af00aeae0d58df36751e0edf0fe13735022e6c6ebb6d7"

# Org: NYU Grossman School of Medicine
ORG_ID = "59cbbcab-ea2b-4dfb-a679-f48df8bd6c64"
# User: liguangyuan798@gmail.com
USER_ID = "3899abde-6661-4bdd-be4d-2c4e3d28e346"
DATA_FORMAT = "FASTA"
DATA_TYPE = "PROTEOME"

# Define a custom file reader class for tracking upload progress
class ProgressFileReader(io.BufferedReader):
    def __init__(self, file_path, *args, **kwargs):
        file = open(file_path, "rb")
        self._total_size = os.path.getsize(file_path)
        self._bytes_read = 0
        self._last_reported_percentage = 0
        super().__init__(file, *args, **kwargs)

    def read(self, size=-1):
        data = super().read(size)
        self._bytes_read += len(data)
        self._report_progress()
        return data

    def _report_progress(self):
        if self._total_size == 0:
            return None

        current_percentage = int((self._bytes_read / self._total_size) * 100)

        # Report progress every 10%
        if current_percentage // 10 > self._last_reported_percentage // 10:
            print(f"Upload progress: {current_percentage}%")
            self._last_reported_percentage = current_percentage


def create_signed_upload_url(file_name):
    url = f"https://api.tesorai.com/v1/orgs/{ORG_ID}/users/{USER_ID}/upload_url"
    headers = {"Authorization": f"ApiKey {API_KEY}", "Content-Type": "application/json"}
    body = {"filename": file_name, "data_type": DATA_TYPE, "data_format": DATA_FORMAT}
    response = requests.post(url, json=body, headers=headers)
    response.raise_for_status()
    return response.json()


def put_file_to_signed_url(signed_url, file_path):
    headers = {"Content-Type": "application/octet-stream"}

    # Use our custom progress tracking file reader
    with ProgressFileReader(file_path) as f:
        response = requests.put(signed_url, data=f, headers=headers)
        response.raise_for_status()

    # Ensure 100% is printed at the end if we didn't reach it during the last chunk
    # This can happen if the last chunk doesn't cross a 10% threshold
    print("Upload progress: 100%")


def create_data_reference(data_reference_name, reference):
    url = f"https://api.tesorai.com/v1/orgs/{ORG_ID}/users/{USER_ID}/data_references"
    headers = {"Authorization": f"ApiKey {API_KEY}", "Content-Type": "application/json"}
    body = {
        "name": data_reference_name,
        "data_type": DATA_TYPE,
        "data_format": DATA_FORMAT,
        "storage_system": "GCS",
        "reference": reference,
    }
    response = requests.post(url, json=body, headers=headers)
    response.raise_for_status()


def upload_single_file(file_path):
    file_name = os.path.basename(file_path)

    print("Requesting Upload URL")
    upload_url_response = create_signed_upload_url(file_name)

    print("Uploading File")
    put_file_to_signed_url(upload_url_response["resumable_upload_session_url"], file_path)

    print("Creating DataReference")
    create_data_reference(file_name, upload_url_response["reference"])

    print("Upload Complete")


# Example usage
if __name__ == "__main__":
    # immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome/osteosarcoma/UCSF_in_house'
    # all_d = subprocess.run('find {} -type d -name "*.d" -exec echo {{}} \;'.format(immuno_dir),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    # for d in all_d:
    #     des_zip = './{}.zip'.format(d.split('/')[-1].split('.')[0])
    #     subprocess.run('zip -r {} {}'.format(des_zip,d),shell=True)
    #     upload_single_file(des_zip)

    # immuno_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/immunopeptidome/bile_duct_cancer/PXD016060'
    # all_raw = subprocess.run('find {} -type f -name "*.raw" -exec echo {{}} \;'.format(immuno_dir),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    # for raw in all_raw:
    #     upload_single_file(raw)

    fasta_dir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/OS/db_fasta_tesorai'
    all_combined_fasta = subprocess.run('find {} -type f -name "*.fasta" -exec echo {{}} \;'.format(fasta_dir),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    all_combined_fasta.append('/gpfs/data/yarmarkovichlab/public/ImmunoVerse/database/contaminants_polish.fasta')
    cmd = 'cat '
    for fasta in all_combined_fasta:
        cmd += '{} '.format(fasta)
    cmd += '> combined_os_pan_cancer.fasta'
    subprocess.run(cmd,shell=True)
    upload_single_file('combined_os_pan_cancer.fasta')