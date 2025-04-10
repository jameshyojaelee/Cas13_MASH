#!/gpfs/commons/home/jameslee/miniconda3/envs/RNA/bin/python
import sys
import subprocess

def download_fastq(srr_id):
    print(f"Downloading fastq for {srr_id}...")
    # Change directory to the fastq output folder before executing fastq-dump
    cmd = f"cd /gpfs/commons/home/jameslee/Cas13/bulk_RNA/fastq && fastq-dump --split-files {srr_id}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except Exception as e:
        print(f"Error downloading {srr_id}: {e}")

def query_and_download(gse_id):
    try:
        from pysradb import SRAweb
    except ImportError:
        sys.exit("pysradb not found. Install via: pip install pysradb")
    db = SRAweb()
    print(f"Querying SRA for {gse_id}...")
    # Changed from gse_to_sra to gse_to_srp to fetch study accessions
    srp_df = db.gse_to_srp(gse_id)
    if srp_df.empty:
        print(f"No study accessions found for {gse_id}")
    else:
        for srp in srp_df['study_accession']:
            # Fetch run accessions from study accession
            metadata = db.sra_metadata(srp, detailed=False)
            if metadata.empty:
                print(f"No run accessions found for study {srp}")
            else:
                for srr in metadata['run_accession']:
                    download_fastq(srr)

if __name__ == '__main__':
    # List of GSE series as described in the README
    gse_ids = ["GSE135251", "GSE130970", "GSE126848"]
    for gse in gse_ids:
        query_and_download(gse)
