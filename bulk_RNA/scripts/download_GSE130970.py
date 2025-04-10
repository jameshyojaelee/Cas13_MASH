#!/gpfs/commons/home/jameslee/miniconda3/envs/RNA/bin/python
import sys
import subprocess

def download_fastq(srr_id):
    print(f"Downloading fastq for {srr_id}...")
    cmd = "cd /gpfs/commons/home/jameslee/Cas13/bulk_RNA/fastq && fastq-dump --split-files " + srr_id
    try:
        subprocess.run(cmd, shell=True, check=True)
    except Exception as e:
        print(f"Error downloading {srr_id}: {e}")

def query_and_download():
    try:
        from pysradb import SRAweb
    except ImportError:
        sys.exit("pysradb not found. Install via: pip install pysradb")
    db = SRAweb()
    gse_id = "GSE130970"  # Fixed GSE series for this script
    print(f"Querying SRA for {gse_id}...")
    srp_df = db.gse_to_srp(gse_id)
    if srp_df.empty:
        print(f"No study accessions found for {gse_id}")
    else:
        for srp in srp_df['study_accession']:
            metadata = db.sra_metadata(srp, detailed=False)
            if metadata.empty:
                print(f"No run accessions found for study {srp}")
            else:
                for srr in metadata['run_accession']:
                    download_fastq(srr)

if __name__ == '__main__':
    query_and_download()
