import pandas as pd
import os
import argparse
from pathlib import Path
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description='Create samplesheet from SRA run table and fastq files')
    
    parser.add_argument('--gse', '-g', required=True,
                       help='GSE accession number (e.g., GSE135251)')
    parser.add_argument('--sra-table', '-s', required=True,
                       help='Path to SRA run table CSV file')
    parser.add_argument('--fastq-dir', '-f', required=True,
                       help='Directory containing fastq files')
    parser.add_argument('--output', '-o',
                       help='Output path for samplesheet CSV. If not provided, will use samplesheets/{gse}_samplesheet.csv')
    parser.add_argument('--fastq-pattern', '-p', default="{srr}_{read}.fastq",
                       help='Pattern for fastq files. Use {srr} for SRR ID and {read} for read number. Default: {srr}_{read}.fastq')
    parser.add_argument('--metadata-cols', '-m', nargs='+',
                       help='Additional metadata columns to include from SRA table')
    
    args = parser.parse_args()
    
    # If output not provided, create default name with GSE in samplesheets directory
    if not args.output:
        args.output = Path("samplesheets") / f"{args.gse}_samplesheet.csv"
    
    return args

def find_fastq_files(fastq_dir, pattern, srr):
    """Find fastq files for a given SRR ID using the specified pattern"""
    r1 = Path(fastq_dir) / pattern.format(srr=srr, read="1")
    r2 = Path(fastq_dir) / pattern.format(srr=srr, read="2")
    
    # Also check for common variations
    patterns = [
        "{srr}_{read}.fastq",
        "{srr}_{read}.fq",
        "{srr}.{read}.fastq",
        "{srr}.{read}.fq",
        "{srr}_R{read}.fastq",
        "{srr}_R{read}.fq",
        "{srr}_{read}.fastq.gz",
        "{srr}_{read}.fq.gz",
        "{srr}.{read}.fastq.gz",
        "{srr}.{read}.fq.gz",
        "{srr}_R{read}.fastq.gz",
        "{srr}_R{read}.fq.gz"
    ]
    
    if not (r1.exists() and r2.exists()):
        for pat in patterns:
            r1 = Path(fastq_dir) / pat.format(srr=srr, read="1")
            r2 = Path(fastq_dir) / pat.format(srr=srr, read="2")
            if r1.exists() and r2.exists():
                break
    
    return r1, r2

def create_samplesheet(sra_table_path, fastq_dir, pattern, metadata_cols=None):
    """Create samplesheet from SRA run table and fastq files"""
    # Read the SRA run table
    try:
        sra_table = pd.read_csv(sra_table_path)
    except Exception as e:
        print(f"Error reading SRA table: {e}")
        sys.exit(1)
    
    # Verify required columns
    required_cols = ['Run']
    if not all(col in sra_table.columns for col in required_cols):
        print(f"Error: SRA table must contain these columns: {required_cols}")
        sys.exit(1)
        
    # Get list of existing fastq files
    fastq_dir = Path(fastq_dir)
    if not fastq_dir.exists():
        print(f"Error: Fastq directory does not exist: {fastq_dir}")
        sys.exit(1)
    
    # Create samplesheet
    samples = []
    for _, row in sra_table.iterrows():
        srr = row['Run']
        r1, r2 = find_fastq_files(fastq_dir, pattern, srr)
        
        if r1.exists() and r2.exists():
            sample = {
                'sample': srr,
                'fastq_1': str(r1.absolute()),
                'fastq_2': str(r2.absolute())
            }
            
            # Add requested metadata columns
            if metadata_cols:
                for col in metadata_cols:
                    if col in row:
                        sample[col] = row[col]
            
            samples.append(sample)
    
    return pd.DataFrame(samples)

def main():
    args = parse_arguments()
    
    # Create samplesheet
    df = create_samplesheet(
        args.sra_table,
        args.fastq_dir,
        args.fastq_pattern,
        args.metadata_cols
    )
    
    # Save samplesheet
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    
    print(f"Created samplesheet with {len(df)} samples at {output_path}")
    
    # Print summary statistics for metadata columns
    if args.metadata_cols:
        print("\nMetadata summaries:")
        for col in args.metadata_cols:
            if col in df.columns:
                print(f"\n{col}:")
                print(df[col].value_counts().sort_index())

if __name__ == '__main__':
    main() 