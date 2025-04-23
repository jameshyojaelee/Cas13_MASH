process DESEQ2 {
    label 'high_memory'
    
    errorStrategy = 'retry'
    maxRetries = 2
    memory = { 128.GB * task.attempt }
    cpus = 16
    
    beforeScript = '''
    module purge
    module load python/3.9.7
    '''
    
    publishDir "${params.outdir}/counts", mode: 'copy'
    
    input:
    path(counts_files)
    val(sample_metadata)
    
    output:
    path("counts_matrix.csv")
    path("sample_metadata.csv")
    
    script:
    """
    #!/usr/bin/env python3
    import os
    import glob
    import pandas as pd
    import numpy as np
    import csv
    
    # Function to extract sample name from filename
    def get_sample_name(filename):
        basename = os.path.basename(filename)
        return basename.replace('.counts.txt', '')
    
    # Parse metadata
    print("Parsing metadata...")
    metadata_str = '${sample_metadata}'.replace('[', '').replace(']', '')
    elements = [x.strip() for x in metadata_str.split(',')]
    
    metadata = []
    for i in range(0, len(elements), 6):
        if i+5 < len(elements):
            sample = elements[i].strip()
            group = elements[i+1].strip()
            disease = elements[i+2].strip()
            try:
                fibrosis = int(elements[i+3].strip())
            except:
                fibrosis = 0
            try:
                nas = int(elements[i+4].strip())
            except:
                nas = 0
            stage = elements[i+5].strip()
            metadata.append([sample, group, disease, fibrosis, nas, stage])
    
    metadata_df = pd.DataFrame(metadata, columns=['sample', 'group', 'disease', 'fibrosis', 'nas', 'stage'])
    
    # Find all count files
    count_files = glob.glob("*.counts.txt")
    print(f"Found {len(count_files)} count files")
    
    if len(count_files) == 0:
        raise Exception("No count files found")
    
    # Process the first file to get gene list
    genes_df = pd.read_csv(count_files[0], sep='\\t', comment='#')
    genes = genes_df['Geneid'].tolist()
    
    # Create empty count matrix
    count_matrix = pd.DataFrame(index=genes)
    
    # Process each count file
    for file in count_files:
        sample_name = get_sample_name(file)
        print(f"Processing {file} for sample {sample_name}")
        
        try:
            df = pd.read_csv(file, sep='\\t', comment='#')
            count_matrix[sample_name] = df.iloc[:, -1].values
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
    
    # Match samples in metadata with count matrix
    common_samples = list(set(metadata_df['sample']).intersection(set(count_matrix.columns)))
    print(f"Found {len(common_samples)} samples in both metadata and count files")
    
    if len(common_samples) == 0:
        print("ERROR: No matching samples")
        print(f"Metadata samples: {', '.join(metadata_df['sample'].head(5).tolist())}")
        print(f"Count file samples: {', '.join(list(count_matrix.columns)[:5])}")
        raise Exception("No matching samples")
    
    # Subset count matrix and metadata to common samples
    filtered_counts = count_matrix[common_samples]
    filtered_metadata = metadata_df[metadata_df['sample'].isin(common_samples)]
    
    # Save outputs
    filtered_counts.to_csv("counts_matrix.csv")
    filtered_metadata.to_csv("sample_metadata.csv", index=False)
    
    print("Successfully generated count matrix and metadata")
    """
} 