import pandas as pd
import os
from multiprocessing import Pool, cpu_count
from functools import partial
import numpy as np

def process_chunk_cds(chunk, cds_values=['DE_NOVO_START_IN_FRAME','DE_NOVO_START_OUT_FRAME','MISSENSE','NONSENSE','START_CODON_SNP'], cds_secondary_values=['INTRON','SILENT']):
    try:
        # Extract the mutations from the chunk, handling None and empty strings
        mutations = chunk.iloc[:, 3].apply(lambda x: x.split('|') if isinstance(x, str) and x.strip() else []).tolist()
        mutations = pd.DataFrame(mutations)
        
        # Debug information
        print(f"Number of rows in mutations DataFrame: {mutations.shape[0]}")
        
        # If no mutations were found or DataFrame is empty, return empty chunk
        if mutations.empty or mutations.shape[1] < 7:  # Need at least 7 fields for valid mutation
            print("Warning: No valid mutations found in chunk")
            return chunk.iloc[[]]  # Return empty DataFrame with same columns
        
        # Select only the needed columns (gene name, primary type, secondary type)
        mutations = mutations.iloc[:, [0, 5, 6]]

        # Create a mask for rows to keep
        keep_mask = []
        for index, row in mutations.iterrows():
            try:
                # Check if we have all required fields and they're not empty
                if row.shape[0] == 3 and not pd.isna(row.iloc[1]) and not pd.isna(row.iloc[2]):
                    condition1 = row.iloc[1] in cds_values  # check if the mutation is in a coding region
                    condition2 = row.iloc[2] not in cds_secondary_values  # check if the secondary SNP type is not an intron or silent mutation

                    print(f"row: {index}")
                    if not condition1:
                        print(f"SNP type: {row.iloc[1]} is not in the list of allowed SNP types")
                    if not condition2:
                        print(f"SNP type: {row.iloc[2]} is in the list of secondary SNP types")
                    
                    keep_mask.append(condition1 and condition2)
                else:
                    print(f"Warning: Invalid mutation format in row {index}")
                    keep_mask.append(False)
            except (IndexError, AttributeError) as e:
                print(f"Error processing row {index}: {str(e)}")
                keep_mask.append(False)

        # Apply the mask to the original chunk
        result = chunk[keep_mask]
        return result if not result.empty else chunk.iloc[[]]  # Return empty DataFrame if no matches found
        
    except Exception as e:
        print(f"Error processing chunk: {str(e)}")
        return chunk.iloc[[]]  # Return empty DataFrame in case of error

def process_vcf_file(filename, vcf_folder_path):
    # Construct the full file path
    vcf_path = os.path.join(vcf_folder_path, filename)

    # Construct the output file path with the modified filename
    output_filename = f"{filename.split('.')[0]}_cds_mutations.csv"
    output_csv_path = os.path.join(vcf_folder_path, output_filename)

    if os.path.exists(output_csv_path):
        return f"File {output_csv_path} already exists"
    
    # Read the CSV file with tab as delimiters
    vcf_df = pd.read_csv(vcf_path, sep='\t', engine='python')
    vcf_df = vcf_df.map(lambda x: x.replace('[', '').replace(']', '') if isinstance(x, str) else x) # Remove square brackets

    # Split the DataFrame into chunks for parallel processing
    num_processes = cpu_count()
    chunk_size = len(vcf_df) // num_processes
    if chunk_size > 0:
        chunks = [vcf_df.iloc[i:i + chunk_size] for i in range(0, len(vcf_df), chunk_size)]
    else:
        chunks = [vcf_df]
    
    # Create a partial function with fixed arguments
    process_chunk_func_cds = partial(process_chunk_cds)
    
    # Process chunks in parallel
    with Pool(processes=num_processes) as pool:
        process_chunk_func_cds = pool.map(process_chunk_func_cds, chunks)
    
    # Combine the processed chunks
    vcf_df = pd.concat(process_chunk_func_cds, ignore_index=True)

    # Save the matching rows to a new CSV file
    vcf_df.to_csv(output_csv_path, index=False)
    return f"Processed {filename}"

def process_chunk_mgl_gene(chunk, genes_df):
    try:
        # Extract the mutations from the chunk, handling None and empty strings
        mutations = chunk.iloc[:, 3].apply(lambda x: x.split('|') if isinstance(x, str) and x.strip() else []).tolist()
        mutations = pd.DataFrame(mutations)
        
        # Debug information
        print(f"Number of rows in mutations DataFrame: {mutations.shape[0]}")
        
        # If no mutations were found or DataFrame is empty, return empty chunk
        if mutations.empty or mutations.shape[1] == 0:
            print("Warning: No mutations found in chunk")
            return chunk.iloc[[]]  # Return empty DataFrame with same columns
        
        # Create a mask for rows to keep
        keep_mask = []
        for index, row in mutations.iterrows():
            try:
                # Check if the first column exists and contains a gene name
                if row.shape[0] > 0 and not pd.isna(row.iloc[0]) and row.iloc[0].strip():
                    condition3 = row.iloc[0] in genes_df.iloc[:, 0].values #check if the mutation is a gene of interest
                    if not condition3:
                        print(f"Gene: {row.iloc[0]} is not in the list of allowed genes")
                    keep_mask.append(condition3)
                else:
                    print(f"Warning: Invalid gene name in row {index}")
                    keep_mask.append(False)
            except (IndexError, AttributeError) as e:
                print(f"Error processing row {index}: {str(e)}")
                keep_mask.append(False)

        # Apply the mask to the original chunk
        result = chunk[keep_mask]
        return result if not result.empty else chunk.iloc[[]]  # Return empty DataFrame if no matches found
        
    except Exception as e:
        print(f"Error processing chunk: {str(e)}")
        return chunk.iloc[[]]  # Return empty DataFrame in case of error

def process_mgl_gene(filename, vcf_folder_path, genes_df):
    # Construct the full file path
    vcf_path = os.path.join(vcf_folder_path, filename)

    # Construct the output file path with the modified filename
    output_filename = f"{filename.split('.')[0]}_mgl_gene_mutations.csv"
    output_csv_path = os.path.join(vcf_folder_path, output_filename)

    if os.path.exists(output_csv_path):
        return f"File {output_csv_path} already exists"
    
    # Read the CSV file with comma as delimiter
    vcf_df = pd.read_csv(vcf_path)
    vcf_df = vcf_df.map(lambda x: x.replace('[', '').replace(']', '') if isinstance(x, str) else x) # Remove square brackets

    # Split the DataFrame into chunks for parallel processing
    num_processes = cpu_count()
    chunk_size = len(vcf_df) // num_processes
    if chunk_size > 0:
        chunks = [vcf_df.iloc[i:i + chunk_size] for i in range(0, len(vcf_df), chunk_size)]
    else:
        chunks = [vcf_df]
    
    # Create a partial function with fixed arguments
    process_chunk_func_mgl_gene = partial(process_chunk_mgl_gene, genes_df=genes_df)
    
    # Process chunks in parallel
    with Pool(processes=num_processes) as pool:
        process_chunk_func_mgl_gene = pool.map(process_chunk_func_mgl_gene, chunks)
    
    # Combine the processed chunks
    vcf_df = pd.concat(process_chunk_func_mgl_gene, ignore_index=True)
    
    # Save the matching rows to a new CSV file
    vcf_df.to_csv(output_csv_path, index=False)
    return f"Processed {filename}"

def main():
    # Define the folder path containing the VCF files
    vcf_folder_path = input("Please enter the VCF folder path: ")
    genes_csv_path = input("Please enter the genes CSV file path: ")

    # Read the genes CSV file
    genes_df = pd.read_csv(genes_csv_path)
    # Get list of all files in the directory
    files_cds = os.listdir(vcf_folder_path)
    
    # Filter files by extension
    files_cds = [f for f in files_cds if f.endswith(".ann.out")]
    
    # Create a partial function with fixed arguments
    process_func_cds = partial(process_vcf_file, 
                         vcf_folder_path=vcf_folder_path)
    
    # Process files sequentially
    for filename in files_cds:
        result = process_func_cds(filename)
        if result:
            print(result)

    files_mgl_gene = os.listdir(vcf_folder_path)
    files_mgl_gene = [f for f in files_mgl_gene if f.endswith("_cds_mutations.csv")]

    process_func_mgl_gene = partial(process_mgl_gene,
                         vcf_folder_path=vcf_folder_path,
                         genes_df=genes_df)
    
    for filename in files_mgl_gene:
        result = process_func_mgl_gene(filename)
        if result:
            print(result)

if __name__ == '__main__':
    main() 