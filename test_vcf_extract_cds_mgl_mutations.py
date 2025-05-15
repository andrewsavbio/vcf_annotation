import pytest
import pandas as pd
import os
import tempfile
import shutil
from vcf_extract_cds_mgl_mutations import (
    process_chunk_cds,
    process_chunk_mgl_gene,
    process_vcf_file,
    process_mgl_gene
)

@pytest.fixture
def test_dir():
    # Create a temporary directory for test files
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Clean up temporary directory after tests
    shutil.rmtree(temp_dir)

@pytest.fixture
def sample_vcf_data():
    return pd.DataFrame({
        'CHROM': ['chr1'] * 18,
        'POS': range(1000, 19000, 1000),
        'TYPE': ['SNP'] * 18,
        'FUNCOTATION': [
            # 7 cds mutations
            # 5 mgl mutations

            # BRCA1 mutations
            'BRCA1|hg38|chr1|1000|1000|DE_NOVO_START_IN_FRAME|MISSENSE',
            'BRCA1|hg38|chr1|2000|2000|MISSENSE|NONSENSE',
            'BRCA1|hg38|chr1|3000|3000|INTRON|MISSENSE',
            'BRCA1|hg38|chr1|4000|4000|MISSENSE|',
            'BRCA1|hg38|chr1|5000|5000|RNA|MISSENSE',
            # TP53 mutations
            'TP53|hg38|chr1|6000|6000|MISSENSE|MISSENSE',
            'TP53|hg38|chr1|7000|7000|NONSENSE|INTRON',
            'TP53|hg38|chr1|8000|8000|DE_NOVO_START_OUT_FRAME|MISSENSE',
            'TP53|hg38|chr1|9000|9000|INTRON|',
            'TP53|hg38|chr1|10000|10000|RNA|',
            # SHROOM3 mutations (should be filtered out because not in genes list)
            'SHROOM3|hg38|chr1|11000|11000|START_CODON_SNP|SILENT',
            'SHROOM3|hg38|chr1|12000|12000|MISSENSE|INTRON',
            'SHROOM3|hg38|chr1|13000|13000|NONSENSE|',
            'SHROOM3|hg38|chr1|14000|14000|INTRON|MISSENSE',
            'SHROOM3|hg38|chr1|15000|15000|RNA|MISSENSE',
            # Additional test cases
            'BRCA1|hg38|chr1|16000|16000|MISSENSE|INTRON',      # Should be filtered out (INTRON in secondary)
            'TP53|hg38|chr1|17000|17000|INTRON|MISSENSE',       # Should be filtered out (INTRON in primary)
            'SHROOM3|hg38|chr1|18000|18000|MISSENSE|MISSENSE'   # Should be filtered out (not in genes list)
        ]
    })

@pytest.fixture
def sample_genes_data():
    return pd.DataFrame({
        'Gene': ['BRCA1', 'TP53']  # SHROOM3 intentionally excluded
    })

@pytest.fixture
def setup_test_files(test_dir, sample_vcf_data, sample_genes_data):
    # Save sample data to temporary files
    vcf_path = os.path.join(test_dir, 'test.ann.out')
    genes_path = os.path.join(test_dir, 'test_genes.csv')
    sample_vcf_data.to_csv(vcf_path, sep='\t', index=False)
    sample_genes_data.to_csv(genes_path, index=False)
    return vcf_path, genes_path

def test_process_chunk_cds(sample_vcf_data):
    # Test with valid CDS mutations
    result = process_chunk_cds(sample_vcf_data)
    assert len(result) == 7  # Should keep rows with valid CDS mutations
    
    # Verify specific mutations are kept
    kept_mutations = result['FUNCOTATION'].tolist()
    # BRCA1 mutations
    assert 'BRCA1|hg38|chr1|1000|1000|DE_NOVO_START_IN_FRAME|MISSENSE' in kept_mutations
    assert 'BRCA1|hg38|chr1|2000|2000|MISSENSE|NONSENSE' in kept_mutations
    assert 'BRCA1|hg38|chr1|4000|4000|MISSENSE|' in kept_mutations
    
    # TP53 mutations
    assert 'TP53|hg38|chr1|6000|6000|MISSENSE|MISSENSE' in kept_mutations
    assert 'TP53|hg38|chr1|8000|8000|DE_NOVO_START_OUT_FRAME|MISSENSE' in kept_mutations
    
    # Verify SHROOM3 mutations
    assert 'SHROOM3|hg38|chr1|11000|11000|START_CODON_SNP|SILENT' not in kept_mutations
    assert 'SHROOM3|hg38|chr1|12000|12000|MISSENSE|INTRON' not in kept_mutations
    assert 'SHROOM3|hg38|chr1|13000|13000|NONSENSE|' in kept_mutations
    assert 'SHROOM3|hg38|chr1|18000|18000|MISSENSE|MISSENSE' in kept_mutations
    
    # Verify filtered mutations
    assert 'BRCA1|hg38|chr1|3000|3000|INTRON|MISSENSE' not in kept_mutations
    assert 'BRCA1|hg38|chr1|5000|5000|RNA|MISSENSE' not in kept_mutations
    assert 'TP53|hg38|chr1|9000|9000|INTRON|' not in kept_mutations
    assert 'TP53|hg38|chr1|10000|10000|RNA|' not in kept_mutations
    assert 'SHROOM3|hg38|chr1|14000|14000|INTRON|MISSENSE' not in kept_mutations
    assert 'SHROOM3|hg38|chr1|15000|15000|RNA|MISSENSE' not in kept_mutations
    
    # Test with empty DataFrame
    empty_df = pd.DataFrame() 
    result = process_chunk_cds(empty_df)
    assert len(result) == 0

def test_process_chunk_mgl_gene(sample_vcf_data, sample_genes_data):
    # Test with valid gene mutations
    result = process_chunk_mgl_gene(sample_vcf_data, sample_genes_data)
    assert len(result) == 12  # Should keep only BRCA1 and TP53 mutations
    
    # Verify all mutations are from valid genes
    kept_genes = [row.split('|')[0] for row in result['FUNCOTATION'].tolist()]
    assert all(gene in ['BRCA1', 'TP53'] for gene in kept_genes)
    assert not any(gene == 'SHROOM3' for gene in kept_genes)
    
    # Test with empty DataFrame
    empty_df = pd.DataFrame()
    result = process_chunk_mgl_gene(empty_df, sample_genes_data)
    assert len(result) == 0

def test_process_vcf_file(test_dir, setup_test_files):
    # Test processing VCF file
    result = process_vcf_file('test.ann.out', test_dir)
    output_path = os.path.join(test_dir, 'test_cds_mutations.csv')
    assert os.path.exists(output_path)
    
    # Verify the output file contains the expected data
    output_df = pd.read_csv(output_path)
    assert len(output_df) == 7

def test_process_mgl_gene(test_dir, setup_test_files, sample_genes_data):
    # First create the CDS mutations file
    vcf_path, genes_path = setup_test_files
    process_vcf_file('test.ann.out', test_dir)
    
    # Verify the CDS mutations file exists and has correct content
    cds_output_path = os.path.join(test_dir, 'test_cds_mutations.csv')
    assert os.path.exists(cds_output_path)
    cds_df = pd.read_csv(cds_output_path)
    assert len(cds_df) == 7
    
    # Test processing MGL gene file
    result = process_mgl_gene('test_cds_mutations.csv', test_dir, sample_genes_data)
    mgl_output_path = os.path.join(test_dir, 'test_cds_mutations_mgl_gene_mutations.csv')
    
    # Verify the MGL gene mutations file exists and has correct content
    assert os.path.exists(mgl_output_path)
    mgl_df = pd.read_csv(mgl_output_path)
    assert len(mgl_df) == 5
    
    # Verify the content of the MGL gene mutations file
    kept_genes = [row.split('|')[0] for row in mgl_df['FUNCOTATION'].tolist()]
    assert all(gene in ['BRCA1', 'TP53'] for gene in kept_genes)
    assert not any(gene == 'SHROOM3' for gene in kept_genes)

def test_process_chunk_cds_edge_cases():
    # Test with None values
    df_with_none = pd.DataFrame({
        'CHROM': ['chr1', 'chr1'],
        'POS': [1000, 2000],
        'TYPE': ['SNP', 'SNP'],
        'FUNCOTATION': ['', 'BRCA1|hg38|chr1|1000|1000|INTRON|SILENT']  # Empty string instead of None
    })
    result = process_chunk_cds(df_with_none)
    assert len(result) == 0

    # Test with malformed FUNCOTATION column
    df_malformed = pd.DataFrame({
        'CHROM': ['chr1'],
        'POS': [1000],
        'TYPE': ['SNP'],
        'FUNCOTATION': ['BRCA1|hg38|chr1|1000|1000|MISSENSE']  # Missing seventh field
    })
    result = process_chunk_cds(df_malformed)
    assert len(result) == 0

    # Test with empty string in FUNCOTATION
    df_empty = pd.DataFrame({
        'CHROM': ['chr1'],
        'POS': [1000],
        'TYPE': ['SNP'],
        'FUNCOTATION': ['']  # Empty string
    })
    result = process_chunk_cds(df_empty)
    assert len(result) == 0

def test_process_chunk_cds_mutation_combinations():
    # Test various mutation combinations across different genes
    test_cases = pd.DataFrame({
        'CHROM': ['chr1'] * 9,
        'POS': range(1000, 10000, 1000),
        'TYPE': ['SNP'] * 9,
        'FUNCOTATION': [
            'BRCA1|hg38|chr1|1000|1000|MISSENSE|INTRON',      # Should be filtered out (INTRON in secondary)
            'TP53|hg38|chr1|2000|2000|INTRON|MISSENSE',       # Should be filtered out (INTRON in primary)
            'SHROOM3|hg38|chr1|3000|3000|MISSENSE|MISSENSE',  # Should be kept
            'BRCA1|hg38|chr1|4000|4000|DE_NOVO_START_IN_FRAME|MISSENSE',  # Should be kept
            'TP53|hg38|chr1|5000|5000|MISSENSE|',             # Should be kept (empty second type)
            'SHROOM3|hg38|chr1|6000|6000|NONSENSE|',          # Should be kept (empty second type)
            'BRCA1|hg38|chr1|7000|7000|INTRON|',              # Should be filtered out (INTRON in primary)
            'TP53|hg38|chr1|8000|8000|RNA|',                  # Should be filtered out (RNA in primary)
            'SHROOM3|hg38|chr1|9000|9000|START_CODON_SNP|INTRON'  # Should be filtered out (not in genes list)
        ]
    })
    
    result = process_chunk_cds(test_cases)
    assert len(result) == 4  # Only 4 combinations should be kept
    
    # Verify specific mutations are kept
    kept_mutations = result['FUNCOTATION'].tolist()
    assert 'SHROOM3|hg38|chr1|3000|3000|MISSENSE|MISSENSE' in kept_mutations
    assert 'BRCA1|hg38|chr1|4000|4000|DE_NOVO_START_IN_FRAME|MISSENSE' in kept_mutations
    assert 'TP53|hg38|chr1|5000|5000|MISSENSE|' in kept_mutations
    assert 'SHROOM3|hg38|chr1|6000|6000|NONSENSE|' in kept_mutations
