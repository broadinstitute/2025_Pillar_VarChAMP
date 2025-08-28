"""
Genomics 2 Proteins Portal
==============================================
Author: Runxi Shen
Created: 2025-08-25
"""

import re
import os
import requests
import glob
import gzip
import time
import json
import warnings
import shutil
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
import pandas as pd
import numpy as np
import polars as pl


def get_uniprot_swissprot_id(protein_name: str) -> str:
    """
    Query UniProt’s REST search API to find the reviewed (Swiss‐Prot) accession
    for a given human gene/protein name. Returns None if not found.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    # Build a query that:
    #  - matches the gene name exactly (using “gene:”)
    #  - restricts to human (organism_id:9606)
    #  - restricts to reviewed (Swiss‐Prot) entries
    query = f"gene:{protein_name} AND organism_id:9606 AND reviewed:true"
    params = {
        "query": query,
        "fields": "accession",
        "format": "json",
        "size": 1,      # only need the top hit
    }
    try:
        resp = requests.get(url, params=params, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        results = data.get("results", [])
        if not results:
            return None
        return results[0]["primaryAccession"]
    except Exception:
        # If the request fails (e.g. no internet), return None
        return None


def extract_g2p_single_endpoint(gene_name: str, uniprot_id: str, endpoint: str, output_dir: str) -> bool:
    """
    Extract data from a single G2P portal API endpoint.
    
    Args:
        gene_name: Gene symbol (e.g., 'AGXT')
        uniprot_id: UniProt accession ID (e.g., 'P21549')
        endpoint: API endpoint path
        output_dir: Directory to save the output
    
    Returns:
        bool: True if successful, False otherwise
    """
    # Construct the G2P API URL
    url = f"https://g2p.broadinstitute.org/api/gene/{gene_name}/protein/{uniprot_id}/{endpoint}"
    
    headers = {
        'accept': 'application/csv'
    }
    
    try:
        print(f"Querying G2P portal endpoint '{endpoint}' for {gene_name} (UniProt: {uniprot_id})...")
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        
        # Check if we got valid data
        if response.content:
            # Extract filename from content-disposition header
            content_disposition = response.headers.get('content-disposition', '')
            if 'filename=' in content_disposition and '"' in content_disposition:
                try:
                    filename = content_disposition.split('filename="')[1].split('"')[0]
                except IndexError:
                    filename = f"G2P_{gene_name}_{uniprot_id}_{endpoint.replace('-', '_')}.tsv"
            else:
                filename = f"G2P_{gene_name}_{uniprot_id}_{endpoint.replace('-', '_')}.tsv"
            
            output_file = Path(output_dir) / filename
            
            # Handle gzipped content
            if response.headers.get('content-encoding') == 'gzip':
                try:
                    # Content is gzipped, decompress it
                    decompressed_content = gzip.decompress(response.content).decode('utf-8')
                    with open(output_file, 'w') as f:
                        f.write(decompressed_content)
                except gzip.BadGzipFile:
                    # Content claims to be gzipped but isn't - save as is
                    with open(output_file, 'wb') as f:
                        f.write(response.content)
            else:
                # Content is not gzipped
                with open(output_file, 'wb') as f:
                    f.write(response.content)
            
            print(f"Successfully saved G2P data to: {output_file}")
            
            # Validate that we got meaningful data
            try:
                # Try to read as TSV first
                df = pd.read_csv(output_file, sep='\t')
                print(f"Retrieved {len(df)} rows with {len(df.columns)} columns")
                if len(df.columns) > 0:
                    print(f"Columns: {list(df.columns[:5])}")  # Show first 5 columns
                return True
            except Exception as e:
                # If TSV parsing fails, check if the file exists and has content
                if output_file.exists() and output_file.stat().st_size > 0:
                    print(f"File saved successfully but could not parse as TSV: {e}")
                    print(f"File size: {output_file.stat().st_size} bytes")
                    return True
                else:
                    print(f"Warning: Could not parse TSV data: {e}")
                    return False
        else:
            print(f"No data returned for {gene_name} (UniProt: {uniprot_id}) from endpoint '{endpoint}'")
            return False
            
    except requests.exceptions.RequestException as e:
        print(f"Error querying G2P portal endpoint '{endpoint}' for {gene_name}: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error for {gene_name} endpoint '{endpoint}': {e}")
        return False


def extract_g2p_protein_features(gene_name: str, uniprot_id: str, output_dir: str = None) -> dict:
    """
    Extract protein features from both G2P portal API endpoints for a given gene and UniProt ID.
    
    Args:
        gene_name: Gene symbol (e.g., 'AGXT')
        uniprot_id: UniProt accession ID (e.g., 'P21549')
        output_dir: Directory to save the output files
    
    Returns:
        dict: Results for each endpoint
    """
    if not output_dir:
        output_dir = "."
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Define the two API endpoints
    endpoints = [
        "gene-transcript-protein-isoform-structure-map",
        "protein-features"
    ]
    
    results = {}
    
    for endpoint in endpoints:
        print(f"\n--- Processing endpoint: {endpoint} ---")
        success = extract_g2p_single_endpoint(gene_name, uniprot_id, endpoint, output_dir)
        results[endpoint] = success
        
        # Be respectful to the API
        time.sleep(0.5)
    
    return results


def process_gene_list(gene_list: list, output_dir: str = None, delay: float = 1.0) -> dict:
    """
    Process a list of genes to extract G2P protein features from both endpoints.
    
    Args:
        gene_list: List of gene symbols
        output_dir: Directory to save outputs
        delay: Delay between genes in seconds to be respectful
    
    Returns:
        dict: Results summary with success/failure counts
    """
    results = {
        "genes_processed": 0, 
        "genes_failed": 0, 
        "no_uniprot": 0, 
        "endpoint_results": {},
        "failed_genes": []
    }
    
    for gene in gene_list:
        print(f"\n{'='*60}")
        print(f"Processing gene: {gene}")
        print(f"{'='*60}")
        
        # Get UniProt ID
        uniprot_id = get_uniprot_swissprot_id(gene)
        if not uniprot_id:
            print(f"Could not find UniProt ID for {gene}")
            results["no_uniprot"] += 1
            results["failed_genes"].append(f"{gene} (no UniProt ID)")
            continue
        
        print(f"Found UniProt ID: {uniprot_id}")
        
        # Extract features from both endpoints
        endpoint_results = extract_g2p_protein_features(gene, uniprot_id, output_dir)
        results["endpoint_results"][gene] = endpoint_results
        
        # Check if at least one endpoint succeeded
        if any(endpoint_results.values()):
            results["genes_processed"] += 1
            print(f"✓ Successfully processed {gene}")
        else:
            results["genes_failed"] += 1
            results["failed_genes"].append(f"{gene} ({uniprot_id}) - all endpoints failed")
            print(f"✗ Failed to process {gene}")
        
        # Be respectful to the API
        if delay > 0:
            time.sleep(delay)
    
    return results


def parse_aa_change(aa_change: str) -> Tuple[str, int, str]:
    """Parse amino acid change string like 'Leu359Pro' into (ref_aa, position, alt_aa)."""
    if pd.isna(aa_change) or not aa_change:
        return None, None, None
    
    match = re.match(r'^([A-Za-z]{3})(\d+)([A-Za-z]{3})$', aa_change)
    if match:
        ref_aa, pos, alt_aa = match.groups()
        return ref_aa, int(pos), alt_aa
    
    return None, None, None


def three_letter_to_one_letter(three_letter: str) -> str:
    """Convert three-letter amino acid code to one-letter code."""
    aa_map = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
    return aa_map.get(three_letter, three_letter)


def validate_variants_against_protein_features(variants_df: pd.DataFrame, protein_features_file: str) -> Dict:
    """Validate that all variants' reference amino acids match the protein features file."""
    if not Path(protein_features_file).exists():
        return {"valid": False, "error": f"File not found: {protein_features_file}"}
    
    try:
        prot_df = pd.read_csv(protein_features_file, sep='\t')
        pos_to_aa = dict(zip(prot_df['residueId'], prot_df['AA']))
        
        mismatches = []
        valid_count = 0
        
        for _, variant_row in variants_df.iterrows():
            aa_change = variant_row['aa_change']
            ref_aa_3letter, position, alt_aa_3letter = parse_aa_change(aa_change)
            
            if ref_aa_3letter is None:
                continue
            
            ref_aa_1letter = three_letter_to_one_letter(ref_aa_3letter)
            
            if position in pos_to_aa:
                protein_aa = pos_to_aa[position]
                if protein_aa == ref_aa_1letter:
                    valid_count += 1
                else:
                    mismatches.append({
                        'variant': aa_change,
                        'position': position,
                        'expected_aa': ref_aa_1letter,
                        'protein_aa': protein_aa
                    })
            else:
                mismatches.append({
                    'variant': aa_change,
                    'position': position,
                    'expected_aa': ref_aa_1letter,
                    'protein_aa': 'POSITION_NOT_FOUND'
                })
        
        return {
            "valid": len(mismatches) == 0,
            "total_variants": len(variants_df),
            "valid_matches": valid_count,
            "mismatches": mismatches,
            "protein_length": len(pos_to_aa)
        }
        
    except Exception as e:
        return {"valid": False, "error": f"Error reading protein features: {e}"}


def extract_protein_features_for_variants(variants_df: pd.DataFrame, protein_features_file: str) -> pd.DataFrame:
    """Extract all protein features for each variant position and append to variants dataframe."""
    enhanced_df = variants_df.copy()
    
    if not Path(protein_features_file).exists():
        print(f"Warning: Protein features file not found: {protein_features_file}")
        return enhanced_df
    
    try:
        # Load protein features
        prot_df = pd.read_csv(protein_features_file, sep='\t')
        
        # Create mapping from position to all features
        protein_features_dict = {}
        for _, row in prot_df.iterrows():
            pos = int(row['residueId'])
            # Convert all feature values to strings to handle mixed types
            features = {col: str(val) if pd.notna(val) else '' for col, val in row.items() if col != 'residueId'}
            protein_features_dict[pos] = features
        
        # Get all possible feature columns
        all_feature_columns = set()
        for features in protein_features_dict.values():
            all_feature_columns.update(features.keys())
        all_feature_columns = sorted(all_feature_columns)
        
        # Initialize new columns with empty strings
        for col in all_feature_columns:
            enhanced_df[f'protein_{col}'] = ''
        
        # Add protein features for each variant
        for idx, row in enhanced_df.iterrows():
            aa_change = row['aa_change']
            ref_aa_3letter, position, alt_aa_3letter = parse_aa_change(aa_change)
            
            if position and position in protein_features_dict:
                features = protein_features_dict[position]
                for col in all_feature_columns:
                    enhanced_df.loc[idx, f'protein_{col}'] = features.get(col, '')
        
        print(f"  ✓ Added {len(all_feature_columns)} protein feature columns")
        return enhanced_df
        
    except Exception as e:
        print(f"  ✗ Error extracting protein features: {e}")
        return enhanced_df
        

def main():
    """
    Main function to extract G2P protein features for example genes.
    """
    # Example gene list - you can expand this
    genes = ["F9"]  # Add more genes as needed
    
    # Process genes
    results = process_gene_list(genes, delay=1.0)
    
    # Print detailed summary
    print("\n" + "="*70)
    print("G2P PROTEIN FEATURES EXTRACTION SUMMARY")
    print("="*70)
    print(f"Genes successfully processed: {results['genes_processed']}")
    print(f"Genes failed: {results['genes_failed']}")
    print(f"Genes without UniProt ID: {results['no_uniprot']}")
    
    # Show endpoint-specific results
    if results['endpoint_results']:
        print(f"\nEndpoint-specific results:")
        for gene, endpoint_results in results['endpoint_results'].items():
            print(f"\n  {gene}:")
            for endpoint, success in endpoint_results.items():
                status = "✓ Success" if success else "✗ Failed"
                print(f"    - {endpoint}: {status}")
    
    if results['failed_genes']:
        print(f"\nFailed genes:")
        for gene in results['failed_genes']:
            print(f"  - {gene}")
    
    print(f"\nOutput files saved to: ./")
    print("="*70)


if __name__ == "__main__":
    main()