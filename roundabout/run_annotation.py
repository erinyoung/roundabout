import logging
import csv
import subprocess
import concurrent.futures
from pathlib import Path

# -----------------------------------------------------------------------------
# AMRFinderPlus
# -----------------------------------------------------------------------------
def _run_single_amrfinder(fasta_path: Path, outdir: Path, db_path: str) -> Path:
    sample_name = fasta_path.stem
    amr_out = outdir / f"{sample_name}_amrfinder.tsv"
    
    cmd = [
        "amrfinder",
        "--nucleotide", str(fasta_path),
        "--database", str(db_path),
        "--output", str(amr_out),
        "--plus"
    ]
    
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return amr_out

def execute_amrfinder_parallel(fasta_paths: list[Path], outdir: str, db_path: str, threads: int) -> dict[str, list[str]]:
    logging.info(f"Running AMRFinderPlus on {len(fasta_paths)} files...")
    amr_outdir = Path(outdir) / "amrfinder_results"
    amr_outdir.mkdir(parents=True, exist_ok=True)
    
    # Replaced 'completed' list with a dictionary
    amr_profiles = {}
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_run_single_amrfinder, f, amr_outdir, db_path): f for f in fasta_paths}
        
        for future in concurrent.futures.as_completed(futures):
            fasta_path = futures[future]
            
            # Use the fasta file name (without extension) as the Isolate ID
            isolate_id = fasta_path.stem 
            
            try:
                # Get the output TSV path from your run function
                tsv_path = future.result()
                
                isolate_genes = []
                
                # Open and parse the TSV immediately
                if tsv_path.exists():
                    with open(tsv_path, mode='r') as file:
                        reader = csv.DictReader(file, delimiter='\t')
                        for row in reader:
                            # Safely extract either 'Element symbol' or 'Gene symbol'
                            gene = row.get('Element symbol') or row.get('Gene symbol')
                            if gene:
                                isolate_genes.append(gene)
                
                # Assign the gene list to the isolate in the dictionary
                amr_profiles[isolate_id] = isolate_genes
                
            except Exception as e:
                logging.error(f"AMRFinderPlus failed for {fasta_path.name}: {e}")

    return amr_profiles
    
# -----------------------------------------------------------------------------
# PlasmidFinder
# -----------------------------------------------------------------------------
def _run_single_plasmidfinder(fasta_path: Path, outdir: Path, db_path: str) -> Path:
    sample_name = fasta_path.stem
    sample_outdir = outdir / sample_name
    sample_outdir.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "plasmidfinder.py",
        "-i", str(fasta_path),
        "-o", str(sample_outdir),
        "-p", str(db_path),
        "-x" 
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return sample_outdir
    except subprocess.CalledProcessError as e:
        logging.error(f"PlasmidFinder crashed on {fasta_path.name}. PF Error:\n{e.stderr}\n{e.stdout}")
        raise

def execute_plasmidfinder_parallel(fasta_paths: list[Path], outdir: str, db_path: str, threads: int) -> dict[str, list[str]]:
    logging.info(f"Running PlasmidFinder on {len(fasta_paths)} files...")
    pf_outdir = Path(outdir) / "plasmidfinder_results"
    pf_outdir.mkdir(parents=True, exist_ok=True)
    
    plasmid_profiles = {}
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_run_single_plasmidfinder, f, pf_outdir, db_path): f for f in fasta_paths}
        
        for future in concurrent.futures.as_completed(futures):
            fasta_path = futures[future]
            isolate_id = fasta_path.stem 
            
            try:
                # future.result() is returning the directory, so let's name it appropriately
                result_dir = Path(future.result()) 
                
                # Point directly to the TSV file inside that directory
                tsv_path = result_dir / "results_tab.tsv"
                
                isolate_plasmids = []
                
                # Open and parse the TSV immediately
                if tsv_path.exists():
                    with open(tsv_path, mode='r') as file:
                        reader = csv.DictReader(file, delimiter='\t')
                        for row in reader:
                            # Safely extract the plasmid column
                            plasmid = row.get('Plasmid')
                            if plasmid:
                                isolate_plasmids.append(plasmid)
                
                # Use set() to remove duplicate hits, then save to dict
                plasmid_profiles[isolate_id] = list(set(isolate_plasmids))
                
            except Exception as e:
                logging.error(f"PlasmidFinder failed for {fasta_path.name}: {e}")
                plasmid_profiles[isolate_id] = []

    return plasmid_profiles

# -----------------------------------------------------------------------------
# Bakta
# -----------------------------------------------------------------------------
def _run_single_bakta(fasta_path: Path, outdir: Path, db_path: str) -> Path:
    sample_name = fasta_path.stem
    sample_outdir = outdir / sample_name
    
    cmd = [
        "bakta",
        "--db", str(db_path),
        "--output", str(sample_outdir),
        "--prefix", sample_name,
        "--force",  
        "--threads", "1",
        "--skip-sorf",  # <-- Bypasses the memory-heavy Diamond step!
        str(fasta_path)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return sample_outdir
    except subprocess.CalledProcessError as e:
        logging.error(f"Bakta crashed on {fasta_path.name}. Bakta Error:\n{e.stderr}\n{e.stdout}")
        raise

def execute_bakta_parallel(fasta_paths: list[Path], outdir: str, db_path: str, threads: int) -> list[Path]:
    logging.info(f"Running Bakta on {len(fasta_paths)} files...")
    bakta_outdir = Path(outdir) / "bakta_results"
    bakta_outdir.mkdir(parents=True, exist_ok=True)
    
    completed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_run_single_bakta, f, bakta_outdir, db_path): f for f in fasta_paths}
        for future in concurrent.futures.as_completed(futures):
            try:
                completed.append(future.result())
            except Exception as e:
                logging.error(f"Bakta failed for {futures[future].name}: {e}")
    return completed