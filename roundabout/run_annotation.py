import logging
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

def execute_amrfinder_parallel(fasta_paths: list[Path], outdir: str, db_path: str, threads: int) -> list[Path]:
    logging.info(f"Running AMRFinderPlus on {len(fasta_paths)} files...")
    amr_outdir = Path(outdir) / "amrfinder_results"
    amr_outdir.mkdir(parents=True, exist_ok=True)
    
    completed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_run_single_amrfinder, f, amr_outdir, db_path): f for f in fasta_paths}
        for future in concurrent.futures.as_completed(futures):
            try:
                completed.append(future.result())
            except Exception as e:
                logging.error(f"AMRFinderPlus failed for {futures[future].name}: {e}")
    return completed

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

def execute_plasmidfinder_parallel(fasta_paths: list[Path], outdir: str, db_path: str, threads: int) -> list[Path]:
    logging.info(f"Running PlasmidFinder on {len(fasta_paths)} files...")
    pf_outdir = Path(outdir) / "plasmidfinder_results"
    pf_outdir.mkdir(parents=True, exist_ok=True)
    
    completed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_run_single_plasmidfinder, f, pf_outdir, db_path): f for f in fasta_paths}
        for future in concurrent.futures.as_completed(futures):
            try:
                completed.append(future.result())
            except Exception as e:
                logging.error(f"PlasmidFinder failed for {futures[future].name}: {e}")
    return completed

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