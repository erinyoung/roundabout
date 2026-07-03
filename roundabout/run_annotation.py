import logging
import csv
import subprocess
import concurrent.futures
from pathlib import Path
import shlex


# -----------------------------------------------------------------------------
# AMRFinderPlus
# -----------------------------------------------------------------------------
def _run_single_amrfinder(fasta_path: Path, outdir: Path, db_path: str, opts: dict) -> Path:
    sample_name = fasta_path.stem
    amr_out = outdir / f"{sample_name}_amrfinder.tsv"
    
    cmd = [
        "amrfinder",
        "--nucleotide", str(fasta_path),
        "--database", str(db_path),
        "--output", str(amr_out),
        "--plus"
    ]
    
    # Inject optional parameters
    if opts.get("organism"):
        cmd.extend(["--organism", opts["organism"]])
        
    if opts.get("ident_min") is not None:
        cmd.extend(["--ident_min", str(opts["ident_min"])])
        
    if opts.get("coverage_min") is not None:
        cmd.extend(["--coverage_min", str(opts["coverage_min"])])
        
    # Inject any extra raw options provided by the user
    if opts.get("options"):
        cmd.extend(shlex.split(opts["options"]))
    
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return amr_out

def execute_amrfinder_parallel(fasta_paths: list[Path], outdir: Path, db_path: str, threads: int, opts: dict) -> dict[str, list[str]]:
    logging.info(f"Running AMRFinderPlus on {len(fasta_paths)} files...")
    
    # Ensure outdir is a Path object, then create the subfolder
    amr_outdir = Path(outdir) / "amrfinder_results"
    amr_outdir.mkdir(parents=True, exist_ok=True)
    
    amr_profiles = {}
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Pass opts into the worker function
        futures = {executor.submit(_run_single_amrfinder, f, amr_outdir, db_path, opts): f for f in fasta_paths}
        
        for future in concurrent.futures.as_completed(futures):
            fasta_path = futures[future]
            isolate_id = fasta_path.stem 
            
            try:
                tsv_path = future.result()
                isolate_genes = []
                
                if tsv_path.exists():
                    with open(tsv_path, mode='r') as file:
                        reader = csv.DictReader(file, delimiter='\t')
                        for row in reader:
                            gene = row.get('Element symbol') or row.get('Gene symbol')
                            
                            if gene:
                                isolate_genes.append(gene)
                
                amr_profiles[isolate_id] = isolate_genes
                
            except Exception as e:
                logging.error(f"AMRFinderPlus failed for {fasta_path.name}: {e}")

    return amr_profiles
    
# -----------------------------------------------------------------------------
# PlasmidFinder
# -----------------------------------------------------------------------------
def _run_single_plasmidfinder(fasta_path: Path, outdir: Path, db_path: str, opts: dict) -> Path:
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
    
    # Inject optional parameters
    if opts.get("mincov") is not None:
        cmd.extend(["-l", str(opts["mincov"])])
        
    if opts.get("threshold") is not None:
        cmd.extend(["-t", str(opts["threshold"])])
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return sample_outdir
    except subprocess.CalledProcessError as e:
        logging.error(f"PlasmidFinder crashed on {fasta_path.name}. PF Error:\n{e.stderr}\n{e.stdout}")
        raise

def execute_plasmidfinder_parallel(fasta_paths: list[Path], outdir: Path, db_path: str, threads: int, opts: dict) -> dict[str, list[str]]:
    logging.info(f"Running PlasmidFinder on {len(fasta_paths)} files...")
    pf_outdir = Path(outdir) / "plasmidfinder_results"
    pf_outdir.mkdir(parents=True, exist_ok=True)
    
    plasmid_profiles = {}
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_run_single_plasmidfinder, f, pf_outdir, db_path, opts): f for f in fasta_paths}
        
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
def _run_single_bakta(fasta_path: Path, outdir: Path, db_path: str, opts: dict) -> Path:
    sample_name = fasta_path.stem
    sample_outdir = outdir / sample_name
    
    cmd = [
        "bakta",
        "--db", str(db_path),
        "--output", str(sample_outdir),
        "--prefix", sample_name,
        "--force",  
        "--threads", "1",
        # plasmids will mostly fail if sorf is not skipped, so this is a temporary workaround
        "--skip-sorf"
    ]
    
    # Map string and integer options to their flags
    arg_mapping = {
        "genus": "--genus",
        "species": "--species",
        "strain": "--strain",
        "plasmid": "--plasmid",
        "prodigal_tf": "--prodigal-tf",
        "translation_table": "--translation-table",
        "gram": "--gram",
        "locus": "--locus",
        "locus_tag": "--locus-tag",
        "locus_tag_increment": "--locus-tag-increment",
        "replicons": "--replicons",
        "regions": "--regions",
        "proteins": "--proteins",
        "hmms": "--hmms"
    }
    
    for key, flag in arg_mapping.items():
        if opts.get(key) is not None:
            # We don't want to pass empty strings if the user didn't specify a flag
            if str(opts[key]).strip(): 
                cmd.extend([flag, str(opts[key])])

    # Map boolean flags (actions)
    bool_mapping = {
        "complete": "--complete",
        "keep_contig_headers": "--keep-contig-headers",
        "compliant": "--compliant",
        "meta": "--meta",
        "partial": "--partial",
        "skip_trna": "--skip-trna",
        "skip_tmrna": "--skip-tmrna",
        "skip_rrna": "--skip-rrna",
        "skip_ncrna": "--skip-ncrna",
        "skip_ncrna_region": "--skip-ncrna-region",
        "skip_crispr": "--skip-crispr",
        "skip_cds": "--skip-cds",
        "skip_pseudo": "--skip-pseudo",
        "skip_gap": "--skip-gap",
        "skip_ori": "--skip-ori",
        "skip_filter": "--skip-filter",
        "skip_plot": "--skip-plot"
    }
    
    for key, flag in bool_mapping.items():
        if opts.get(key):
            cmd.append(flag)

    # Inject any extra raw options provided by the user
    if opts.get("options"):
        cmd.extend(shlex.split(opts["options"]))
        
    # Append the input fasta path as the positional argument at the very end
    cmd.append(str(fasta_path))
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return sample_outdir
    except subprocess.CalledProcessError as e:
        logging.error(f"Bakta crashed on {fasta_path.name}. Bakta Error:\n{e.stderr}\n{e.stdout}")
        raise

def execute_bakta_parallel(fasta_paths: list[Path], outdir: Path, db_path: str, threads: int, opts: dict) -> list[Path]:
    logging.info(f"Running Bakta on {len(fasta_paths)} files...")
    bakta_outdir = Path(outdir) / "bakta_results"
    bakta_outdir.mkdir(parents=True, exist_ok=True)
    
    completed = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_run_single_bakta, f, bakta_outdir, db_path, opts): f for f in fasta_paths}
        
        for future in concurrent.futures.as_completed(futures):
            try:
                completed.append(future.result())
            except Exception as e:
                logging.error(f"Bakta failed for {futures[future].name}: {e}")
                
    return completed