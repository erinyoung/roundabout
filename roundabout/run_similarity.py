import subprocess
import logging
import csv
from pathlib import Path
from collections import defaultdict

def execute_skani(input_fastas: list[Path], refseq_fasta: Path, outdir: Path, threads: int) -> Path:
    """
    Executes 'skani dist' to calculate all-vs-all ANI and alignment fractions.
    Runs inputs against themselves AND the global RefSeq database in a single pass.
    """
    skani_outdir = Path(outdir) / "skani_results"
    skani_outdir.mkdir(exist_ok=True)
    
    query_list_path = skani_outdir / "skani_queries.txt"
    ref_list_path = skani_outdir / "skani_refs.txt"
    skani_tsv = skani_outdir / "skani_matrix.tsv"
    
    # Queries are always just the input sequences
    with open(query_list_path, 'w') as qf:
        for f in input_fastas:
            qf.write(f"{f.absolute()}\n")
            
    # References include inputs (for local grouping) + RefSeq (for global grouping)
    with open(ref_list_path, 'w') as rf:
        for f in input_fastas:
            rf.write(f"{f.absolute()}\n")
        if refseq_fasta and refseq_fasta.exists():
            rf.write(f"{refseq_fasta.absolute()}\n")
            
    cmd = [
        "skani", "dist",
        "--ql", str(query_list_path),
        "--rl", str(ref_list_path),
        "--ri",  # Treats multi-fasta records (RefSeq) as independent reference genomes
        "-o", str(skani_tsv),
        "-t", str(threads)
    ]
    
    logging.info("Running skani identity profiling (Local + Global)...")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return skani_tsv
    except subprocess.CalledProcessError as e:
        logging.error(f"skani execution failed:\n{e.stderr}")
        raise


def parse_skani_matrix(skani_tsv: Path) -> dict:
    """
    Reads the single skani TSV output into a nested matrix dictionary.
    """
    matrix = {}
    with open(skani_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            q_name = Path(row['Query_file']).stem
            
            # If the hit came from the RefSeq multi-fasta, we want the specific sequence header.
            # If it came from our staged inputs, we want the file stem to match our pipeline logic.
            if 'refseq_plasmids_dl' in row['Ref_file']:
                r_name = row['Ref_name'].strip()
            else:
                r_name = Path(row['Ref_file']).stem 
            
            if q_name not in matrix:
                matrix[q_name] = {}
                
            matrix[q_name][r_name] = {
                'ani': float(row['ANI']),
                'q_cov': float(row['Align_fraction_query']),
                'r_cov': float(row['Align_fraction_ref'])
            }
    return matrix


def build_sentinel_groups(similarity_matrix: dict, input_sequences: list[str], refseq_sequences: list[str], min_identity: float, min_coverage: float):
    """
    Evaluates sequence similarity to build ego-networks for each input sentinel.
    Categorizes relationships into Strict (symmetric) and Contained (asymmetric)
    across both Local (inputs) and Global (inputs + refseq) scopes.
    """
    local_strict = {}
    local_contained = {}
    global_strict = {}
    global_contained = {}
    
    for sentinel in input_sequences:
        ls_members = [sentinel]
        lc_members = [sentinel]
        
        # 1. Local Scope (Sentinel vs Input Sequences)
        for target in input_sequences:
            if target == sentinel:
                continue
                
            metrics = similarity_matrix.get(sentinel, {}).get(target, {})
            ani = metrics.get('ani', 0.0)
            q_cov = metrics.get('q_cov', 0.0)
            r_cov = metrics.get('r_cov', 0.0)
            
            if ani >= min_identity:
                if q_cov >= min_coverage or r_cov >= min_coverage:
                    lc_members.append(target)
                if q_cov >= min_coverage and r_cov >= min_coverage:
                    ls_members.append(target)
                    
        local_strict[f"{sentinel}"] = ls_members
        local_contained[f"{sentinel}"] = lc_members
        
        # 2. Global Scope (Local Scope + RefSeq Sequences)
        gs_members = list(ls_members)
        gc_members = list(lc_members)
        
        for ref_target in refseq_sequences:
            metrics = similarity_matrix.get(sentinel, {}).get(ref_target, {})
            ani = metrics.get('ani', 0.0)
            q_cov = metrics.get('q_cov', 0.0)
            r_cov = metrics.get('r_cov', 0.0)
            
            if ani >= min_identity:
                if q_cov >= min_coverage or r_cov >= min_coverage:
                    gc_members.append(ref_target)
                if q_cov >= min_coverage and r_cov >= min_coverage:
                    gs_members.append(ref_target)
                    
        global_strict[f"{sentinel}"] = gs_members
        global_contained[f"{sentinel}"] = gc_members
        
    return local_strict, local_contained, global_strict, global_contained