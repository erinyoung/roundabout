import subprocess
import logging
from pathlib import Path
import pandas as pd

def execute_skani(input_fastas: list[Path], refseq_fasta: Path | None, outdir: Path, threads: int, opts: dict) -> Path:
    """
    Executes 'skani dist' to calculate all-vs-all ANI and alignment fractions.
    Runs inputs against themselves AND the global RefSeq database in a single pass.
    """
    skani_outdir = Path(outdir) / "skani_results"
    skani_outdir.mkdir(exist_ok=True, parents=True)
    
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
        # Runs cleanly even if refseq_fasta is missing
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
    
    # ---------------------------------------------------------
    # Inject Skani CLI Options
    # ---------------------------------------------------------
    
    # Map values/integers
    val_mapping = {
        "min_af": "--min-af",
        "both_min_af": "--both-min-af",
        "n": "-n",
        "c": "-c",
        "m": "-m",
        "s": "-s"
    }
    for key, flag in val_mapping.items():
        if opts.get(key) is not None:
            cmd.extend([flag, str(opts[key])])

    # Map boolean flags
    bool_mapping = {
        "ci": "--ci",
        "detailed": "--detailed",
        "short_header": "--short-header",
        "fast": "--fast",
        "medium": "--medium",
        "slow": "--slow",
        "small_genomes": "--small-genomes",
        "faster_small": "--faster-small",
        "median": "--median",
        "no_learned_ani": "--no-learned-ani",
        "no_marker_index": "--no-marker-index",
        "robust": "--robust"
    }
    for key, flag in bool_mapping.items():
        if opts.get(key):
            cmd.append(flag)
            
    # ---------------------------------------------------------

    logging.info("Running skani identity profiling...")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return skani_tsv 
    except subprocess.CalledProcessError as e:
        logging.error(f"skani execution failed:\n{e.stderr}")
        raise

def parse_skani_results(skani_tsv: Path) -> pd.DataFrame:
    """
    Reads the single skani TSV output into a pandas dataframe.
    """
    try:
        # skani default columns often include: Query file, Reference file, ANI, etc.
        # Ensure your skani version matches these column names or rename them after loading.
        df = pd.read_csv(skani_tsv, sep='\t')
        return df
    except Exception as e:
        logging.error(f"Failed to parse skani TSV: {e}")
        raise

def create_local_ani_matrix(skani_df: pd.DataFrame, input_fastas: list[Path]) -> pd.DataFrame:
    """
    Filters out RefSeq references and creates a square ANI matrix for input sequences only.
    """
    input_filenames = [f.name for f in input_fastas]
    
    # Check what columns skani actually generated to avoid KeyErrors
    cols = skani_df.columns.tolist()
    
    # Dynamically map the columns based on standard skani outputs or your previous manual parsing
    q_col = 'Query_file' if 'Query_file' in cols else ('query' if 'query' in cols else None)
    r_col = 'Ref_file' if 'Ref_file' in cols else ('target' if 'target' in cols else None)
    ani_col = 'ANI' if 'ANI' in cols else ('ani' if 'ani' in cols else None)
    
    if not all([q_col, r_col, ani_col]):
        logging.error(f"Skani columns mismatch. Available columns are: {cols}")
        raise KeyError(f"Missing required skani columns. Found: {cols}")

    # Filter for local vs local
    local_df = skani_df[
        skani_df[q_col].apply(lambda x: Path(x).name in input_filenames) & 
        skani_df[r_col].apply(lambda x: Path(x).name in input_filenames)
    ].copy()

    # Clean up paths to just show filenames for the heatmap
    local_df['query_name'] = local_df[q_col].apply(lambda x: Path(x).stem)
    local_df['target_name'] = local_df[r_col].apply(lambda x: Path(x).stem)
    
    # Pivot to create a square matrix
    matrix_df = local_df.pivot(index='query_name', columns='target_name', values=ani_col)
    
    # skani only reports hits above a certain threshold (usually ~80% ANI).
    # Fill missing comparisons with 0 for the heatmap.
    matrix_df = matrix_df.fillna(0)
    
    return matrix_df

def extract_global_hits(skani_df: pd.DataFrame, input_fastas: list[Path], refseq_meta_csv: Path | None = None) -> pd.DataFrame:
    """
    Filters the raw skani dataframe to extract only the comparisons 
    between local input sequences and the global RefSeq database, 
    and merges biological metadata if available.
    """
    input_filenames = [f.name for f in input_fastas]
    cols = skani_df.columns.tolist()
    
    q_col = 'Query_file' if 'Query_file' in cols else ('query' if 'query' in cols else cols[0])
    r_col = 'Ref_file' if 'Ref_file' in cols else ('target' if 'target' in cols else cols[1])
    ani_col = 'ANI' if 'ANI' in cols else ('ani' if 'ani' in cols else 'ANI')
    
    # Identify where the actual FASTA header lives (Skani uses Ref_name when --ri is passed)
    ref_name_col = 'Ref_name' if 'Ref_name' in cols else r_col
    
    # Filter: Query FILE is in our inputs, Reference FILE is NOT in our inputs (meaning it's RefSeq)
    global_df = skani_df[
        skani_df[q_col].apply(lambda x: Path(x).name in input_filenames) & 
        ~skani_df[r_col].apply(lambda x: Path(x).name in input_filenames)
    ].copy()
    
    if not global_df.empty:
        # Clean up the query paths to just show the names
        global_df['Query_Name'] = global_df[q_col].apply(lambda x: Path(x).stem)
        
        # Safely extract the RefSeq accession from the Ref_name column!
        def extract_accession(val):
            seq_id = str(val).split()[0]  # Grab just the accession, drop the description
            base = Path(seq_id).name
            # Strip standard fasta extensions if it somehow treated it as a file
            for ext in ['.fasta', '.fa', '.fna', '.gz']:
                if base.endswith(ext):
                    base = base[:-len(ext)]
            return base

        global_df['RefSeq_Hit'] = global_df[ref_name_col].apply(extract_accession)
        
        # Merge the metadata if the CSV is provided
        if refseq_meta_csv and refseq_meta_csv.exists():
            try:
                meta_df = pd.read_csv(refseq_meta_csv)
                
                # Ensure types match for the merge
                meta_df['accession'] = meta_df['accession'].astype(str)
                global_df['RefSeq_Hit'] = global_df['RefSeq_Hit'].astype(str)
                
                # Perform a Left Join on the accession number
                global_df = global_df.merge(
                    meta_df, 
                    how='left', 
                    left_on='RefSeq_Hit', 
                    right_on='accession'
                )
                
                # Drop the duplicate 'accession' column created by the merge
                if 'accession' in global_df.columns:
                    global_df.drop(columns=['accession'], inplace=True)
                    
            except Exception as e:
                logging.error(f"Failed to merge RefSeq metadata: {e}")

        # Sort by Query sequence, then by highest ANI match
        global_df = global_df.sort_values(by=['Query_Name', ani_col], ascending=[True, False])
        
        # Reorder columns to put the cleaned names first
        cols_to_keep = ['Query_Name', 'RefSeq_Hit'] + [c for c in global_df.columns if c not in ['Query_Name', 'RefSeq_Hit', q_col, r_col, ref_name_col]]
        global_df = global_df[cols_to_keep]
        
    return global_df


