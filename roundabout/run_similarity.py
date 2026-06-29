import subprocess
import logging
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def execute_skani(input_fastas: list[Path], refseq_fasta: Path | None, outdir: Path, threads: int) -> Path:
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
        # Handle TODO: Runs cleanly even if refseq_fasta is missing
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
    
    logging.info("Running skani identity profiling...")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return skani_tsv # Return just the path, let the pipeline orchestrate the next steps
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

def extract_global_hits(skani_df: pd.DataFrame, input_fastas: list[Path]) -> pd.DataFrame:
    """
    Filters the raw skani dataframe to extract only the comparisons 
    between local input sequences and the global RefSeq database.
    """
    input_filenames = [f.name for f in input_fastas]
    cols = skani_df.columns.tolist()
    
    q_col = 'Query_file' if 'Query_file' in cols else ('query' if 'query' in cols else cols[0])
    r_col = 'Ref_file' if 'Ref_file' in cols else ('target' if 'target' in cols else cols[1])
    ani_col = 'ANI' if 'ANI' in cols else ('ani' if 'ani' in cols else 'ANI')
    
    # Filter: Query is in our inputs, Reference is NOT in our inputs (meaning it's RefSeq)
    global_df = skani_df[
        skani_df[q_col].apply(lambda x: Path(x).name in input_filenames) & 
        ~skani_df[r_col].apply(lambda x: Path(x).name in input_filenames)
    ].copy()
    
    if not global_df.empty:
        # Clean up the paths to just show the names
        global_df['Query_Name'] = global_df[q_col].apply(lambda x: Path(x).stem)
        global_df['RefSeq_Hit'] = global_df[r_col].apply(lambda x: Path(x).stem)
        
        # Sort by Query sequence, then by highest ANI match
        global_df = global_df.sort_values(by=['Query_Name', ani_col], ascending=[True, False])
        
        # Reorder columns to put the cleaned names first
        cols_to_keep = ['Query_Name', 'RefSeq_Hit'] + [c for c in global_df.columns if c not in ['Query_Name', 'RefSeq_Hit', q_col, r_col]]
        global_df = global_df[cols_to_keep]
        
    return global_df

import warnings
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
import logging

def visualize_global_matches_scatter(global_df: pd.DataFrame, out_path: Path):
    """
    Visualizes global hits using a strip plot.
    Highly efficient and readable for massive reference databases (e.g., 10k+ hits).
    """
    if global_df.empty:
        logging.getLogger(__name__).warning("Empty global dataframe. Skipping scatter visualization.")
        return
        
    # Silence Matplotlib's internal logger for categorical units
    logging.getLogger('matplotlib.category').setLevel(logging.WARNING)
    logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

    # Dynamically find the ANI column name
    ani_col = 'ANI' if 'ANI' in global_df.columns else 'ani'

    # Clean data
    global_df = global_df.copy()
    global_df[ani_col] = pd.to_numeric(global_df[ani_col], errors='coerce')
    global_df['Query_Name'] = global_df['Query_Name'].astype(str)

    num_queries = global_df['Query_Name'].nunique()
    fig_height = max(6, num_queries * 0.5)

    plt.figure(figsize=(10, fig_height))

    # Plot all RefSeq hits as semi-transparent dots
    sns.stripplot(
        data=global_df, 
        x=ani_col, 
        y='Query_Name', 
        jitter=0.2, 
        alpha=0.4, 
        size=4,
        color='teal'
    )

    # Add a visual threshold line for high-confidence matches (95%)
    plt.axvline(x=95.0, color='black', linestyle='--', alpha=0.5, label='95% ANI Threshold')

    plt.title('Global Database (RefSeq) Matches per Input', pad=15)
    plt.xlabel('ANI (%)')
    plt.ylabel('Local Sequences')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def visualize_skani_matrix(matrix_df: pd.DataFrame, out_path: Path):
    """
    Visualizes the square skani matrix. Dynamically scales figure size 
    and label formatting to prevent unreadable, cluttered outputs.
    """
    num_seqs = len(matrix_df)
    if num_seqs == 0:
        logging.warning("Empty matrix provided for visualization. Skipping heatmap.")
        return

    # Dynamically scale the figure size based on the number of sequences (minimum 10x8)
    fig_width = max(10, num_seqs * 0.4)
    fig_height = max(8, num_seqs * 0.35)
    plt.figure(figsize=(fig_width, fig_height))
    
    # If there are too many sequences, the text numbers in the boxes will overlap. 
    # Turn them off if n > 20.
    show_annot = True if num_seqs <= 20 else False
    
    heatmap = sns.heatmap(
        matrix_df, 
        cmap='viridis', 
        annot=show_annot, 
        fmt=".2f" if show_annot else "", 
        cbar_kws={'label': 'ANI (%)'},
        linewidths=0.5 if num_seqs <= 50 else 0 # Add gridlines only for smaller matrices
    )
    
    # Rotate axis labels to fit long sequence names
    plt.xticks(rotation=45, ha='right', fontsize=max(6, 12 - (num_seqs // 10)))
    plt.yticks(rotation=0, fontsize=max(6, 12 - (num_seqs // 10)))
    
    plt.title('Local Isolates All-vs-All ANI', pad=20)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()



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