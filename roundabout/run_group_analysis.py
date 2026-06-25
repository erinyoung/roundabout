import csv
import logging
import subprocess
from pathlib import Path

from .run_daisyblast import execute_daisyblast

def run_nucmer_pairwise(query: Path, ref: Path, out_prefix: Path, threads: int) -> float:
    """
    Runs nucmer between two sequences and calculates a custom distance score.
    Returns: 100.0 - % Identity of aligned regions (or 100.0 if no alignment).
    Compatible with both MUMmer 3 (single-threaded) and MUMmer 4 engines.
    """
    # Remove the '-t' option entirely to remain compatible with older system-installed nucmer engines
    nucmer_cmd = [
        "nucmer",
        "--maxmatch",
        "-p", str(out_prefix),
        str(ref),
        str(query)
    ]
    
    res = subprocess.run(nucmer_cmd, capture_output=True, text=True)
    if res.returncode != 0:
        logging.error(f"REAL NUCMER ERROR OUTPUT:\nSTDOUT: {res.stdout}\nSTDERR: {res.stderr}")
        raise subprocess.CalledProcessError(res.returncode, nucmer_cmd)
    
    delta_file = f"{out_prefix}.delta"
    coords_cmd = ["show-coords", "-r", "-c", "-l", delta_file]
    result = subprocess.run(coords_cmd, check=True, capture_output=True, text=True)
    
    for ext in [".delta", ".ntref", ".mgaps"]:
        f_path = Path(f"{out_prefix}{ext}")
        if f_path.exists():
            f_path.unlink()

    lines = result.stdout.strip().split('\n')
    total_aligned_len = 0
    weighted_identity_sum = 0.0
    
    data_started = False
    for line in lines:
        if "====" in line:
            data_started = True
            continue
        if not data_started or not line.strip():
            continue
            
        parts = line.split()
        if len(parts) >= 10:
            try:
                idy_val = None
                for token in parts:
                    if '.' in token and token.replace('.', '', 1).isdigit():
                        val = float(token)
                        if 0.0 <= val <= 100.0:
                            idy_val = val
                            break
                
                align_len = int(parts[4]) 
                
                if idy_val is not None and align_len > 0:
                    total_aligned_len += align_len
                    weighted_identity_sum += (idy_val * align_len)
            except (ValueError, IndexError):
                continue

    if total_aligned_len == 0:
        return 100.0  
        
    avg_identity = weighted_identity_sum / total_aligned_len
    return 100.0 - avg_identity


def generate_amr_highlights(ref_name: str, amr_dir: Path, output_csv: Path):
    """
    Parses the AMRFinderPlus TSV for the current reference sequence and generates
    a MinkeMap-compatible highlights CSV for any detected AMR genes.
    """
    tsv_file = amr_dir / f"{ref_name}_amrfinder.tsv"
    if not tsv_file.exists():
        return False
        
    highlights = []
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('Element symbol') or row.get('Gene symbol') or row.get('Sequence name')
            start = row.get('Start')
            stop = row.get('Stop')
            
            if gene and start and stop and row.get('Element type') == 'AMR':
                highlights.append({
                    'start': start,
                    'end': stop,
                    'color': '#ffcccc',  # Light transparent red wedge
                    'label': gene
                })
                
    if highlights:
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['start', 'end', 'color', 'label'])
            writer.writeheader()
            writer.writerows(highlights)
        return True
        
    return False


# Change the parameter from 'base_outdir' to 'group_outdir'
def analyze_cohort_group(group_name: str, sample_names: list[str], staging_dir: Path, group_outdir: Path, args):
    """
    Runs DaisyBlast, Nucmer All-vs-All, HeatCluster, and MinkeMap sequentially for a cohort.
    """
    logging.info(f"=== Processing Cohort Group: {group_name} ({len(sample_names)} sequences) ===")
    
    # We retrieve the base_outdir by going up one level from group_outdir
    # so we can still find the amrfinder_results directory for MinkeMap highlights
    base_outdir = group_outdir.parent
    amr_dir = base_outdir / "amrfinder_results"
    
    group_fastas = []
    for name in sample_names:
        fasta_path = staging_dir / f"{name}.fasta"
        if fasta_path.exists():
            group_fastas.append(fasta_path)

    if len(group_fastas) < 2:
        logging.warning(f"Not enough valid FASTAs to compare in group {group_name}. Skipping.")
        return

    threads = getattr(args, 'threads', 4)

    # -------------------------------------------------------------------------
    # A. DaisyBlast
    # -------------------------------------------------------------------------
    execute_daisyblast(group_fastas, str(group_outdir), args)

    # TODO replace nucmer with skani for all-vs-all distance matrix generation, as skani is faster and more efficient for large datasets.
    # TODO: put skani results into heatcluster    

    # # -------------------------------------------------------------------------
    # # B. Nucmer All-vs-All Distance Matrix
    # # -------------------------------------------------------------------------
    # logging.info(f"Running pairwise Nucmer alignments for {group_name}...")
    # melted_matrix_tsv = group_outdir / "nucmer_melted_matrix.txt"
    
    # with open(melted_matrix_tsv, 'w') as out_f:
    #     for i, query_f in enumerate(group_fastas):
    #         for j, ref_f in enumerate(group_fastas):
    #             q_name = query_f.stem
    #             r_name = ref_f.stem
                
    #             if i == j:
    #                 out_f.write(f"{q_name}\t{r_name}\t0.0\n")
    #             elif i > j:
    #                 continue 
    #             else:
    #                 prefix = group_outdir / f"temp_{q_name}_vs_{r_name}"
    #                 #try:
    #                     #distance = run_nucmer_pairwise(query_f, ref_f, prefix, threads)
    #                 #except Exception as e:
    #                     #logging.error(f"Nucmer alignment failed between {q_name} and {r_name}: {e}")
    #                     #distance = 100.0
                    
    #                 #out_f.write(f"{q_name}\t{r_name}\t{distance:.4f}\n")
    #                 #out_f.write(f"{r_name}\t{q_name}\t{distance:.4f}\n")

    # # -------------------------------------------------------------------------
    # # C. HeatCluster
    # # -------------------------------------------------------------------------
    # logging.info(f"Running HeatCluster for {group_name}...")
    # # Updated to just use group_name (e.g., group_1_heatcluster.png)
    # heat_out_png = group_outdir / f"{group_name}_heatcluster.png"
    
    # if melted_matrix_tsv.exists():
    #     heatcluster_cmd = [
    #         "heatcluster",
    #         "-i", str(melted_matrix_tsv),
    #         "--format", "melted",  
    #         "-o", str(heat_out_png),
    #         "--title", f"Sequence Distance (Nucmer): {group_name}"
    #     ]
    #     try:
    #         subprocess.run(heatcluster_cmd, check=True, capture_output=True, text=True)
    #         logging.info(f"Saved HeatCluster plot to {heat_out_png.name}")
    #     except subprocess.CalledProcessError as e:
    #         logging.error(f"HeatCluster failed for {group_name}. Error:\n{e.stderr}")

    # -------------------------------------------------------------------------
    # D. MinkeMap (Concentric Reference Permutations with AMR Highlights)
    # -------------------------------------------------------------------------
    logging.info(f"Running MinkeMap permutations for {group_name}...")
    minke_outdir = group_outdir / "minkemap_results"
    minke_outdir.mkdir(exist_ok=True)
    
    for ref_fasta in group_fastas:
        ref_name = ref_fasta.stem
        query_fastas = [str(f.absolute()) for f in group_fastas if f != ref_fasta]
        output_filename = f"minkemap_base_{ref_name}.png"
        
        minkemap_cmd = [
            "minkemap",
            "-r", str(ref_fasta.absolute()),
            "-i"
        ] + query_fastas + [
            "-o", output_filename,
            "--outdir", str(minke_outdir),
            "--title", f"MinkeMap Reference Base: {ref_name}",
            "--gc-skew"
        ]
        
        highlights_csv = minke_outdir / f"highlights_{ref_name}.csv"
        if generate_amr_highlights(ref_name, amr_dir, highlights_csv):
            logging.info(f"   Found AMR genes for {ref_name}. Adding visual wedges to plot.")
            minkemap_cmd += ["--highlights", str(highlights_csv)]
        
        try:
            logging.info(f"   -> Mapping cohort against reference base: {ref_name}...")
            subprocess.run(minkemap_cmd, check=True, capture_output=True, text=True)
            
            if highlights_csv.exists():
                highlights_csv.unlink()
                
        except subprocess.CalledProcessError as e:
            logging.error(f"MinkeMap failed for reference base {ref_name}. Error:\n{e.stderr}")