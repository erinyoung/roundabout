import os
import csv
import logging
import subprocess
from pathlib import Path

def get_symmetric_ani(sim_matrix: dict, a: str, b: str) -> float:
    """Safely retrieves the highest ANI between two sequences from the skani matrix."""
    ani_ab = sim_matrix.get(a, {}).get(b, {}).get('ani')
    ani_ba = sim_matrix.get(b, {}).get(a, {}).get('ani')
    
    val_ab = ani_ab if ani_ab is not None else 0.0
    val_ba = ani_ba if ani_ba is not None else 0.0
    return max(val_ab, val_ba)

def order_sequences_by_similarity(samples: list[str], sim_matrix: dict) -> list[str]:
    """
    Uses a greedy chain algorithm to optimally order sequences so highly 
    similar neighbors are adjacent (crucial for PyGenomeViz synteny ribbons).
    """
    if len(samples) <= 2:
        return samples
        
    unplaced = set(samples)
    
    best_ani = -1.0
    seed_pair = (samples[0], samples[1])
    
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            ani = get_symmetric_ani(sim_matrix, samples[i], samples[j])
            if ani > best_ani:
                best_ani = ani
                seed_pair = (samples[i], samples[j])
                
    ordered_chain = [seed_pair[0], seed_pair[1]]
    unplaced.remove(seed_pair[0])
    unplaced.remove(seed_pair[1])
    
    while unplaced:
        best_ani_overall = -1.0
        best_candidate = None
        best_position = None 
        
        left_end = ordered_chain[0]
        right_end = ordered_chain[-1]
        
        for candidate in unplaced:
            ani_left = get_symmetric_ani(sim_matrix, candidate, left_end)
            if ani_left > best_ani_overall:
                best_ani_overall = ani_left
                best_candidate = candidate
                best_position = 'left'
                
            ani_right = get_symmetric_ani(sim_matrix, candidate, right_end)
            if ani_right > best_ani_overall:
                best_ani_overall = ani_right
                best_candidate = candidate
                best_position = 'right'
                
        if best_position == 'left':
            ordered_chain.insert(0, best_candidate)
        else:
            ordered_chain.append(best_candidate)
            
        unplaced.remove(best_candidate)
        
    return ordered_chain

def generate_amr_highlights(ref_name: str, amr_dir: Path, output_csv: Path):
    """Parses AMRFinder TSV to generate MinkeMap AMR highlights."""
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
            el_type = str(row.get('Element type') or "").strip().lower()
            
            if gene and start and stop and el_type == 'amr':
                highlights.append({
                    'start': start,
                    'end': stop,
                    'color': '#ffcccc',
                    'label': gene
                })
                
    if highlights:
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['start', 'end', 'color', 'label'])
            writer.writeheader()
            writer.writerows(highlights)
        return True
        
    return False

def extract_fasta_from_refseq(name: str, refseq_fasta: Path, outdir: Path):
    """Extracts a specific sequence from the multi-FASTA RefSeq database on the fly."""
    if not refseq_fasta or not refseq_fasta.exists():
        return None
        
    out_path = outdir / f"{name}.fasta"
    if out_path.exists():
        return out_path 
        
    with open(refseq_fasta, 'r') as f:
        lines = []
        capturing = False
        for line in f:
            if line.startswith(">"):
                if capturing:
                    break  # We got our sequence, stop reading
                if name in line:
                    capturing = True
                    lines.append(line)
            elif capturing:
                lines.append(line)
                
    if capturing:
        outdir.mkdir(parents=True, exist_ok=True)
        with open(out_path, 'w') as f:
            f.writelines(lines)
        return out_path
        
    return None

def analyze_cohort_group(group_name: str, sample_names: list[str], staging_dir: Path, group_outdir: Path, sim_matrix: dict, args):
    """
    Runs downstream analysis for a specific cohort.
    Uses the pre-calculated skani matrix to bypass redundant alignments.
    """
    logging.info(f"=== Processing Cohort Group: {group_name} ({len(sample_names)} sequences) ===")
    
    base_outdir = group_outdir.parent
    amr_dir = base_outdir / "amrfinder_results"
    
    # -------------------------------------------------------------------------
    # A. Write Cohort Distance Matrix & Run HeatCluster
    # -------------------------------------------------------------------------
    logging.info(f"Generating HeatCluster plot for {group_name}...")
    melted_matrix_tsv = group_outdir / f"{group_name}_distance_matrix.txt"
    heat_out_png = group_outdir / f"{group_name}_heatcluster.png"
    
    with open(melted_matrix_tsv, 'w') as out_f:
        for i, q_name in enumerate(sample_names):
            out_f.write(f"{q_name}\t{q_name}\t0.0000\n")
            for j in range(i + 1, len(sample_names)):
                r_name = sample_names[j]
                ani = get_symmetric_ani(sim_matrix, q_name, r_name)
                dist = 100.0 - ani if ani > 0 else 100.0
                
                out_f.write(f"{q_name}\t{r_name}\t{dist:.4f}\n")
                out_f.write(f"{r_name}\t{q_name}\t{dist:.4f}\n")

    heatcluster_cmd = [
        "heatcluster",
        "-i", str(melted_matrix_tsv),
        "--format", "melted",  
        "-o", str(heat_out_png),
        "--title", f"Sequence Distance: {group_name}"
    ]
    try:
        subprocess.run(heatcluster_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"HeatCluster failed for {group_name}. Error:\n{e.stderr}")

    # -------------------------------------------------------------------------
    # B. MinkeMap 
    # -------------------------------------------------------------------------
    logging.info(f"Running MinkeMap permutations for {group_name}...")
    minke_outdir = group_outdir / "minkemap_results"
    minke_outdir.mkdir(exist_ok=True)
    
    # Identify the RefSeq Database Path
    refseq_db_str = getattr(args, 'refseq_plasmid_dl_db', None)
    refseq_fasta = Path(str(refseq_db_str)) / "refseq_plasmids_dl.fasta" if refseq_db_str else None
    refseq_extract_dir = group_outdir / "refseq_fastas"
    
    # Gather all FASTAs (Inputs + On-the-fly extracted RefSeq)
    group_fastas = []
    for name in sample_names:
        staging_path = staging_dir / f"{name}.fasta"
        if staging_path.exists():
            group_fastas.append(staging_path)
        else:
            extracted_path = extract_fasta_from_refseq(name, refseq_fasta, refseq_extract_dir)
            if extracted_path:
                group_fastas.append(extracted_path)
    
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
            minkemap_cmd += ["--highlights", str(highlights_csv)]
        
        try:
            subprocess.run(minkemap_cmd, check=True, capture_output=True, text=True)
            if highlights_csv.exists(): highlights_csv.unlink()
        except subprocess.CalledProcessError as e:
            logging.error(f"MinkeMap failed for {ref_name}. Error:\n{e.stderr}")

    # -------------------------------------------------------------------------
    # C. PyGenomeViz (Optimally Ordered via MMseqs2, MUMmer, and BLAST)
    # -------------------------------------------------------------------------
    logging.info(f"Running PyGenomeViz synteny plots for {group_name}...")
    pgv_outdir = group_outdir / "pygenomeviz_results"
    pgv_outdir.mkdir(exist_ok=True)
    
    # Optimally order the sequences
    ordered_names = order_sequences_by_similarity(sample_names, sim_matrix)
    
    # Gather files (Prefer Bakta GenBank for feature arrows, fallback to raw FASTA)
    bakta_dir = base_outdir / "bakta_results"
    ordered_pgv_paths = []
    
    for name in ordered_names:
        gbff_path = bakta_dir / name / f"{name}.gbff"
        gbk_path = bakta_dir / name / f"{name}.gbk"
        fasta_path = staging_dir / f"{name}.fasta"
        refseq_path = refseq_extract_dir / f"{name}.fasta"
        
        if gbff_path.exists():
            ordered_pgv_paths.append(str(gbff_path))
        elif gbk_path.exists():
            ordered_pgv_paths.append(str(gbk_path))
        elif fasta_path.exists():
            ordered_pgv_paths.append(str(fasta_path))
        elif refseq_path.exists():
            # If it's a RefSeq sequence, it will use the extracted FASTA we made for MinkeMap!
            ordered_pgv_paths.append(str(refseq_path))
        else:
            if " " not in name:
                logging.warning(f"   -> Missing both GenBank and FASTA for {name}. Dropping from PyGenomeViz.")

    if len(ordered_pgv_paths) >= 2:
        # 1. MMSeqs2 (Removed the -t flag)
        mmseqs_out = pgv_outdir / f"{group_name}_mmseqs.png"
        pgv_mmseqs_cmd = ["pgv-mmseqs", *ordered_pgv_paths, "-o", str(mmseqs_out)]
        try:
            subprocess.run(pgv_mmseqs_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"pgv-mmseqs failed for {group_name}:\n{e.stderr}")

        # 2. MUMmer (Removed the -t flag)
        mummer_out = pgv_outdir / f"{group_name}_mummer.png"
        pgv_mummer_cmd = ["pgv-mummer", *ordered_pgv_paths, "-o", str(mummer_out)]
        try:
            subprocess.run(pgv_mummer_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"pgv-mummer failed for {group_name}:\n{e.stderr}")

        # 3. BLAST (Removed the -t flag)
        blast_out = pgv_outdir / f"{group_name}_blast.png"
        pgv_blast_cmd = ["pgv-blast", *ordered_pgv_paths, "-o", str(blast_out)]
        try:
            subprocess.run(pgv_blast_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"pgv-blast failed for {group_name}:\n{e.stderr}")