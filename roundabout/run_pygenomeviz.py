import logging
from pathlib import Path
import pandas as pd
import numpy as np

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

from pygenomeviz import GenomeViz
from pygenomeviz.align import Blast, MUMmer, ProgressiveMauve, MMseqs, AlignCoord

def get_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """
    Parses a FASTA file and returns a dictionary.
    FIXED: Overrides the internal fasta header ID with the filename stem
    to prevent SegmentNotFoundErrors across multi-aligners.
    """
    length = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                length += len(line.strip())
                
    # Use the filename stem as the segment name, making it completely bulletproof
    return {fasta_path.stem: length}

def sort_fastas_by_ani(fasta_paths: list[Path], local_matrix_df: pd.DataFrame) -> list[Path]:
    """
    Sorts a list of FASTA paths mathematically using hierarchical clustering 
    derived from the local ANI matrix. Highly similar sequences will be adjacent.
    """
    if len(fasta_paths) < 3:
        return fasta_paths

    def clean_id(val):
        return Path(str(val)).stem.split('.')[0].strip()

    valid_paths = []
    valid_ids = []
    
    for path in fasta_paths:
        cid = clean_id(path)
        if cid in local_matrix_df.index:
            valid_paths.append(path)
            valid_ids.append(cid)

    if len(valid_paths) < 3:
        return fasta_paths

    # Extract sub-matrix and convert similarity to distance
    sub_matrix = local_matrix_df.loc[valid_ids, valid_ids].copy().fillna(80.0)
    distance_matrix = 100.0 - sub_matrix.values

    # Force strict symmetry for SciPy
    distance_matrix = (distance_matrix + distance_matrix.T) / 2.0
    np.fill_diagonal(distance_matrix, 0.0)

    # Perform clustering and extract optimal order
    condensed_dist = squareform(distance_matrix)
    Z = linkage(condensed_dist, method='average') 
    leaf_order = leaves_list(Z)

    return [valid_paths[i] for i in leaf_order]

def run_pygenomeviz_blast(
    fasta_paths: list[Path], 
    local_matrix_df: pd.DataFrame,
    out_path: Path, 
    min_length: int = 500, 
    min_identity: float = 80.0
):
    """
    Takes a list of FASTA paths, extracts their lengths to build tracks,
    runs nucleotide BLAST, and plots the synteny links.
    """
    if len(fasta_paths) < 2:
        logging.warning("Need at least 2 FASTA files to run BLAST synteny. Skipping.")
        return

    # --- ADDED: Sort the sequences using the ANI matrix before plotting ---
    sorted_fastas = sort_fastas_by_ani(fasta_paths, local_matrix_df)
    # ----------------------------------------------------------------------

    logging.info(f"Running pygenomeviz BLAST on {len(sorted_fastas)} sequences...")

    # 1. Initialize GenomeViz Canvas
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar()

    # 2. Build Tracks from FASTA lengths in the newly sorted order
    for fasta in sorted_fastas:
        seq_lengths = get_fasta_lengths(fasta)
        # Use the filename stem as the track label
        gv.add_feature_track(fasta.stem, seq_lengths)

    # 3. Run Nucleotide BLAST directly on the sorted FASTAs
    try:
        align_coords = Blast(
            [str(p) for p in sorted_fastas], 
            seqtype="nucleotide"
        ).run()
    except Exception as e:
        logging.error(f"pygenomeviz BLAST failed: {e}")
        return

    # Filter out tiny or low-quality hits to keep the plot clean
    align_coords = AlignCoord.filter(
        align_coords, 
        length_thr=min_length, 
        identity_thr=min_identity
    )

    # 4. Draw the Alignment Links
    if align_coords:
        # Dynamically set the colorbar floor based on the worst passing hit
        min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
        color, inverted_color = "grey", "red"
        
        for ac in align_coords:
            gv.add_link(
                ac.query_link, 
                ac.ref_link, 
                color=color, 
                inverted_color=inverted_color, 
                v=ac.identity, 
                vmin=min_ident
            )
            
        gv.set_colorbar([color, inverted_color], vmin=min_ident)
    else:
        logging.warning("No BLAST alignments passed the thresholds for these FASTAs.")

    # 5. Save out the figure
    gv.savefig(out_path, dpi=300)

def run_pygenomeviz_mummer(
    fasta_paths: list[Path], 
    local_matrix_df: pd.DataFrame, 
    out_path: Path, 
    min_length: int = 500, 
    min_identity: float = 80.0
):
    """
    Takes a list of FASTA paths, extracts their lengths to build tracks,
    runs MUMmer (nucmer), and plots the synteny links.
    """
    if len(fasta_paths) < 2:
        logging.warning("Need at least 2 FASTA files to run MUMmer synteny. Skipping.")
        return

    # Sort the sequences using the ANI matrix before plotting
    sorted_fastas = sort_fastas_by_ani(fasta_paths, local_matrix_df)

    logging.info(f"Running pygenomeviz MUMmer on {len(sorted_fastas)} sequences...")

    # Initialize GenomeViz Canvas
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar()

    # Build Tracks from FASTA lengths in the sorted order
    for fasta in sorted_fastas:
        seq_lengths = get_fasta_lengths(fasta)
        gv.add_feature_track(fasta.stem, seq_lengths)

    # Run Nucleotide MUMmer (nucmer) directly on the sorted FASTAs
    try:
        align_coords = MUMmer(
            [str(p) for p in sorted_fastas], 
            seqtype="nucleotide"
        ).run()
    except Exception as e:
        logging.error(f"pygenomeviz MUMmer failed: {e}")
        return

    # Filter out tiny or low-quality hits
    align_coords = AlignCoord.filter(
        align_coords, 
        length_thr=min_length, 
        identity_thr=min_identity
    )

    # Draw the Alignment Links
    if align_coords:
        min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
        color, inverted_color = "grey", "red"
        
        for ac in align_coords:
            gv.add_link(
                ac.query_link, 
                ac.ref_link, 
                color=color, 
                inverted_color=inverted_color, 
                v=ac.identity, 
                vmin=min_ident
            )
            
        gv.set_colorbar([color, inverted_color], vmin=min_ident)
    else:
        logging.warning("No MUMmer alignments passed the thresholds for these FASTAs.")

    # Save out the figure
    gv.savefig(out_path, dpi=300)


def run_pygenomeviz_pmauve(
    fasta_paths: list[Path], 
    local_matrix_df: pd.DataFrame, 
    out_path: Path, 
    min_length: int = 500, 
    min_identity: float = 80.0
):
    """
    Takes a list of FASTA paths, extracts their lengths to build tracks,
    runs progressiveMauve, and plots the synteny links.
    """
    if len(fasta_paths) < 2:
        logging.warning("Need at least 2 FASTA files to run progressiveMauve synteny. Skipping.")
        return

    # Sort the sequences using the ANI matrix before plotting
    sorted_fastas = sort_fastas_by_ani(fasta_paths, local_matrix_df)

    logging.info(f"Running pygenomeviz progressiveMauve on {len(sorted_fastas)} sequences...")

    # Initialize GenomeViz Canvas
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar()

    # Build Tracks from FASTA lengths in the sorted order
    for fasta in sorted_fastas:
        seq_lengths = get_fasta_lengths(fasta)
        # Using the clean structural positional binding
        gv.add_feature_track(fasta.stem, seq_lengths)

    # Run progressiveMauve directly on the sorted FASTAs
    try:
        align_coords = ProgressiveMauve(
            [str(p) for p in sorted_fastas]
        ).run()
    except Exception as e:
        logging.error(f"pygenomeviz progressiveMauve failed: {e}")
        return

    # Filter out tiny hits
    align_coords = AlignCoord.filter(
        align_coords, 
        length_thr=min_length
    )

    # Draw the Alignment Links
    if align_coords:
        color, inverted_color = "grey", "red"
        
        for ac in align_coords:
            gv.add_link(
                ac.query_link, 
                ac.ref_link, 
                color=color, 
                inverted_color=inverted_color
            )
    else:
        logging.warning("No progressiveMauve alignments passed the length threshold.")

    # Save out the figure
    gv.savefig(out_path, dpi=300)
def run_pygenomeviz_mmseqs(
    fasta_paths: list[Path], 
    local_matrix_df: pd.DataFrame, 
    out_path: Path, 
    min_length: int = 500, 
    min_identity: float = 80.0
):
    """
    Takes a list of FASTA paths, extracts their lengths to build tracks,
    runs MMseqs (translated protein alignment), and plots the synteny links.
    """
    if len(fasta_paths) < 2:
        logging.warning("Need at least 2 FASTA files to run MMseqs synteny. Skipping.")
        return

    # Sort the sequences using the ANI matrix before plotting
    sorted_fastas = sort_fastas_by_ani(fasta_paths, local_matrix_df)

    logging.info(f"Running pygenomeviz MMseqs on {len(sorted_fastas)} sequences...")

    # Initialize GenomeViz Canvas
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar()

    # Build Tracks from FASTA lengths in the sorted order
    for fasta in sorted_fastas:
        seq_lengths = get_fasta_lengths(fasta)
        gv.add_feature_track(fasta.stem, seq_lengths)

    # Run MMseqs directly on the sorted FASTAs (seqtype="protein" runs 6-frame translation)
    try:
        align_coords = MMseqs(
            [str(p) for p in sorted_fastas]
        ).run()
    except Exception as e:
        logging.error(f"pygenomeviz MMseqs failed: {e}")
        return

    # Filter out tiny or low-quality hits
    align_coords = AlignCoord.filter(
        align_coords, 
        length_thr=min_length, 
        identity_thr=min_identity
    )

    # Draw the Alignment Links
    if align_coords:
        min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
        color, inverted_color = "grey", "red"
        
        for ac in align_coords:
            gv.add_link(
                ac.query_link, 
                ac.ref_link, 
                color=color, 
                inverted_color=inverted_color, 
                v=ac.identity, 
                vmin=min_ident
            )
            
        gv.set_colorbar([color, inverted_color], vmin=min_ident)
    else:
        logging.warning("No MMseqs alignments passed the thresholds for these FASTAs.")

    # Save out the figure
    gv.savefig(out_path, dpi=300)