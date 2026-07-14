import itertools
import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


from Bio import SeqIO


def get_cds_coords_from_gbk(gbk_path: str) -> list:
    """Parses a GenBank file and returns a list of (start, end) tuples for all CDS features."""
    coords = []

    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                # Biopython uses 0-based indexing, so we cast to standard 1-based coordinates
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                coords.append((start, end))

    return coords


def visualize_nucmer_heatmap(
    matrix_df: pd.DataFrame, out_path: Path, metric_type: str = "identity"
):
    """
    Visualizes pairwise nucmer metrics using hierarchical clustering.
    Adapts scale and labels based on the metric_type: 'identity', 'coverage', or 'snps'.
    """
    num_seqs = len(matrix_df)
    if num_seqs == 0:
        logging.warning("Empty matrix provided for visualization. Skipping heatmap.")
        return

    show_annot = False
    fig_width = max(8, num_seqs * 0.4)
    fig_height = max(8, num_seqs * 0.4)

    # Configure heatmap dynamics based on what we are plotting
    if metric_type == "identity":
        cbar_label = "Average Identity (%)"
        vmin, vmax = 80, 100  # Tailored for high sequence similarity
        cmap = "viridis"
    elif metric_type == "coverage":
        cbar_label = "Alignment Coverage (%)"
        vmin, vmax = 0, 100  # Plasmids can vary wildly in coverage length
        cmap = "magma"
    elif metric_type == "snps":
        cbar_label = "SNP Count"
        vmin, vmax = None, None  # Automatically scale based on the data minimum/maximum
        cmap = "rocket_r"  # Reverted so more SNPs = hotter color
    else:
        cbar_label = metric_type
        vmin, vmax = None, None
        cmap = "viridis"

    try:
        cg = sns.clustermap(
            matrix_df,
            cmap=cmap,
            annot=show_annot,
            figsize=(fig_width, fig_height),
            cbar_kws={"label": cbar_label},
            linewidths=0.5 if num_seqs <= 30 else 0,
            vmin=vmin,
            vmax=vmax,
            dendrogram_ratio=(0.01, 0.01),
            cbar_pos=(1.02, 0.15, 0.03, 0.7),
        )
    except (FloatingPointError, ValueError) as e:
        logging.warning(f"Clustering failed ({e}). Plotting unclustered fallback.")
        cg = sns.clustermap(
            matrix_df,
            cmap=cmap,
            annot=show_annot,
            figsize=(fig_width, fig_height),
            cbar_kws={"label": cbar_label},
            linewidths=0.5 if num_seqs <= 30 else 0,
            vmin=vmin,
            vmax=vmax,
            dendrogram_ratio=(0.01, 0.01),
            cbar_pos=(1.02, 0.15, 0.03, 0.7),
            row_cluster=False,
            col_cluster=False,
        )

    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)
    cg.ax_heatmap.set_xlabel("")
    cg.ax_heatmap.set_ylabel("")

    plt.setp(
        cg.ax_heatmap.get_xticklabels(),
        rotation=45,
        ha="right",
        fontsize=max(8, 12 - (num_seqs // 10)),
    )
    plt.setp(
        cg.ax_heatmap.get_yticklabels(),
        rotation=0,
        fontsize=max(8, 12 - (num_seqs // 10)),
    )

    cg.figure.suptitle(
        f"Cohort Pairwise Comparison: {metric_type.title()}", y=1.02, fontsize=14
    )

    cg.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def parse_delta_with_snps(
    delta_path: str, cds_coords: Optional[List[tuple]] = None
) -> dict:
    """
    Runs show-snps on a delta file and returns both the total number of SNPs
    and the number of SNPs falling within the provided CDS coordinates.
    """
    cmd = ["show-snps", "-T", "-H", delta_path]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout_clean = result.stdout.strip()
        if not stdout_clean:
            return {"total_snps": 0, "coding_snps": 0}

        lines = stdout_clean.split("\n")
        total_snps = len(lines)
        coding_snps = 0

        # Filter for coding SNPs if coordinates are provided
        if cds_coords:
            for line in lines:
                parts = line.split("\t")
                if not parts:
                    continue

                snp_pos = int(parts[0])
                # Check if SNP falls in any CDS region
                if any(start <= snp_pos <= end for start, end in cds_coords):
                    coding_snps += 1

        return {"total_snps": total_snps, "coding_snps": coding_snps}

    except subprocess.CalledProcessError as e:
        logging.error(f"show-snps failed for {delta_path}: {e.stderr}")
        return {"total_snps": 0, "coding_snps": 0}
    except FileNotFoundError:
        logging.error("show-snps command not found. Ensure MUMmer is in your PATH.")
        return {"total_snps": 0, "coding_snps": 0}


def build_matrices_for_group(
    group_id: str,
    sequence_ids: list,
    group_outdir: str,
    cds_map: Optional[Dict[str, List[tuple]]] = None,
):
    """
    Parses delta files and constructs square DataFrames for Identity, Coverage,
    Total SNPs, and Coding SNPs.
    """
    df_identity = pd.DataFrame(100.0, index=sequence_ids, columns=sequence_ids)
    df_coverage = pd.DataFrame(100.0, index=sequence_ids, columns=sequence_ids)
    df_total_snps = pd.DataFrame(0, index=sequence_ids, columns=sequence_ids)
    df_coding_snps = pd.DataFrame(0, index=sequence_ids, columns=sequence_ids)

    pairs = list(itertools.combinations(sequence_ids, 2))

    for seq_a, seq_b in pairs:
        delta_file = os.path.join(group_outdir, f"{seq_a}_vs_{seq_b}.delta")

        # Nucmer is directional. The sequence acting as the reference determines the coordinates.
        reference_seq = seq_a

        if not os.path.exists(delta_file):
            delta_file = os.path.join(group_outdir, f"{seq_b}_vs_{seq_a}.delta")
            reference_seq = seq_b

        if os.path.exists(delta_file):
            # Fetch coordinates for whichever sequence acted as the reference
            ref_cds_coords = cds_map.get(reference_seq) if cds_map else None

            coords_metrics = parse_delta_with_coords(delta_file)
            snp_metrics = parse_delta_with_snps(delta_file, cds_coords=ref_cds_coords)

            # Populate Identity & Coverage
            df_identity.loc[seq_a, seq_b] = df_identity.loc[seq_b, seq_a] = (
                coords_metrics["identity"]
            )
            avg_cov = (
                coords_metrics["coverage_ref"] + coords_metrics["coverage_que"]
            ) / 2
            df_coverage.loc[seq_a, seq_b] = df_coverage.loc[seq_b, seq_a] = avg_cov

            # Populate Total and Coding SNPs
            df_total_snps.loc[seq_a, seq_b] = df_total_snps.loc[seq_b, seq_a] = (
                snp_metrics["total_snps"]
            )
            df_coding_snps.loc[seq_a, seq_b] = df_coding_snps.loc[seq_b, seq_a] = (
                snp_metrics["coding_snps"]
            )
        else:
            df_identity.loc[seq_a, seq_b] = df_identity.loc[seq_b, seq_a] = 0.0
            df_coverage.loc[seq_a, seq_b] = df_coverage.loc[seq_b, seq_a] = 0.0
            df_total_snps.loc[seq_a, seq_b] = df_total_snps.loc[seq_b, seq_a] = 0
            df_coding_snps.loc[seq_a, seq_b] = df_coding_snps.loc[seq_b, seq_a] = 0

    return df_identity, df_coverage, df_total_snps, df_coding_snps


def parse_delta_with_coords(delta_path: str) -> dict:
    """
    Runs show-coords on a delta file and calculates total alignment coverage
    and average identity between the reference and query.
    """
    # Run show-coords with tab-delimited output (-T)
    # -r sorts by ref, -l includes sequence lengths, -c includes coverage
    cmd = ["show-coords", "-T", "-r", "-l", "-c", delta_path]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"show-coords failed for {delta_path}: {e.stderr}")
        return {"identity": 0.0, "coverage_ref": 0.0, "coverage_que": 0.0}
    except FileNotFoundError:
        logging.error("show-coords command not found. Ensure MUMmer is in your PATH.")
        return {"identity": 0.0, "coverage_ref": 0.0, "coverage_que": 0.0}

    lines = result.stdout.strip().split("\n")

    # Find where the data actually starts (skipping MUMmer headers)
    data_start_idx = 0
    for i, line in enumerate(lines):
        if line.startswith("NUCMER") or line.startswith("[S1]"):
            data_start_idx = i + 1
            continue

    data_lines = lines[data_start_idx:]
    if not data_lines:
        return {"identity": 0.0, "coverage_ref": 0.0, "coverage_que": 0.0}

    total_iden = 0.0
    total_aln_len = 0

    # We track covered positions using sets to accurately handle multiple
    # overlapping alignment blocks without double-counting coverage.
    ref_covered_bases = set()
    que_covered_bases = set()
    ref_len = 1
    que_len = 1

    for line in data_lines:
        parts = line.split("\t")
        if len(parts) < 13:
            continue

        # show-coords -T layout:
        # [S1] [E1] [S2] [E2] [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [% COV R] [% COV Q] [TAG R] [TAG Q]
        s1, e1 = sorted([int(parts[0]), int(parts[1])])
        s2, e2 = sorted([int(parts[2]), int(parts[3])])
        aln_len_ref = int(parts[4])
        identity = float(parts[6])

        ref_len = max(ref_len, int(parts[7]))
        que_len = max(que_len, int(parts[8]))

        # Calculate weighted identity based on alignment block length
        total_iden += identity * aln_len_ref
        total_aln_len += aln_len_ref

        # Log specific bases covered
        ref_covered_bases.update(range(s1, e1 + 1))
        que_covered_bases.update(range(s2, e2 + 1))

    if total_aln_len == 0:
        return {"identity": 0.0, "coverage_ref": 0.0, "coverage_que": 0.0}

    weighted_identity = total_iden / total_aln_len
    coverage_ref = (len(ref_covered_bases) / ref_len) * 100
    coverage_que = (len(que_covered_bases) / que_len) * 100

    return {
        "identity": weighted_identity,
        "coverage_ref": coverage_ref,
        "coverage_que": coverage_que,
    }


def run_nucmer_cohorts(
    pipeline_groups: Dict[str, List[str]],
    fasta_map: Dict[str, str],
    outdir: Path,
    nucmer_opts: Optional[str] = "",
) -> None:
    """Runs MUMmer's nucmer tool, parses metrics, and builds pairwise heatmaps."""

    # Process nucmer options safely
    opts_list = [opt for opt in (nucmer_opts or "").split(" ") if opt]

    for group_id, sequence_ids in pipeline_groups.items():
        logging.info(
            f"Processing Nucmer visualizations for {group_id} ({len(sequence_ids)} members)..."
        )

        local_analysis_dir = Path(outdir) / group_id / "nucmer_results"
        local_analysis_dir.mkdir(parents=True, exist_ok=True)

        # Generate all unique pairs within the group (n choose 2)
        pairs = list(itertools.combinations(sequence_ids, 2))

        for seq_a, seq_b in pairs:
            path_a = fasta_map.get(seq_a)
            path_b = fasta_map.get(seq_b)

            if not path_a or not path_b:
                missing = seq_a if not path_a else seq_b
                logging.warning(
                    f"[Warning] Skipping pair {seq_a} <-> {seq_b}: '{missing}' not found in fasta_map."
                )
                continue

            out_prefix = os.path.join(local_analysis_dir, f"{seq_a}_vs_{seq_b}")

            # Added the opts_list into the command stream
            command = ["nucmer"] + opts_list + ["--prefix", out_prefix, path_a, path_b]

            try:
                subprocess.run(
                    command,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            except subprocess.CalledProcessError as e:
                logging.error(f"[Error] Nucmer failed for {seq_a} vs {seq_b}!")
                logging.error(f"Stderr: {e.stderr.strip()}")
            except FileNotFoundError:
                logging.warning(
                    "[Critical] 'nucmer' command not found. Is MUMmer installed and in your PATH?"
                )
                return

        # --- Visualization & Export Step (Inside the main group loop) ---
        logging.info(f"Building matrices, CSVs, and heatmaps for {group_id}...")

        # Unpack all three dataframes now
        # TODO : limit SNPs to CDS
        df_id, df_cov, df_snps, df_cds_snps = build_matrices_for_group(
            group_id, sequence_ids, str(local_analysis_dir)
        )

        # Define output file destinations
        id_csv_path = local_analysis_dir / f"{group_id}_identity_matrix.csv"
        cov_csv_path = local_analysis_dir / f"{group_id}_coverage_matrix.csv"
        snps_csv_path = local_analysis_dir / f"{group_id}_snps_matrix.csv"

        id_heatmap_path = local_analysis_dir / f"{group_id}_identity_heatmap.png"
        cov_heatmap_path = local_analysis_dir / f"{group_id}_coverage_heatmap.png"
        snps_heatmap_path = local_analysis_dir / f"{group_id}_snps_heatmap.png"

        # 1. EXPORT ALL CSV FILES
        df_id.to_csv(id_csv_path)
        df_cov.to_csv(cov_csv_path)
        df_snps.to_csv(snps_csv_path)

        # 2. GENERATE ALL HEATMAPS
        visualize_nucmer_heatmap(df_id, id_heatmap_path, metric_type="identity")
        visualize_nucmer_heatmap(df_cov, cov_heatmap_path, metric_type="coverage")
        visualize_nucmer_heatmap(df_snps, snps_heatmap_path, metric_type="snps")

        logging.info(
            f"Identity, Coverage, and SNP data fully processed for {group_id}!\n"
        )
