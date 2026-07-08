import logging
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.colors as mcolors

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.align import Blast, MUMmer, ProgressiveMauve, MMseqs, AlignCoord


def get_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """
    Parses a FASTA file and returns a dictionary of {sequence_id: length}.
    This ensures segment names match the IDs that BLAST uses.
    """
    seq_lengths = {}
    seq_id = None
    length = 0

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save the previous sequence if it exists
                if seq_id is not None:
                    seq_lengths[seq_id] = length
                # Extract the BLAST-compatible sequence ID (first word after '>')
                seq_id = line[1:].split()[0]
                length = 0
            else:
                length += len(line)

    # Add the final sequence
    if seq_id is not None:
        seq_lengths[seq_id] = length

    return seq_lengths


def get_input_order(members: list[str], local_matrix_df: pd.DataFrame) -> list[Path]:
    """
    Sorts a list of FASTA paths mathematically using hierarchical clustering
    derived from the local ANI matrix. Highly similar sequences will be adjacent.
    """
    if len(members) < 3:
        return members

    # Extract sub-matrix and convert similarity to distance
    sub_matrix = local_matrix_df.loc[members, members].copy().fillna(80.0)
    distance_matrix = 100.0 - sub_matrix.values

    # Force strict symmetry for SciPy
    distance_matrix = (distance_matrix + distance_matrix.T) / 2.0
    np.fill_diagonal(distance_matrix, 0.0)

    # Perform clustering and extract optimal order
    condensed_dist = squareform(distance_matrix)
    Z = linkage(condensed_dist, method="average")
    leaf_order = leaves_list(Z)

    return [members[i] for i in leaf_order]


def initial_gv_canvas(pgv_opts: dict) -> GenomeViz:
    """
    Initializes the pyGenomeViz canvas using the shared layout options.
    """
    gv = GenomeViz(
        fig_width=pgv_opts.get("fig_width", 15.0),
        fig_track_height=pgv_opts.get("fig_track_height", 1.0),
        feature_track_ratio=pgv_opts.get("feature_track_ratio", 0.25),
        track_align_type=pgv_opts.get("track_align_type", "center"),
    )

    return gv


def render_links_and_save(
    gv: GenomeViz, align_coords: list, out_path: Path, pgv_opts: dict, method: str
):
    """
    Generalized renderer. Applies different link styles and colorbars based
    on the alignment method (e.g., identity vs blocks).
    """
    if not align_coords:
        logging.warning(f"No alignment links passed filtering thresholds for {method}.")
        gv.savefig(out_path, dpi=pgv_opts.get("dpi", 300))
        return

    # 1. Grab link styles
    color = pgv_opts.get("normal_link_color", "grey")
    inverted_color = pgv_opts.get("inverted_link_color", "red")
    curve = pgv_opts.get("curve", False)

    # 2. Check alignment data types
    # 2. Check alignment data types
    has_identity = any(getattr(ac, "identity", None) is not None for ac in align_coords)

    # Safely check for any block identifier (handles different pygenomeviz versions)
    def get_block_id(ac):
        return (
            getattr(ac, "group", None)
            or getattr(ac, "lcb_id", None)
            or getattr(ac, "group_id", None)
        )

    has_group = any(get_block_id(ac) is not None for ac in align_coords)

    # Route A: Identity Gradient (BLAST, MUMmer)
    if has_identity and method != "pmauve":
        vmin = float(pgv_opts.get("identity_thr", 30.0))
        for ac in align_coords:
            gv.add_link(
                ac.query_link,
                ac.ref_link,
                color=color,
                inverted_color=inverted_color,
                v=ac.identity,
                vmin=vmin,
                vmax=100,
                curve=curve,
            )
        gv.set_colorbar(
            [color, inverted_color],
            vmin=vmin,
            vmax=100,
            bar_width=pgv_opts.get("cbar_width", 0.01),
            bar_height=pgv_opts.get("cbar_height", 0.2),
        )

    # Route B: Block Clustering (progressiveMauve)
    # Route B: Block Clustering (progressiveMauve)
    elif has_group or method == "pmauve":
        block_cmap = pgv_opts.get("block_cmap", "gist_rainbow")

        # 1. Collect block group IDs (natively 'group' in pygenomeviz for Mauve)
        group_ids = [
            ac.group for ac in align_coords if getattr(ac, "group", None) is not None
        ]

        # Failsafe: if the parser stripped them, fall back to a uniform single group
        vmin_group = min(group_ids) if group_ids else 1
        vmax_group = max(group_ids) if group_ids else 1

        for ac in align_coords:
            group_val = getattr(ac, "group", None) or 1

            # 2. Let pygenomeviz natively map the group ID to both the track blocks
            # AND the linking ribbons using the internal colormap logic.
            gv.add_link(
                ac.query_link,
                ac.ref_link,
                v=group_val,
                vmin=vmin_group,
                vmax=vmax_group,
                cmap=block_cmap,
                curve=curve,
            )

    # Route C: Fallback (Solid colors, no identity gradient)
    else:
        for ac in align_coords:
            gv.add_link(
                ac.query_link,
                ac.ref_link,
                color=color,
                inverted_color=inverted_color,
                curve=curve,
            )

    # 3. Handle global canvas elements
    scale_size = pgv_opts.get("scale_labelsize", 15.0)
    if pgv_opts.get("show_scale_bar", False):
        gv.set_scale_bar(labelsize=scale_size)
    elif pgv_opts.get("show_scale_xticks", False):
        gv.set_scale_xticks(labelsize=scale_size)

    # 4. Save
    gv.savefig(out_path, dpi=pgv_opts.get("dpi", 300))


def run_pygenomeviz_blast(
    input_paths: list[Path],
    out_path: Path,
    file_type: str = "fasta",
    pgv_opts: dict = None,
):
    """
    Builds tracks from sequences, runs BLAST, and plots synteny.
    Automatically handles nucleotide BLAST for FASTAs and protein BLAST for GBFFs.
    """
    if pgv_opts is None:
        pgv_opts = {}

    logging.info(
        f"Running pygenomeviz BLAST on {len(input_paths)} {file_type} sequences..."
    )
    gv = initial_gv_canvas(pgv_opts)

    track_labelsize = pgv_opts.get("track_labelsize", 20.0)

    # Parse Feature Colors
    feature_colors = {}
    for item in pgv_opts.get("feature_type2color", ["CDS:black"]):
        if ":" in item:
            ftype, fcolor = item.split(":", 1)
            feature_colors[ftype] = fcolor

    label_type = pgv_opts.get("feature_labeltype", "None")
    label_type = None if label_type == "None" else label_type

    # Initialize lists to hold the specific inputs BLAST needs for each type
    blast_inputs = []
    blast_seqtype = "nucleotide"

    # 1. Build Tracks & Collect BLAST Inputs
    if file_type == "fasta":
        blast_seqtype = "nucleotide"

        for fasta in input_paths:
            seq_lengths_dict = get_fasta_lengths(fasta)
            gv.add_feature_track(
                fasta.stem, seq_lengths_dict, labelsize=track_labelsize
            )
            # FASTAs just need string paths passed to BLAST
            blast_inputs.append(str(fasta))

    elif file_type == "gbff":
        blast_seqtype = "protein"  # Switch to protein alignment

        for gbk_path in input_paths:
            gbk = Genbank(gbk_path)
            # GenBank files should have the parsed objects passed to BLAST
            blast_inputs.append(gbk)

            track = gv.add_feature_track(
                gbk.name, gbk.get_seqid2size(), labelsize=track_labelsize
            )

            for ftype, fcolor in feature_colors.items():
                for seqid, features in gbk.get_seqid2features(ftype).items():
                    segment = track.get_segment(seqid)
                    segment.add_features(
                        features,
                        plotstyle=pgv_opts.get("feature_plotstyle", "arrow"),
                        fc=fcolor,
                        lw=pgv_opts.get("feature_linewidth", 0.0),
                        label_type=label_type,
                        text_kws={"size": pgv_opts.get("feature_labelsize", 8.0)},
                    )

    # 2. Run Alignment
    try:
        align_coords = Blast(blast_inputs, seqtype=blast_seqtype).run()
    except Exception as e:
        logging.error(f"pygenomeviz BLAST failed: {e}")
        return

    # 3. Filter and Render
    align_coords = AlignCoord.filter(
        align_coords,
        length_thr=pgv_opts.get("length_thr", 100),
        identity_thr=pgv_opts.get("identity_thr", 30.0),
    )

    render_links_and_save(gv, align_coords, out_path, pgv_opts, method="blast")


def run_pygenomeviz_mummer(
    input_paths: list[Path],
    out_path: Path,
    file_type: str = "fasta",
    pgv_opts: dict = None,
):
    """
    Builds tracks from sequences, runs MUMmer, and plots synteny.
    Automatically handles nucleotide MUMmer (nucmer) for FASTAs
    and protein MUMmer (promer) for GBFFs.
    """
    if pgv_opts is None:
        pgv_opts = {}

    logging.info(
        f"Running pygenomeviz MUMmer on {len(input_paths)} {file_type} sequences..."
    )
    gv = initial_gv_canvas(pgv_opts)

    track_labelsize = pgv_opts.get("track_labelsize", 20.0)

    # Parse Feature Colors
    feature_colors = {}
    for item in pgv_opts.get("feature_type2color", ["CDS:black"]):
        if ":" in item:
            ftype, fcolor = item.split(":", 1)
            feature_colors[ftype] = fcolor

    label_type = pgv_opts.get("feature_labeltype", "None")
    label_type = None if label_type == "None" else label_type

    # Initialize lists to hold the specific inputs MUMmer needs
    mummer_inputs = []
    mummer_seqtype = "nucleotide"

    # 1. Build Tracks & Collect MUMmer Inputs
    if file_type == "fasta":
        mummer_seqtype = "nucleotide"

        for fasta in input_paths:
            seq_lengths_dict = get_fasta_lengths(fasta)
            gv.add_feature_track(
                fasta.stem, seq_lengths_dict, labelsize=track_labelsize
            )
            mummer_inputs.append(str(fasta))

    elif file_type == "gbff":
        mummer_seqtype = "protein"

        for gbk_path in input_paths:
            gbk = Genbank(gbk_path)
            mummer_inputs.append(gbk)

            track = gv.add_feature_track(
                gbk.name, gbk.get_seqid2size(), labelsize=track_labelsize
            )

            for ftype, fcolor in feature_colors.items():
                for seqid, features in gbk.get_seqid2features(ftype).items():
                    segment = track.get_segment(seqid)
                    segment.add_features(
                        features,
                        plotstyle=pgv_opts.get("feature_plotstyle", "arrow"),
                        fc=fcolor,
                        lw=pgv_opts.get("feature_linewidth", 0.0),
                        label_type=label_type,
                        text_kws={"size": pgv_opts.get("feature_labelsize", 8.0)},
                    )

    # 2. Run MUMmer Alignment
    try:
        # This acts exactly like the Blast class!
        align_coords = MUMmer(mummer_inputs, seqtype=mummer_seqtype).run()
    except Exception as e:
        logging.error(f"pygenomeviz MUMmer failed: {e}")
        return

    # 3. Filter Coordinates
    align_coords = AlignCoord.filter(
        align_coords,
        length_thr=pgv_opts.get("length_thr", 100),
        identity_thr=pgv_opts.get("identity_thr", 30.0),
    )
    logging.info(f"MUMmer produced {len(align_coords)} links")

    # 4. Render (Using the exact same generalized renderer we built!)
    render_links_and_save(gv, align_coords, out_path, pgv_opts, method="mummer")


def run_pygenomeviz_pmauve(
    input_paths: list[Path],
    out_path: Path,
    file_type: str = "fasta",
    pgv_opts: dict = None,
):
    """
    Builds tracks from sequences, runs progressiveMauve, and plots matching
    rainbow synteny blocks and links, fixing the native ColorCycler bug.
    """
    if pgv_opts is None:
        pgv_opts = {}

    if len(input_paths) < 2:
        logging.warning("Need at least 2 files to run progressiveMauve. Skipping.")
        return

    logging.info(
        f"Running pygenomeviz progressiveMauve on {len(input_paths)} {file_type} sequences..."
    )

    # 1. Initialize Canvas
    gv = initial_gv_canvas(pgv_opts)
    track_labelsize = pgv_opts.get("track_labelsize", 20.0)
    block_plotstyle = pgv_opts.get("block_plotstyle", "box")
    block_cmap = pgv_opts.get("block_cmap", "gist_rainbow")

    # Override link colors to match the block color if requested
    normal_link_color = pgv_opts.get("normal_link_color", "#323e4f")
    inverted_link_color = pgv_opts.get("inverted_link_color", "#3d6e70")

    pmauve_inputs = []
    name2seqlen = {}

    # 2. Build Flattened Tracks & Collect Inputs
    if file_type == "fasta":
        for fasta in input_paths:
            seq_lengths_dict = get_fasta_lengths(fasta)
            total_length = sum(seq_lengths_dict.values())
            name2seqlen[fasta.stem] = total_length
            pmauve_inputs.append(str(fasta))

    elif file_type == "gbff":
        for gbk_path in input_paths:
            gbk = Genbank(gbk_path)
            total_length = sum(gbk.get_seqid2size().values())
            name2seqlen[gbk.name] = total_length
            pmauve_inputs.append(gbk)

    # 3. Run progressiveMauve Alignment
    try:
        align_coords = ProgressiveMauve(
            pmauve_inputs, refid=pgv_opts.get("refid", 0)
        ).run()
    except Exception as e:
        logging.error(f"pygenomeviz progressiveMauve failed: {e}")
        return

    # Filter out tiny blocks based on length threshold
    align_coords = AlignCoord.filter(
        align_coords, length_thr=pgv_opts.get("length_thr", 500)
    )

    if not align_coords:
        logging.warning("No progressiveMauve alignments passed the thresholds.")
        gv.savefig(out_path, dpi=pgv_opts.get("dpi", 300))
        return

    # 4. Parse Synteny Blocks exactly like the official CLI script
    from collections import defaultdict

    name2blocks = defaultdict(list)
    for ac in align_coords:
        if ac.query_block not in name2blocks[ac.query_name]:
            name2blocks[ac.query_name].append(ac.query_block)
        if ac.ref_block not in name2blocks[ac.ref_name]:
            name2blocks[ac.ref_name].append(ac.ref_block)

    # 5. Build Global Coordinates-to-Color Map
    # Gather all unique blocks across the entire genome to establish color groups
    all_unique_blocks = set()
    for blocks in name2blocks.values():
        for b in blocks:
            all_unique_blocks.add((b[0], b[1]))  # Group by start and end positions

    sorted_unique_blocks = sorted(list(all_unique_blocks))
    num_unique_blocks = len(sorted_unique_blocks)

    # Map each physical block boundary to a consistent rainbow color
    cmap = mpl.colormaps[block_cmap]
    block2color = {}
    for idx, b_coords in enumerate(sorted_unique_blocks):
        norm_val = idx / (num_unique_blocks - 1) if num_unique_blocks > 1 else 0.0
        block2color[b_coords] = mcolors.to_hex(cmap(norm_val))

    # 6. Add Tracks and Features using the Unified Color Map
    for name, seqlen in name2seqlen.items():
        track = gv.add_feature_track(
            name=name,
            segments={name: seqlen},
            labelsize=track_labelsize,
        )
        for block in name2blocks[name]:
            # Fetch the color tied strictly to the coordinates of this block
            hex_color = block2color[(block[0], block[1])]
            track.add_feature(*block, plotstyle=block_plotstyle, fc=hex_color)

    # 7. Plot Connecting Links
    for ac in align_coords:
        # BONUS: Dynamically color links to match the block color they connect!
        hex_color = block2color.get(
            (ac.query_block[0], ac.query_block[1]), normal_link_color
        )

        gv.add_link(
            ac.query_link,
            ac.ref_link,
            curve=pgv_opts.get("curve", False),
            color=hex_color,
            inverted_color=hex_color,
        )

    # 8. Global Canvas Configurations & Save
    scale_size = pgv_opts.get("scale_labelsize", 15.0)
    if pgv_opts.get("show_scale_bar", False):
        gv.set_scale_bar(labelsize=scale_size)
    elif pgv_opts.get("show_scale_xticks", False):
        gv.set_scale_xticks(labelsize=scale_size)

    gv.savefig(out_path, dpi=pgv_opts.get("dpi", 300))


def run_pygenomeviz_mmseqs(
    input_paths: list[Path],
    out_path: Path,
    file_type: str = "gbff",
    pgv_opts: dict = None,
):
    """
    Builds tracks from GenBank sequences, runs MMseqs (translated protein alignment),
    and plots the synteny links with an identity colorbar. Only supports GBFF files.
    """
    if pgv_opts is None:
        pgv_opts = {}

    # Strict guard rail: Only allow execution on GBFF files
    if file_type != "gbff":
        logging.warning(
            f"MMseqs workflow in this pipeline is configured strictly for 'gbff' files. "
            f"Skipping execution for file_type='{file_type}'."
        )
        return

    if len(input_paths) < 2:
        logging.warning("Need at least 2 GBFF files to run MMseqs synteny. Skipping.")
        return

    logging.info(f"Running pygenomeviz MMseqs on {len(input_paths)} sequences...")

    # 1. Initialize Canvas
    gv = initial_gv_canvas(pgv_opts)
    track_labelsize = pgv_opts.get("track_labelsize", 20.0)

    # Parse Feature Colors
    feature_colors = {}
    for item in pgv_opts.get("feature_type2color", ["CDS:black"]):
        if ":" in item:
            ftype, fcolor = item.split(":", 1)
            feature_colors[ftype] = fcolor

    label_type = pgv_opts.get("feature_labeltype", "None")
    label_type = None if label_type == "None" else label_type

    mmseqs_inputs = []

    # 2. Build Tracks from Contig Sizes & Add GenBank Features
    for gbk_path in input_paths:
        gbk = Genbank(gbk_path)
        mmseqs_inputs.append(gbk)

        # MMseqs is contig-aware, so we pass the full parsed {seqid: size} dictionary
        track = gv.add_feature_track(
            gbk.name, gbk.get_seqid2size(), labelsize=track_labelsize
        )

        # Populate track features (CDS, tRNA, etc.)
        for ftype, fcolor in feature_colors.items():
            for seqid, features in gbk.get_seqid2features(ftype).items():
                segment = track.get_segment(seqid)
                segment.add_features(
                    features,
                    plotstyle=pgv_opts.get("feature_plotstyle", "arrow"),
                    fc=fcolor,
                    lw=pgv_opts.get("feature_linewidth", 0.0),
                    label_type=label_type,
                    text_kws={"size": pgv_opts.get("feature_labelsize", 8.0)},
                )

    # 3. Run MMseqs Alignment
    try:
        # MMseqs wrapper handles 6-frame translated alignments natively for Genbank objects
        align_coords = MMseqs(mmseqs_inputs).run()
    except Exception as e:
        logging.error(f"pygenomeviz MMseqs failed: {e}")
        return

    # 4. Filter Coordinates based on identity and length thresholds
    align_coords = AlignCoord.filter(
        align_coords,
        length_thr=pgv_opts.get("length_thr", 100),
        identity_thr=pgv_opts.get("identity_thr", 30.0),
    )
    logging.info(f"MMseqs produced {len(align_coords)} links")

    # 5. Render using our shared identity-gradient renderer (Route A)
    render_links_and_save(gv, align_coords, out_path, pgv_opts, method="mmseqs")


def run_pygenomeviz(
    pipeline_groups: dict[str, list[str]],
    fasta_map: dict[str, Path],
    gbff_map: dict[str, Path],
    local_matrix_df: pd.DataFrame,
    outdir: Path,
    method: str = "blast",
    pgv_opts: dict = None,
):
    """
    Main function to run pygenomeviz synteny visualizations based on the specified method.
    """

    for group_id, members in pipeline_groups.items():
        logging.info(f"Running pyGenomeViz for {group_id} ({len(members)} members)")

        group_out_dir = Path(outdir) / group_id / "pygenomeviz_results" / method
        group_out_dir.mkdir(parents=True, exist_ok=True)

        sorted_members = get_input_order(members, local_matrix_df)
        # members = ['4051899_2', '4051901_2', '4051902_2', '4051903_2', '4051904_2', '4051905_2', '4051906_2', '4051907_2', '4051908_2']
        # sorted_members = ['4051906_2', '4051904_2', '4051901_2', '4051903_2', '4051908_2', '4051907_2', '4051905_2', '4051899_2', '4051902_2']

        # 2. Map group member IDs to their physical staged FASTA paths
        fasta_paths = [fasta_map[m_id] for m_id in sorted_members if m_id in fasta_map]
        # group_fasta_paths = [PosixPath('results/staging_fastas/4051908_2.fasta'), PosixPath('results/staging_fastas/4051905_2.fasta'), PosixPath('results/staging_fastas/4051904_2.fasta'), PosixPath('results/staging_fastas/4051903_2.fasta'), PosixPath('results/staging_fastas/4051901_2.fasta'), PosixPath('results/staging_fastas/4051902_2.fasta')]

        if not fasta_paths:
            logging.warning(
                f"No valid FASTA paths found for group {group_id}, skipping."
            )
            continue

        if len(fasta_paths) < 2:
            logging.warning("Need at least 2 FASTA files to run synteny. Skipping.")
            continue

        gbff_paths = [gbff_map[m_id] for m_id in sorted_members if m_id in gbff_map]
        # gbff_paths = [PosixPath('results/bakta_results/4051908_2/4051908_2.gbff'), PosixPath('results/bakta_results/4051905_2/4051905_2.gbff'), PosixPath('results/bakta_results/4051904_2/4051904_2.gbff'), PosixPath('results/bakta_results/4051903_2/4051903_2.gbff'), PosixPath('results/bakta_results/4051901_2/4051901_2.gbff'), PosixPath('results/bakta_results/4051902_2/4051902_2.gbff')]

        if not gbff_paths:
            logging.warning(
                f"No valid FASTA paths found for group {group_id}, skipping."
            )
            continue

        if len(gbff_paths) < 2:
            logging.warning("Need at least 2 FASTA files to run synteny. Skipping.")
            continue

        # Call the appropriate method
        if method == "blast":
            if len(fasta_paths) >= 2:
                run_pygenomeviz_blast(
                    fasta_paths,
                    group_out_dir
                    / f"pygenomeviz_synteny_{group_id}_{method}_fasta.png",
                    "fasta",
                    pgv_opts,
                )
            if len(gbff_paths) >= 2:
                run_pygenomeviz_blast(
                    gbff_paths,
                    group_out_dir / f"pygenomeviz_synteny_{group_id}_{method}_gbff.png",
                    "gbff",
                    pgv_opts,
                )
        elif method == "mummer":
            if len(fasta_paths) >= 2:
                run_pygenomeviz_mummer(
                    fasta_paths,
                    group_out_dir
                    / f"pygenomeviz_synteny_{group_id}_{method}_fasta.png",
                    "fasta",
                    pgv_opts,
                )
            if len(gbff_paths) >= 2:
                run_pygenomeviz_mummer(
                    gbff_paths,
                    group_out_dir / f"pygenomeviz_synteny_{group_id}_{method}_gbff.png",
                    "gbff",
                    pgv_opts,
                )
        elif method == "pmauve":
            if len(fasta_paths) >= 2:
                run_pygenomeviz_pmauve(
                    fasta_paths,
                    group_out_dir
                    / f"pygenomeviz_synteny_{group_id}_{method}_fasta.png",
                    "fasta",
                    pgv_opts,
                )
        elif method == "mmseqs":
            if len(gbff_paths) >= 2:
                run_pygenomeviz_mmseqs(
                    gbff_paths,
                    group_out_dir / f"pygenomeviz_synteny_{group_id}_{method}_gbff.png",
                    "gbff",
                    pgv_opts,
                )
        else:
            logging.error(f"Unknown pygenomeviz method: {method}")
