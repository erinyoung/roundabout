import os
import logging
from pathlib import Path

from .database import run_setup
from .run_annotation import (
    execute_amrfinder_parallel,
    execute_plasmidfinder_parallel,
    execute_bakta_parallel,
)
from .run_grouping_by_similarity import define_groups_by_similarity
from .run_grouping_by_amr import define_groups_by_amr
from .run_grouping_by_pf import define_groups_by_plasmidfinder
from .run_grouping_summary import run_grouping_summary
from .run_extract_references import get_top_refseq_ids_for_sample, extract_refseq_fastas

from .run_similarity import (
    execute_skani,
    parse_skani_results,
    create_local_ani_matrix,
    extract_global_hits,
)
from .run_visualize_similarity import (
    visualize_as_heatmap,
    visualize_global_matches_scatter,
)

from .run_minkemap import run_minkemap_cohorts
from .run_daisyblast import run_daisyblast_cohorts
from .run_cohort_heatmaps import run_cohort_heatmaps
from .run_pygenomeviz import run_pygenomeviz


def stage_and_split_fastas(
    input_dir: Path, staging_dir: Path, min_length: int = 0, max_length: int = None
) -> list[dict]:
    """
    Reads FASTA files, splits multi-FASTAs into single sequences,
    filters them by length, and saves them as individual files.
    Returns a list of dictionaries containing file metadata.
    """
    logging.info(
        f"Staging and splitting FASTA files from {input_dir} into {staging_dir}"
    )
    staging_dir.mkdir(parents=True, exist_ok=True)

    staged_data = []

    for ext in ("*.fa", "*.fasta", "*.fna"):
        for file_path in input_dir.glob(ext):
            with open(file_path, "r") as fasta_file:
                content = fasta_file.read()

            records = content.split(">")

            for rec in records[1:]:
                lines = rec.strip().split("\n")
                if not lines:
                    continue

                # Calculate the actual sequence length (ignore the header)
                sequence = "".join(lines[1:])
                seq_length = len(sequence)

                # Apply length filters
                if min_length and seq_length < min_length:
                    continue
                if max_length and seq_length > max_length:
                    continue

                header_id = lines[0].split()[0]
                # Sanitize the ID for the filesystem
                contig_id = "".join(
                    c for c in header_id if c.isalnum() or c in ("_", "-")
                )

                # Capture the precise naming metadata
                parent_name = file_path.stem
                combined_id = f"{parent_name}_{contig_id}"
                out_path = staging_dir / f"{combined_id}.fasta"

                with open(out_path, "w") as out_f:
                    out_f.write(f">{rec.strip()}\n")

                # Store it all in a dictionary
                staged_data.append(
                    {
                        "path": out_path,
                        "parent_name": parent_name,
                        "contig_id": contig_id,
                        "combined_id": combined_id,
                        "length": seq_length,  # Added length to metadata just in case!
                    }
                )

    return staged_data


def run_pipeline(args):
    """Orchestrates the Roundabout horizontal gene transfer pipeline."""

    # -------------------------------------------------------------------------
    # Extract and Organize Input Arguments
    # -------------------------------------------------------------------------

    # System & Grouping Options
    threads = args.threads
    outdir = Path(args.outdir)
    min_ani = args.min_ani
    min_num_seqs = args.min_num_seqs
    min_ani_fraction_query = args.min_ani_align_fraction_query
    min_ani_fraction_ref = args.min_ani_align_fraction_ref
    num_references = args.num_ref
    min_contig_length = args.min_contig_length
    max_contig_length = args.max_contig_length

    # Bakta Options
    bakta_opts = {
        "db": args.bakta_db,
        "genus": args.bakta_genus,
        "species": args.bakta_species,
        "strain": args.bakta_strain,
        "plasmid": args.bakta_plasmid,
        "complete": args.bakta_complete,
        "prodigal_tf": args.bakta_prodigal_tf,
        "translation_table": args.bakta_translation_table,
        "gram": args.bakta_gram,
        "locus": args.bakta_locus,
        "locus_tag": args.bakta_locus_tag,
        "locus_tag_increment": args.bakta_locus_tag_increment,
        "keep_contig_headers": args.bakta_keep_contig_headers,
        "compliant": args.bakta_compliant,
        "replicons": args.bakta_replicons,
        "regions": args.bakta_regions,
        "proteins": args.bakta_proteins,
        "hmms": args.bakta_hmms,
        "meta": args.bakta_meta,
        "partial": args.bakta_partial,
        "skip_trna": args.bakta_skip_trna,
        "skip_tmrna": args.bakta_skip_tmrna,
        "skip_rrna": args.bakta_skip_rrna,
        "skip_ncrna": args.bakta_skip_ncrna,
        "skip_ncrna_region": args.bakta_skip_ncrna_region,
        "skip_crispr": args.bakta_skip_crispr,
        "skip_cds": args.bakta_skip_cds,
        "skip_pseudo": args.bakta_skip_pseudo,
        "skip_gap": args.bakta_skip_gap,
        "skip_ori": args.bakta_skip_ori,
        "skip_filter": args.bakta_skip_filter,
        "skip_plot": args.bakta_skip_plot,
        "options": args.bakta_options,
    }

    # AMRFinderPlus Options
    amrfinder_opts = {
        "db": args.amrfinder_db,
        "organism": args.amr_organism,
        "gene": args.amr_gene,
        "ident_min": args.amr_ident_min,
        "coverage_min": args.amr_coverage_min,
        "options": args.amrfinder_options,
    }

    # PlasmidFinder Options
    plasmidfinder_opts = {
        "db": args.plasmidfinder_db,
        "filter_string": args.plasmidfinder_string,
        "mincov": args.plasmidfinder_mincov,
        "threshold": args.plasmidfinder_threshold,
    }

    # refseq-plasmid-dl Options
    # TODO : add these to database download function
    refseq_opts = {
        "db": args.refseq_plasmid_dl_db,
        "organism": args.refseq_organism,
        "taxid": args.refseq_taxid,
        "strain": args.refseq_strain,
        "isolate": args.refseq_isolate,
        "host": args.refseq_host,
        "plasmid_name": args.refseq_plasmid_name,
        "geo_loc_name": args.refseq_geo_loc_name,
        "isolation_source": args.refseq_isolation_source,
        "min_length": args.refseq_min_length,
        "max_length": args.refseq_max_length,
        "topology": args.refseq_topology,
        "min_date": args.refseq_min_date,
        "max_date": args.refseq_max_date,
        "min_collection_date": args.refseq_min_collection_date,
        "max_collection_date": args.refseq_max_collection_date,
    }

    # Skani Options
    skani_opts = {
        "min_af": args.skani_min_af,
        "both_min_af": args.skani_both_min_af,
        "ci": args.skani_ci,
        "detailed": args.skani_detailed,
        "n": args.skani_n,
        "short_header": args.skani_short_header,
        "fast": args.skani_presets == "fast",
        "medium": args.skani_presets == "medium",
        "slow": args.skani_presets == "slow",
        "small_genomes": args.skani_presets == "small-genomes",
        "c": args.skani_c,
        "faster_small": args.skani_faster_small,
        "m": args.skani_m,
        "median": args.skani_median,
        "no_learned_ani": args.skani_no_learned_ani,
        "no_marker_index": args.skani_no_marker_index,
        "robust": args.skani_robust,
        "s": args.skani_s,
    }

    # Daisyblast Options
    daisyblast_opts = {
        "evalue": args.daisyblast_evalue,
        "min_pident": args.daisyblast_min_pident,
        "min_length": args.daisyblast_min_length,
        "num_groups": args.daisyblast_num_groups,
    }

    # MinkeMap Options
    minkemap_opts = {
        "palette": args.minkemap_palette,
        "track_width": args.minkemap_track_width,
        "track_gap": args.minkemap_track_gap,
        "dpi": args.minkemap_dpi,
        "no_backbone": args.minkemap_no_backbone,
        "no_legend": args.minkemap_no_legend,
        "label_size": args.minkemap_label_size,
        "title": args.minkemap_title,
        "gc_skew": args.minkemap_gc_skew,
        "annotations": args.minkemap_annotations,
        "highlights": args.minkemap_highlights,
        "exclude_genes": args.minkemap_exclude_genes,
        "min_identity": args.minkemap_min_identity,
        "min_coverage": args.minkemap_min_coverage,
    }

    # Heatcluster Options
    # heatcluster_opts = {
    #     "out": args.heatcluster_out,
    #     "k": args.heatcluster_k,
    #     "t": args.heatcluster_t,
    #     "auto_k": args.heatcluster_auto_k,
    #     "pca": args.heatcluster_pca,
    #     "pca_out": args.heatcluster_pca_out,
    #     "title": args.heatcluster_title,
    #     "cmap": args.heatcluster_cmap,
    #     "dpi": args.heatcluster_dpi,
    #     "no_annot": args.heatcluster_no_annot,
    #     "no_plot": args.heatcluster_no_plot,
    #     "width": args.heatcluster_width,
    #     "height": args.heatcluster_height,
    #     "font_scale": args.heatcluster_font_scale,
    #     "vmin": args.heatcluster_vmin,
    #     "vmax": args.heatcluster_vmax,
    #     "hide_below": args.heatcluster_hide_below,
    #     "hide_above": args.heatcluster_hide_above,
    #     "no_cluster": args.heatcluster_no_cluster,
    #     "dendrogram": args.heatcluster_dendrogram
    # }

    # PyGenomeViz Options
    pgv_opts = {
        "identity_thr": args.pgv_identity_thr,
        "evalue_thr": args.pgv_evalue_thr,
        "length_thr": args.pgv_length_thr,
        "cmap": args.pgv_cmap,
        "fig_width": args.pgv_fig_width,
        "fig_track_height": args.pgv_fig_track_height,
        "track_align_type": args.pgv_track_align_type,
        "feature_track_ratio": args.pgv_feature_track_ratio,
        "show_scale_bar": args.pgv_show_scale_bar,
        "show_scale_xticks": args.pgv_show_scale_xticks,
        "curve": args.pgv_curve,
        "formats": args.pgv_formats,
        "dpi": args.pgv_dpi,
        "track_labelsize": args.pgv_track_labelsize,
        "scale_labelsize": args.pgv_scale_labelsize,
        "normal_link_color": args.pgv_normal_link_color,
        "inverted_link_color": args.pgv_inverted_link_color,
        "segment_space": args.pgv_segment_space,
        "feature_type2color": args.pgv_feature_type2color,
        "pseudo_color": args.pgv_pseudo_color,
        "feature_plotstyle": args.pgv_feature_plotstyle,
        "feature_linewidth": args.pgv_feature_linewidth,
        "feature_labeltrack": args.pgv_feature_labeltrack,
        "feature_labeltype": args.pgv_feature_labeltype,
        "feature_labelsize": args.pgv_feature_labelsize,
        "cbar_width": args.pgv_cbar_width,
        "cbar_height": args.pgv_cbar_height,
        "block_plotstyle": args.pgv_block_plotstyle,
        "block_cmap": args.pgv_block_cmap,
        "skip_blast": args.pgv_skip_blast,
        "skip_mummer": args.pgv_skip_mummer,
        "skip_pmauve": args.pgv_skip_pmauve,
        "skip_mmseqs": args.pgv_skip_mmseqs,
    }

    logging.info("Checking databases...")
    db_dict = run_setup(args)

    logging.info("Starting Roundabout Pipeline.")
    logging.info(f"Outputs will be saved to: {outdir}")
    os.makedirs(outdir, exist_ok=True)

    # -------------------------------------------------------------------------
    # STEP 1: Stage and Split FASTAs
    # -------------------------------------------------------------------------
    if not args.fastas:
        logging.error("No input FASTA directory provided.")
        raise ValueError("Missing input FASTA directory.")

    fasta_dir = Path(args.fastas)
    if not fasta_dir.exists() or not fasta_dir.is_dir():
        logging.error(f"Invalid input directory: {fasta_dir}")
        raise FileNotFoundError(f"Invalid input directory: {fasta_dir}")

    staging_dir = outdir / "staging_fastas"
    staged_fasta_data = stage_and_split_fastas(
        fasta_dir,
        staging_dir,
        min_length=min_contig_length,
        max_length=max_contig_length,
    )

    # print(staged_fasta_data)
    # staged_fasta_data = [{'path': PosixPath('results/staging_fastas/4051899_2.fasta'), 'parent_name': '4051899', 'contig_id': '2', 'combined_id': '4051899_2', 'length': 133243}, {'path': PosixPath('results/staging_fastas/4051900_3.fasta'), 'parent_name': '4051900', 'contig_id': '3', 'combined_id': '4051900_3', 'length': 67221}, {'path': PosixPath('results/staging_fastas/4051900_4.fasta'), 'parent_name': '4051900', 'contig_id': '4', 'combined_id': '4051900_4', 'length': 17248}, {'path': PosixPath('results/staging_fastas/4051901_2.fasta'), 'parent_name': '4051901', 'contig_id': '2', 'combined_id': '4051901_2', 'length': 133243}, {'path': PosixPath('results/staging_fastas/4051902_2.fasta'), 'parent_name': '4051902', 'contig_id': '2', 'combined_id': '4051902_2', 'length': 133243}, {'path': PosixPath('results/staging_fastas/4051903_2.fasta'), 'parent_name': '4051903', 'contig_id': '2', 'combined_id': '4051903_2', 'length': 133243}, {'path': PosixPath('results/staging_fastas/4051904_2.fasta'), 'parent_name': '4051904', 'contig_id': '2', 'combined_id': '4051904_2', 'length': 133243}, {'path': PosixPath('results/staging_fastas/4051904_3.fasta'), 'parent_name': '4051904', 'contig_id': '3', 'combined_id': '4051904_3', 'length': 67221}, {'path': PosixPath('results/staging_fastas/4051905_2.fasta'), 'parent_name': '4051905', 'contig_id': '2', 'combined_id': '4051905_2', 'length': 133243}, {'path': PosixPath('results/staging_fastas/4051905_3.fasta'), 'parent_name': '4051905', 'contig_id': '3', 'combined_id': '4051905_3', 'length': 67221}, {'path': PosixPath('results/staging_fastas/4051906_2.fasta'), 'parent_name': '4051906', 'contig_id': '2', 'combined_id': '4051906_2', 'length': 133243}, {'path': PosixPath('results/staging_fastas/4051906_3.fasta'), 'parent_name': '4051906', 'contig_id': '3', 'combined_id': '4051906_3', 'length': 67221}, {'path': PosixPath('results/staging_fastas/4051907_2.fasta'), 'parent_name': '4051907', 'contig_id': '2', 'combined_id': '4051907_2', 'length': 134557}, {'path': PosixPath('results/staging_fastas/4051908_2.fasta'), 'parent_name': '4051908', 'contig_id': '2', 'combined_id': '4051908_2', 'length': 133243}]

    if not staged_fasta_data:
        logging.error("No valid FASTA sequences found to process.")
        raise FileNotFoundError("No valid FASTA sequences found to process.")

    logging.info(
        f"Successfully staged {len(staged_fasta_data)} individual FASTA sequences."
    )

    # If your execution functions just need a flat list of Path objects:
    staged_fasta_paths = [item["path"] for item in staged_fasta_data]

    # -------------------------------------------------------------------------
    # STEP 2: Annotation
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # AMRFinderPlus Execution
    # -------------------------------------------------------------------------
    amr_dict = {}
    if db_dict.get("amrfinder"):
        amr_dict = execute_amrfinder_parallel(
            staged_fasta_paths, outdir, db_dict["amrfinder"], threads, amrfinder_opts
        )
    else:
        logging.warning("AMRFinderPlus database missing; skipping AMRFinderPlus.")

    total_isolates = len(amr_dict)
    isolates_with_genes = sum(1 for genes in amr_dict.values() if genes)
    empty_isolates = total_isolates - isolates_with_genes

    logging.info(
        f"AMRFinderPlus parsing complete: {total_isolates} total isolates processed. "
        f"{isolates_with_genes} carried AMR genes, {empty_isolates} had no hits."
    )

    # print(amr_dict)
    # line below can be used for debugging to skip the amrfinder step
    # amr_dict = {
    #     "4051900_4": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-5",
    #         "mph(A)",
    #         "mrx(A)",
    #     ],
    #     "4051900_3": [],
    #     "4051899_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-5",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    #     "4051901_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-98",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    #     "4051902_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-98",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    #     "4051903_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-98",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    #     "4051904_3": [],
    #     "4051904_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-98",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    #     "4051905_3": [],
    #     "4051905_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-98",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    #     "4051906_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-5",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    #     "4051906_3": [],
    #     "4051907_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-5",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    #     "4051908_2": [
    #         "dfrA12",
    #         "aadA2",
    #         "qacEdelta1",
    #         "sul1",
    #         "ble",
    #         "blaNDM-98",
    #         "mph(A)",
    #         "mrx(A)",
    #         "tet(A)",
    #         "traT",
    #     ],
    # }

    # -------------------------------------------------------------------------
    # PlasmidFinder Execution
    # -------------------------------------------------------------------------
    pf_dict = {}
    if db_dict.get("plasmidfinder"):
        pf_dict = execute_plasmidfinder_parallel(
            staged_fasta_paths,
            outdir,
            db_dict["plasmidfinder"],
            threads,
            plasmidfinder_opts,
        )
    else:
        logging.warning("PlasmidFinder database missing; skipping PlasmidFinder.")

    total_isolates = len(pf_dict)
    isolates_with_plasmids = sum(1 for plasmids in pf_dict.values() if plasmids)
    empty_isolates = total_isolates - isolates_with_plasmids

    logging.info(
        f"PlasmidFinder parsing complete: {total_isolates} isolates processed. "
        f"{isolates_with_plasmids} carried plasmids, {empty_isolates} had no hits."
    )

    # print(pf_dict)  # <-- Debugging line to inspect the PlasmidFinder results
    # line below can be used for debugging to skip the plasmidfinder step for dev purposes
    # pf_dict = {
    #     "4051901_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    #     "4051900_4": [],
    #     "4051900_3": ["IncM1"],
    #     "4051899_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    #     "4051902_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    #     "4051903_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    #     "4051904_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    #     "4051904_3": ["IncM1"],
    #     "4051905_3": ["IncM1"],
    #     "4051905_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    #     "4051906_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    #     "4051906_3": ["IncM1"],
    #     "4051907_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    #     "4051908_2": ["IncFII", "IncFIA", "IncFIB(AP001918)"],
    # }

    # TODO: Insertion sequences
    # TODO: Prophages

    sample_ids = sorted(list(set(amr_dict.keys()) | set(pf_dict.keys())))

    # -------------------------------------------------------------------------
    # Bakta Execution
    # -------------------------------------------------------------------------
    staged_gbk_data = []
    if db_dict.get("bakta"):
        staged_gbk_data = execute_bakta_parallel(
            staged_fasta_paths, outdir, db_dict["bakta"], threads, bakta_opts
        )
    else:
        logging.warning("Bakta database missing; skipping Bakta.")

    # print(staged_gbk_data)
    # staged_gbk_data = [PosixPath('results/bakta_results/4051899_2/4051899_2.gbff'), PosixPath('results/bakta_results/4051900_3/4051900_3.gbff'), PosixPath('results/bakta_results/4051900_4/4051900_4.gbff'), PosixPath('results/bakta_results/4051901_2/4051901_2.gbff'), PosixPath('results/bakta_results/4051902_2/4051902_2.gbff'), PosixPath('results/bakta_results/4051903_2/4051903_2.gbff'), PosixPath('results/bakta_results/4051904_2/4051904_2.gbff'), PosixPath('results/bakta_results/4051904_3/4051904_3.gbff'), PosixPath('results/bakta_results/4051905_2/4051905_2.gbff'), PosixPath('results/bakta_results/4051905_3/4051905_3.gbff'), PosixPath('results/bakta_results/4051906_2/4051906_2.gbff'), PosixPath('results/bakta_results/4051906_3/4051906_3.gbff'), PosixPath('results/bakta_results/4051907_2/4051907_2.gbff'), PosixPath('results/bakta_results/4051908_2/4051908_2.gbff')]

    # -------------------------------------------------------------------------
    # STEP 3: Sequence Comparison with Skani
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Skani Execution
    # -------------------------------------------------------------------------
    if db_dict.get("refseq_plasmid_dl"):
        refseq_db_dir = Path(db_dict["refseq_plasmid_dl"])
        refseq_fasta = refseq_db_dir / "refseq_plasmids_dl.fasta"
        refseq_meta_csv = refseq_db_dir / "refseq_plasmids_dl_metadata.csv"

        if not refseq_fasta.exists():
            refseq_fasta = None
            logging.warning(
                "RefSeq multi-FASTA not found. Running SKANI on local inputs only."
            )

        if not refseq_meta_csv.exists():
            refseq_meta_csv = None
            logging.warning(
                "RefSeq metadata CSV not found. Global hits will lack context."
            )
    else:
        refseq_fasta = None
        refseq_meta_csv = None
        logging.warning(
            "RefSeq database not provided. Running SKANI on local inputs only."
        )

    # Execute Skani
    skani_tsv = execute_skani(
        staged_fasta_paths, refseq_fasta, Path(args.outdir), threads, skani_opts
    )

    # print(skani_tsv)
    # Use the value below for dev
    # skani_tsv = 'results/skani_results/skani_matrix.tsv'

    # Parse results into a raw pandas DataFrame
    raw_skani_df = parse_skani_results(skani_tsv)

    # Extract matrices and hits
    local_matrix_df = create_local_ani_matrix(raw_skani_df, staged_fasta_paths)
    global_hits_df = extract_global_hits(
        raw_skani_df, staged_fasta_paths, refseq_meta_csv
    )

    # -------------------------------------------------------------------------
    # Visualizations & Clustering
    # -------------------------------------------------------------------------
    logging.info("Generating clustered ANI heatmap...")

    # 1. Local Clustermap
    logging.info("Visualizing ANI of local input sequences as heatmap")
    skani_outdir = outdir / "skani_results"
    basic_heatmap_path = skani_outdir / "local_ani_clustermap.png"
    visualize_as_heatmap(local_matrix_df, basic_heatmap_path)
    logging.info(f"Saved ANI clustermap to {basic_heatmap_path}")

    # 2. Global Scatter Plot
    logging.info("Generating global RefSeq matches scatter plot")
    global_scatter_path = skani_outdir / "global_ani_scatter.png"
    visualize_global_matches_scatter(global_hits_df, global_scatter_path)
    logging.info(f"Saved global scatter plot to {global_scatter_path}")

    # -------------------------------------------------------------------------
    # STEP 4: Grouping
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Group Assignment & Reporting
    # -------------------------------------------------------------------------

    logging.info("Creating groups based on ANI similarity")
    similarity_groups = define_groups_by_similarity(
        raw_skani_df=raw_skani_df,
        local_matrix_df=local_matrix_df,
        min_ani=min_ani,
        min_ani_align_fraction_ref=min_ani_fraction_ref,
        min_ani_align_fraction_query=min_ani_fraction_query,
        outdir=outdir,
    )

    # print(similarity_groups)
    # similarity_groups = [['4051899_2', '4051901_2', '4051902_2', '4051903_2', '4051904_2', '4051905_2', '4051906_2', '4051907_2', '4051908_2'], ['4051900_3', '4051904_3', '4051905_3', '4051906_3'], ['4051900_4']]
    # similarity_groups = []

    logging.info("Creating groups based on AMR gene detection")
    amr_groups = define_groups_by_amr(
        amr_dict, target_gene=args.amr_gene, outdir=outdir
    )
    # print(amr_groups)
    # amr_groups = [['4051901_2', '4051902_2', '4051903_2', '4051904_2', '4051905_2', '4051908_2'], ['4051899_2', '4051900_4', '4051906_2', '4051907_2']]
    # amr_groups = []

    logging.info("Creating groups based on plasmidfinder markers")
    plasmidfinder_groups = define_groups_by_plasmidfinder(
        pf_dict=pf_dict, target_replicon=args.plasmidfinder_string, outdir=outdir
    )

    # print(plasmidfinder_groups)
    # plasmidfinder_groups = [['4051899_2', '4051901_2', '4051902_2', '4051903_2', '4051904_2', '4051905_2', '4051906_2', '4051907_2', '4051908_2'], ['4051900_3', '4051904_3', '4051905_3', '4051906_3']]
    # plasmidfinder_groups = []

    # -------------------------------------------------------------------------
    # Combine and Deduplicate Groups
    # -------------------------------------------------------------------------

    logging.info("Combining and Deduplicating groups")
    all_groups = amr_groups + plasmidfinder_groups + similarity_groups

    # Deduplicate by converting to sorted tuples
    unique_groups_set = {tuple(sorted(g)) for g in all_groups}

    # Convert back to lists, filtering out any groups smaller than min_num_seqs
    unique_groups = [list(g) for g in unique_groups_set if len(g) >= min_num_seqs]

    # Sort them by size (largest groups first) just to keep the output tidy
    unique_groups.sort(key=lambda g: (len(g), g), reverse=True)

    # print(unique_groups)
    # unique_groups = [['4051899_2', '4051901_2', '4051902_2', '4051903_2', '4051904_2', '4051905_2', '4051906_2', '4051907_2', '4051908_2'], ['4051901_2', '4051902_2', '4051903_2', '4051904_2', '4051905_2', '4051908_2'], ['4051900_3', '4051904_3', '4051905_3', '4051906_3'], ['4051899_2', '4051900_4', '4051906_2', '4051907_2'], ['4051900_4']]

    # -------------------------------------------------------------------------
    # Run Grouping Summary
    # -------------------------------------------------------------------------
    logging.info("Creating summaries for groups")
    pipeline_groups = run_grouping_summary(
        groups=unique_groups,
        staged_data=staged_fasta_data,
        amr_dict=amr_dict,
        plasmidfinder_dict=pf_dict,
        skani_df=local_matrix_df,
        outdir=outdir,
        min_num_seqs=args.min_num_seqs,
        target_amr_gene=args.amr_gene,
        target_pf_string=args.plasmidfinder_string,
    )

    # print(pipeline_groups)
    # pipeline_groups = {'group_0001': ['4051899_2', '4051900_4', '4051901_2', '4051902_2', '4051903_2', '4051904_2', '4051905_2', '4051906_2', '4051907_2', '4051908_2'], 'group_0002': ['4051901_2', '4051902_2', '4051903_2', '4051904_2', '4051905_2', '4051908_2'], 'group_0003': ['4051900_3', '4051904_3', '4051905_3', '4051906_3'], 'group_0004': ['4051899_2', '4051900_4', '4051906_2', '4051907_2']}

    # -------------------------------------------------------------------------
    # STEP 5: Process the Cohorts
    # -------------------------------------------------------------------------
    logging.info("Starting cohort comparative visualizations")

    fasta_map = {item["combined_id"]: item["path"] for item in staged_fasta_data}
    gbff_map = {path.stem: path for path in staged_gbk_data}
    amr_outdir = Path(args.outdir) / "amrfinder_results"

    # Default if no RefSeq database is available
    sample_to_ref_paths = {}

    # Identify the RefSeq database path if available
    if db_dict.get("refseq_plasmid_dl"):
        logging.info("Finding matches from refseq_plasmid_dl")
        refseq_db_dir = Path(db_dict["refseq_plasmid_dl"])
        refseq_fasta = refseq_db_dir / "refseq_plasmids_dl.fasta"

        sample_to_refs = {}
        all_needed_refs = set()

        # 1. Collect all required RefSeq IDs across ALL samples first
        for sample_id in sample_ids:
            logging.info(f"Finding reference matches for {sample_id}")
            refs = get_top_refseq_ids_for_sample(
                global_hits_df,
                sample_id,
                num_refs=num_references,
                min_ani=min_ani,
                min_af_query=min_ani_fraction_query,
                min_af_ref=min_ani_fraction_ref,
            )
            sample_to_refs[sample_id] = refs
            all_needed_refs.update(refs)

        # 2. Extract everything from the massive multi-FASTA in ONE single pass
        logging.info(f"Extracting matches from {refseq_fasta}")
        staged_ref_paths = extract_refseq_fastas(all_needed_refs, refseq_fasta, outdir)

        # Helper to mirror the filename cleaning logic used during extraction
        def clean_id(val):
            return Path(str(val)).stem.split(".")[0].strip()

        # 3. Map each sample back to its specific staged RefSeq FASTA paths
        logging.info("Mapping samples to references")
        for sample_id, refs in sample_to_refs.items():
            sample_to_ref_paths[sample_id] = []
            for ref in refs:
                cleaned_ref = clean_id(ref)
                if cleaned_ref in staged_ref_paths:
                    sample_to_ref_paths[sample_id].append(staged_ref_paths[cleaned_ref])

    # TODO : run all these in parallel with concurrent
    logging.info("Running MinkeMap")
    run_minkemap_cohorts(
        pipeline_groups=pipeline_groups,
        fasta_map=fasta_map,
        amr_outdir=amr_outdir,
        outdir=outdir,
        minkemap_opts=minkemap_opts,
        sample_to_ref_paths=sample_to_ref_paths,
        target_amr_gene=args.amr_gene,
    )

    logging.info("Running DaisyBlast")
    run_daisyblast_cohorts(
        pipeline_groups=pipeline_groups,
        fasta_map=fasta_map,
        outdir=outdir,
        daisyblast_opts=daisyblast_opts,
    )

    pgv_methods = []

    if not args.pgv_skip_blast:
        pgv_methods.append("blast")
    if not args.pgv_skip_mummer:
        pgv_methods.append("mummer")
    if not args.pgv_skip_mmseqs:
        pgv_methods.append("mmseqs")
    if not args.pgv_skip_pmauve:
        pgv_methods.append("pmauve")

    # pgv_methods = ['mmseqs']
    for method in pgv_methods:
        logging.info(f"Running pyGenomeViz with {method} alignment method")
        run_pygenomeviz(
            pipeline_groups=pipeline_groups,
            fasta_map=fasta_map,
            gbff_map=gbff_map,
            local_matrix_df=local_matrix_df,
            outdir=outdir,
            method=method,
            pgv_opts=pgv_opts,
        )

    # TODO: create nucmer comparisions of each cohort with and without references to generate visualizations
    # TODO: get heatmap of nucmer visualization

    run_cohort_heatmaps(
        pipeline_groups=pipeline_groups,
        local_matrix_df=local_matrix_df,
        global_hits_df=global_hits_df,
        outdir=outdir,
        sample_to_ref_paths=sample_to_ref_paths,
    )

    logging.info("Comparative genomic visualizations complete.")

    logging.info("Roundabout pipeline execution finished successfully!")
