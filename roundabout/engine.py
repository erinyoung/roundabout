import os
import logging
from pathlib import Path

from .database import run_setup
from .run_annotation import (
    execute_amrfinder_parallel,
    execute_plasmidfinder_parallel,
    execute_bakta_parallel
)
from .run_grouping import (
    define_groups_by_similarity,
    define_groups_by_amr,
    define_groups_by_plasmidfinder,
    build_consensus_groups,
    write_group_summary
)
from .run_similarity import (
    execute_skani, 
    parse_skani_results,
    create_local_ani_matrix,
    extract_global_hits
)
from .run_visualize_similarity import (
    visualize_as_heatmap,
    visualize_global_matches_scatter
)


def stage_and_split_fastas(input_dir: Path, staging_dir: Path, min_length: int = 0, max_length: int = None) -> list[dict]:
    """
    Reads FASTA files, splits multi-FASTAs into single sequences, 
    filters them by length, and saves them as individual files.
    Returns a list of dictionaries containing file metadata.
    """
    logging.info(f"Staging and splitting FASTA files from {input_dir} into {staging_dir}")
    staging_dir.mkdir(parents=True, exist_ok=True)
    
    staged_data = []
    
    for ext in ("*.fa", "*.fasta", "*.fna"):
        for file_path in input_dir.glob(ext):
            with open(file_path, 'r') as fasta_file:
                content = fasta_file.read()
                
            records = content.split('>')
            
            for rec in records[1:]: 
                lines = rec.strip().split('\n')
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
                contig_id = "".join(c for c in header_id if c.isalnum() or c in ('_', '-'))
                
                # Capture the precise naming metadata
                parent_name = file_path.stem
                combined_id = f"{parent_name}_{contig_id}"
                out_path = staging_dir / f"{combined_id}.fasta"
                
                with open(out_path, 'w') as out_f:
                    out_f.write(f">{rec.strip()}\n")
                    
                # Store it all in a dictionary
                staged_data.append({
                    "path": out_path,              
                    "parent_name": parent_name,    
                    "contig_id": contig_id,        
                    "combined_id": combined_id,
                    "length": seq_length           # Added length to metadata just in case!
                })
                
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
        "options": args.bakta_options
    }

    # AMRFinderPlus Options
    amrfinder_opts = {
        "db": args.amrfinder_db,
        "organism": args.amr_organism,
        "gene": args.amr_gene,
        "ident_min": args.amr_ident_min,
        "coverage_min": args.amr_coverage_min,
        "options": args.amrfinder_options
    }

    # PlasmidFinder Options
    plasmidfinder_opts = {
        "db": args.plasmidfinder_db,
        "filter_string": args.plasmidfinder_string,
        "mincov": args.plasmidfinder_mincov,
        "threshold": args.plasmidfinder_threshold
    }

    # refseq-plasmid-dl Options
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
        "max_collection_date": args.refseq_max_collection_date
    }

    # Skani Options
    skani_opts = {
        "min_af": args.skani_min_af,
        "both_min_af": args.skani_both_min_af,
        "ci": args.skani_ci,
        "detailed": args.skani_detailed,
        "n": args.skani_n,
        "short_header": args.skani_short_header,
        "fast": args.skani_fast,
        "medium": args.skani_medium,
        "slow": args.skani_slow,
        "small_genomes": args.skani_small_genomes,
        "c": args.skani_c,
        "faster_small": args.skani_faster_small,
        "m": args.skani_m,
        "median": args.skani_median,
        "no_learned_ani": args.skani_no_learned_ani,
        "no_marker_index": args.skani_no_marker_index,
        "robust": args.skani_robust,
        "s": args.skani_s
    }

    # Daisyblast Options
    daisyblast_opts = {
        "evalue": args.daisyblast_evalue,
        "min_pident": args.daisyblast_min_pident,
        "min_length": args.daisyblast_min_length,
        "num_groups": args.daisyblast_num_groups
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
        "min_coverage": args.minkemap_min_coverage
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
        "skip_blast": args.pgv_skip_blast,
        "skip_mummer": args.pgv_skip_mummer,
        "skip_mmseqs": args.pgv_skip_mmseqs,
        "skip_pmauve": args.pgv_skip_pmauve,
        "min_identity": args.pgv_min_identity,
        "min_length": args.pgv_min_length
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
        max_length=max_contig_length
    )
    
    if not staged_fasta_data:
        logging.error("No valid FASTA sequences found to process.")
        raise FileNotFoundError("No valid FASTA sequences found to process.")
        
    logging.info(f"Successfully staged {len(staged_fasta_data)} individual FASTA sequences.")

    # If your execution functions just need a flat list of Path objects:
    staged_fasta_paths = [item["path"] for item in staged_fasta_data]

    # -------------------------------------------------------------------------
    # STEP 2: Annotation
    # -------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------
    # AMRFinderPlus Execution
    # -------------------------------------------------------------------------
    # amr_dict = {}
    # TODO filter genes OUT HERE so the full list still makes it to the final report, but only the filtered ones are used for grouping
    # if db_dict.get("amrfinder"):
    #     amr_dict = execute_amrfinder_parallel(
    #         staged_fasta_paths, 
    #         outdir, 
    #         db_dict["amrfinder"], 
    #         threads, 
    #         amrfinder_opts
    #     )
    # else:
    #     logging.warning("AMRFinderPlus database missing; skipping AMRFinderPlus.")

    # total_isolates = len(amr_dict)
    # isolates_with_genes = sum(1 for genes in amr_dict.values() if genes)
    # empty_isolates = total_isolates - isolates_with_genes

    # logging.info(f"AMRFinderPlus parsing complete: {total_isolates} total isolates processed. "
    #              f"{isolates_with_genes} carried AMR genes, {empty_isolates} had no hits.")
    
    # print(amr_dict)  # <-- Debugging line to inspect the AMRFinderPlus results
    # line below can be used for debugging to skip the amrfinder step
    amr_dict = {'4051900_4': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-5', 'mph(A)', 'mrx(A)'], '4051900_3': [], '4051899_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-5', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT'], '4051901_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-98', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT'], '4051902_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-98', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT'], '4051903_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-98', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT'], '4051904_3': [], '4051904_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-98', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT'], '4051905_3': [], '4051905_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-98', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT'], '4051906_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-5', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT'], '4051906_3': [], '4051907_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-5', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT'], '4051908_2': ['dfrA12', 'aadA2', 'qacEdelta1', 'sul1', 'ble', 'blaNDM-98', 'mph(A)', 'mrx(A)', 'tet(A)', 'traT']}


    # -------------------------------------------------------------------------
    # PlasmidFinder Execution
    # -------------------------------------------------------------------------
    pf_dict = {}
    # TODO filter plasmid identifiers OUT HERE so the full list still makes it to the final report, but only the filtered ones are used for grouping
    # if db_dict.get("plasmidfinder"):
    #     pf_dict = execute_plasmidfinder_parallel(
    #         staged_fasta_paths, 
    #         outdir, 
    #         db_dict["plasmidfinder"], 
    #         threads,
    #         plasmidfinder_opts
    #     )
    # else:
    #     logging.warning("PlasmidFinder database missing; skipping PlasmidFinder.")

    # total_isolates = len(pf_dict)
    # isolates_with_plasmids = sum(1 for plasmids in pf_dict.values() if plasmids)
    # empty_isolates = total_isolates - isolates_with_plasmids

    # logging.info(f"PlasmidFinder parsing complete: {total_isolates} isolates processed. "
    #              f"{isolates_with_plasmids} carried plasmids, {empty_isolates} had no hits.")

    # print(pf_dict)  # <-- Debugging line to inspect the PlasmidFinder results
    # line below can be used for debugging to skip the plasmidfinder step for dev purposes
    pf_dict = {'4051901_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)'], '4051900_4': [], '4051900_3': ['IncM1'], '4051899_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)'], '4051902_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)'], '4051903_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)'], '4051904_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)'], '4051904_3': ['IncM1'], '4051905_3': ['IncM1'], '4051905_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)'], '4051906_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)'], '4051906_3': ['IncM1'], '4051907_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)'], '4051908_2': ['IncFII', 'IncFIA', 'IncFIB(AP001918)']}

    # TODO: Insertion sequences
    # TODO: Prophages

    # -------------------------------------------------------------------------
    # Bakta Execution
    # -------------------------------------------------------------------------
    # if db_dict.get("bakta"):
    #     execute_bakta_parallel(
    #         staged_fasta_paths, 
    #         outdir, 
    #         db_dict["bakta"], 
    #         threads,
    #         bakta_opts
    #     )
    # else:
    #     logging.warning("Bakta database missing; skipping Bakta.")

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
            logging.warning("RefSeq multi-FASTA not found. Running SKANI on local inputs only.")
            
        if not refseq_meta_csv.exists():
            refseq_meta_csv = None
            logging.warning("RefSeq metadata CSV not found. Global hits will lack context.")
    else:
        refseq_fasta = None
        refseq_meta_csv = None
        logging.warning("RefSeq database not provided. Running SKANI on local inputs only.")
            
    # Execute Skani
    # skani_tsv = execute_skani(staged_fasta_paths, refseq_fasta, Path(args.outdir), threads, skani_opts)

    # print(skani_tsv)
    skani_tsv = 'results/skani_results/skani_matrix.tsv'

    # Parse results into a raw pandas DataFrame
    raw_skani_df = parse_skani_results(skani_tsv)
    
    # Extract matrices and hits
    local_matrix_df = create_local_ani_matrix(raw_skani_df, staged_fasta_paths)
    global_hits_df = extract_global_hits(raw_skani_df, staged_fasta_paths, refseq_meta_csv)

    # -------------------------------------------------------------------------
    # Visualizations & Clustering
    # -------------------------------------------------------------------------
    logging.info("Generating clustered ANI heatmap...")
    
    # We already defined outdir = Path(args.outdir) at the top of run_pipeline
    skani_outdir = outdir / "skani_results"
    
    # 1. Local Clustermap
    basic_heatmap_path = skani_outdir / "local_ani_clustermap.png"
    visualize_as_heatmap(local_matrix_df, basic_heatmap_path)
    logging.info(f"Saved ANI clustermap to {basic_heatmap_path}")

    # 2. Global Scatter Plot
    logging.info("Generating global RefSeq matches scatter plot...")
    global_scatter_path = skani_outdir / "global_ani_scatter.png"
    visualize_global_matches_scatter(global_hits_df, global_scatter_path)
    logging.info(f"Saved global scatter plot to {global_scatter_path}")

    # -------------------------------------------------------------------------
    # STEP 4: Grouping
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Group Assignment & Reporting
    # -------------------------------------------------------------------------
    seq_to_similarity_group, unique_similarity_groups = define_groups_by_similarity(
        local_matrix_df, 
        global_hits_df, 
        min_ani,
        min_ani_fraction_query,
        min_ani_fraction_ref,
        num_references
    )
    
    seq_to_amr_group, unique_amr_groups = define_groups_by_amr(amr_dict, target_gene=args.amr_gene)

    seq_to_plasmidfinder_group, unique_plasmidfinder_groups = define_groups_by_plasmidfinder(
        pf_dict=pf_dict, 
        target_replicon=args.plasmidfinder_string
    )

    seq_to_consensus, unique_consensus_groups = build_consensus_groups(
        seq_to_sim=seq_to_similarity_group,
        seq_to_amr=seq_to_amr_group,
        seq_to_pf=seq_to_plasmidfinder_group
    )

    # Generate and save the sheet using the outdir path
    group_report_csv = write_group_summary(
        master_groups=unique_consensus_groups,
        unique_similarity_groups=unique_similarity_groups,
        amr_dict=amr_dict,
        pf_dict=pf_dict,
        local_matrix_df=local_matrix_df,
        global_hits_df=global_hits_df,
        outdir=Path(args.outdir), # or however your outdir is defined
        target_amr_gene=args.amr_gene,
        target_pf_replicon=args.plasmidfinder_string
    )

    exit(0)
    
    # -------------------------------------------------------------------------
    # STEP 5: Process the Cohorts
    # -------------------------------------------------------------------------

    # TODO: run minkemap on each cohort to generate a visual representation of the groupings with and without references
    # TODO: run daisyblast on each cohort to generate visualizations without references
    # TODO: run pygenomeviz on bakta out of inputs to create visualizations without references
    # TODO: create nucmer comparisions of each cohort with and without references to generate visualizations
    # TODO: pass nucmer output to heatcluster for visualization
    # TODO: create a heatmap of the ANI matrix for each cohort with and without references to generate visualizations

    bakta_dir = Path(args.outdir) / "bakta_results"

    logging.info("Roundabout pipeline execution finished successfully!")