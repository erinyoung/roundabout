import os
import logging
import itertools
import statistics
from pathlib import Path

from .database import run_setup
from .run_annotation import (
    execute_amrfinder_parallel,
    execute_plasmidfinder_parallel,
    execute_bakta_parallel
)
from .run_parsing import (
    parse_amrfinder_results, 
    parse_plasmidfinder_results
)
from .run_similarity import (
    execute_skani, 
    parse_skani_results,
    create_local_ani_matrix,
    visualize_skani_matrix,
    extract_global_hits,
    visualize_global_matches_scatter,  # Updated import
    build_sentinel_groups
)
from .run_group_analysis import analyze_cohort_group
from .run_heatcluster import run_heatcluster

def stage_and_split_fastas(input_dir: Path, staging_dir: Path) -> list[Path]:
    """
    Reads FASTA files, splits multi-FASTAs into single sequences, 
    and saves them as individual files in the staging directory.
    """
    logging.info(f"Staging and splitting FASTA files from {input_dir} into {staging_dir}")
    staging_dir.mkdir(parents=True, exist_ok=True)
    
    staged_paths = []
    
    for ext in ("*.fa", "*.fasta", "*.fna"):
        for file_path in input_dir.glob(ext):
            with open(file_path, 'r') as fasta_file:
                content = fasta_file.read()
                
            records = content.split('>')
            
            for rec in records[1:]: 
                lines = rec.strip().split('\n')
                if not lines:
                    continue
                    
                header_id = lines[0].split()[0]
                # Sanitize the ID for the filesystem
                safe_id = "".join(c for c in header_id if c.isalnum() or c in ('_', '-'))
                out_path = staging_dir / f"{file_path.stem}_{safe_id}.fasta"
                
                with open(out_path, 'w') as out_f:
                    out_f.write(f">{rec.strip()}\n")
                    
                staged_paths.append(out_path)
                
    return staged_paths

def run_pipeline(args):
    """Orchestrates the Roundabout horizontal gene transfer pipeline."""
    
    logging.info("Checking databases...")
    db_dict = run_setup(args)

    logging.info("Starting Roundabout Pipeline.")
    logging.info(f"Outputs will be saved to: {args.outdir}")
    os.makedirs(args.outdir, exist_ok=True)
    
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

    staging_dir = Path(args.outdir) / "staging_fastas"
    staged_fasta_paths = stage_and_split_fastas(fasta_dir, staging_dir)
    
    if not staged_fasta_paths:
        logging.error("No valid FASTA sequences found to process.")
        raise FileNotFoundError("No valid FASTA sequences found to process.")
        
    logging.info(f"Successfully staged {len(staged_fasta_paths)} individual FASTA sequences.")

    # -------------------------------------------------------------------------
    # STEP 2: Annotation
    # -------------------------------------------------------------------------
    
    # # AMRFinderPlus
    # amr_dict = {}
    # if db_dict.get("amrfinder"):
    #     amr_dict = execute_amrfinder_parallel(staged_fasta_paths, args.outdir, db_dict["amrfinder"], args.threads)
    # else:
    #     logging.warning("AMRFinderPlus database missing; skipping AMRFinderPlus.")

    # total_isolates = len(amr_dict)
    # isolates_with_genes = sum(1 for genes in amr_dict.values() if genes)
    # empty_isolates = total_isolates - isolates_with_genes

    # logging.info(f"AMRFinderPlus parsing complete: {total_isolates} total isolates processed. "
    #             f"{isolates_with_genes} carried AMR genes, {empty_isolates} had no hits.")

    # # PlasmidFinder
    # pf_dict = {}
    # if db_dict.get("plasmidfinder"):
    #     pf_dict = execute_plasmidfinder_parallel(staged_fasta_paths, args.outdir, db_dict["plasmidfinder"], args.threads)
    # else:
    #     logging.warning("PlasmidFinder database missing; skipping PlasmidFinder.")

    # total_isolates = len(pf_dict)
    # isolates_with_plasmids = sum(1 for plasmids in pf_dict.values() if plasmids)
    # empty_isolates = total_isolates - isolates_with_plasmids

    # logging.info(f"PlasmidFinder parsing complete: {total_isolates} isolates processed. "
    #              f"{isolates_with_plasmids} carried plasmids, {empty_isolates} had no hits.")

    # TODO: Insertion sequences
    # TODO: Prophages

    # Bakta
    # if db_dict.get("bakta"):
    #    execute_bakta_parallel(staged_fasta_paths, args.outdir, db_dict["bakta"], args.threads)
    # else:
    #    logging.warning("Bakta database missing; skipping Bakta.")

    # -------------------------------------------------------------------------
    # STEP 3: Sequence Comparison with Skani
    # -------------------------------------------------------------------------

    if db_dict.get("refseq_plasmid_dl"):
        db_path = Path(db_dict["refseq_plasmid_dl"]) / "refseq_plasmids_dl.fasta"
        refseq_fasta = db_path if db_path.exists() else None
        if not refseq_fasta:
            logging.warning("RefSeq multi-FASTA not found. Running SKANI on local inputs only.")
    else:
        refseq_fasta = None
        logging.warning("RefSeq database not provided. Running SKANI on local inputs only.")
            
    # Execute Skani
    # Execute Skani
    skani_tsv = execute_skani(staged_fasta_paths, refseq_fasta, Path(args.outdir), args.threads)

    # Parse results into a raw pandas DataFrame
    raw_skani_df = parse_skani_results(skani_tsv)
    
    # ---> ADD THIS LINE <---
    skani_outdir = Path(args.outdir) / "skani_results"
    
    # 1. Handle Local Matrix and Visualization
    local_matrix_df = create_local_ani_matrix(raw_skani_df, staged_fasta_paths)
    heatmap_out = skani_outdir / "local_ani_heatmap.png"
    visualize_skani_matrix(local_matrix_df, heatmap_out)
    
    # Save the local matrix to CSV for your records
    local_matrix_df.to_csv(skani_outdir / "local_ani_matrix.csv")

    # 2. Handle Global (RefSeq) Hits
    global_hits_df = extract_global_hits(raw_skani_df, staged_fasta_paths)
    
    if not global_hits_df.empty:
        # Save CSV 
        global_csv_out = skani_outdir / "global_refseq_hits.csv"
        global_hits_df.to_csv(global_csv_out, index=False)
        logging.info(f"Saved global RefSeq hits to {global_csv_out}")
        
        # Generate Global Scatter Plot
        global_scatter_out = skani_outdir / "global_ani_scatter.png"
        visualize_global_matches_scatter(global_hits_df, global_scatter_out)
        logging.info(f"Generated global ANI scatter plot at {global_scatter_out}")
    else:
        logging.info("No global RefSeq hits were found or RefSeq was not provided.")

    exit(0)

    # -------------------------------------------------------------------------
    # STEP 4: Grouping
    # -------------------------------------------------------------------------

    # Similarity grouping
    input_names = [f.stem for f in staged_fasta_paths]
    
    # Dynamically extract any RefSeq hits from the matrix
    refseq_names = set()
    for targets in sim_matrix.values():
        for r_name in targets.keys():
            if r_name not in input_names:
                refseq_names.add(r_name)
    
    min_identity = getattr(args, 'min_identity', 95.0)
    min_coverage = getattr(args, 'min_coverage', 80.0) 
    
    l_strict, l_contain, g_strict, g_contain = build_sentinel_groups(
        sim_matrix, input_names, list(refseq_names), min_identity, min_coverage
    )

    # Grouping based on AMRFinderPlus
    amr_dir = Path(args.outdir) / "amrfinder_results"
    amr_groups = parse_amrfinder_results(amr_dir, substring=args.amr_gene)

    # Grouping based on PlasmidFinder results
    
    pf_dir = Path(args.outdir) / "plasmidfinder_results"   
    inc_groups = parse_plasmidfinder_results(pf_dir, substring=args.plasmidfinder_string)
    
    # 3. Create discrete groups for exactly what the parsers found
    unique_groups = {}
    
    def register_cohort(samples: list, category: str, detail: str):
        if len(samples) > 1:
            cohort = frozenset(samples)
            if cohort not in unique_groups:
                unique_groups[cohort] = {'amr': [], 'inc': [], 'identity': []}
            unique_groups[cohort][category].append(detail)

    for amr_gene, samples in amr_groups.items(): register_cohort(samples, 'amr', amr_gene)
    for inc_gene, samples in inc_groups.items(): register_cohort(samples, 'inc', inc_gene)
    for sentinel, samples in l_strict.items(): register_cohort(samples, 'identity', f"Local Strict (Sentinel: {sentinel})")
    for sentinel, samples in l_contain.items(): register_cohort(samples, 'identity', f"Local Contained (Sentinel: {sentinel})")
    for sentinel, samples in g_strict.items(): register_cohort(samples, 'identity', f"Global Strict (Sentinel: {sentinel})")
    for sentinel, samples in g_contain.items(): register_cohort(samples, 'identity', f"Global Contained (Sentinel: {sentinel})")
    
    logging.info("=" * 60)
    logging.info("ROUNDABOUT PIPELINE GROUPING SUMMARY")
    logging.info("=" * 60)
    logging.info(f"Consolidated into {len(unique_groups)} distinct sequence cohorts.")
    logging.info("=" * 60)
    
    # -------------------------------------------------------------------------
    # STEP 5: Process the Cohorts
    # -------------------------------------------------------------------------

    # TODO: run minkemap on each cohort to generate a visual representation of the groupings with and without references
    # TODO: run daisyblast on each cohort to generate visualizations without references
    # TODO: run pygenomeviz on bakta out of inputs to create visualizations without references
    # TODO: create nucmer comparisions of each cohort with and without references to generate visualizations
    # TODO: pass nucmer output to heatcluster for visualization

    bakta_dir = Path(args.outdir) / "bakta_results"

    for idx, (samples_set, traits) in enumerate(unique_groups.items(), start=1):
        group_name = f"group_{idx}"
        samples_list = list(samples_set)
        
        # Calculate pairwise similarities for this specific cohort using the skani matrix
        pairwise_anis = []
        for a, b in itertools.combinations(samples_list, 2):
            ani = sim_matrix.get(a, {}).get(b, {}).get('ani')
            if ani is None:
                ani = sim_matrix.get(b, {}).get(a, {}).get('ani')
            pairwise_anis.append(ani if ani is not None else 0.0)

        if pairwise_anis:
            min_sim = f"{min(pairwise_anis):.2f}%"
            max_sim = f"{max(pairwise_anis):.2f}%"
            avg_sim = f"{statistics.mean(pairwise_anis):.2f}%"
        else:
            min_sim = max_sim = avg_sim = "N/A"
            
        group_outdir = Path(args.outdir) / group_name
        group_outdir.mkdir(parents=True, exist_ok=True)
        
        # Determine ALL shared unfiltered traits for this cohort
        all_shared_amrs = [gene for gene, members in unfiltered_amr.items() if set(samples_list).issubset(set(members))]
        all_shared_incs = [inc for inc, members in unfiltered_inc.items() if set(samples_list).issubset(set(members))]
        
        # The info.txt now states exactly WHY this specific group was formed, plus ALL shared traits
        info_file = group_outdir / "info.txt"
        with open(info_file, 'w') as f:
            f.write(f"Group Name: {group_name}\n")
            f.write(f"Input Sequences ({len(samples_list)}): {', '.join(samples_list)}\n")
            f.write(f"Grouped By Shared AMR Trigger: {', '.join(traits['amr']) if traits['amr'] else 'N/A'}\n")
            f.write(f"Grouped By Shared Plasmid Trigger: {', '.join(traits['inc']) if traits['inc'] else 'N/A'}\n")
            f.write(f"Grouped By Identity Trigger: {', '.join(traits['identity']) if traits['identity'] else 'N/A'}\n")
            f.write("-" * 40 + "\n")
            f.write(f"ALL Shared AMR Genes (Unfiltered): {', '.join(all_shared_amrs) if all_shared_amrs else 'None'}\n")
            f.write(f"ALL Shared Plasmid Replicons (Unfiltered): {', '.join(all_shared_incs) if all_shared_incs else 'None'}\n")
            f.write("-" * 40 + "\n")
            f.write(f"Min Similarity: {min_sim}\n")
            f.write(f"Max Similarity: {max_sim}\n")
            f.write(f"Average Similarity: {avg_sim}\n")
            
        # Execute downstream analysis (Passing sim_matrix and threads!)
        analyze_cohort_group(group_name, samples_list, staging_dir, group_outdir, sim_matrix, args)

    logging.info("Roundabout pipeline execution finished successfully!")