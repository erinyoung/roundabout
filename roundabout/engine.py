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
    parse_skani_matrix, 
    build_sentinel_groups
)
from .run_group_analysis import analyze_cohort_group

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
    # STEP 2: Parallel Annotation
    # -------------------------------------------------------------------------
    
    # AMRFinderPlus
    if db_dict.get("amrfinder"):
        execute_amrfinder_parallel(staged_fasta_paths, args.outdir, db_dict["amrfinder"], args.threads)
    else:
        logging.warning("AMRFinderPlus database missing; skipping AMRFinderPlus.")

    # PlasmidFinder
    if db_dict.get("plasmidfinder"):
        execute_plasmidfinder_parallel(staged_fasta_paths, args.outdir, db_dict["plasmidfinder"], args.threads)
    else:
        logging.warning("PlasmidFinder database missing; skipping PlasmidFinder.")

    # Bakta
    #if db_dict.get("bakta"):
    #    execute_bakta_parallel(staged_fasta_paths, args.outdir, db_dict["bakta"], args.threads)
    #else:
    #    logging.warning("Bakta database missing; skipping Bakta.")

    # -------------------------------------------------------------------------
    # STEP 3: Sequence Comparison & Grouping
    # -------------------------------------------------------------------------
 
    amr_dir = Path(args.outdir) / "amrfinder_results"
    pf_dir = Path(args.outdir) / "plasmidfinder_results"
    
    # 1. Execute the biological marker parsers
    amr_groups = parse_amrfinder_results(amr_dir, substring=args.amr_gene)
    inc_groups = parse_plasmidfinder_results(pf_dir, substring=args.plasmidfinder_string)
    
    # 2. Execute skani for Sequence Identity Grouping
    refseq_fasta = None
    if db_dict.get("refseq_plasmid_dl"):
        db_path = Path(db_dict["refseq_plasmid_dl"]) / "refseq_plasmids_dl.fasta"
        if db_path.exists():
            refseq_fasta = db_path
        else:
            logging.warning("RefSeq multi-FASTA not found. Skipping global identity grouping.")
            
    # Execute Skani and parse matrix
    skani_tsv = execute_skani(staged_fasta_paths, refseq_fasta, Path(args.outdir), args.threads)
    sim_matrix = parse_skani_matrix(skani_tsv)
    
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
    
    # 3. Create discrete groups for exactly what the parsers found
    unique_groups = {}
    
    def register_cohort(samples: list, category: str, detail: str):
        """Registers a group exactly as defined. Only merges if the exact same sequences share another trait."""
        if len(samples) > 1:
            cohort = frozenset(samples)
            if cohort not in unique_groups:
                unique_groups[cohort] = {'amr': [], 'inc': [], 'identity': []}
            unique_groups[cohort][category].append(detail)

    # Rule 1: Group by shared AMR (divergent or not)
    for amr_gene, samples in amr_groups.items():
        register_cohort(samples, 'amr', amr_gene)
        
    # Rule 2: Group by shared Plasmid targets
    for inc_gene, samples in inc_groups.items():
        register_cohort(samples, 'inc', inc_gene)
        
    # Rule 3: Group by sequence identity (Sentinel Networks)
    for sentinel, samples in l_strict.items():
        register_cohort(samples, 'identity', f"Local Strict (Sentinel: {sentinel})")
    for sentinel, samples in l_contain.items():
        register_cohort(samples, 'identity', f"Local Contained (Sentinel: {sentinel})")
    for sentinel, samples in g_strict.items():
        register_cohort(samples, 'identity', f"Global Strict (Sentinel: {sentinel})")
    for sentinel, samples in g_contain.items():
        register_cohort(samples, 'identity', f"Global Contained (Sentinel: {sentinel})")
    
    logging.info("=" * 60)
    logging.info("ROUNDABOUT PIPELINE GROUPING SUMMARY")
    logging.info("=" * 60)
    logging.info(f"Consolidated into {len(unique_groups)} distinct sequence cohorts.")
    logging.info("=" * 60)
    
    # -------------------------------------------------------------------------
    # STEP 4: Process the Cohorts and Write Output
    # -------------------------------------------------------------------------

    # TODO : list all AMR genes that are shared within a cohort, and all plasmid replicons that are shared within a cohort, in the info.txt file for each group, even if they were filtered out.
    # TODO : for input files that share plasmids, or AMR genes, or are similar: run daisyblast, heatcluster (with subset of skani results), and minkemap
    # TODO : for input files that are similar to the RefSeq plasmids, run minkemap, run heatcluster, but with a subset of the matrix that drops sequences that are too divergent
    # TODO : run heatcluster on the entire skani matrix for all input sequences
    # TODO : replace nucmer with skani. there should be no nucmer references
    # TODO : ensure that threads information is being passed to all relevant functions

    for idx, (samples_set, traits) in enumerate(unique_groups.items(), start=1):
        group_name = f"group_{idx}"
        samples_list = list(samples_set)
        
        # Calculate pairwise similarities for this specific cohort using the skani matrix
        pairwise_anis = []
        for a, b in itertools.combinations(samples_list, 2):
            # Check a -> b
            ani = sim_matrix.get(a, {}).get(b, {}).get('ani')
            if ani is None:
                # Check b -> a (since matrix might be directional)
                ani = sim_matrix.get(b, {}).get(a, {}).get('ani')
                
            # If skani didn't report an alignment, they are too divergent (0%)
            pairwise_anis.append(ani if ani is not None else 0.0)

        if pairwise_anis:
            min_sim = f"{min(pairwise_anis):.2f}%"
            max_sim = f"{max(pairwise_anis):.2f}%"
            avg_sim = f"{statistics.mean(pairwise_anis):.2f}%"
        else:
            min_sim = max_sim = avg_sim = "N/A"
            
        group_outdir = Path(args.outdir) / group_name
        group_outdir.mkdir(parents=True, exist_ok=True)
        
        # The info.txt now states exactly WHY this specific group was formed
        info_file = group_outdir / "info.txt"
        with open(info_file, 'w') as f:
            f.write(f"Group Name: {group_name}\n")
            f.write(f"Input Sequences ({len(samples_list)}): {', '.join(samples_list)}\n")
            f.write(f"Grouped By Shared AMR: {', '.join(traits['amr']) if traits['amr'] else 'N/A'}\n")
            f.write(f"Grouped By Shared Plasmid: {', '.join(traits['inc']) if traits['inc'] else 'N/A'}\n")
            f.write(f"Grouped By Identity: {', '.join(traits['identity']) if traits['identity'] else 'N/A'}\n")
            f.write(f"Min Similarity: {min_sim}\n")
            f.write(f"Max Similarity: {max_sim}\n")
            f.write(f"Average Similarity: {avg_sim}\n")
            
        # Execute downstream analysis
        analyze_cohort_group(group_name, samples_list, staging_dir, group_outdir, args)

    logging.info("Roundabout pipeline execution finished successfully!")