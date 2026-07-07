import logging
from pathlib import Path

# Import the daisyblast modules directly to bypass the CLI
from daisyblast.blast_runner import perform_self_blast
from daisyblast.parse_blast_and_shatter import parse_blast_and_shatter
from daisyblast.trim_blast import trim_blast
from daisyblast.grouping import hit_grouping
from daisyblast.visualize_groups import visualize_groups
from daisyblast.dotplot import generate_dotplots
from daisyblast.summary import generate_summary_plots

def run_daisyblast_cohorts(pipeline_groups: dict, fasta_map: dict, outdir: Path, args):
    """
    Cycles through each cohort group and runs the DaisyBlast analysis exactly once
    per group, isolating the output into its own group-specific directory.
    """
    daisy_base_dir = Path(outdir) / "daisyblast_results"
    
    # Safely extract your tool-prefixed CLI parameters
    evalue = getattr(args, 'daisyblast_evalue', 1e-10)
    min_pident = getattr(args, 'daisyblast_min_pident', 90.0)
    min_len = getattr(args, 'daisyblast_min_length', 200)
    num_groups = getattr(args, 'daisyblast_num_groups', 20)

    for group_id, members in pipeline_groups.items():
        logging.info(f"Running DaisyBlast natively for {group_id} ({len(members)} members)...")
        
        # 1. Create isolated subdirectory: results/daisyblast_results/group_0001
        group_out_dir = daisy_base_dir / group_id
        group_out_dir.mkdir(parents=True, exist_ok=True)
        
        # 2. Map group member IDs to their physical staged FASTA paths
        group_fasta_paths = [fasta_map[m_id] for m_id in members if m_id in fasta_map]
        
        if not group_fasta_paths:
            logging.warning(f"No valid FASTA paths found for group {group_id}, skipping.")
            continue
            
        try:
            # 3. Run Self-BLAST
            blast_result_file, header_map = perform_self_blast(
                input_files=group_fasta_paths, 
                output_dir=group_out_dir, 
                evalue=evalue, 
                min_pident=min_pident
            )

            # 4. Parse and Shatter
            shattered_file = parse_blast_and_shatter(
                blast_file=blast_result_file,
                outdir=group_out_dir,
                min_len=min_len,
                identity_cutoff=min_pident
            )

            # 5. Trim
            trimmed_file = trim_blast(blast_result_file, shattered_file, group_out_dir)

            # 6. Grouping
            grouped_file = hit_grouping(
                bed_file=shattered_file, 
                blast_file=trimmed_file, 
                outdir=group_out_dir, 
                max_groups=num_groups
            )

            # 7. Visualizations
            visualize_groups(grouped_file, group_out_dir)
            generate_dotplots(blast_result_file, group_out_dir)
            generate_summary_plots(blast_result_file, group_out_dir)

            logging.info(f"DaisyBlast execution for {group_id} completed successfully.")
            
        except Exception as e:
            logging.error(f"DaisyBlast failed during execution for {group_id}: {e}")
            # Keep moving to process the next cohort if one fails
            continue