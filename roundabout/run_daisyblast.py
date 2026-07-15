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


def run_daisyblast_cohorts(
    pipeline_groups: dict, fasta_map: dict, outdir: Path, daisyblast_opts: dict
):
    """
    Cycles through each cohort group and runs the DaisyBlast analysis exactly once
    per group, isolating the output into its own group-specific directory.
    """

    # Safely extract your CLI parameters from the passed dictionary
    evalue = daisyblast_opts.get("evalue", 1e-10)
    min_pident = daisyblast_opts.get("min_pident", 90.0)
    min_len = daisyblast_opts.get("min_length", 200)
    num_groups = daisyblast_opts.get("num_groups", 20)

    for group_id, members in pipeline_groups.items():
        logging.info(f"Running DaisyBlast for {group_id} ({len(members)} members)")

        # 1. Create isolated subdirectory: results/[group_id]/daisyblast_results
        group_out_dir = Path(outdir) / group_id / "daisyblast_results"
        group_out_dir.mkdir(parents=True, exist_ok=True)

        # 2. Map group member IDs to their physical staged FASTA paths
        group_fasta_paths = [fasta_map[m_id] for m_id in members if m_id in fasta_map]

        if not group_fasta_paths:
            logging.warning(
                f"No valid FASTA paths found for group {group_id}, skipping."
            )
            continue

        try:
            # 3. Run Self-BLAST
            logging.info(f"Running Self-BLAST in DaisyBlast for {group_id}")
            blast_result_file, header_map = perform_self_blast(
                input_files=group_fasta_paths,
                output_dir=group_out_dir,
                evalue=evalue,
                min_pident=min_pident,
            )

            # 4. Parse and Shatter
            logging.info(f"Running Parse and Shatter in DaisyBlast for {group_id}")
            shattered_file = parse_blast_and_shatter(
                blast_file=blast_result_file,
                outdir=group_out_dir,
                min_len=min_len,
                identity_cutoff=min_pident,
            )

            # 5. Trim
            logging.info(f"Running Trim in DaisyBlast for {group_id}")
            trimmed_file = trim_blast(blast_result_file, shattered_file, group_out_dir)

            # 6. Grouping
            logging.info(f"Running Grouping in DaisyBlast for {group_id}")
            grouped_file = hit_grouping(
                bed_file=shattered_file,
                blast_file=trimmed_file,
                outdir=group_out_dir,
                max_groups=num_groups,
            )

            # 7. Visualizations
            logging.info(f"Running Visualizations in DaisyBlast for {group_id}")
            visualize_groups(grouped_file, group_out_dir)
            generate_dotplots(blast_result_file, group_out_dir)
            generate_summary_plots(blast_result_file, group_out_dir)

            logging.info(f"DaisyBlast execution for {group_id} completed successfully.")

        except Exception as e:
            logging.error(f"DaisyBlast failed during execution for {group_id}: {e}")
            # Keep moving to process the next cohort if one fails
            continue
