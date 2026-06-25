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

def execute_daisyblast(staged_fasta_paths: list[Path], outdir: str, args) -> Path:
    """
    Wrapper to execute the DaisyBlast pipeline natively via its Python API.
    """
    logging.info("Running daisyblast natively for sequence similarity and alignments...")
    
    daisy_out = Path(outdir) / "daisyblast_results"
    daisy_out.mkdir(parents=True, exist_ok=True)
    
    try:
        # 1. Run Self-BLAST
        blast_result_file, header_map = perform_self_blast(
            input_files=staged_fasta_paths, 
            output_dir=daisy_out, 
            evalue=getattr(args, 'evalue', 1e-10), 
            min_pident=getattr(args, 'min_pident', 90.0)
        )

        # 2. Parse and Shatter
        shattered_file = parse_blast_and_shatter(
            blast_file=blast_result_file,
            outdir=daisy_out,
            min_len=getattr(args, 'min_length', 500),
            identity_cutoff=getattr(args, 'min_pident', 90.0)
        )

        # 3. Trim
        trimmed_file = trim_blast(blast_result_file, shattered_file, daisy_out)

        # 4. Grouping
        grouped_file = hit_grouping(
            bed_file=shattered_file, 
            blast_file=trimmed_file, 
            outdir=daisy_out, 
            max_groups=getattr(args, 'num_groups', 20)
        )

        # 5. Visualizations
        visualize_groups(grouped_file, daisy_out)
        generate_dotplots(blast_result_file, daisy_out)
        generate_summary_plots(blast_result_file, daisy_out)

        logging.info("daisyblast execution completed successfully.")
        return daisy_out
        
    except Exception as e:
        logging.error(f"daisyblast failed during execution: {e}")
        raise