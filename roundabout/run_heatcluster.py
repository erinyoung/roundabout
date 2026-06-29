import subprocess
import logging
from pathlib import Path

# -----------------------------------------------------------------------------
# Heatcluster
# -----------------------------------------------------------------------------
def run_heatcluster(args):
    """
    Runs heatcluster to generate an ANI heatmap from a skani TSV matrix.
    """

    skani_file = Path(f"{args.outdir}/skani_results/skani_matrix.tsv")
    output_plot = Path(f"{args.outdir}/skani_plot.png")
    
    cmd = [
        "heatcluster",
        "-i", skani_file,
        "--format", "skani",
        "--cmap", "viridis",
        "-o", output_plot
    ]
    

    if skani_file.is_file():
        logging.info("Generating heatmap for skani output using heatcluster...")
    
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.info(f"Heatmap successfully saved to {output_plot}")
            return output_plot
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Heatcluster crashed on {skani_file.name}. Heatcluster Error:\n{e.stderr}\n{e.stdout}")
            raise