import argparse
import shutil
import subprocess
import importlib.metadata
import logging
import sys
import re
from roundabout.engine import run_pipeline

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

try:
    __version__ = importlib.metadata.version("roundabout")
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.1.0-dev"


def extract_version_from_text(text: str, name: str) -> str:
    """Isolates and strips down text to just the version footprint."""
    lines = [line.strip() for line in text.split("\n") if line.strip()]

    # Look for a line that contains numbers (the version footprint)
    for line in lines:
        # Skip lines that are just repeating the tool name/command
        if line.lower() == name.lower():
            continue

        # Common bio-info patterns: 'version 3.1', 'v0.1.0', 'blastn: 2.12.0+'
        # Look for words like version, v, or simply numbers separated by dots
        match = re.search(r"(?:version\s+|v)?\d+(?:\.\d+)+[-a-zA-Z0-9]*", line, re.IGNORECASE)
        if match:
            return match.group(0).strip()
            
        # Fallback: If line contains a digit and isn't just the tool name, grab it
        if any(char.isdigit() for char in line):
            return line

    return "Unknown Version"


def get_python_module_version(module_name: str, version_flag: str = "-v") -> str:
    """Runs a Python module and extracts its isolated version string."""
    try:
        cmd = [sys.executable, "-m", module_name, version_flag]
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output = result.stdout + result.stderr
        if not output.strip():
            return "Not found"
            
        # Clean up submodule prefix if present (e.g., pygenomeviz.gv_blast -> gv_blast)
        short_name = module_name.split(".")[-1]
        return extract_version_from_text(output, short_name)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "Not found"


def get_cli_version(cli_tool: str, version_flag: str = "-v") -> str:
    """Checks a PATH binary and extracts its isolated version string."""
    if shutil.which(cli_tool) is None:
        return "Not found"

    try:
        cmd = [cli_tool, version_flag]
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output = result.stdout + result.stderr
        if not output.strip():
            return "Not found"

        return extract_version_from_text(output, cli_tool)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "Not found"


def run_dependency_checks(check=True):
    """Checks the environment for required dependencies and their versions."""

    python_modules = {
        "plasmidfinder": "-v",
    }

    cli_binaries = {
        "daisyblast": "-v",
        # TODO : update heatcluster so that it accepts a pandas dataframe
        # "heatcluster": "-v",
        "minkemap": "-v",
        "refseq-plasmid-dl": "-v",
        "git": "--version",
        "kma": "-v",
        "bakta": "--version",
        "skani": "-V",
        "amrfinder": "-v",
        "minimap2": "--version",
        "blastn": "-version",
        "nucmer": "--version",
        "mmseqs": "",
        # TODO: check version of pygenomeviz
        #"pgv-mmseqs": "-v",
        #"pgv-mummer": "-v",
        #"progressiveMauve": "--version", 
        #"pgv-blast": "-v",
    }

    all_passed = True
    logging.info("Checking Tool Dependencies:")
    logging.info(f"[OK] roundabout ({__version__})")

    # Check Python Modules
    for mod, flag in python_modules.items():
        version = get_python_module_version(mod, flag)
        if version != "Not found":
            logging.info(f"[OK] {mod} ({version})")
        else:
            logging.warning(f"[FAIL] {mod} (Missing)")
            all_passed = False

    # Check Standard CLI Binaries
    for bin_name, flag in cli_binaries.items():
        version = get_cli_version(bin_name, flag)
        if version != "Not found":
            logging.info(f"[OK] {bin_name} ({version})")
        else:
            logging.warning(f"[FAIL] {bin_name} (Missing)")
            all_passed = False

    if all_passed:
        logging.info("Success: All dependencies are satisfied.")
        if check:
            sys.exit(0)
    else:
        logging.error("Error: Some dependencies are missing. Please install them before running.")
        sys.exit(1)


def main(args=None):

    cli_parser = argparse.ArgumentParser(
        description="Roundabout: Plasmid clustering and visualization.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # ---------------------------------------------------------
    # Standard Options
    # ---------------------------------------------------------
    cli_parser.add_argument('-f', '--fastas', type=str, help="Input directory containing plasmid fasta files")
    cli_parser.add_argument('-o', '--outdir', type=str, default="results", help="Pipeline output directory")
    cli_parser.add_argument('-t', '--threads', type=int, default=4, help="Number of concurrent threads to use")
    cli_parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}', help="Show the application version and exit")
    cli_parser.add_argument('-c', '--check', action='store_true', help="Show the application dependencies and exit")
    cli_parser.add_argument('--min-contig-length', type=int, default=0, help="Minimum sequence length to keep during staging")
    cli_parser.add_argument('--max-contig-length', type=int, default=1000000, help="Maximum sequence length to keep during staging (default: 1000000)")

    # ---------------------------------------------------------
    # Grouping Options
    # ---------------------------------------------------------
    grouping_parser = cli_parser.add_argument_group("Grouping Options")
    grouping_parser.add_argument('--min-ani', type=float, default=95.0, help="Minimum ANI for ANI-based grouping (default: 95.0)")
    grouping_parser.add_argument('--min-ani-align-fraction-ref', type=float, default=95.0, help="Minimum aligned fraction for ANI-based grouping (default: 95.0)")
    grouping_parser.add_argument('--min-ani-align-fraction-query', type=float, default=95.0, help="Minimum aligned fraction for ANI-based grouping (default: 95.0)")
    grouping_parser.add_argument('--num-ref', type=int, default=5, help="Number of reference sequences to use for identity-based grouping (default: 5)")

    # ---------------------------------------------------------
    # Database Options
    # ---------------------------------------------------------
    db_parser = cli_parser.add_argument_group("Database Options")
    db_parser.add_argument('-d', '--db-check', action='store_true', help="Check databases, then exit")
    db_parser.add_argument('--setup-db', action='store_true', help="Download/configure databases, then exit")
    db_parser.add_argument('--force', action='store_true', help="Force update of existing databases")
    db_parser.add_argument('--skip-db', action='store_true', help="Skip database setup and use existing databases if found and ignore if not.")

    # ---------------------------------------------------------
    # Bakta Options
    # ---------------------------------------------------------
    bakta_parser = cli_parser.add_argument_group("Bakta Options")
    bakta_parser.add_argument('--bakta-db', type=str, help="Specific path to Bakta database")
    bakta_parser.add_argument('--bakta-db-type', type=str, default="light", help="Type of Bakta database to download (e.g., light, full)")
    bakta_parser.add_argument('--bakta-options', type=str, help="Additional options for Bakta")
    
    # Bakta Organism
    bakta_parser.add_argument('--bakta-genus', type=str, help="Genus name")
    bakta_parser.add_argument('--bakta-species', type=str, help="Species name")
    bakta_parser.add_argument('--bakta-strain', type=str, help="Strain name")
    bakta_parser.add_argument('--bakta-plasmid', type=str, help="Plasmid name")
    
    # Bakta Annotation
    bakta_parser.add_argument('--bakta-complete', action='store_true', help="All sequences are complete replicons")
    bakta_parser.add_argument('--bakta-prodigal-tf', type=str, help="Path to existing Prodigal training file")
    bakta_parser.add_argument('--bakta-translation-table', type=int, choices=[11, 4, 25], default=11, help="Translation table (default = 11)")
    bakta_parser.add_argument('--bakta-gram', type=str, choices=['+', '-', '?'], default='?', help="Gram type for signal peptide predictions (default = ?)")
    bakta_parser.add_argument('--bakta-locus', type=str, default="contig", help="Locus prefix (default = 'contig')")
    bakta_parser.add_argument('--bakta-locus-tag', type=str, help="Locus tag prefix (default = autogenerated)")
    bakta_parser.add_argument('--bakta-locus-tag-increment', type=int, choices=[1, 5, 10], default=1, help="Locus tag increment (default = 1)")
    bakta_parser.add_argument('--bakta-keep-contig-headers', action='store_true', help="Keep original contig/sequence headers")
    bakta_parser.add_argument('--bakta-compliant', action='store_true', help="Force Genbank/ENA/DDJB compliance")
    bakta_parser.add_argument('--bakta-replicons', '-r', type=str, help="Replicon information table (tsv/csv)")
    bakta_parser.add_argument('--bakta-regions', type=str, help="Path to pre-annotated regions in GFF3 or Genbank format")
    bakta_parser.add_argument('--bakta-proteins', type=str, help="Fasta file of trusted protein sequences for CDS annotation")
    bakta_parser.add_argument('--bakta-hmms', type=str, help="HMM file of trusted hidden markov models in HMMER format")
    bakta_parser.add_argument('--bakta-meta', action='store_true', help="Run in metagenome mode (affects CDS prediction)")
    bakta_parser.add_argument('--bakta-partial', action='store_true', help="Predict partial (truncated) genes spanning linear sequence ends")
    
    # Bakta Workflow (Skips)
    bakta_parser.add_argument('--bakta-skip-trna', action='store_true', help="Skip tRNA detection & annotation")
    bakta_parser.add_argument('--bakta-skip-tmrna', action='store_true', help="Skip tmRNA detection & annotation")
    bakta_parser.add_argument('--bakta-skip-rrna', action='store_true', help="Skip rRNA detection & annotation")
    bakta_parser.add_argument('--bakta-skip-ncrna', action='store_true', help="Skip ncRNA detection & annotation")
    bakta_parser.add_argument('--bakta-skip-ncrna-region', action='store_true', help="Skip ncRNA region detection & annotation")
    bakta_parser.add_argument('--bakta-skip-crispr', action='store_true', help="Skip CRISPR array detection & annotation")
    bakta_parser.add_argument('--bakta-skip-cds', action='store_true', help="Skip CDS detection & annotation")
    bakta_parser.add_argument('--bakta-skip-pseudo', action='store_true', help="Skip pseudogene detection & annotation")
    bakta_parser.add_argument('--bakta-skip-gap', action='store_true', help="Skip gap detection & annotation")
    bakta_parser.add_argument('--bakta-skip-ori', action='store_true', help="Skip oriC/oriT detection & annotation")
    bakta_parser.add_argument('--bakta-skip-filter', action='store_true', help="Skip feature overlap filters")
    bakta_parser.add_argument('--bakta-skip-plot', action='store_true', help="Skip generation of circular genome plots")

    # ---------------------------------------------------------
    # NCBI AMRFinderPlus Options
    # ---------------------------------------------------------
    amr_parser = cli_parser.add_argument_group("NCBI AMRFinderPlus Options")
    amr_parser.add_argument('--amrfinder-db', type=str, help="Specific path to AMRFinderPlus database")
    amr_parser.add_argument('--amr-organism', type=str, help="Organism for AMRFinderPlus")
    amr_parser.add_argument('--amr-gene', type=str, help="String used to filter AMRFinderPlus results by gene name (case-insensitive)")
    amr_parser.add_argument('--amr-ident-min', type=float, default=0.9, help="Minimum proportion of identical amino acids in alignment for hit (0..1). -1 uses curated threshold.")
    amr_parser.add_argument('--amr-coverage-min', type=float, default=0.5, help="Minimum coverage of the reference protein (0..1)")
    amr_parser.add_argument('--amrfinder-options', type=str, help="Additional options for AMRFinderPlus")

    # ---------------------------------------------------------
    # PlasmidFinder Options
    # ---------------------------------------------------------
    plasmidfinder_parser = cli_parser.add_argument_group("PlasmidFinder Options")
    plasmidfinder_parser.add_argument('--plasmidfinder-db', type=str, help="Specific path to PlasmidFinder database")
    plasmidfinder_parser.add_argument('--plasmidfinder-string', type=str, help="Filter PlasmidFinder results by string (case-insensitive)")
    plasmidfinder_parser.add_argument('--plasmidfinder-mincov', type=float, help="Minimum coverage")
    plasmidfinder_parser.add_argument('--plasmidfinder-threshold', type=float, help="Minimum threshold for identity")

    # ---------------------------------------------------------
    # refseq-plasmid-dl Options
    # ---------------------------------------------------------
    refseq_parser = cli_parser.add_argument_group("refseq-plasmid-dl Options")
    refseq_parser.add_argument('--refseq-plasmid-dl-db', type=str, help="Specific path to refseq-plasmid-dl database")
    refseq_parser.add_argument('--refseq-organism', type=str, help="Filter by species/organism (substring match)")
    refseq_parser.add_argument('--refseq-taxid', type=str, help="Filter by NCBI Taxonomy ID")
    refseq_parser.add_argument('--refseq-strain', type=str, help="Filter by strain")
    refseq_parser.add_argument('--refseq-isolate', type=str, help="Filter by isolate")
    refseq_parser.add_argument('--refseq-host', type=str, help="Filter by host organism")
    refseq_parser.add_argument('--refseq-plasmid-name', type=str, help="Filter by plasmid name")
    refseq_parser.add_argument('--refseq-geo-loc-name', type=str, help="Filter by geographic location")
    refseq_parser.add_argument('--refseq-isolation-source', type=str, help="Filter by isolation source")
    refseq_parser.add_argument('--refseq-min-length', type=int, help="Minimum sequence length (bp)")
    refseq_parser.add_argument('--refseq-max-length', type=int, help="Maximum sequence length (bp)")
    refseq_parser.add_argument('--refseq-topology', type=str, choices=['circular', 'linear', 'all'], help="Filter by topology")
    refseq_parser.add_argument('--refseq-min-date', type=str, help="Include only records updated after YYYY-MM-DD")
    refseq_parser.add_argument('--refseq-max-date', type=str, help="Include only records updated before YYYY-MM-DD")
    refseq_parser.add_argument('--refseq-min-collection-date', type=str, help="Include only records collected after YYYY-MM-DD")
    refseq_parser.add_argument('--refseq-max-collection-date', type=str, help="Include only records collected before YYYY-MM-DD")

    # ---------------------------------------------------------
    # Skani Options
    # ---------------------------------------------------------
    skani_parser = cli_parser.add_argument_group("Skani Options")
    
    # Output Filters & Formatting
    skani_parser.add_argument('--skani-min-af', type=float, default=15.0, help="Only output ANI values where one genome has aligned fraction > this value (default: 15)")
    skani_parser.add_argument('--skani-both-min-af', type=float, help="Only output ANI values where both genomes have aligned fraction > this value")
    skani_parser.add_argument('--skani-ci', action='store_true', help="Output [5%%,95%%] ANI confidence intervals using percentile bootstrap")
    skani_parser.add_argument('--skani-detailed', action='store_true', help="Print additional info including contig N50s and more")
    skani_parser.add_argument('--skani-n', type=int, help="Max number of results to show for each query (default: unlimited)")
    skani_parser.add_argument('--skani-short-header', action='store_true', help="Only display the first part of contig names (before first whitespace)")

    # Presets
    skani_parser.add_argument('--skani-fast', action='store_true', help="Faster skani mode (Alias for -c 200)")
    skani_parser.add_argument('--skani-medium', action='store_true', help="Medium skani mode (Alias for -c 70)")
    skani_parser.add_argument('--skani-slow', action='store_true', help="Slower skani mode (Alias for -c 30)")
    skani_parser.add_argument('--skani-small-genomes', action='store_true', help="Mode for small genomes such as viruses or plasmids (Alias for: -c 30 -m 200 --faster-small)")

    # Algorithm Parameters
    skani_parser.add_argument('--skani-c', type=int, default=125, help="Compression factor / k-mer subsampling rate (default: 125)")
    skani_parser.add_argument('--skani-faster-small', action='store_true', help="Filter genomes with < 20 marker k-mers more aggressively")
    skani_parser.add_argument('--skani-m', type=int, default=1000, help="Marker k-mer compression factor (default: 1000)")
    skani_parser.add_argument('--skani-median', action='store_true', help="Estimate median identity instead of average (mean) identity")
    skani_parser.add_argument('--skani-no-learned-ani', action='store_true', help="Disable regression model for ANI prediction")
    skani_parser.add_argument('--skani-no-marker-index', action='store_true', help="Do not use hash-table inverted index for faster ANI filtering")
    skani_parser.add_argument('--skani-robust', action='store_true', help="Estimate mean after trimming off 10%%/90%% quantiles")
    skani_parser.add_argument('--skani-s', type=float, default=80.0, help="Screen out pairs with approximately < %% identity using k-mer sketching (default: 80)")

    # ---------------------------------------------------------
    # Daisyblast Options
    # ---------------------------------------------------------
    daisyblast_parser = cli_parser.add_argument_group("Daisyblast Options")
    daisyblast_parser.add_argument('--daisyblast-evalue', type=float, default=1e-10, help="E-value threshold for Daisyblast (default: 1e-10)")
    daisyblast_parser.add_argument('--daisyblast-min_pident', type=float, default=90.0, help="Minimum percent identity for Daisyblast (default: 90.0)")
    daisyblast_parser.add_argument('--daisyblast-min_length', type=int, default=200, help="Minimum length for Daisyblast (default: 200)")
    daisyblast_parser.add_argument('--daisyblast-num_groups', type=int, default=20, help="Number of groups for Daisyblast (default: 20)")

    # ---------------------------------------------------------
    # MinkeMap Options
    # ---------------------------------------------------------
    minkemap_parser = cli_parser.add_argument_group("MinkeMap Options")
    minkemap_parser.add_argument('--minkemap-palette', type=str, help="Color scheme presets (whale, viridis, etc.) or comma-separated hex codes")
    minkemap_parser.add_argument('--minkemap-track-width', type=float, default=6.0, help="Width of each ring (default: 6)")
    minkemap_parser.add_argument('--minkemap-track-gap', type=float, default=4.0, help="Gap between rings (default: 4)")
    minkemap_parser.add_argument('--minkemap-dpi', type=int, default=300, help="Resolution for images (default: 300)")
    minkemap_parser.add_argument('--minkemap-no-backbone', action='store_true', help="Hides the central black reference axis")
    minkemap_parser.add_argument('--minkemap-no-legend', action='store_true', help="Hides the legends")
    minkemap_parser.add_argument('--minkemap-label-size', type=float, default=6.0, help="Font size for gene labels (default: 6)")
    minkemap_parser.add_argument('--minkemap-title', type=str, help="Add a title to the top of the plot")
    minkemap_parser.add_argument('--minkemap-gc-skew', action='store_true', help="Adds a GC Skew track (Blue+/Orange-) in the center")
    minkemap_parser.add_argument('--minkemap-annotations', type=str, help="Path to custom annotations CSV")
    minkemap_parser.add_argument('--minkemap-highlights', type=str, help="Path to background highlights CSV")
    minkemap_parser.add_argument('--minkemap-exclude-genes', type=str, help="Comma-separated list of keywords to hide from the gene track")
    minkemap_parser.add_argument('--minkemap-min-identity', type=float, help="Minimum %% identity (0-100) to display an alignment block")
    minkemap_parser.add_argument('--minkemap-min-coverage', type=float, help="Minimum %% query coverage (0-100) required to include a read/contig")

    # ---------------------------------------------------------
    # Heatcluster Options
    # ---------------------------------------------------------
    # TODO : add back in when heatcluster is updated to accept a pandas dataframe instead of a file path
    # heatcluster_parser = cli_parser.add_argument_group("Heatcluster Options")
    # heatcluster_parser.add_argument('--heatcluster-out', type=str, help="Output filename for cluster assignments (e.g., clusters.csv)")
    # heatcluster_parser.add_argument('--heatcluster-k', type=int, help="Split tree into K groups. Overrides --cluster-t")
    # heatcluster_parser.add_argument('--heatcluster-t', type=float, help="Split tree by distance threshold (e.g., 10 SNPs or 0.05 ANI)")
    # heatcluster_parser.add_argument('--heatcluster-auto-k', action='store_true', help="Automatically detect optimal K using Silhouette Scores")
    # heatcluster_parser.add_argument('--heatcluster-pca', action='store_true', help="Generate a PCA scatter plot")
    # heatcluster_parser.add_argument('--heatcluster-pca-out', type=str, help="Filename for PCA plot")
    # heatcluster_parser.add_argument('--heatcluster-title', type=str, help="Title of the plot")
    # heatcluster_parser.add_argument('--heatcluster-cmap', type=str, help="Matplotlib colormap (e.g., Reds_r, viridis)")
    # heatcluster_parser.add_argument('--heatcluster-dpi', type=int, default=300, help="DPI for output image (default: 300)")
    # heatcluster_parser.add_argument('--heatcluster-no-annot', action='store_true', help="Do not show numbers inside the heatmap cells")
    # heatcluster_parser.add_argument('--heatcluster-no-plot', action='store_true', help="Skip generating the heatmap image (Computation only)")
    # heatcluster_parser.add_argument('--heatcluster-width', type=float, help="Force figure width in inches")
    # heatcluster_parser.add_argument('--heatcluster-height', type=float, help="Force figure height in inches")
    # heatcluster_parser.add_argument('--heatcluster-font-scale', type=float, help="Scale all font sizes by this factor")
    # heatcluster_parser.add_argument('--heatcluster-vmin', type=float, help="Minimum value for color scale")
    # heatcluster_parser.add_argument('--heatcluster-vmax', type=float, help="Maximum value for color scale")
    # heatcluster_parser.add_argument('--heatcluster-hide-below', type=float, help="Mask values lower than this")
    # heatcluster_parser.add_argument('--heatcluster-hide-above', type=float, help="Mask values higher than this")
    # heatcluster_parser.add_argument('--heatcluster-no-cluster', action='store_true', help="Disable hierarchical clustering")
    # heatcluster_parser.add_argument('--heatcluster-dendrogram', action='store_true', help="Show the dendrogram (tree) and borders")

    # ---------------------------------------------------------
    # PyGenomeViz Options
    # ---------------------------------------------------------
    pgv_parser = cli_parser.add_argument_group("PyGenomeViz Options")
    pgv_parser.add_argument('--pgv-skip-blast', action='store_true', help="Skip PyGenomeViz alignment/visualization using BLAST")
    pgv_parser.add_argument('--pgv-skip-mummer', action='store_true', help="Skip PyGenomeViz alignment/visualization using MUMmer")
    pgv_parser.add_argument('--pgv-skip-mmseqs', action='store_true', help="Skip PyGenomeViz alignment/visualization using MMseqs")
    pgv_parser.add_argument('--pgv-skip-pmauve', action='store_true', help="Skip PyGenomeViz alignment/visualization using progressiveMauve")
    pgv_parser.add_argument('--pgv-min-identity', type=float, default=30.0, help="Minimum identity threshold for alignment rendering (default: 30)")
    pgv_parser.add_argument('--pgv-min-length', type=int, default=100, help="Minimum length threshold for alignment rendering (default: 100)")

    args = cli_parser.parse_args(args)

    run_dependency_checks(args.check)

    if not args.fastas and not args.check and not args.db_check and not args.setup_db:
        logging.error("Error: No directory with input FASTA files provided.")
        cli_parser.print_help()
        sys.exit(1)

    # Only run pipeline if we aren't just doing DB setup/checks
    if args.setup_db or args.db_check:
        logging.info("Database operations completed.")
        sys.exit(0)

    run_pipeline(args)


if __name__ == "__main__":
    main()