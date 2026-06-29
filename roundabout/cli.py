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
    """ Checks the environment for required dependencies and their versions."""

    python_modules = {
        "plasmidfinder": "-v",
    }

    cli_binaries = {
        "daisyblast": "-v",
        # TODO: add heatcluster
        #"heatcluster": "-v",
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
        "pgv-mmseqs": "-v",
        "pgv-mummer": "-v",
        # TODO: add progressiveMauve
        "pgv-blast": "-v",
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

    # Standard flags
    cli_parser = argparse.ArgumentParser(
        description="Roundabout: Plasmid clustering and visualization.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    cli_parser.add_argument(
        '-f', '--fastas', 
        type=str, 
        help="Input directory containing plasmid fasta files"
    )

    cli_parser.add_argument(
        '-o', '--outdir', 
        type=str, 
        default="results", 
        help="Pipeline output directory"
    )

    cli_parser.add_argument(
        '-t', '--threads', 
        type=int, 
        default=4, 
        help="Number of concurrent threads to use"
    )

    cli_parser.add_argument(
        '-v', '--version', 
        action='version', 
        version=f'%(prog)s {__version__}',
        help="Show the application version and exit"
    )

    cli_parser.add_argument(
        '-c', '--check', 
        action='store_true', 
        help="Show the application dependencies and exit"
    )

    # ---------------------------------------------------------
    # Grouping Options
    # ---------------------------------------------------------
    grouping_parser = cli_parser.add_argument_group("Grouping Options")
    grouping_parser.add_argument(
        '--min-identity',
        type=float,
        default=95.0,
        help="Minimum sequence identity percentage for identity-based grouping (default: 95.0)"
    )
    grouping_parser.add_argument(
        '--min-coverage',
        type=float,
        default=80.0,
        help="Minimum alignment fraction (coverage) for identity-based grouping (default: 80.0)"
    )

    db_parser = cli_parser.add_argument_group("Database Options")
    db_parser.add_argument(
        '-d', '--db-check', 
        action='store_true', 
        help="Check databases, then exit"
    )
    db_parser.add_argument(
        '--setup-db', 
        action='store_true', 
        help="Download/configure databases, then exit"
    )
    db_parser.add_argument(
        '--force', 
        action='store_true', 
        help="Force update of existing databases"
    )
    db_parser.add_argument(
        '--skip-db', 
        action='store_true', 
        help="Skip database setup and use existing databases if found and ignore if not."
    )
    
    # Using Bakta
    bakta_parser = cli_parser.add_argument_group("Bakta Options")
    bakta_parser.add_argument(
        '--bakta-db', 
        type=str,
        help="Specific path to Bakta database"
    )
    bakta_parser.add_argument(
        '--bakta-db-type', 
        type=str,
        default="light",
        help="If using roundabout to download Bakta database, specify the type (e.g., light, full)"
    )
    bakta_parser.add_argument(
        '--bakta-options', 
        type=str, 
        help="Additional options for Bakta"
    )

    # Using AMRFinderPlus
    amr_parser = cli_parser.add_argument_group("NCBI AMRFinderPlus Options")
    amr_parser.add_argument(
        '--amrfinder-db', 
        type=str, 
        help="Specific path to AMRFinderPlus database"
    )
    amr_parser.add_argument(
        '--amr-organism',
        type=str,
        help="Organism for AMRFinderPlus"
    )
    amr_parser.add_argument(
        '--amr-gene',
        type=str,
        help="String used to filter AMRFinderPlus results by gene name (case-insensitive)"
    )
    amr_parser.add_argument(
        '--amrfinder-options', 
        type=str, 
        help="Additional options for AMRFinderPlus"
    )

    # Using PlasmidFinder
    plasmidfinder_parser = cli_parser.add_argument_group("PlasmidFinder Options")
    plasmidfinder_parser.add_argument(
        '--plasmidfinder-db', 
        type=str, 
        help="Specific path to PlasmidFinder database"
    )
    plasmidfinder_parser.add_argument(
        '--plasmidfinder-string',
        type=str,
        help="String used to filter PlasmidFinder results by string (case-insensitive)"
    )
    plasmidfinder_parser.add_argument(
        '--plasmidfinder-options', 
        type=str, 
        help="Additional options for PlasmidFinder"
    )

    # Using refseq-plasmid-dl
    refseq_parser = cli_parser.add_argument_group("refseq-plasmid-dl Options")
    refseq_parser.add_argument(
        '--refseq-plasmid-dl-db', 
        type=str, 
        help="Specific path to refseq-plasmid-dl database"
    )
    refseq_parser.add_argument(
        '--query', 
        type=str, 
        help="Query string to pull down plasmids using pathofetch"
    )

    args = cli_parser.parse_args(args)

    run_dependency_checks(args.check)

    if not args.fastas:
        logging.error("Error: No directory with input FASTA files provided.")
        sys.exit(1)

    run_pipeline(args)


if __name__ == "__main__":
    main()