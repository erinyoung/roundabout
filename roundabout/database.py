import os
import sys
import shutil
import subprocess
import logging
import json

def check_amrfinder_db(db_dir: str):
    """Checks if the AMRFinderPlus database exists and is properly configured."""

    logging.info("Evaluating AMRFinderPlus database")

    cmd = ["amrfinder", "--database_version"]

    if db_dir:
        cmd.extend(["--database", db_dir])

    try:
        # Run the command and capture the output streams
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        output = result.stdout + result.stderr
        
        db_version = "Not found"
        
        # Parse the output line by line
        for line in output.splitlines():
            line = line.strip()
            if line.startswith("Database directory:"):
                # Split at the first colon, strip whitespace, and remove the single quotes
                db_dir = line.split(":", 1)[1].strip().strip("'")
            elif line.startswith("Database version:"):
                # Split at the first colon and strip whitespace
                db_version = line.split(":", 1)[1].strip()

        logging.info(f"AMRFinderPlus database directory: {db_dir}")
        logging.info(f"AMRFinderPlus database version: {db_version}")

                
        return db_dir

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logging.error(f"Failed to run amrfinder --database_version: {e}")
        return None

def check_bakta_db(db_dir: str):
    """Checks if the Bakta database exists and is properly configured."""
    if not db_dir:
        return None

    logging.info(f"Evaluating Bakta database at: {db_dir}")
    bakta_db_file = os.path.join(db_dir, "bakta.db")
    version_file = os.path.join(db_dir, "version.json")
    
    if not os.path.exists(bakta_db_file) or not os.path.exists(version_file):
        logging.error("Bakta database check failed: Missing required files.")
        return None

    try:
        with open(version_file, 'r') as f:
            data = json.load(f)
            logging.info(f"Success: Found Bakta '{data.get('type', 'unknown')}' database v{data.get('major', '?')}.{data.get('minor', '?')}")
        return db_dir
    except Exception as e:
        logging.error(f"Failed to read Bakta version.json: {e}")
        return None


def check_plasmidfinder_db(db_dir: str):
    """Checks if the PlasmidFinder database exists and is indexed."""
    if not db_dir:
        return None

    logging.info(f"Evaluating PlasmidFinder database at: {db_dir}")
    version_file = os.path.join(db_dir, "VERSION")
    config_file = os.path.join(db_dir, "config")
    
    if not os.path.exists(db_dir) or not os.path.exists(config_file) or not os.path.exists(version_file):
        logging.error(f"Missing core PlasmidFinder files in: {db_dir}")
        return None
        
    kma_indexed = any(f.endswith(".seq.b") for f in os.listdir(db_dir))
    if not kma_indexed:
        logging.error(f"PlasmidFinder database at {db_dir} is not indexed with KMA.")
        return None

    try:
        with open(version_file, 'r') as f:
            logging.info(f"Success: Found PlasmidFinder database v{f.read().strip()}")
        return db_dir
    except Exception as e:
        logging.error(f"Failed to read PlasmidFinder VERSION: {e}")
        return None


def check_refseq_plasmid_dl_db(db_dir: str):
    """Checks if the refseq plasmid database exists."""
    if not db_dir:
        return None

    logging.info(f"Evaluating refseq Plasmid database at: {db_dir}")
    final_fasta = os.path.join(db_dir, "refseq_plasmids_dl.fasta")
    metadata_csv = os.path.join(db_dir, "refseq_plasmids_dl_metadata.csv")
    report_csv = os.path.join(db_dir, "refseq_plasmids_dl_report.csv")
    
    missing_files = [os.path.basename(f) for f in [final_fasta, metadata_csv, report_csv] if not os.path.exists(f)]
    if missing_files:
        logging.error(f"Missing refseq database files in {db_dir}: {', '.join(missing_files)}")
        return None

    try:
        with open(report_csv, 'r') as f:
            for line in f:
                if line.startswith("Sequences Written to FASTA"):
                    seq_count = line.strip().split(',')[1].strip()
                    logging.info(f"Success: Found refseq Plasmid database ({seq_count} sequences)")
                    break
        return db_dir
    except Exception as e:
        logging.error(f"Failed to read refseq report: {e}")
        return None

def setup_amrfinder_database(target_dir: str = "", force: bool = False):
    """Download the AMRFinderPlus database."""

    logging.info(f"Downloading AMRFinderPlus database to {target_dir or 'default location'}")
    
    cmd = ["amrfinder", "-u"]
    
    # Append flags based on conditions
    if force:
        cmd.append("--force_update")
        
    if target_dir:
        cmd.extend(["--database", target_dir])
        
    # Execute the final command list
    subprocess.run(cmd, check=True)

def setup_bakta_database(target_dir: str, db_type: str = "light", force: bool = False):
    """Downloads and configures the Bakta database."""
    
    if target_dir:
        db_dir = target_dir
    else:
        if db_type == "light":
            db_dir = os.path.join(os.getcwd(), "db-light")
        else: # defaults to full
            db_dir = os.path.join(os.getcwd(), "db")

    logging.info(f"Downloading Bakta database to db_dir: {db_dir} with type: {db_type}")
    
    bakta_db_file = os.path.join(db_dir, "bakta.db")
    version_file = os.path.join(db_dir, "version.json")

    # Check if the files exist on the system
    if os.path.exists(bakta_db_file) and os.path.exists(version_file):
        logging.info("Bakta database already exists.")

        if not force:
            logging.info("Skipping download since --force flag is not set.")
            return db_dir
        else:
            logging.info("Force flag detected. Removing existing Bakta DB for update.")
            # Remove db_dir, not target_dir (which might be None)
            shutil.rmtree(db_dir, ignore_errors=True)
    
    # Start with the base command
    cmd = ["bakta_db", "download", "--type", db_type]
        
    if target_dir:
        cmd.extend(["--database", target_dir])
        
    # Execute the final command list
    subprocess.run(cmd, check=True)

    return db_dir


def setup_plasmidfinder_database(target_dir, force: bool = False):
    """Downloads and indexes the PlasmidFinder database via Git and KMA."""
    
    logging.info(f"Setting up PlasmidFinder database at {target_dir or 'default location'}")

    # Resolve the target directory
    db_dir = target_dir if target_dir else os.path.join(os.getcwd(), "plasmidfinder_db")
    
    # Check if the database already exists
    if os.path.exists(db_dir) and os.listdir(db_dir):
        logging.info("PlasmidFinder database directory already exists.")
        
        if not force:
            logging.info("Skipping download and indexing since --force flag is not set.")
            return db_dir
        else:
            logging.info("Force flag detected. Removing existing PlasmidFinder DB for update.")
            shutil.rmtree(db_dir, ignore_errors=True)
            
    # Define the Git clone command
    repo_url = "https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git"
    clone_cmd = ["git", "clone", repo_url, db_dir]
    
    # Execute download and indexing
    try:
        # Step A: Clone the database
        logging.info(f"Cloning repository: {repo_url}")
        subprocess.run(clone_cmd, check=True, capture_output=True, text=True)
        logging.info("PlasmidFinder database downloaded successfully.")
        
        # Step B: Index the database using the provided INSTALL.py script
        logging.info("Indexing PlasmidFinder database using KMA...")
        index_cmd = ["python3", "INSTALL.py"]
        
        # The cwd=db_dir argument ensures the script runs *inside* the cloned folder
        subprocess.run(index_cmd, cwd=db_dir, check=True, capture_output=True, text=True)
        logging.info("PlasmidFinder database indexed successfully.")
        
    except subprocess.CalledProcessError as e:
        # e.stderr captures the specific command line error output from Git or Python/KMA
        logging.error(f"Failed to setup PlasmidFinder database. Error: {e.stderr}")
        raise
        
    return db_dir


def setup_refseq_plasmid_dl_database(target_dir: str, force: bool = False) -> str:
    """Downloads the refseq plasmid database using refseq-plasmid-dl."""
    
    logging.info(f"Setting up refseq Plasmid database at {target_dir or 'default location'}")

    # Resolve the target directory
    db_dir = target_dir if target_dir else os.path.join(os.getcwd(), "refseq_plasmids_db")
    final_fasta = os.path.join(db_dir, "refseq_plasmids_dl.fasta")
    
    # Check if the database already exists
    if os.path.exists(db_dir) and os.path.exists(final_fasta):
        logging.info("refseq Plasmid database already exists.")
        
        if not force:
            logging.info("Skipping download since --force flag is not set.")
            return db_dir
        else:
            logging.info("Force flag detected. Removing existing refseq DB for update.")
            shutil.rmtree(db_dir, ignore_errors=True)
            
    cmd = ["refseq-plasmid-dl", "--outdir", db_dir, "--topology", "circular"]

    if force:
        cmd.append("--force")
    
    # Execute download
    try:
        logging.info("Running refseq-plasmid-dl (this may take a while depending on NCBI servers)...")
        # I did not capture_output here so the user can see your tool's native progress/logging
        subprocess.run(cmd, check=True) 
        logging.info("refseq Plasmid database downloaded and curated successfully.")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to setup refseq Plasmid database. Error code: {e.returncode}")
        raise

    return db_dir


def run_setup(args):
    """Main handler for setting up databases."""

    # AMRFinder
    if not args.skip_db:
        logging.info("Downloading AMRFinderPlus database")
        setup_amrfinder_database(args.amrfinder_db, args.force)
    amr_db = check_amrfinder_db(args.amrfinder_db)

    # Bakta
    bakta_db = args.bakta_db
    if not args.skip_db:
        logging.info("Downloading Bakta database")
        bakta_db = setup_bakta_database(args.bakta_db, args.bakta_db_type, args.force)
    bakta_db = check_bakta_db(bakta_db)
    
    # PlasmidFinder
    plasmidfinder_db = args.plasmidfinder_db
    if not args.skip_db:
        logging.info("Downloading PlasmidFinder database")
        plasmidfinder_db = setup_plasmidfinder_database(args.plasmidfinder_db, args.force)
    plasmidfinder_db = check_plasmidfinder_db(plasmidfinder_db)
    
    # refseq (Fixed Typo Here)
    refseq_plasmid_dl_db = args.refseq_plasmid_dl_db 
    if not args.skip_db:
        logging.info("Downloading refseq-plasmid-dl database")
        refseq_plasmid_dl_db = setup_refseq_plasmid_dl_database(args.refseq_plasmid_dl_db, args.force)
    refseq_plasmid_dl_db = check_refseq_plasmid_dl_db(refseq_plasmid_dl_db)

    # Final Validation
    if not args.skip_db and not all([amr_db, bakta_db, plasmidfinder_db, refseq_plasmid_dl_db]):
        logging.error("One or more databases failed to set up properly. Please check the logs for details.")
        sys.exit(1)

    if args.db_check:
        logging.info("Database check complete. Exiting as per --db-check flag.")
        sys.exit(0)

    db_paths = {
        "amrfinder": amr_db,
        "bakta": bakta_db,
        "plasmidfinder": plasmidfinder_db,
        "refseq_plasmid_dl": refseq_plasmid_dl_db
    }

    return db_paths