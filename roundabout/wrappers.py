import subprocess
import os

def run_pathofetch(query, outdir):
    """
    Executes pathofetch to download plasmid FASTA files based on a query string.
    """
    cmd = [
        "pathofetch",
        "--query", query,
        "--outdir", outdir
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_pling(input_dir, outdir):
    """
    Executes pling to calculate containment/rearrangement distances 
    and group plasmids into outbreak clusters.
    """
    cmd = [
        "pling",
        "--input", input_dir,
        "--output", outdir
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    # Optional parsing logic: You would read pling's output files here
    # to return a structured cluster map back to engine.py if needed.
    return outdir


def run_bakta(fasta_path, sample_outdir, bakta_db):
    """
    Executes Bakta to generate comprehensive structural and functional annotations.
    Replaces the old prokka.nf workflow step.
    """
    cmd = [
        "bakta",
        "--db", bakta_db,
        "--output", sample_outdir,
        fasta_path
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_amrfinder(fasta_path, sample_outdir, amr_db):
    """
    Executes ncbi-amrfinderplus to locate antimicrobial resistance genes.
    Replaces amrfinder.nf.
    """
    # Create an explicit isolated folder inside the sample directory
    amr_dir = os.path.join(sample_outdir, "ncbi-AMRFinderplus")
    os.makedirs(amr_dir, exist_ok=True)
    
    sample_name = os.path.splitext(os.path.basename(fasta_path))[0]
    output_txt = os.path.join(amr_dir, f"{sample_name}_amrfinder.txt")

    cmd = [
        "amrfinder",
        "--nucleotide", fasta_path,
        "--database", amr_db,
        "--name", sample_name,
        "--output", output_txt
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_plasmidfinder(fasta_path, sample_outdir, pf_db):
    """
    Executes plasmidfinder.py to type and identify plasmid replicons.
    Replaces plasmidfinder.nf.
    """
    pf_dir = os.path.join(sample_outdir, "plasmidfinder")
    os.makedirs(pf_dir, exist_ok=True)

    cmd = [
        "plasmidfinder.py",
        "-i", fasta_path,
        "-o", pf_dir,
        "-p", pf_db
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_daisyblast(input_dir, outdir):
    """
    Executes daisyblast to run rapid all-to-all sequence similarity alignments.
    Replaces blastn.nf and the old python divider logic.
    """
    cmd = [
        "daisyblast",
        "--input", input_dir,
        "--output", outdir
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_heatcluster(daisyblast_dir, outdir):
    """
    Executes heatcluster to generate distance matrix heatmaps 
    from the daisyblast alignment data.
    """
    cmd = [
        "heatcluster",
        "--input", daisyblast_dir,
        "--output", outdir
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_minkemap(staged_fastas, pipeline_outdir, outdir):
    """
    Executes MinkeMap to produce final comparative circular plasmid layouts.
    Replaces the entire circos.nf processing loop.
    """
    os.makedirs(outdir, exist_ok=True)
    manifest_csv = os.path.join(outdir, "minkemap_manifest.csv")
    
    # Dynamic manifest creation: Tell MinkeMap where to look for raw data tracks
    with open(manifest_csv, "w") as f:
        f.write("sample,read1,read2,type\n")
        for fasta in staged_fastas:
            sample_name = os.path.splitext(os.path.basename(fasta))[0]
            f.write(f"{sample_name},{fasta},,fasta\n")

    # Pick a reference file from annotations to serve as the structural anchor ring
    # (e.g., pulling the first sample's Bakta GenBank output as an initial mapping target)
    samples_dir = os.path.join(pipeline_outdir, "samples")
    first_sample = os.listdir(samples_dir)[0]
    reference_gbk = os.path.join(samples_dir, first_sample, f"{first_sample}.gbk")

    cmd = [
        "minkemap",
        "-r", reference_gbk,
        "-f", manifest_csv,
        "--gc-skew",
        "--output-dir", outdir
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)