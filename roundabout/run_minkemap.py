import os
import csv
import logging
import subprocess
from pathlib import Path

def generate_amr_highlights_from_tsv(reference_id: str, target_amr_gene: str, amr_outdir: Path, out_csv: Path) -> bool:
    """
    Parses the AMRFinderPlus TSV to extract AMR gene coordinates.
    Writes matches to a single MinkeMap-compatible highlights.csv file.
    """
    amr_tsv = amr_outdir / f"{reference_id}_amrfinder.tsv"
    
    if not amr_tsv.exists():
        logging.warning(f"AMRFinder TSV not found at {amr_tsv}. Cannot generate highlights.")
        return False
        
    found_target = False
    target = target_amr_gene.lower() if target_amr_gene else None
    
    try:
        with open(amr_tsv, 'r') as tsvfile, open(out_csv, 'w', newline='') as csvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            writer = csv.writer(csvfile)
            
            writer.writerow(["start", "end", "color", "label"])
            
            for row in reader:
                element_symbol = row.get('Element symbol', '').strip()
                element_name = row.get('Element name', '').strip()
                
                is_match = False
                if target:
                    if target in element_symbol.lower() or target in element_name.lower():
                        is_match = True
                else:
                    is_match = True 
                    
                if is_match:
                    start = row.get('Start')
                    end = row.get('Stop') 
                    label = element_symbol if element_symbol else element_name
                    
                    if start and end:
                        writer.writerow([start, end, "#FF0000", label])
                        found_target = True

        if not found_target:
            os.remove(out_csv)
            
    except Exception as e:
        logging.error(f"Failed to parse AMRFinder TSV for {reference_id}: {e}")
        return False
        
    return found_target


def execute_minkemap_run(central_fasta: Path, query_fastas: list, out_dir: Path, out_prefix: str, minkemap_opts: dict, highlight_csv: Path = None):
    """
    Executes a singular MinkeMap command given an absolute output flavor directory.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    output_image_name = f"{out_prefix}.png"

    cmd = [
        "minkemap",
        "-r", str(central_fasta),
        "-o", output_image_name,
        "--outdir", str(out_dir)
    ]
    
    cmd.append("-i")
    cmd.extend([str(q) for q in query_fastas])
    
    # Aesthetic mappings
    if minkemap_opts.get("palette"):
        cmd.extend(["--palette", minkemap_opts["palette"]])
    if minkemap_opts.get("gc_skew"):
        cmd.append("--gc-skew")
    if minkemap_opts.get("no_backbone"):
        cmd.append("--no-backbone")
        
    if highlight_csv and highlight_csv.exists():
        cmd.extend(["--highlights", str(highlight_csv)])

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"MinkeMap execution error for prefix {out_prefix}: {e.stderr}")


def run_minkemap_cohorts(
    pipeline_groups,
    fasta_map,
    amr_outdir,
    outdir,
    minkemap_opts,
    sample_to_ref_paths=None,
    target_amr_gene=None,
):
    """
    Cycles through groups and sets each local sequence as a central reference ring.
    Runs flavors across completely isolated directories.
    """

    minkemap_base_dir = Path(outdir) / "minkemap_results"
    
    for group_id, members in pipeline_groups.items():
        logging.info(f"Processing MinkeMap visualizations for {group_id} ({len(members)} members)...")
        
        for central_id in members:
            central_fasta = fasta_map.get(central_id)
            if not central_fasta or not os.path.exists(central_fasta):
                continue
                
            # Local queries (everything inside the cohort group except the center)
            local_queries = [fasta_map[q_id] for q_id in members if q_id != central_id and q_id in fasta_map]
            
            if not local_queries:
                continue

            # Establish safe isolated environment directories
            sample_root = minkemap_base_dir / group_id / central_id
            local_analysis_dir = sample_root / "local"
            refseq_analysis_dir = sample_root / "with_refseq"

            # -----------------------------------------------------------------
            # FLAVOR 1: Local Cohort Sequence Analyses Only
            # -----------------------------------------------------------------
            logging.info(f"Processing MinkeMap visualizations for {central_id} in {group_id}")
            # standard flavor
            execute_minkemap_run(
                central_fasta=central_fasta,
                query_fastas=local_queries,
                out_dir=local_analysis_dir,
                out_prefix=f"{group_id}_center_{central_id}_local_standard",
                minkemap_opts=minkemap_opts
            )
            
            # highlighted flavor
            local_highlight_csv = local_analysis_dir / "highlights.csv"
            if generate_amr_highlights_from_tsv(central_id, target_amr_gene, amr_outdir, local_highlight_csv):
                logging.info(f"Adding AMR highlights to MinkeMap visualizations for {central_id} in {group_id}")
                execute_minkemap_run(
                    central_fasta=central_fasta,
                    query_fastas=local_queries,
                    out_dir=local_analysis_dir,
                    out_prefix=f"{group_id}_center_{central_id}_local_highlighted",
                    minkemap_opts=minkemap_opts,
                    highlight_csv=local_highlight_csv
                )

            # -----------------------------------------------------------------
            # FLAVOR 2: Local Cohort Sequences + staged RefSeq references
            # -----------------------------------------------------------------
            staged_refs = []

            if sample_to_ref_paths:
                staged_refs = sample_to_ref_paths.get(central_id, [])

            if staged_refs:
                logging.info(
                    f"Processing MinkeMap visualizations for {central_id} "
                    f"in {group_id} with RefSeq context references"
                )

                combined_queries = local_queries + staged_refs

                execute_minkemap_run(
                    central_fasta=central_fasta,
                    query_fastas=combined_queries,
                    out_dir=refseq_analysis_dir,
                    out_prefix=f"{group_id}_center_{central_id}_refseq_standard",
                    minkemap_opts=minkemap_opts
                )
                    # This takes too long to run and is not necessary for the current analysis, so it has been commented out.
                    # # highlighted with refseq
                    # refseq_highlight_csv = refseq_analysis_dir / "highlights.csv"
                    # if generate_amr_highlights_from_tsv(central_id, target_amr_gene, amr_outdir, refseq_highlight_csv):
                    #     logging.info(f"Adding AMR highlights to MinkeMap visualizations for {central_id} in {group_id} with RefSeq context references")
                    #     execute_minkemap_run(
                    #         central_fasta=central_fasta,
                    #         query_fastas=combined_queries,
                    #         out_dir=refseq_analysis_dir,
                    #         out_prefix=f"{group_id}_center_{central_id}_refseq_highlighted",
                    #         minkemap_opts=minkemap_opts,
                    #         highlight_csv=refseq_highlight_csv
                    #     )
                        
    logging.info("MinkeMap analysis variants organized and completed successfully.")