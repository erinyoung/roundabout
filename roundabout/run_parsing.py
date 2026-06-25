import csv
import logging
import re
from pathlib import Path
from collections import defaultdict

def parse_amrfinder_results(amr_dir: Path, substring: str = None) -> dict:
    """Scans AMRFinderPlus TSVs and groups samples by the AMR genes they carry."""
    logging.info("Parsing AMRFinderPlus results...")
    amr_groups = defaultdict(list)
    
    if substring:
        # A true, unconstrained substring match
        pattern = re.compile(re.escape(substring), re.IGNORECASE)
    
    for tsv_file in amr_dir.glob("*_amrfinder.tsv"):
        sample_name = tsv_file.name.replace("_amrfinder.tsv", "")
        
        with open(tsv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene = row.get('Element symbol') or row.get('Gene symbol') or row.get('Sequence name')
                el_type = str(row.get('Element type') or row.get('Type') or row.get('Class') or "").strip().lower()
                
                if gene and (el_type in ['amr', 'virulence', 'stress'] or el_type == ""):
                    if substring and not pattern.search(gene):
                        continue
                    
                    cleaned_gene = gene.strip()
                    if sample_name not in amr_groups[cleaned_gene]:
                        amr_groups[cleaned_gene].append(sample_name)
                        
    return dict(amr_groups)


def parse_plasmidfinder_results(pf_dir: Path, substring: str = None) -> dict:
    """Scans PlasmidFinder results and groups samples by their Inc groups."""
    logging.info("Parsing PlasmidFinder results...")
    inc_groups = defaultdict(list)
    
    if substring:
        # A true, unconstrained substring match for plasmid replicons
        pattern = re.compile(re.escape(substring), re.IGNORECASE)
    
    for sample_dir in pf_dir.iterdir():
        if not sample_dir.is_dir():
            continue
            
        sample_name = sample_dir.name
        results_tsv = sample_dir / "results_tab.tsv"
        
        if not results_tsv.exists():
            continue
            
        with open(results_tsv, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                inc_gene = row.get('Plasmid')
                if inc_gene:
                    if substring and not pattern.search(inc_gene):
                        continue
                        
                    if sample_name not in inc_groups[inc_gene]:
                        inc_groups[inc_gene].append(sample_name)
                        
    return dict(inc_groups)

def build_sentinel_groups(similarity_matrix: dict, input_sequences: list, refseq_sequences: list, min_identity: float, min_coverage: float) -> tuple[dict, dict, dict, dict]:
    """
    Evaluates sequence similarity to build ego-networks for each input sentinel.
    Categorizes relationships into Strict (symmetric) and Contained (asymmetric)
    across both Local (inputs) and Global (inputs + refseq) scopes.
    """
    local_strict = {}
    local_contained = {}
    global_strict = {}
    global_contained = {}
    
    for sentinel in input_sequences:
        ls_members = [sentinel]
        lc_members = [sentinel]
        
        # ---------------------------------------------------------
        # 1. Local Scope (Sentinel vs Input Sequences)
        # ---------------------------------------------------------
        for target in input_sequences:
            if target == sentinel:
                continue
                
            metrics = similarity_matrix.get(sentinel, {}).get(target, {})
            ani = metrics.get('ani', 0.0)
            q_cov = metrics.get('q_cov', 0.0)
            r_cov = metrics.get('r_cov', 0.0)
            
            if ani >= min_identity:
                # Metric 1: Containment (At least one sequence overlaps > threshold)
                if q_cov >= min_coverage or r_cov >= min_coverage:
                    lc_members.append(target)
                    
                # Metric 2: Strict (Both sequences overlap > threshold)
                if q_cov >= min_coverage and r_cov >= min_coverage:
                    ls_members.append(target)
                    
        local_strict[f"{sentinel}"] = ls_members
        local_contained[f"{sentinel}"] = lc_members
        
        # ---------------------------------------------------------
        # 2. Global Scope (Local Scope + RefSeq Sequences)
        # ---------------------------------------------------------
        gs_members = list(ls_members)
        gc_members = list(lc_members)
        
        for ref_target in refseq_sequences:
            metrics = similarity_matrix.get(sentinel, {}).get(ref_target, {})
            ani = metrics.get('ani', 0.0)
            q_cov = metrics.get('q_cov', 0.0)
            r_cov = metrics.get('r_cov', 0.0)
            
            if ani >= min_identity:
                if q_cov >= min_coverage or r_cov >= min_coverage:
                    gc_members.append(ref_target)
                if q_cov >= min_coverage and r_cov >= min_coverage:
                    gs_members.append(ref_target)
                    
        global_strict[f"{sentinel}"] = gs_members
        global_contained[f"{sentinel}"] = gc_members
        
    return local_strict, local_contained, global_strict, global_contained