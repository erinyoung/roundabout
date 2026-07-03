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
