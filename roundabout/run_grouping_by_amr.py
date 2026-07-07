import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)


def define_groups_by_amr(
    amr_dict: dict,
    outdir: Path,
    target_gene: str = None,
) -> list[list[str]]:
    """
    Groups sequences by identical AMR profiles.

    If target_gene is provided, only matching genes are used
    when defining the grouping profile. Sequences with no matching
    genes (or no AMR genes at all) are excluded from the groups.
    """

    if target_gene:
        logging.info(f"Grouping by AMR genes matching '{target_gene}'...")
        search = target_gene.lower()
    else:
        logging.info("Grouping by complete AMR profiles...")

    profile_to_group = {}

    for seq, genes in amr_dict.items():

        if target_gene:
            profile = tuple(sorted(
                g for g in genes
                if search in g.lower()
            ))
        else:
            profile = tuple(sorted(genes))

        # Skip sequences that don't have any matching genes
        if not profile:
            continue

        profile_to_group.setdefault(profile, []).append(seq)

    # Standard list of lists, sorted by size then alphabetically
    groups = [sorted(group) for group in profile_to_group.values()]
    groups.sort(key=lambda g: (len(g), g), reverse=True)

    # ---------------------------------------------------------
    # Build detailed metadata for the JSON output
    # ---------------------------------------------------------
    json_payload = []
    
    for i, group_seqs in enumerate(groups, start=1):
        
        # Pull original full gene profiles for every sequence in this group
        group_gene_sets = [set(amr_dict.get(seq, [])) for seq in group_seqs]
        
        if group_gene_sets:
            shared_genes = sorted(set.intersection(*group_gene_sets))
            all_genes = sorted(set.union(*group_gene_sets))
        else:
            shared_genes = []
            all_genes = []

        # Re-derive the filtered profile that caused this group to form
        first_seq = group_seqs[0]
        if target_gene:
            filtered_profile = sorted([
                g for g in amr_dict.get(first_seq, []) 
                if search in g.lower()
            ])
        else:
            filtered_profile = None

        json_payload.append({
            "group_id": f"amr_group_{i:04d}",
            "sequences": group_seqs,
            "filtered_amr_profile": filtered_profile,
            "shared_amr_genes": shared_genes,
            "all_amr_genes": all_genes
        })

    # Ensure the subdirectory exists
    outfile = outdir / "amrfinder_results" / "amr_groups.json"
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w") as f:
        json.dump(json_payload, f, indent=2)

    logging.info(f"Wrote {len(groups)} AMR groups to {outfile}")

    # Return the clean list of lists for downstream processing
    return groups