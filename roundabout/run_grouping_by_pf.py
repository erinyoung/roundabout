import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)


def define_groups_by_plasmidfinder(
    pf_dict: dict,
    outdir: Path,
    target_replicon: str = None,
) -> list[list[str]]:
    """
    Groups sequences by identical PlasmidFinder replicon profiles.

    If target_replicon is provided, only matching replicons are used
    when defining the grouping profile. Sequences with no matching
    replicons are excluded from the groups.
    """

    if target_replicon:
        logging.info(f"Grouping by replicons matching '{target_replicon}'...")
        search = target_replicon.lower()
    else:
        logging.info("Grouping by complete replicon profiles...")

    profile_to_group = {}

    for seq, replicons in pf_dict.items():

        if target_replicon:
            profile = tuple(sorted(
                r for r in replicons
                if search in r.lower()
            ))
        else:
            profile = tuple(sorted(replicons))

        # Skip sequences that don't have any matching replicons
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
        
        # Pull original full replicon profiles for every sequence in this group
        group_replicon_sets = [set(pf_dict.get(seq, [])) for seq in group_seqs]
        
        if group_replicon_sets:
            shared_replicons = sorted(set.intersection(*group_replicon_sets))
            all_replicons = sorted(set.union(*group_replicon_sets))
        else:
            shared_replicons = []
            all_replicons = []

        # Re-derive the filtered profile that caused this group to form
        first_seq = group_seqs[0]
        if target_replicon:
            filtered_profile = sorted([
                r for r in pf_dict.get(first_seq, []) 
                if search in r.lower()
            ])
        else:
            filtered_profile = None

        json_payload.append({
            "group_id": f"plasmidfinder_group_{i:04d}",
            "sequences": group_seqs,
            "filtered_replicon_profile": filtered_profile,
            "shared_replicons": shared_replicons,
            "all_replicons": all_replicons
        })

    # Ensure the subdirectory exists
    outfile = outdir / "plasmidfinder_results" / "plasmidfinder_groups.json"
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w") as f:
        json.dump(json_payload, f, indent=2)

    logging.info(f"Wrote {len(groups)} PlasmidFinder groups to {outfile}")

    # Return the clean list of lists for downstream processing
    return groups