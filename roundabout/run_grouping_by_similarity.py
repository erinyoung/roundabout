import json
import logging
from pathlib import Path
from itertools import combinations
import pandas as pd


def define_groups_by_similarity(
    raw_skani_df: pd.DataFrame,
    local_matrix_df: pd.DataFrame,
    min_ani: float,
    min_ani_align_fraction_query: float,
    min_ani_align_fraction_ref: float,
    outdir: Path,
) -> list[list[str]]:
    """
    Groups sequences that share the same local ANI neighborhood based on
    stringent ANI and alignment fraction thresholds.
    """
    logging.info(
        f"Defining similarity groups using boundaries: "
        f">={min_ani}% ANI, >={min_ani_align_fraction_query}% Query AF, >={min_ani_align_fraction_ref}% Ref AF"
    )

    # 1. Standardize column names to lowercase to ensure lookup safety
    df_clean = raw_skani_df.copy()
    df_clean.columns = [c.lower() for c in df_clean.columns]

    # Map standard Skani matrix headers dynamically
    query_col = next(
        (
            c
            for c in df_clean.columns
            if "query" in c and ("name" in c or "id" in c or "file" in c)
        ),
        "query_name",
    )
    ref_col = next(
        (
            c
            for c in df_clean.columns
            if "ref" in c and ("name" in c or "id" in c or "file" in c)
        ),
        "refseq_hit",
    )
    ani_col = "ani" if "ani" in df_clean.columns else "ani"
    af_query_col = "align_fraction_query"
    af_ref_col = "align_fraction_ref"

    # 2. Build a mapping of valid structural neighbors
    # For every isolate, it must know exactly which other isolates passed ALL three filters
    valid_neighbors = {str(seq): {str(seq)} for seq in local_matrix_df.index}

    for _, row in df_clean.iterrows():
        q = (
            Path(str(row[query_col])).header_id
            if hasattr(row[query_col], "header_id")
            else Path(str(row[query_col])).stem.split(".")[0]
        )
        r = (
            Path(str(row[ref_col])).header_id
            if hasattr(row[ref_col], "header_id")
            else Path(str(row[ref_col])).stem.split(".")[0]
        )

        # Ensure we are looking strictly at local vs local comparisons inside the matrix index
        if q not in valid_neighbors or r not in valid_neighbors:
            continue

        try:
            ani_val = float(row[ani_col])
            af_q_val = float(row[af_query_col]) if af_query_col in row else 100.0
            af_r_val = float(row[af_ref_col]) if af_ref_col in row else 100.0

            # Enforce the complete boundary checklist
            if (
                ani_val >= min_ani
                and af_q_val >= min_ani_align_fraction_query
                and af_r_val >= min_ani_align_fraction_ref
            ):

                valid_neighbors[q].add(r)
                valid_neighbors[r].add(q)
        except (ValueError, TypeError):
            continue

    # 3. Cluster profiles based on shared neighborhoods
    profile_to_group = {}
    for seq, neighbors in valid_neighbors.items():
        profile = tuple(sorted(list(neighbors)))
        profile_to_group.setdefault(profile, []).append(seq)

    groups = [sorted(group) for group in profile_to_group.values()]
    groups.sort(key=lambda g: (len(g), g), reverse=True)

    # ---------------------------------------------------------
    # Build detailed metadata summary payloads
    # ---------------------------------------------------------
    # TODO : redo summary file into something more informative
    json_payload = []

    for i, group_seqs in enumerate(groups, start=1):
        ani_values = []
        if len(group_seqs) > 1:
            for a, b in combinations(group_seqs, 2):
                try:
                    val = local_matrix_df.loc[a, b]
                    if pd.notna(val) and float(val) > 0:
                        ani_values.append(float(val))
                except KeyError:
                    continue

        min_val = min(ani_values) if ani_values else None
        max_val = max(ani_values) if ani_values else None
        avg_val = sum(ani_values) / len(ani_values) if ani_values else None

        json_payload.append(
            {
                "group_id": f"similarity_group_{i:04d}",
                "sequences": group_seqs,
                "min_ani": min_val,
                "max_ani": max_val,
                "avg_ani": avg_val,
            }
        )

    outfile = outdir / "skani_results" / "similarity_groups.json"
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w") as f:
        json.dump(json_payload, f, indent=2)

    logging.info(f"Wrote {len(groups)} similarity groups to {outfile}")
    return groups

import logging
from pathlib import Path
import pandas as pd

def define_groups_by_ani(
    raw_skani_df: pd.DataFrame,
    local_matrix_df: pd.DataFrame,
    min_ani: float,
    min_ani_align_fraction_query: float,
    min_ani_align_fraction_ref: float,
    outdir: Path,
) -> list[list[str]]:
    """
    For each local sample, finds all passing skani matches. 
    Returns a sorted, deduplicated list of these match lists.
    """
    logging.info(
        f"Defining ANI similarity groups using boundaries: "
        f">={min_ani}% ANI, >={min_ani_align_fraction_query}% Query AF, >={min_ani_align_fraction_ref}% Ref AF"
    )

    # 1. Get list of local samples
    local_samples = set(local_matrix_df.columns) - {"query_name"}

    # 2. Extract clean stems from file paths
    skani_filtered = raw_skani_df.copy()
    skani_filtered["parsed_query"] = skani_filtered["Query_file"].apply(lambda x: Path(x).stem)
    skani_filtered["parsed_ref"] = skani_filtered["Ref_file"].apply(lambda x: Path(x).stem)

    # 3. Keep only rows where both query and ref are local samples
    skani_filtered = skani_filtered[
        skani_filtered["parsed_query"].isin(local_samples) & 
        skani_filtered["parsed_ref"].isin(local_samples)
    ]

    # 4. Apply metric filters (handling the 0-1 scale of Align_fraction_query)
    af_query_cutoff = min_ani_align_fraction_query / 100.0 if min_ani_align_fraction_query > 1.0 else min_ani_align_fraction_query

    passing_df = skani_filtered[
        (skani_filtered["ANI"] >= min_ani) &
        (skani_filtered["Align_fraction_ref"] >= min_ani_align_fraction_ref) &
        (skani_filtered["Align_fraction_query"] >= af_query_cutoff)
    ]

    # 5. Build raw match lists for each sample
    raw_groups = []
    for sample in local_samples:
        # Find all references where this sample is the query
        matches = passing_df[passing_df["parsed_query"] == sample]["parsed_ref"].tolist()
        
        # Combine the sample itself with all its distinct matches
        full_match_set = {sample}.union(matches)
        raw_groups.append(tuple(sorted(full_match_set)))

    # 6. Deduplicate lists at the end and sort by size
    unique_groups = [list(g) for g in set(raw_groups)]
    unique_groups.sort(key=lambda g: (len(g), g), reverse=True)

    logging.info(f"Found {len(unique_groups)} ANI similarity groups")
    return unique_groups