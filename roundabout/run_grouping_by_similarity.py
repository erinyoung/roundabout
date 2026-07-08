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
