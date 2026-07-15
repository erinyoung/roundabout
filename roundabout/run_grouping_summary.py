import json
import logging
from pathlib import Path
from itertools import combinations

import pandas as pd

logging.basicConfig(level=logging.INFO)


def calculate_jaccard(set_a: set, set_b: set) -> float:
    """
    Calculate Jaccard similarity between two sets.

    Returns:
        0.0 if both sets are empty or no overlap exists.
    """

    if not set_a and not set_b:
        return 1.0

    union = set_a | set_b

    if not union:
        return 0.0

    return len(set_a & set_b) / len(union)


def get_group_jaccard_statistics(
    seq_ids: list[str],
    feature_dict: dict,
) -> dict:
    """
    Calculate pairwise Jaccard similarity statistics
    for features within a group.
    """

    values = []

    for a, b in combinations(seq_ids, 2):

        features_a = set(feature_dict.get(a, []))
        features_b = set(feature_dict.get(b, []))

        values.append(
            calculate_jaccard(
                features_a,
                features_b,
            )
        )

    return summarize_values(values)


def calculate_gc_content(fasta_path: Path) -> float:
    """
    Calculate GC percentage from a FASTA file.
    """

    sequence = []

    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                sequence.append(line.strip())

    sequence = "".join(sequence).upper()

    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")

    return (gc_count / len(sequence)) * 100


def build_sequence_metadata(staged_data: list[dict]) -> dict:
    """
    Build metadata lookup from staged FASTA information.

    Returns
    -------
    dict

    {
        sequence_id: {
            length: int,
            gc: float,
            path: Path
        }
    }
    """

    metadata = {}

    logging.info("Calculating sequence GC content...")

    for record in staged_data:
        seq_id = record["combined_id"]

        metadata[seq_id] = {
            "length": record["length"],
            "gc": calculate_gc_content(record["path"]),
            "path": record["path"],
        }

    return metadata


def summarize_values(values: list[float]) -> dict:
    """
    Calculate summary statistics.
    """

    if not values:
        return {
            "mean": None,
            "min": None,
            "max": None,
            "std": None,
        }

    series = pd.Series(values)

    return {
        "mean": float(series.mean()),
        "min": float(series.min()),
        "max": float(series.max()),
        "std": float(series.std()),
    }


def shared_features(
    seq_ids: list[str],
    feature_dict: dict,
) -> list[str]:
    """
    Return features present in every sequence.
    """

    if not seq_ids:
        return []

    profiles = [set(feature_dict.get(seq, [])) for seq in seq_ids]

    shared = set.intersection(*profiles)

    return sorted(shared)


def unique_feature_count(
    seq_ids: list[str],
    feature_dict: dict,
) -> int:
    """
    Count unique features across all sequences.
    """

    features = set()

    for seq in seq_ids:
        features.update(feature_dict.get(seq, []))

    return len(features)


def get_group_ani_statistics(
    seq_ids: list[str],
    skani_df: pd.DataFrame,
) -> dict:
    """
    Calculate within-group pairwise ANI statistics.

    Self comparisons are excluded.
    """

    ani_values = []

    for a, b in combinations(seq_ids, 2):

        try:
            ani = skani_df.loc[a, b]
        except KeyError:
            continue

        if pd.notna(ani):
            ani_values.append(float(ani))

    return summarize_values(ani_values)


def run_grouping_summary(
    groups: list[list[str]],
    staged_data: list[dict],
    amr_dict: dict,
    plasmidfinder_dict: dict,
    skani_df: pd.DataFrame,
    outdir: Path,
    min_num_seqs: int,
    target_amr_gene: str | None = None,
    target_pf_string: str | None = None,
) -> dict:
    """
    Generate summary statistics for sequence groups and individual samples.
    Saves as a combined JSON and separate TSVs.
    """

    outdir.mkdir(parents=True, exist_ok=True)
    metadata = build_sequence_metadata(staged_data)

    # -------------------------------------------------------------------------
    # Generate Filtered Dictionaries Internally
    # -------------------------------------------------------------------------
    filtered_amr_dict = None
    if target_amr_gene:
        filtered_amr_dict = {
            seq: [g for g in genes if target_amr_gene.lower() in g.lower()]
            for seq, genes in amr_dict.items()
        }

    filtered_plasmidfinder_dict = None
    if target_pf_string:
        filtered_plasmidfinder_dict = {
            seq: [r for r in replicons if target_pf_string.lower() in r.lower()]
            for seq, replicons in plasmidfinder_dict.items()
        }

    logging.info(
        f"Summarizing {len(groups)} unique groups and {len(metadata)} samples..."
    )

    group_summaries = []
    seq_to_groups = {}  # Dictionary to track which groups each sequence ends up in

    # -------------------------------------------------------------------------
    # 1. Build Group Summaries
    # -------------------------------------------------------------------------
    for i, seq_ids in enumerate(groups, start=1):

        group_id = f"group_{i:04d}"

        # Map each sequence to this group ID for the sample summary later
        for seq in seq_ids:
            seq_to_groups.setdefault(seq, []).append(group_id)

        lengths = [metadata[s]["length"] for s in seq_ids if s in metadata]
        gc_values = [metadata[s]["gc"] for s in seq_ids if s in metadata]

        result = {
            "group_id": group_id,
            "sequence_ids": sorted(seq_ids),
            "number_of_sequences": len(seq_ids),
            "seq_num_below_threshold": len(seq_ids) < min_num_seqs,
            "length": summarize_values(lengths),
            "gc_content": summarize_values(gc_values),
            "ani": get_group_ani_statistics(seq_ids, skani_df),
            "shared_amr": shared_features(seq_ids, amr_dict),
            "shared_plasmidfinder": shared_features(seq_ids, plasmidfinder_dict),
            "unique_amr_count": unique_feature_count(seq_ids, amr_dict),
            "unique_plasmidfinder_count": unique_feature_count(
                seq_ids, plasmidfinder_dict
            ),
            "amr_jaccard": get_group_jaccard_statistics(seq_ids, amr_dict),
            "plasmidfinder_jaccard": get_group_jaccard_statistics(
                seq_ids, plasmidfinder_dict
            ),
        }

        if filtered_amr_dict:
            result["shared_filtered_amr"] = shared_features(seq_ids, filtered_amr_dict)
        if filtered_plasmidfinder_dict:
            result["shared_filtered_plasmidfinder"] = shared_features(
                seq_ids, filtered_plasmidfinder_dict
            )

        group_summaries.append(result)

    # -------------------------------------------------------------------------
    # 2. Build Sample Summaries
    # -------------------------------------------------------------------------
    sample_summaries = []

    # Iterate over every sequence we staged, even if it wasn't placed in a group
    for seq_id, meta in metadata.items():
        sample_summaries.append(
            {
                "sequence_id": seq_id,
                "length": meta["length"],
                "gc_content": meta["gc"],
                "assigned_groups": sorted(seq_to_groups.get(seq_id, [])),
                "amr_genes": sorted(amr_dict.get(seq_id, [])),
                "plasmidfinder_replicons": sorted(plasmidfinder_dict.get(seq_id, [])),
            }
        )

    # ---------------------------------------------------------
    # Output 1: Combined JSON
    # ---------------------------------------------------------
    final_summary = {"groups": group_summaries, "samples": sample_summaries}

    json_outfile = outdir / "group_summary.json"
    with open(json_outfile, "w") as f:
        json.dump(final_summary, f, indent=2)
    logging.info(f"Wrote complete JSON summary to {json_outfile}")

    # ---------------------------------------------------------
    # Output 2: Group TSV
    # ---------------------------------------------------------
    df_group = pd.json_normalize(group_summaries)
    list_columns_group = [
        col for col in df_group.columns if df_group[col].apply(type).eq(list).any()
    ]
    for col in list_columns_group:
        df_group[col] = df_group[col].apply(
            lambda x: ",".join(x) if isinstance(x, list) else x
        )

    group_tsv = outdir / "group_summary.tsv"
    df_group.to_csv(group_tsv, sep="\t", index=False)
    logging.info(f"Wrote Group TSV to {group_tsv}")

    # ---------------------------------------------------------
    # Output 3: Sample TSV
    # ---------------------------------------------------------
    df_sample = pd.json_normalize(sample_summaries)
    list_columns_sample = [
        col for col in df_sample.columns if df_sample[col].apply(type).eq(list).any()
    ]
    for col in list_columns_sample:
        df_sample[col] = df_sample[col].apply(
            lambda x: ",".join(x) if isinstance(x, list) else x
        )

    sample_tsv = outdir / "sample_summary.tsv"
    df_sample.to_csv(sample_tsv, sep="\t", index=False)
    logging.info(f"Wrote Sample TSV to {sample_tsv}")

    # ---------------------------------------------------------
    # Return Group Mapping
    # ---------------------------------------------------------
    group_mapping = {g["group_id"]: g["sequence_ids"] for g in group_summaries}

    # Create a quick summary of group sizes to keep the log clean
    mapping_summary = {group: len(seqs) for group, seqs in group_mapping.items()}
    logging.info(
        f"Final group assignments (Group: Sequence Count) -> {mapping_summary}"
    )

    return group_mapping
