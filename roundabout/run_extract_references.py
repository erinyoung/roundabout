import os
import logging
from pathlib import Path
import pandas as pd


def sort_dataframe_for_refseq_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sorts a DataFrame based on a compound homology score derived from
    alignment fractions and ANI.
    """
    required_metrics = ["align_fraction_query", "align_fraction_ref", "ani"]
    if all(m in df.columns for m in required_metrics):
        # Calculate a compound fraction score (0.0 to 1.0)
        df["homology_score"] = (
            (df["align_fraction_query"] / 100.0)
            * (df["align_fraction_ref"] / 100.0)
            * (df["ani"] / 100.0)
        )
        # Sort primarily by this combined score
        df = df.sort_values(by="homology_score", ascending=False)
    else:
        # Fallback to standard sorting if columns are missing
        available_sort_cols = [c for c in required_metrics if c in df.columns]
        df = df.sort_values(by=available_sort_cols, ascending=False)

    return df


def get_top_refseq_ids_for_sample(
    global_hits_df: pd.DataFrame,
    sample_id: str,
    num_refs: int = 5,
    min_ani: float = 95.0,
    min_af_query: float = 15.0,
    min_af_ref: float = 15.0,
) -> list[str]:
    """
    Identifies the top N RefSeq reference identifiers for a single sample
    based on homology and structural alignment thresholds.
    """
    if global_hits_df is None or global_hits_df.empty:
        return []

    df_copy = global_hits_df.copy()
    df_copy.columns = [c.lower() for c in df_copy.columns]
    df_copy = df_copy.loc[:, ~df_copy.columns.duplicated()]

    def clean_id(val):
        return Path(str(val)).stem.split(".")[0].strip()

    # Bypasses index alignment to prevent duplicate label crashes
    df_copy["clean_query"] = [clean_id(x) for x in df_copy["query_name"]]
    clean_sample = clean_id(sample_id)

    sample_hits = df_copy[df_copy["clean_query"] == clean_sample].copy()

    if sample_hits.empty:
        return []

    if "ani" in sample_hits.columns:
        sample_hits = sample_hits[sample_hits["ani"] >= min_ani]

    if "align_fraction_query" in sample_hits.columns:
        sample_hits = sample_hits[sample_hits["align_fraction_query"] >= min_af_query]

    if "align_fraction_ref" in sample_hits.columns:
        sample_hits = sample_hits[sample_hits["align_fraction_ref"] >= min_af_ref]

    if sample_hits.empty:
        return []

    sample_hits = sort_dataframe_for_refseq_hits(sample_hits)

    # Use actual reference headers matching the multi-FASTA formatting requirements
    return sample_hits["refseq_hit"].drop_duplicates().head(num_refs).tolist()


def extract_refseq_fastas(
    refseq_ids: set[str],
    refseq_fasta_path: Path,
    outdir: Path,
) -> dict[str, Path]:
    """
    Extracts a batch collection of unique RefSeq records from the master multi-FASTA
    in a single streaming pass.

    Returns
    -------
    dict
        Mapping of cleaned RefSeq ID -> staged FASTA Path.
    """
    if not refseq_ids:
        return {}

    if not refseq_fasta_path or not refseq_fasta_path.exists():
        logging.warning(f"RefSeq master FASTA not found: {refseq_fasta_path}")
        return {}

    def clean_id(val):
        return Path(str(val)).stem.split(".")[0].strip()

    needed_headers = {str(r).strip() for r in refseq_ids}
    needed_clean = {clean_id(r) for r in refseq_ids}

    stage_dir = Path(outdir) / "refseq_plasmids_dl_staged"
    stage_dir.mkdir(parents=True, exist_ok=True)

    staged_paths = {}

    try:
        with open(refseq_fasta_path, "r") as f:
            current_id = None
            current_clean = None
            current_seq = []

            for line in f:
                if line.startswith(">"):
                    # Write out the completed sequence block if it matches our targets
                    if current_id and (
                        current_id in needed_headers or current_clean in needed_clean
                    ):
                        out_path = stage_dir / f"{current_clean}.fasta"
                        if not out_path.exists():
                            with open(out_path, "w") as out_f:
                                out_f.write("".join(current_seq))
                        staged_paths[current_clean] = out_path

                    raw_header = line[1:].strip()
                    current_id = raw_header.split()[0]
                    current_clean = clean_id(current_id)
                    current_seq = [line]
                else:
                    if current_id:
                        current_seq.append(line)

            # Write out the final record in the file if it matches
            if current_id and (
                current_id in needed_headers or current_clean in needed_clean
            ):
                out_path = stage_dir / f"{current_clean}.fasta"
                if not out_path.exists():
                    with open(out_path, "w") as out_f:
                        out_f.write("".join(current_seq))
                staged_paths[current_clean] = out_path

    except Exception as e:
        logging.error(f"Error extracting RefSeq FASTAs in single-pass stream: {e}")
        return {}

    return staged_paths
