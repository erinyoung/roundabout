import logging
from pathlib import Path
import pandas as pd

from .run_visualize_similarity import visualize_as_heatmap

def run_cohort_heatmaps(
    pipeline_groups: dict,
    local_matrix_df: pd.DataFrame,
    global_hits_df: pd.DataFrame,
    outdir: Path,
    sample_to_ref_paths: dict,
):
    """
    Generates two distinct heatmaps per cohort:
      1. A square Local-vs-Local matrix using local_matrix_df.
      2. A rectangular Input-vs-RefSeq matrix using global_hits_df.
    """

    # -----------------------------------------------------------------
    # Prep Local Matrix
    # -----------------------------------------------------------------
    local_matrix_df = local_matrix_df.copy()
    local_matrix_df.index = local_matrix_df.index.astype(str)
    local_matrix_df.columns = local_matrix_df.columns.astype(str)

    # -----------------------------------------------------------------
    # Prep Global Hits
    # -----------------------------------------------------------------
    if global_hits_df is not None and not global_hits_df.empty:
        df_clean = global_hits_df.copy()
        df_clean.columns = [c.lower() for c in df_clean.columns]
        df_clean = df_clean.loc[:, ~df_clean.columns.duplicated()]

        def clean_id(val):
            return Path(str(val)).stem.split(".")[0].strip()

        df_clean["query_clean"] = df_clean["query_name"].apply(clean_id)
        df_clean["ref_clean"] = df_clean["refseq_hit"].apply(clean_id)
    else:
        df_clean = None

    for group_id, members in pipeline_groups.items():
        logging.info(f"Generating heatmaps for {group_id}...")

        # Identify group members present in our local matrix
        local_ids = [str(m) for m in members if str(m) in local_matrix_df.index]
        clean_local_ids = [clean_id(x) for x in local_ids]

        # ========================================================
        # 1. LOCAL VS LOCAL HEATMAP (Square)
        # ========================================================
        if len(local_ids) >= 2:
            # Extract the square block from the local matrix
            local_group_matrix = local_matrix_df.loc[local_ids, local_ids].copy()

            # Apply clean IDs for tidy plot labels
            local_group_matrix.index = clean_local_ids
            local_group_matrix.columns = clean_local_ids

            local_out_dir = (
                Path(outdir) / group_id / "ani_heatmap_results" / "local_only"
            )
            local_out_dir.mkdir(parents=True, exist_ok=True)
            local_heatmap_path = local_out_dir / f"{group_id}_local_ani_heatmap.png"

            visualize_as_heatmap(
                matrix_df=local_group_matrix, out_path=local_heatmap_path
            )
        else:
            logging.warning(f"{group_id}: Not enough members for a local heatmap.")

        # ========================================================
        # 2. INPUT VS REFSEQ HEATMAP (Rectangular)
        # ========================================================
        if df_clean is None:
            continue

        # Collect top references for this specific cohort
        ref_ids = set()
        for member in members:
            ref_paths = sample_to_ref_paths.get(member) or sample_to_ref_paths.get(
                clean_id(member)
            )
            if ref_paths:
                for path in ref_paths:
                    if path:
                        ref_ids.add(clean_id(path.stem))

        cohort_refs = list(ref_ids)
        if not cohort_refs or not clean_local_ids:
            continue

        # Filter global hits to ONLY Rows = Inputs, Columns = Refs
        cohort_rows = df_clean[
            df_clean["query_clean"].isin(clean_local_ids)
            & df_clean["ref_clean"].isin(cohort_refs)
        ].copy()

        if cohort_rows.empty:
            logging.warning(
                f"{group_id}: No matching RefSeq alignments found in global data."
            )
            continue

        # Pivot directly into a rectangular format
        rect_matrix = cohort_rows.pivot_table(
            index="query_clean", columns="ref_clean", values="ani", aggfunc="max"
        )

        # Ensure all inputs and refs are present in the layout, even if some lacked matches
        rect_matrix = rect_matrix.reindex(index=clean_local_ids, columns=cohort_refs)

        # Fill missing intersection gaps with the background floor
        rect_matrix = rect_matrix.fillna(80.0)

        # Save the rectangular heatmap
        refseq_out_dir = Path(outdir) / group_id / "ani_heatmap_results" / "with_refseq"
        refseq_out_dir.mkdir(parents=True, exist_ok=True)
        refseq_heatmap_path = refseq_out_dir / f"{group_id}_input_vs_refseq_heatmap.png"

        visualize_as_heatmap(matrix_df=rect_matrix, out_path=refseq_heatmap_path)

    logging.info("Cohort matrix heatmap generation finished completely.")
