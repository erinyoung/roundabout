import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import logging


def visualize_as_heatmap(matrix_df: pd.DataFrame, out_path: Path):
    """
    Visualizes the square skani matrix using hierarchical clustering.
    Sorts similar plasmids together, hides dendrogram lines, and neatly places the legend.
    Falls back to unclustered heatmap layout if variance errors occur.
    """
    num_seqs = len(matrix_df)
    if num_seqs == 0:
        logging.warning("Empty matrix provided for visualization. Skipping heatmap.")
        return

    show_annot = False

    fig_width = max(8, num_seqs * 0.4)
    fig_height = max(8, num_seqs * 0.4)

    # --- FIXED: Try clustering first, drop clustering parameters if zero-variance / NaN errors trip ---
    try:
        cg = sns.clustermap(
            matrix_df,
            cmap="viridis",
            annot=show_annot,
            figsize=(fig_width, fig_height),
            cbar_kws={"label": "ANI (%)"},
            linewidths=0.5 if num_seqs <= 30 else 0,
            vmin=80,
            vmax=100,
            dendrogram_ratio=(0.01, 0.01),
            cbar_pos=(1.02, 0.15, 0.03, 0.7),
        )
    except (FloatingPointError, ValueError) as e:
        logging.warning(
            f"Clustering failed for cohort map ({e}). Plotting unclustered symmetric grid fallback."
        )
        # Re-run clustermap but disable the fastcluster linkage algorithm entirely
        cg = sns.clustermap(
            matrix_df,
            cmap="viridis",
            annot=show_annot,
            figsize=(fig_width, fig_height),
            cbar_kws={"label": "ANI (%)"},
            linewidths=0.5 if num_seqs <= 30 else 0,
            vmin=80,
            vmax=100,
            dendrogram_ratio=(0.01, 0.01),
            cbar_pos=(1.02, 0.15, 0.03, 0.7),
            row_cluster=False,  # <-- Disables row dendrogram math
            col_cluster=False,  # <-- Disables column dendrogram math
        )
    # --------------------------------------------------------------------------------------------------

    # Hide the dendrogram tree lines (keeps the sorting, drops the ugly visuals)
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)

    # Strip off the raw "query_name" and "target_name" axis titles
    cg.ax_heatmap.set_xlabel("")
    cg.ax_heatmap.set_ylabel("")

    # Rotate axis labels to fit the sequence names neatly
    plt.setp(
        cg.ax_heatmap.get_xticklabels(),
        rotation=45,
        ha="right",
        fontsize=max(8, 12 - (num_seqs // 10)),
    )
    plt.setp(
        cg.ax_heatmap.get_yticklabels(),
        rotation=0,
        fontsize=max(8, 12 - (num_seqs // 10)),
    )

    # Add a clean title at the top
    cg.figure.suptitle("Local Isolates Clustered ANI", y=1.02, fontsize=14)

    # bbox_inches='tight' ensures the colorbar on the right doesn't get cut off
    cg.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def visualize_global_matches_scatter(global_df: pd.DataFrame, out_path: Path):
    """
    Visualizes global hits using a strip plot.
    Dots are colored by the source organism, capping at the Top 10 most frequent
    and pooling the rest into an 'Other Species' category.
    """
    if global_df.empty:
        logging.getLogger(__name__).warning(
            "Empty global dataframe. Skipping scatter visualization."
        )
        return

    logging.getLogger("matplotlib.category").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)

    ani_col = "ANI" if "ANI" in global_df.columns else "ani"

    # Clean data
    global_df = global_df.copy()
    global_df[ani_col] = pd.to_numeric(global_df[ani_col], errors="coerce")
    global_df["Query_Name"] = global_df["Query_Name"].astype(str)

    # 1. Collapse the long list into Top 10 + Other
    if "organism" in global_df.columns:
        global_df["organism"] = global_df["organism"].fillna("Unknown")

        global_df["organism"] = global_df["organism"].apply(
            lambda x: " ".join(str(x).split()[:2]) if x != "Unknown" else x
        )

        # Determine the top 10 most frequent species across the whole dataset
        top_10_organisms = (
            global_df["organism"].value_counts().nlargest(10).index.tolist()
        )

        # If 'Unknown' is in the top 10, remove it so it doesn't mask true taxa
        if "Unknown" in top_10_organisms:
            top_10_organisms.remove("Unknown")
            # Grab one more to make it a true 10 if available
            extra_top = global_df["organism"].value_counts().nlargest(11).index.tolist()
            if len(extra_top) > 10:
                top_10_organisms.append(extra_top[-1])

        # Replace any organism outside the top 10 with "Other Species"
        global_df["Plot_Color_Group"] = global_df["organism"].apply(
            lambda x: x if x in top_10_organisms or x == "Unknown" else "Other Species"
        )
        has_taxa = True
    else:
        has_taxa = False
        logging.getLogger(__name__).warning(
            "Organism column missing. Defaulting to single color."
        )

    num_queries = global_df["Query_Name"].nunique()
    fig_height = max(6, num_queries * 0.6)

    fig, ax = plt.subplots(figsize=(11, fig_height))

    # 2. Plot the data using our collapsed category column
    sns.stripplot(
        data=global_df,
        x=ani_col,
        y="Query_Name",
        hue="Plot_Color_Group" if has_taxa else None,
        palette=(
            "tab10" if has_taxa else None
        ),  # 'tab10' perfectly matches our 10 discrete categories
        color="teal" if not has_taxa else None,
        jitter=0.2,
        alpha=0.6,
        size=4,
        ax=ax,
    )

    ax.axvline(
        x=95.0, color="black", linestyle="--", alpha=0.5, label="95% ANI Threshold"
    )

    ax.set_title("Global Database (RefSeq) Matches per Input", pad=20, fontsize=14)
    ax.set_xlabel("ANI (%)", labelpad=10)
    ax.set_ylabel("Local Sequences", labelpad=10)

    # Grab the clean legend elements
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0.0,
        fontsize=9,
        title="Top Organisms",
    )

    # Manual constraints to completely accommodate the layout margins
    plt.subplots_adjust(left=0.15, right=0.75, top=0.90, bottom=0.10)

    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
