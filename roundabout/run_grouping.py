import pandas as pd
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)

def define_groups_by_similarity(local_matrix_df: pd.DataFrame, global_hits_df: pd.DataFrame, min_ani: float , min_ani_fraction_query: float, min_ani_fraction_ref: float, num_references: int) -> tuple[dict, dict]:

    """
    1. Finds similar local inputs for each query.
    2. Expands to RefSeq hits above the identity threshold.
    3. Combines them into groups and removes exact duplicates.
    4. Assigns alphanumeric IDs to unique groups.
    """
    logging.info(f"Defining plasmid groups using {min_ani}% ANI threshold...")
    
    raw_groups = {}
    input_seqs = local_matrix_df.index.tolist()
    
    # Pre-filter the global hits to speed up lookups
    if not global_hits_df.empty:
        # sort by 1) ANI, then Align_fraction_query, then Align_fraction_ref. All descending.
        # TODO : add in min_ani_fraction_query and min_ani_fraction_ref filtering
        valid_globals = global_hits_df[global_hits_df['ANI'] >= min_ani]
        
    else:
        valid_globals = pd.DataFrame()
        
    for seq in input_seqs:
        # 1. Local similarities: Get row, filter by threshold, extract names
        local_row = local_matrix_df.loc[seq]
        local_matches = local_row[local_row >= min_ani].index.tolist()
        
        # 2. Global similarities: Find RefSeq hits for this specific query
        global_matches = []
        if not valid_globals.empty:
            # TODO : keep only the top N references based on num_references
            hits = valid_globals[valid_globals['Query_Name'] == seq]['RefSeq_Hit'].tolist()
            global_matches.extend(hits)
            
        # Combine and sort to create a consistent, hashable tuple
        group_members = tuple(sorted(set(local_matches + global_matches)))
        raw_groups[seq] = group_members
        
    # 3 & 4. Deduplicate groups and assign alphanumeric names
    # TODO: fix naming scheme from 001 to a random alphanumeric string to avoid confusion with other groups
    unique_groups = {}
    seq_to_group = {}
    group_counter = 1
    
    for seq, members in raw_groups.items():
        # Check if this exact set of members already has a group ID assigned
        assigned_id = None
        for g_id, g_members in unique_groups.items():
            if members == g_members:
                assigned_id = g_id
                break
        
        # If it's a completely new unique group, generate a new ID
        if not assigned_id:
            assigned_id = f"Group_{group_counter:03d}"
            unique_groups[assigned_id] = members
            group_counter += 1
            
        seq_to_group[seq] = assigned_id
        
    logging.info(f"Collapsed {len(input_seqs)} inputs into {len(unique_groups)} unique groups.")
    
    # seq_to_group maps: {"4051900_3": "Group_001"}
    # unique_groups maps: {"Group_001": ("4051900_3", "4051904_3", "NZ_JAJIZO010000007.1")}
    return seq_to_group, unique_groups

def write_group_summary(unique_groups: dict, local_matrix_df: pd.DataFrame, global_hits_df: pd.DataFrame, outdir: Path):
    """
    Generates a structured CSV summary report for each unique plasmid group.
    """
    logging.info("Compiling plasmid group summary report...")
    
    summary_rows = []
    
    for group_id, members in unique_groups.items():
        # Separate local inputs from global RefSeq references
        local_members = [m for m in members if m in local_matrix_df.index]
        ref_members = [m for m in members if m not in local_matrix_df.index]
        
        total_seqs = len(members)
        num_local = len(local_members)
        num_ref = len(ref_members)
        
        # Collect all pairwise similarity scores within this group to find Min/Max/Avg
        group_anis = []
        
        # 1. Gather from local-vs-local matrix
        if num_local > 1:
            for i in range(len(local_members)):
                for j in range(i + 1, len(local_members)):
                    ani = local_matrix_df.loc[local_members[i], local_members[j]]
                    if ani > 0: # Skani fills missing with 0, don't let it tank the minimum
                        group_anis.append(ani)
                        
        # 2. Gather from global hits dataframe
        if num_local > 0 and num_ref > 0:
            group_global_hits = global_hits_df[
                global_hits_df['Query_Name'].isin(local_members) & 
                global_hits_df['RefSeq_Hit'].isin(ref_members)
            ]
            if not group_global_hits.empty:
                group_anis.extend(pd.to_numeric(group_global_hits['ANI'], errors='coerce').dropna().tolist())
        
        # Calculate similarity bounds safely
        min_sim = round(min(group_anis), 2) if group_anis else 100.00
        max_sim = round(max(group_anis), 2) if group_anis else 100.00
        avg_sim = round(sum(group_anis) / len(group_anis), 2) if group_anis else 100.00
        
        # Pull the dominant organism representation for this group from global hits metadata
        top_organism = "N/A (No Global References)"
        if num_ref > 0 and 'organism' in global_hits_df.columns:
            group_global_hits = global_hits_df[
                global_hits_df['Query_Name'].isin(local_members) & 
                global_hits_df['RefSeq_Hit'].isin(ref_members)
            ]
            if not group_global_hits.empty:
                top_organism = group_global_hits['organism'].value_counts().idxmax()

        # Build row
        summary_rows.append({
            "Group_ID": group_id,
            "Total_Sequences": total_seqs,
            "Num_Local_Sequences": num_local,
            "Num_Reference_Sequences": num_ref,
            "Min_Similarity_ANI": min_sim,
            "Max_Similarity_ANI": max_sim,
            "Avg_Similarity_ANI": avg_sim,
            "Top_Global_Species": top_organism,
            "Local_Isolate_IDs": ", ".join(local_members)
        })
        
    # Convert to DataFrame and save
    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values(by="Group_ID")
    
    report_path = outdir / "skani_results" / "plasmid_groups_summary.csv"
    summary_df.to_csv(report_path, index=False)
    logging.info(f"Group summary report successfully saved to: {report_path}")
    
    return report_path