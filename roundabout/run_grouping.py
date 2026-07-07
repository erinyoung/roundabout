import pandas as pd
import logging
from pathlib import Path
import hashlib

logging.basicConfig(level=logging.INFO)

def define_groups_by_similarity(
    local_matrix_df: pd.DataFrame, 
    global_hits_df: pd.DataFrame, 
    min_ani: float, 
    min_ani_fraction_query: float, 
    min_ani_fraction_ref: float, 
    num_references: int
) -> tuple[dict, dict]:
    """
    1. Finds similar local inputs for each query based on ANI thresholds.
    2. Groups queries *solely* by their local similarity profile.
    3. Assigns a unique alphanumeric ID (MD5 hash) to each local group.
    4. Expands to RefSeq hits above the thresholds, keeping top N.
    5. Returns a structured dictionary preserving distinct references per local sequence.
    """
    logging.info(f"Defining plasmid groups using {min_ani}% ANI threshold...")
    
    input_seqs = local_matrix_df.index.tolist()
    
    # --- 1. PRE-FILTER & SORT GLOBAL HITS ---
    if not global_hits_df.empty:
        # Filter by ANI and both query/ref alignment fractions
        valid_globals = global_hits_df[
            (global_hits_df['ANI'] >= min_ani) &
            (global_hits_df['Align_fraction_query'] >= min_ani_fraction_query) &
            (global_hits_df['Align_fraction_ref'] >= min_ani_fraction_ref)
        ].copy()
        
        # Sort by ANI, then AF query, then AF ref (all descending) to prep for Top N
        valid_globals.sort_values(
            by=['ANI', 'Align_fraction_query', 'Align_fraction_ref'],
            ascending=[False, False, False],
            inplace=True
        )
    else:
        valid_globals = pd.DataFrame()

    # --- 2. EXTRACT LOCAL GROUPS AND GLOBAL REFERENCES ---
    raw_local_groups = {}
    seq_to_refs = {}
        
    for seq in input_seqs:
        # Local similarities: Get row, filter by threshold, extract names
        local_row = local_matrix_df.loc[seq]
        local_matches = tuple(sorted(local_row[local_row >= min_ani].index.tolist()))
        raw_local_groups[seq] = local_matches
        
        # Global similarities: Find top N RefSeq hits for this specific query
        if not valid_globals.empty:
            hits = valid_globals[valid_globals['Query_Name'] == seq]['RefSeq_Hit'].head(num_references).tolist()
        else:
            hits = []
        
        seq_to_refs[seq] = hits
        
    # --- 3. DEDUPLICATE GROUPS & ASSIGN ALPHANUMERIC IDs ---
    unique_groups = {}
    seq_to_group = {}
    seen_local_groups = {} # Maps the local_matches tuple to an assigned group ID
    
    for seq, local_members in raw_local_groups.items():
        if local_members not in seen_local_groups:
            # Generate a short alphanumeric hash for the group ID based on its members
            group_hash_str = "_".join(local_members).encode('utf-8')
            assigned_id = f"Group_{hashlib.md5(group_hash_str).hexdigest()[:8].upper()}"
            
            seen_local_groups[local_members] = assigned_id
            
            # Initialize the structured payload for this unique group
            unique_groups[assigned_id] = {
                "local_members": list(local_members),
                "references": {}
            }
        else:
            assigned_id = seen_local_groups[local_members]
            
        seq_to_group[seq] = assigned_id
        
        # Attach this specific sequence's references to the group dictionary
        unique_groups[assigned_id]["references"][seq] = seq_to_refs[seq]
        
    logging.info(f"Collapsed {len(input_seqs)} inputs into {len(unique_groups)} unique groups.")
    
    # seq_to_group maps: {"4051900_3": "Group_8A3F9B2C"}
    # unique_groups maps: {
    #     "Group_8A3F9B2C": {
    #         "local_members": ["4051900_3", "4051904_3"],
    #         "references": {
    #             "4051900_3": ["NZ_JAJIZO010000007.1", "NZ_..."],
    #             "4051904_3": ["NC_123456.1"]
    #         }
    #     }
    # }
   
    return seq_to_group, unique_groups


def define_groups_by_amr(amr_dict: dict, target_gene: str = None) -> tuple[dict, dict]:
    """
    Groups isolates based on their complete profile of AMR genes.
    If target_gene is provided, sequences are grouped by the exact combination 
    of matching genes they possess.
    Drops any groups that do not have > 1 members.
    Returns a mapping of seq->group_id, and group_id->{members, shared_genes, filtered_genes}.
    """
    if target_gene:
        logging.info(f"Defining plasmid groups for AMR genes matching '{target_gene}'...")
        search_term = target_gene.lower()
    else:
        logging.info("Defining plasmid groups by exact AMR gene profiles...")

    # --- 1. DETERMINE PROFILES FOR EACH SEQUENCE ---
    profile_to_seqs = {}
    seq_to_full_genes = {seq: sorted(genes) for seq, genes in amr_dict.items()}

    for seq, genes in amr_dict.items():
        if target_gene:
            # Keep only genes that match the target string
            active_profile = tuple(sorted(g for g in genes if search_term in g.lower()))
        else:
            # Use the entire gene profile
            active_profile = tuple(sorted(genes))
        
        if active_profile not in profile_to_seqs:
            profile_to_seqs[active_profile] = []
        profile_to_seqs[active_profile].append(seq)

    # --- 2. BUILD THE DICTIONARIES (Filtering out singletons) ---
    seq_to_group = {}
    unique_groups = {}

    for active_profile, seqs in profile_to_seqs.items():
        # DROP groups that have no matching AMR genes at all
        if not active_profile:
            continue
            
        # DROP groups that do not have strictly > 1 members
        if len(seqs) <= 1:
            continue
            
        # Sort members to ensure consistent hashing
        members = tuple(sorted(seqs))
        
        # Generate MD5 Hash ID
        cohort_str = "_".join(members)
        short_hash = hashlib.md5(cohort_str.encode('utf-8')).hexdigest()[:8].upper()
        group_id = f"Group_{short_hash}"
        
        # Find the universal "Shared genes" for this group (Intersection of all full gene lists)
        shared_genes = set(seq_to_full_genes[members[0]])
        for seq in members[1:]:
            shared_genes.intersection_update(seq_to_full_genes[seq])
            
        # Construct the payload
        group_payload = {
            "local_members": list(members),
            "Shared genes": sorted(list(shared_genes))
        }
        
        # Explicitly append the filtered genes that caused them to group together
        if target_gene:
            group_payload["Shared filtered genes"] = list(active_profile)
            
        unique_groups[group_id] = group_payload
        
        # Map each sequence to its group (Singletons will NOT be in this dictionary)
        for seq in members:
            seq_to_group[seq] = group_id

    logging.info(f"Grouped inputs into {len(unique_groups)} distinct AMR profiles (excluding singletons).")


    return seq_to_group, unique_groups

def define_groups_by_plasmidfinder(pf_dict: dict, target_replicon: str = None) -> tuple[dict, dict]:
    """
    Groups isolates based on their complete profile of PlasmidFinder replicons.
    If target_replicon is provided, sequences are grouped by the exact combination 
    of matching replicons they possess.
    Drops any empty profiles and groups that do not have > 1 members.
    Returns a mapping of seq->group_id, and group_id->{members, shared_reps, filtered_reps}.
    """
    if target_replicon:
        logging.info(f"Defining groups for replicons matching '{target_replicon}'...")
        search_term = target_replicon.lower()
    else:
        logging.info("Defining groups by exact replicon profiles...")

    # --- 1. DETERMINE PROFILES FOR EACH SEQUENCE ---
    profile_to_seqs = {}
    seq_to_full_reps = {seq: sorted(reps) for seq, reps in pf_dict.items()}

    for seq, reps in pf_dict.items():
        if target_replicon:
            # Keep only replicons that match the target string
            active_profile = tuple(sorted(r for r in reps if search_term in r.lower()))
        else:
            # Use the entire replicon profile
            active_profile = tuple(sorted(reps))
        
        if active_profile not in profile_to_seqs:
            profile_to_seqs[active_profile] = []
        profile_to_seqs[active_profile].append(seq)

    # --- 2. BUILD THE DICTIONARIES (Filtering out singletons and empties) ---
    seq_to_group = {}
    unique_groups = {}

    for active_profile, seqs in profile_to_seqs.items():
        # DROP groups that have no matching replicons at all
        if not active_profile:
            continue
            
        # DROP groups that do not have strictly > 1 members
        if len(seqs) <= 1:
            continue
            
        # Sort members to ensure consistent hashing
        members = tuple(sorted(seqs))
        
        # Generate MD5 Hash ID
        cohort_str = "_".join(members)
        short_hash = hashlib.md5(cohort_str.encode('utf-8')).hexdigest()[:8].upper()
        group_id = f"Group_{short_hash}"
        
        # Find the universal "Shared replicons" for this group 
        shared_reps = set(seq_to_full_reps[members[0]])
        for seq in members[1:]:
            shared_reps.intersection_update(seq_to_full_reps[seq])
            
        # Construct the payload with a consistent schema
        group_payload = {
            "local_members": list(members),
            "Shared replicons": sorted(list(shared_reps)),
            "Shared filtered replicons": list(active_profile) 
        }
            
        unique_groups[group_id] = group_payload
        
        # Map each sequence to its group 
        for seq in members:
            seq_to_group[seq] = group_id

    logging.info(f"Grouped inputs into {len(unique_groups)} distinct PlasmidFinder profiles (excluding singletons/empties).")

    return seq_to_group, unique_groups

def build_consensus_groups(seq_to_sim: dict, seq_to_amr: dict, seq_to_pf: dict) -> tuple[dict, dict]:
    """
    Combines the outputs of Similarity, AMR, and PlasmidFinder groupings.
    Isolates are grouped together ONLY if they share the exact same profile across all three methods.
    """
    logging.info("Building consensus groups across Similarity, AMR, and PlasmidFinder...")
    
    # 1. Gather every unique sequence that appears in ANY of the three dictionaries
    all_seqs = set(seq_to_sim.keys()) | set(seq_to_amr.keys()) | set(seq_to_pf.keys())
    
    # 2. Map each profile tuple to its matching sequences
    profile_to_seqs = {}
    
    for seq in all_seqs:
        # Create a tuple of their 3 group IDs. Use "None" if they were dropped/filtered out.
        master_profile = (
            seq_to_sim.get(seq, None),
            seq_to_amr.get(seq, None),
            seq_to_pf.get(seq, None)
        )
        
        if master_profile not in profile_to_seqs:
            profile_to_seqs[master_profile] = []
        profile_to_seqs[master_profile].append(seq)
        
    # 3. Build the final downstream dictionaries
    seq_to_master = {}
    master_groups = {}
    
    for profile, seqs in profile_to_seqs.items():
        # Optional: if you want to drop singletons here too, add `if len(seqs) <= 1: continue`
        
        members = tuple(sorted(seqs))
        cohort_str = "_".join(members)
        master_id = f"Consensus_{hashlib.md5(cohort_str.encode('utf-8')).hexdigest()[:8].upper()}"
        
        master_groups[master_id] = {
            "local_members": list(members),
            "sim_group": profile[0],
            "amr_group": profile[1],
            "pf_group": profile[2]
        }
        
        for seq in members:
            seq_to_master[seq] = master_id
            
    logging.info(f"Generated {len(master_groups)} highly specific consensus groups.")

    print(seq_to_master)
    print(master_groups)
    
    return seq_to_master, master_groups

import pandas as pd
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)

def write_group_summary(
    master_groups: dict, 
    unique_similarity_groups: dict,
    amr_dict: dict,
    pf_dict: dict,
    local_matrix_df: pd.DataFrame, 
    global_hits_df: pd.DataFrame, 
    outdir: Path,
    target_amr_gene: str = None,
    target_pf_replicon: str = None
) -> Path:
    """
    Generates a structured CSV summary report for each Consensus Plasmid Group,
    detailing similarities, references, and exact AMR/PlasmidFinder compositions.
    """
    logging.info("Compiling consensus plasmid group summary report...")
    
    summary_rows = []
    
    # Pre-process search terms
    amr_search = target_amr_gene.lower() if target_amr_gene else None
    pf_search = target_pf_replicon.lower() if target_pf_replicon else None
    
    for group_id, group_data in master_groups.items():
        local_members = group_data['local_members']
        sim_group_id = group_data['sim_group']
        num_local = len(local_members)
        
        # --- 1. GATHER REFERENCES ---
        # Pull references from the original similarity dictionary if a sim group exists
        references_dict = {}
        all_ref_members = set()
        
        if sim_group_id and sim_group_id in unique_similarity_groups:
            references_dict = unique_similarity_groups[sim_group_id].get('references', {})
            for refs in references_dict.values():
                all_ref_members.update(refs)
                
        num_ref = len(all_ref_members)
        
        # Format a clean string for the CSV: "Seq1: [RefA, RefB] | Seq2: [RefC]"
        ref_strings = []
        for seq, refs in references_dict.items():
            if seq in local_members and refs:
                ref_strings.append(f"{seq}: [{', '.join(refs)}]")
        formatted_references = " | ".join(ref_strings) if ref_strings else "None"

        # --- 2. CALCULATE SIMILARITY BOUNDS ---
        group_anis = []
        
        # Local similarities
        if num_local > 1:
            for i in range(len(local_members)):
                for j in range(i + 1, len(local_members)):
                    ani = local_matrix_df.loc[local_members[i], local_members[j]]
                    if ani > 0: 
                        group_anis.append(ani)
                        
        # Global similarities
        if num_local > 0 and num_ref > 0:
            group_global_hits = global_hits_df[
                global_hits_df['Query_Name'].isin(local_members) & 
                global_hits_df['RefSeq_Hit'].isin(all_ref_members)
            ]
            if not group_global_hits.empty:
                group_anis.extend(pd.to_numeric(group_global_hits['ANI'], errors='coerce').dropna().tolist())
        
        min_sim = round(min(group_anis), 2) if group_anis else 100.00
        max_sim = round(max(group_anis), 2) if group_anis else 100.00
        avg_sim = round(sum(group_anis) / len(group_anis), 2) if group_anis else 100.00

        # --- 3. CALCULATE AMR METADATA ---
        amr_lists = [amr_dict.get(seq, []) for seq in local_members]
        
        amr_found = set().union(*amr_lists) if amr_lists else set()
        amr_shared = set.intersection(*map(set, amr_lists)) if amr_lists else set()
        
        filtered_amr_found = {g for g in amr_found if amr_search in g.lower()} if amr_search else set()
        filtered_amr_shared = {g for g in amr_shared if amr_search in g.lower()} if amr_search else set()

        # --- 4. CALCULATE PLASMIDFINDER METADATA ---
        pf_lists = [pf_dict.get(seq, []) for seq in local_members]
        
        pf_found = set().union(*pf_lists) if pf_lists else set()
        pf_shared = set.intersection(*map(set, pf_lists)) if pf_lists else set()
        
        filtered_pf_found = {r for r in pf_found if pf_search in r.lower()} if pf_search else set()
        filtered_pf_shared = {r for r in pf_shared if pf_search in r.lower()} if pf_search else set()

        # --- 5. COMPILE ROW ---
        summary_rows.append({
            "Consensus_Group_ID": group_id,
            "Num_Local_Sequences": num_local,
            "Min_Similarity_ANI": min_sim,
            "Max_Similarity_ANI": max_sim,
            "Avg_Similarity_ANI": avg_sim,
            "Local_Isolate_IDs": ", ".join(local_members),
            "Mapped_References": formatted_references,
            "AMR_Genes_Found": ", ".join(sorted(amr_found)),
            "AMR_Genes_Shared": ", ".join(sorted(amr_shared)),
            "Filtered_AMR_Genes_Found": ", ".join(sorted(filtered_amr_found)) if amr_search else "N/A",
            "Filtered_AMR_Genes_Shared": ", ".join(sorted(filtered_amr_shared)) if amr_search else "N/A",
            "PF_Replicons_Found": ", ".join(sorted(pf_found)),
            "PF_Replicons_Shared": ", ".join(sorted(pf_shared)),
            "Filtered_PF_Replicons_Found": ", ".join(sorted(filtered_pf_found)) if pf_search else "N/A",
            "Filtered_PF_Replicons_Shared": ", ".join(sorted(filtered_pf_shared)) if pf_search else "N/A"
        })
        
    # Convert to DataFrame and save
    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values(by="Consensus_Group_ID")
    
    report_path = outdir / "consensus_plasmid_groups_summary.csv"
    summary_df.to_csv(report_path, index=False)
    logging.info(f"Consensus group summary report successfully saved to: {report_path}")
    
    return report_path