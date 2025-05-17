import streamlit as st
import numpy as np

from molstorage import (
    get_compound_safety_data,
    get_name_and_smiles,
    classify_acid_base,
    get_mp_bp,
    compound_state,
    prioritize_pictograms,
    initialize_storage_groups,
    chemsort_multiple_order_3
)

# Page config
st.set_page_config(
    page_title="MolStorage - Chemical Storage Manager",
    page_icon=None,
    layout="wide"
)

# Session state init
if 'stored_compounds' not in st.session_state:
    st.session_state.stored_compounds = []
if 'compound_input' not in st.session_state:
    st.session_state.compound_input = ''
if 'processed_compounds' not in st.session_state:
    st.session_state.processed_compounds = []
if 'storage_groups' not in st.session_state:
    st.session_state.storage_groups = initialize_storage_groups()
if 'displayed_compounds' not in st.session_state:
    st.session_state.displayed_compounds = []

# CSS Styling for compactness
st.markdown("""
<style>
    .storage-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 10px;
    }
    .compound-button {
        margin: 2px 0;
        width: 100%;
    }
    .state-section {
        margin-top: 8px;
    }
    .main-header {
        color: #FFFFFF;
        text-align: center;
        margin-bottom: 0.5rem;
    }
    .storage-group-header {
        background-color: #000000;
        padding: 5px 8px;
        border-radius: 5px;
        margin-bottom: 6px;
    }
    .compound-list {
        max-height: 180px;
        overflow-y: auto;
        padding-left: 10px;
    }
    .detail-box {
        border: 1px solid #ddd;
        border-radius: 5px;
        padding: 8px 12px;
        margin-bottom: 10px;
        font-size: 0.9rem;
        max-width: 300px;
        overflow-wrap: break-word;
    }
    .remove-button {
        float: right;
        font-weight: bold;
        color: red;
        cursor: pointer;
        border: none;
        background: transparent;
        font-size: 1.1rem;
        line-height: 1rem;
        padding: 0;
    }
    .sidebar-compounds {
        max-height: 250px;
        overflow-y: auto;
        margin-bottom: 10px;
        padding-left: 6px;
    }
    .compound-name {
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)

# Functions
def add_compound():
    compound = st.session_state.compound_input.strip()
    if not compound:
        st.warning("Please enter a compound name.")
        return
    if compound in st.session_state.stored_compounds:
        st.info(f"Compound '{compound}' is already in the list.")
        return
    st.session_state.stored_compounds.append(compound)
    st.session_state.compound_input = ''  # Clear input after successful add
    st.success(f"Added compound: {compound}")

def remove_compound(compound):
    if compound in st.session_state.stored_compounds:
        st.session_state.stored_compounds.remove(compound)
    # Remove also from processed and displayed
    st.session_state.processed_compounds = [c for c in st.session_state.processed_compounds if c['name'] != compound]
    st.session_state.displayed_compounds = [c for c in st.session_state.displayed_compounds if c['name'] != compound]
    # Re-run to update UI
    st.rerun()

def display_pictogram(picto_name):
    """
    Convert a pictogram name to a short text code.
    
    Args:
        picto_name (str): The name of the pictogram
        
    Returns:
        str: A short text code for the pictogram
    """
    # Dictionary mapping pictogram names to short text codes
    picto_codes = {
        "Explosive": "[Explosive]",
        "Flammable": "[Flammable]",
        "Oxidizer": "[Oxidizer]",
        "Compressed Gas": "[Compressed Gas]",
        "Corrosive": "[Corrosive]",
        "Acute Toxic": "[Acute Toxic]",
        "Health Hazard": "[Health Hazard]",
        "Irritant": "[Irritant]",
        "Environmental Hazard": "[Environmental Hazard]"
    }
    
    # Return the text code or unknown
    return picto_codes.get(picto_name, "[?]")

def process_compounds(compounds: list):
    if not compounds:
        st.warning("No compounds to process.")
        return False

    processed_names = {c['name'] for c in st.session_state.processed_compounds}
    new_compounds = [c for c in compounds if c not in processed_names]
    if not new_compounds:
        st.info("No new compounds to process.")
        return True

    results = []
    with st.spinner(f"Processing {len(new_compounds)} new compounds..."):
        for compound in new_compounds:
            try:
                cid, pictos, hazards = get_compound_safety_data(compound)
                if not cid:
                    st.error(f"Could not find data for '{compound}'. Skipping.")
                    continue

                name, iupac, smiles = get_name_and_smiles(cid)
                
                # Extra safety for SMILES - some compounds might have invalid SMILES
                try:
                    acid_base_class = classify_acid_base(name, iupac, smiles, hazards)
                except Exception as e:
                    st.warning(f"Could not classify acid/base for '{compound}'. Setting as unknown. Error: {str(e)}")
                    acid_base_class = "unknown"
                
                # Ensure acid_base_class is a string for compatibility with chemsort_multiple_order_3
                if isinstance(acid_base_class, (list, tuple)):
                    acid_base_class = ", ".join(acid_base_class)
                elif acid_base_class is None:
                    acid_base_class = "unknown"
                
                try:
                    mp_c, bp_c, mp_f, bp_f = get_mp_bp(compound)
                except Exception as e:
                    st.warning(f"Could not retrieve melting/boiling points for '{compound}'. Error: {str(e)}")
                    mp_c, bp_c, mp_f, bp_f = None, None, None, None
                    
                try:
                    state = compound_state(mp_c, bp_c, mp_f, bp_f)
                except Exception as e:
                    st.warning(f"Could not determine state for '{compound}'. Setting as unknown. Error: {str(e)}")
                    state = "state unknown"
                    
                try:
                    sorted_picto = prioritize_pictograms(pictos)
                except Exception as e:
                    st.warning(f"Could not prioritize pictograms for '{compound}'. Error: {str(e)}")
                    sorted_picto = []

                compound_info = {
                    "name": compound,
                    "iupac": iupac,
                    "smiles": smiles,
                    "sorted_pictograms": sorted_picto,
                    "hazard_statements": hazards,
                    "acid_base_class": acid_base_class,  # Now guaranteed to be a string
                    "state_room_temp": state,
                    "melting_point_c": mp_c,
                    "boiling_point_c": bp_c
                }
                results.append(compound_info)
            except Exception as e:
                st.error(f"Error processing '{compound}': {str(e)}")
                continue

    if results:
        st.session_state.processed_compounds.extend(results)
        # Update storage groups with all processed compounds, do NOT reset
        try:
            # Ensure all state keys exist in all storage groups
            for group_name, group_data in st.session_state.storage_groups.items():
                for state_key in ["solid", "liquid", "gas", "state unknown"]:
                    if state_key not in group_data:
                        group_data[state_key] = []
                        
            updated_storage = chemsort_multiple_order_3(st.session_state.processed_compounds, st.session_state.storage_groups)
            st.session_state.storage_groups = updated_storage
        except Exception as e:
            st.error(f"Error during compound sorting: {str(e)}")
            # Since we already added the compounds to processed_compounds, return True
            # to indicate partial success even though sorting failed
            return True
        return True
    else:
        st.error("Failed to process new compounds.")
        return False

def add_to_displayed(compound):
    # Avoid duplicates
    if compound not in st.session_state.displayed_compounds:
        st.session_state.displayed_compounds.append(compound)

def remove_from_displayed(compound_name):
    st.session_state.displayed_compounds = [c for c in st.session_state.displayed_compounds if c['name'] != compound_name]

def display_compound_details(compound):
    # Show minimal compact details in expander
    with st.expander(f"{compound['name']} Details"):
        cols = st.columns(2)
        with cols[0]:
            st.markdown(f"**IUPAC:** {compound.get('iupac', 'N/A')}")
            st.markdown(f"**SMILES:** {compound.get('smiles', 'N/A')}")
            st.markdown(f"**State:** {compound.get('state_room_temp', 'N/A')}")
            mp = compound.get('melting_point_c')
            bp = compound.get('boiling_point_c')
            st.markdown(f"**Melting Point:** {mp if mp is not None else 'N/A'} °C")
            st.markdown(f"**Boiling Point:** {bp if bp is not None else 'N/A'} °C")
        with cols[1]:
            abc = compound.get('acid_base_class', 'N/A')
            # Ensure the acid_base_class is displayed properly
            if isinstance(abc, (list, tuple)):
                abc_str = ", ".join(abc)
            elif abc is None:
                abc_str = "Unknown"
            else:
                abc_str = str(abc)
            st.markdown(f"**Acid/Base Class:** {abc_str}")
            if compound.get('sorted_pictograms', []):
                picto_codes = [display_pictogram(p) for p in compound['sorted_pictograms']]
                st.markdown(f"**Pictograms:** {' '.join(picto_codes)}")
            if compound.get('hazard_statements', []):
                st.markdown("**Hazard Statements:**")
                for h in compound['hazard_statements']:
                    st.markdown(f"- {h}")

# --- Layout ---

st.markdown("<h1 class='main-header'>MolStorage - Chemical Storage Manager</h1>", unsafe_allow_html=True)

with st.sidebar:
    st.header("Input Compounds")
    # Using form to handle input + add compound button better
    with st.form("compound_form", clear_on_submit=True):
        compound_input = st.text_input("Enter Compound Name", key='compound_input')
        submitted = st.form_submit_button("Add Compound")
        if submitted:
            if compound_input.strip() == "":
                st.warning("Please enter a compound name.")
            elif compound_input.strip() in st.session_state.stored_compounds:
                st.info(f"Compound '{compound_input.strip()}' is already in the list.")
            else:
                st.session_state.stored_compounds.append(compound_input.strip())
                st.success(f"Added compound: {compound_input.strip()}")

    st.markdown("### Stored Compounds")
    if st.session_state.stored_compounds:
        for compound in st.session_state.stored_compounds:
            col1, col2 = st.columns([0.85, 0.15])
            with col1:
                st.markdown(f"• {compound}", unsafe_allow_html=True)
            with col2:
                if st.button("❌", key=f"remove_{compound}"):
                    remove_compound(compound)
    else:
        st.write("No compounds added yet.")

    if st.button("Process Compounds"):
        if process_compounds(st.session_state.stored_compounds):
            st.success("Compounds processed successfully!")

    if st.button("Clear All"):
        st.session_state.stored_compounds = []
        st.session_state.processed_compounds = []
        st.session_state.storage_groups = initialize_storage_groups()
        st.session_state.displayed_compounds = []
        st.success("All data cleared!")

if st.session_state.processed_compounds:
    st.header("Storage Groups")
    # Filter out empty groups
    non_empty_groups = {g: s for g, s in st.session_state.storage_groups.items()
                        if any(c for state_compounds in s.values() for c in state_compounds)}
    if not non_empty_groups:
        st.info("No compounds have been sorted into storage groups yet.")
    else:
        cols_per_row = 3
        group_names = list(non_empty_groups.keys())
        rows = np.ceil(len(group_names) / cols_per_row).astype(int)

        for row in range(rows):
            cols = st.columns(cols_per_row)
            for col_idx in range(cols_per_row):
                group_idx = row * cols_per_row + col_idx
                if group_idx < len(group_names):
                    group_name = group_names[group_idx]
                    states = non_empty_groups[group_name]
                    with cols[col_idx]:
                        st.markdown(f"<div class='storage-group-header'><h3>{group_name.replace('_', ' ').title()}</h3></div>", unsafe_allow_html=True)
                        for state_key, compounds in states.items():
                            if compounds:
                                with st.expander(f"{state_key.capitalize()} ({len(compounds)})", expanded=False):
                                    for i, compound in enumerate(compounds):
                                        # Create a unique key with index to avoid duplicates
                                        button_key = f"btn_{group_name}_{state_key}_{compound['name']}_{i}"
                                        
                                        # Get pictograms as emojis
                                        pictos = ""
                                        if "sorted_pictograms" in compound and compound["sorted_pictograms"]:
                                            pictos = "".join(display_pictogram(p) for p in compound["sorted_pictograms"][:2])
                                        
                                        # Create button with name and pictogram emojis
                                        btn_label = f"{pictos} {compound['name']}"
                                        if st.button(btn_label, key=button_key):
                                            add_to_displayed(compound)

# Compound details below storage groups
if st.session_state.displayed_compounds:
    st.header("Selected Compound Details")
    cols = st.columns(len(st.session_state.displayed_compounds))
    for i, compound in enumerate(st.session_state.displayed_compounds):
        with cols[i]:
            # Safely get pictograms
            pictos_str = ""
            if "sorted_pictograms" in compound and compound["sorted_pictograms"]:
                pictos_str = "".join(display_pictogram(p) for p in compound["sorted_pictograms"][:2])
                
            remove_key = f"remove_display_{compound['name']}_{i}"  # Added index for uniqueness

            # Safely get pictograms as emojis
            pictos_str = ""
            if "sorted_pictograms" in compound and compound["sorted_pictograms"]:
                pictos_str = "".join(display_pictogram(p) for p in compound["sorted_pictograms"][:2])
                
            remove_key = f"remove_display_{compound['name']}_{i}"  # Added index for uniqueness

            # Full styled box with compound name and details
            st.markdown(f"""
            <div class='detail-box'>
                <div style='display: flex; justify-content: space-between; align-items: center;'>
                    <p style='font-weight: bold; font-size: 16px; margin: 0;'>{compound['name']} {pictos_str}</p>
                    <form action='' method='post'>
                        <button name='{remove_key}' style='border: none; background: none; font-size: 16px; cursor: pointer;'>❌</button>
                    </form>
                </div>
            """, unsafe_allow_html=True)

            display_compound_details(compound)

            st.markdown("</div>", unsafe_allow_html=True)

            # Handle removal after render (Streamlit limitation workaround)
            if st.session_state.get(remove_key):
                remove_from_displayed(compound['name'])
                st.session_state.pop(remove_key)
                st.rerun()