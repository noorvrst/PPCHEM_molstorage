import streamlit as st
import numpy as np
import time

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

# Set page configuration
st.set_page_config(
    page_title="ChemStorM - Chemical Storage Manager",
    page_icon=None,
    layout="wide",
)

# Initialize session state
if 'compounds_list' not in st.session_state:
    st.session_state.compounds_list = []
if 'storage' not in st.session_state:
    st.session_state.storage = initialize_storage_groups()
if 'compound_data' not in st.session_state:
    st.session_state.compound_data = {}
if 'selected_compounds' not in st.session_state:
    st.session_state.selected_compounds = []

# Functions for managing compounds
def add_compound():
    if st.session_state.new_compound and st.session_state.new_compound not in st.session_state.compounds_list:
        st.session_state.compounds_list.append(st.session_state.new_compound)
        st.session_state.new_compound = ""

def remove_compound(compound_name):
    if compound_name in st.session_state.compounds_list:
        st.session_state.compounds_list.remove(compound_name)

def process_compounds():
    if not st.session_state.compounds_list:
        st.warning("Please add at least one compound to the list")
        return

    progress_bar = st.progress(0)
    status_text = st.empty()

    compounds_to_sort = []
    total = len(st.session_state.compounds_list)

    # Reset storage groups
    st.session_state.storage = initialize_storage_groups()

    # Process each compound
    for i, compound in enumerate(st.session_state.compounds_list):
        status_text.text(f"Processing {compound}...")

        # Step 1: Get safety data
        cid, pictos, hazards = get_compound_safety_data(compound)

        # Step 2: Get name and SMILES
        name, iupac, smiles = get_name_and_smiles(cid)

        # Step 3: Classify acid/base
        acid_base_class = classify_acid_base(name, iupac, smiles, hazards)

        # Step 4: Get melting and boiling points
        mp_c, bp_c, mp_f, bp_f = get_mp_bp(compound)

        # Step 5: Determine state at room temperature
        state = compound_state(mp_c, bp_c, mp_f, bp_f)

        # Step 6: Prioritize pictograms
        sorted_pictos = prioritize_pictograms(pictos)

        # Create compound info dictionary
        compound_info = {
            'name': name if name != "Unknown" else compound,
            'cid': cid,
            'iupac': iupac,
            'smiles': smiles,
            'sorted_pictograms': sorted_pictos,
            'hazard_statements': hazards,
            'acid_base_class': acid_base_class,
            'state_room_temp': state,
            'mp_celsius': mp_c,
            'bp_celsius': bp_c,
            'mp_fahrenheit': mp_f,
            'bp_fahrenheit': bp_f
        }

        # Store compound data
        st.session_state.compound_data[name if name != "Unknown" else compound] = compound_info
        compounds_to_sort.append(compound_info)

        # Update progress
        progress = (i + 1) / total
        progress_bar.progress(progress)
        time.sleep(0.1)

    # Sort compounds into storage groups
    status_text.text("Sorting compounds into storage groups...")
    st.session_state.storage = chemsort_multiple_order_3(compounds_to_sort, st.session_state.storage)

    progress_bar.empty()
    status_text.text("Processing complete!")
    time.sleep(0.5)
    status_text.empty()

def add_to_selected_compounds(compound_name):
    if compound_name not in st.session_state.selected_compounds:
        st.session_state.selected_compounds.append(compound_name)

def remove_from_selected_compounds(compound_name):
    if compound_name in st.session_state.selected_compounds:
        st.session_state.selected_compounds.remove(compound_name)

def clear_selected_compounds():
    st.session_state.selected_compounds = []

# Sidebar
with st.sidebar:
    st.title("ChemStorM - Chemical Storage Manager")
    st.subheader("Add Compounds")

    st.text_input("Enter compound name:", key="new_compound", on_change=add_compound)

    st.subheader("Current Compounds")
    if st.session_state.compounds_list:
        for compound in st.session_state.compounds_list:
            col1, col2 = st.columns([3, 1])
            with col1:
                st.write(compound)
            with col2:
                if st.button("❌", key=f"remove_{compound}"):
                    remove_compound(compound)
    else:
        st.write("No compounds added yet")

    st.button("Store Compounds", on_click=process_compounds, type="primary")

# Main area
st.title("Chemical Storage Groups")

# Display storage groups
tab_overview, tab_details = st.tabs(["Storage Overview", "Compound Details"])

with tab_overview:
    any_compounds_sorted = any(
        len(compounds) > 0
        for group in st.session_state.storage.values()
        for state in group.values()
        for compounds in [state]
    )

    if not any_compounds_sorted:
        st.info("Add compounds and click 'Store Compounds' to see storage groups")
    else:
        col1, col2 = st.columns(2)

        active_groups = []
        for group_name, group_data in st.session_state.storage.items():
            has_compounds = any(len(compounds) > 0 for state, compounds in group_data.items())
            if has_compounds:
                active_groups.append(group_name)

        for i, group_name in enumerate(active_groups):
            display_col = col1 if i % 2 == 0 else col2
            group_data = st.session_state.storage[group_name]

            with display_col:
                st.markdown(
                    f"<h4 style='border:2px solid white; padding: 0.4em; border-radius: 8px;'>{group_name.replace('_', ' ').title()}</h4>",
                    unsafe_allow_html=True
                )

                for state in ["solid", "liquid", "gas", "unknown"]:
                    if state in group_data and group_data[state]:
                        with st.expander(f"{state.title()} Compounds ({len(group_data[state])})"):
                            for compound in group_data[state]:
                                name = compound["name"]
                                col_a, col_b, col_c = st.columns([3, 1, 1])

                                with col_a:
                                    st.write(name)
                                with col_b:
                                    if st.button("Details", key=f"details_{group_name}_{state}_{name}"):
                                        add_to_selected_compounds(name)

with tab_details:
    if not st.session_state.selected_compounds:
        st.info("Select compounds from the Storage Overview tab to view their details here")
    else:
        st.button("Clear selection", on_click=clear_selected_compounds)

        compound_tabs = st.tabs(st.session_state.selected_compounds)

        for i, compound_name in enumerate(st.session_state.selected_compounds):
            with compound_tabs[i]:
                if compound_name in st.session_state.compound_data:
                    data = st.session_state.compound_data[compound_name]

                    col1, col2 = st.columns(2)

                    with col1:
                        st.write(f"**Name:** {data['name']}")
                        st.write(f"**IUPAC:** {data['iupac']}")
                        st.write(f"**CID:** {data['cid']}")
                        st.write(f"**SMILES:** {data['smiles']}")
                        st.write(f"**Physical State:** {data['state_room_temp'].title()}")
                        st.write(f"**Acid/Base Classification:** {data['acid_base_class'].title()}")

                        temp_data = []
                        if data['mp_celsius'] is not None:
                            temp_data.append(f"Melting Point [°C]: {data['mp_celsius']}")
                        if data['bp_celsius'] is not None:
                            temp_data.append(f"Boiling Point [°C]: {data['bp_celsius']}")
                        if data['mp_fahrenheit'] is not None:
                            temp_data.append(f"Melting Point [°F]: {data['mp_fahrenheit']}")
                        if data['bp_fahrenheit'] is not None:
                            temp_data.append(f"Boiling Point [°F]: {data['bp_fahrenheit']}")

                        if temp_data:
                            st.write("**Temperature Data:**")
                            for temp in temp_data:
                                st.write(f"- {temp}")
                        else:
                            st.write("**Temperature Data:** Not available")

                    with col2:

                        if data['sorted_pictograms']:
                            st.write("**GHS Pictograms:**")
                            for picto in data['sorted_pictograms']:
                                st.write(f"- {picto}")
                        else:
                            st.write("**GHS Pictograms:** None found")

                        if data['hazard_statements']:
                            st.write("**Hazard Statements:**")
                            for hazard in data['hazard_statements']:
                                st.write(f"- {hazard}")
                        else:
                            st.write("**Hazard Statements:** None found")
                else:
                    st.error(f"No data available for {compound_name}")