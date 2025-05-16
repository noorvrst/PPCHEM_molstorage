import streamlit as st
from typing import List, Dict

# Function to get compound safety data (CID, pictograms, hazard statements)
def get_compound_safety_data(compound_name: str):
    # Dummy values for now (you would replace this with actual logic to fetch data)
    return "dummy_cid", ["flammable", "toxic"], ["H225", "H301"]

# Function to get compound name details (IUPAC name, SMILES, etc.)
def get_name_and_smiles(cid: str) -> Dict:
    # Dummy logic (replace with actual fetching of IUPAC, SMILES, etc.)
    if cid == "dummy_cid":
        return {"record_title": "Dummy Record", "iupac_name": "Acetone", "smiles": "CC(=O)C"}
    else:
        # Handle missing or incorrect CID
        return {"record_title": "Unknown", "iupac_name": "Unknown", "smiles": "Unknown"}

# Function to classify the chemical as acid/base
def classify_acid_base(compound_name: str, iupac_name: str, smiles: str, hazard_statements: List[str]) -> List[str]:
    # Simple classification logic (replace with actual logic)
    if "acid" in compound_name.lower():
        return ["acid"]
    elif "base" in compound_name.lower():
        return ["base"]
    else:
        return ["unknown"]

# Function to prioritize hazard pictograms based on severity
def prioritize_pictograms(pictograms: List[str]) -> List[str]:
    # Simple prioritization logic (replace with actual logic)
    return sorted(pictograms)

# Function to process a list of chemical names
def process_chemicals(compounds: List[str]):
    processed_data = []

    for compound_name in compounds:
        # Step 1: Get the PubChem CID, GHS pictograms, and hazard statements
        cid, pictograms, hazard_statements = get_compound_safety_data(compound_name)

        # Step 2: Get the compound's name details (IUPAC, SMILES, etc.)
        name_and_smiles = get_name_and_smiles(cid)
        
        record_title = name_and_smiles.get("record_title", "Unknown")
        iupac_name = name_and_smiles.get("iupac_name", "Unknown")
        smiles = name_and_smiles.get("smiles", "Unknown")

        # Step 3: Classify the compound as acid, base, or unknown
        acid_base_class = classify_acid_base(compound_name, iupac_name, smiles, hazard_statements)

        # Step 4: Prioritize the pictograms
        sorted_pictograms = prioritize_pictograms(pictograms)

        # Step 5: Store the compound data
        compound_data = {
            "name": compound_name,
            "cid": cid,
            "record_title": record_title,
            "iupac_name": iupac_name,
            "smiles": smiles,
            "hazard_statements": hazard_statements,
            "acid_base_class": acid_base_class,
            "sorted_pictograms": sorted_pictograms
        }

        processed_data.append(compound_data)

    return processed_data

# Function to display the results
def display_results(processed_data):
    # Here we assume that `chemsort_multiple_order` is some function that sorts chemicals into storage groups.
    # For now, it's just a dummy function.
    def chemsort_multiple_order(data):
        return {"Group 1": data}  # Dummy group for demonstration

    # Sort and categorize chemicals into storage groups
    storage_groups = chemsort_multiple_order(processed_data)

    st.header("Chemical Safety and Storage Information")

    for group, chemicals in storage_groups.items():
        st.subheader(f"Storage Group: {group}")
        for chem in chemicals:
            st.write(f"**Chemical Name**: {chem['name']}")
            st.write(f"**Record Title**: {chem.get('record_title', 'Not Available')}")
            st.write(f"**IUPAC Name**: {chem.get('iupac_name', 'Not Available')}")
            st.write(f"**SMILES**: {chem.get('smiles', 'Not Available')}")
            st.write(f"**Hazard Statements**: {', '.join(chem['hazard_statements']) if chem.get('hazard_statements') else 'Not Available'}")
            st.write(f"**Acid/Base Class**: {', '.join(chem['acid_base_class']) if chem.get('acid_base_class') else 'Not Available'}")
            st.write(f"**Sorted Pictograms**: {', '.join(chem['sorted_pictograms']) if chem.get('sorted_pictograms') else 'Not Available'}")
            st.write("---")

# Main function to run the app
def main():
    # Example chemical list to process
    compounds = ["Acetone", "Ethanol", "Methanol", "UnknownCompound"]

    # Process the chemicals
    processed_data = process_chemicals(compounds)
    
    # Display the results in the Streamlit app
    display_results(processed_data)

if __name__ == "__main__":
    main()