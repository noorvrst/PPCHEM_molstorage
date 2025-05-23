import re
import json
import requests
import pubchempy as pcp # type: ignore
from rdkit import Chem
from typing import Optional, List, Dict, Any, Union, Tuple
import math

def get_compound_safety_data(compound_name: str, debug: bool = False) -> Tuple[str, List[str], List[str]]:
    """
    Fetches the PubChem CID, GHS pictograms, and hazard statements for a chemical compound.

    Args:
        compound_name (str): The common name of the compound (e.g., "acetone").
        debug (bool): If True, prints debug messages.

    Returns:
        Tuple[str, List[str], List[str]]: A tuple containing:
            - CID as a string (or empty string on failure)
            - List of pictogram names
            - List of hazard statements
    """
    try:
        # Step 1: Get CID (Compound ID)
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/cids/JSON"
        search_response = requests.get(search_url)
        if debug:
            print(f"Search URL: {search_url}")
            print(f"Search Response Status: {search_response.status_code}")
        search_response.raise_for_status()
        search_data = search_response.json()
        if debug:
            print(f"Search Data: {json.dumps(search_data, indent=2)[:500]}...")
        if "IdentifierList" not in search_data or "CID" not in search_data["IdentifierList"]:
            if debug:
                print(f"Compound '{compound_name}' not found.")
            return "", [], []
        cid = str(search_data["IdentifierList"]["CID"][0])
        if debug:
            print(f"Found CID: {cid}")

        # Step 2: Try various safety data headings
        headings = [
            "GHS+Classification",
            "Safety+and+Hazards",
            ""
        ]

        for heading in headings:
            safety_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
            if heading:
                safety_url += f"?heading={heading}"
            if debug:
                print(f"Trying Safety URL: {safety_url}")
            safety_response = requests.get(safety_url)
            safety_response.raise_for_status()
            safety_data = safety_response.json()

            pictograms = set()
            hazard_statements = set()
            hazard_codes = set()

            def process_section(section):
                if "Information" in section:
                    for info in section["Information"]:
                        if info.get("Name") == "Pictogram(s)" and "Value" in info:
                            for markup in info["Value"].get("StringWithMarkup", []):
                                for mark in markup.get("Markup", []):
                                    if "Extra" in mark:
                                        pictograms.add(mark["Extra"].strip())

                        if info.get("Name") == "GHS Hazard Statements" and "Value" in info:
                            for markup in info["Value"].get("StringWithMarkup", []):
                                statement = markup.get("String", "").strip()
                                if statement and "%" in statement:
                                    match = re.match(r"^(H\d{3})", statement)
                                    if match:
                                        code = match.group(1)
                                        if code not in hazard_codes:
                                            hazard_statements.add(statement)
                                            hazard_codes.add(code)

                for subsection in section.get("Section", []):
                    process_section(subsection)

            if "Record" in safety_data and "Section" in safety_data["Record"]:
                for section in safety_data["Record"]["Section"]:
                    process_section(section)

                if pictograms or hazard_statements:
                    return cid, list(pictograms), list(hazard_statements)

        if debug:
            print(f"No safety data found for '{compound_name}'")
        return cid, [], []

    except requests.exceptions.RequestException as e:
        if debug:
            print(f"Network error: {str(e)}")
        return "", [], []
    except json.JSONDecodeError as e:
        if debug:
            print(f"JSON parse error: {str(e)}")
        return "", [], []
    except Exception as e:
        if debug:
            print(f"Unexpected error: {str(e)}")
        return "", [], []
    
    
def get_name_and_smiles(cid: str) -> Tuple[str, str, str]:
    """Returns the Record Title (generic name), IUPAC name, and SMILES from a given CID using PubChemPy.

    Args:
        cid (str): The PubChem Compound ID (CID) as a string.

    Returns:
        tuple[str, str, str]: A tuple containing the Record Title, IUPAC name, and SMILES string.
    """
    if cid is None:
        print("Warning: CID is None")
        return "Unknown", "Unknown", "Unknown"

    try:
        compound = pcp.Compound.from_cid(cid)
        record_title = compound.synonyms[0] if compound.synonyms else "Unknown"
        iupac_name = compound.iupac_name or "Unknown"
        smiles = compound.isomeric_smiles or compound.canonical_smiles or "Unknown"
        return record_title, iupac_name, smiles
    except Exception as e:
        print(f"Error retrieving compound for CID '{cid}': {e}")
        return "Unknown", "Unknown", "Unknown"


def classify_acid_base(name: str, iupac_name: str, smiles: str, ghs_statements: List[str]) -> Union[str, Tuple[str, ...]]:
    """Classify a compound as 'acid', 'base', or 'neutral' based on name, IUPAC, SMILES structure,
    and GHS hazard statements.

    Args:
        name (str): Generic name of the compound.
        iupac_name (str): The IUPAC name of the compound.
        smiles (str): The SMILES string representing the compound's structure.
        ghs_statements (List[str]): GHS hazard statements related to the compound.

    Returns:
        str: A classification such as 'acid', 'base', 'neutral', or 'amphoteric'.
    """
    name = name.lower()
    iupac_name = iupac_name.lower()
    full_name = name + " " + iupac_name
    result = []


    if "acid" in full_name:
        result.append("acid")
    if any(base_word in full_name for base_word in ["hydroxide", "amine", "ammonia", "amide"]):
        result.append("base")

    # 2. Structure-based rules using SMARTS
    acid_smarts = {
        "Carboxylic acid": "[CX3](=O)[OX2H1]",
        "Sulfonic acid": "S(=O)(=O)[OH]",
        "Phosphonic acid": "P(=O)([OH])[OH]",
        "Phenol": "c[OH]",
        "Imide": "[NX3](C(=O))[CX3](=O)",
        "Thiol": "[SH]"
    }

    base_smarts = {
        "Amide": "[NX3][CX3](=O)[#6]",
        "Primary amine": "[NX3;H2][CX4]",
        "Secondary amine": "[NX3;H1][CX4][CX4]",
        "Tertiary amine": "[NX3]([CX4])([CX4])",
        "Aniline": "c1ccc(cc1)[NH2]",
        "Guanidine": "[NX3;H2][C](=[NX3;H1])[NX3;H2]",
        "Imidazole": "c1cnc[nH]1",
        "Pyridine": "n1ccccc1",
        "Alkoxide": "[O-][#6]",
        "Ammonia": "[NX3;H3]"
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"

    found_acid = [label for label, smarts in acid_smarts.items() if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))]
    found_base = [label for label, smarts in base_smarts.items() if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))]

    if found_acid:
        result.append("acid")
    if found_base:
        result.append("base")

    # Final classification
    if "acid" in result and "base" not in result:
        acid_base_class = "acid"
    elif "base" in result and "acid" not in result:
        acid_base_class = "base"
    elif "acid" in result and "base" in result:
        acid_base_class = "amphoteric"
    else:
        acid_base_class = "neutral"

    return acid_base_class


def get_mp_bp(compound_name: str) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Given the name of a chemical compound, fetch its melting and boiling points 
    from PubChem in both Celsius and Fahrenheit, and return their average values.

    Args:
        compound_name (str): The name of the compound (e.g., "water", "ethanol").

    Returns:
        Tuple containing:
            - avg_melting_point_celsius (float or None): Average melting point in °C, if available.
            - avg_boiling_point_celsius (float or None): Average boiling point in °C, if available.
            - avg_melting_point_fahrenheit (float or None): Average melting point in °F, if available.
            - avg_boiling_point_fahrenheit (float or None): Average boiling point in °F, if available.

        Returns None values for any temperatures that could not be extracted.
    """
    
    try:
        compound = pcp.get_compounds(compound_name, 'name')[0]
        cid = compound.cid
    except IndexError:
        return None, None, None, None

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
    except:
        return None, None, None, None

    # Lists to store temperatures in Celsius and Fahrenheit
    mp_celsius = []
    bp_celsius = []
    mp_fahrenheit = []
    bp_fahrenheit = []

    def search_sections(sections_list):
        nonlocal mp_celsius, bp_celsius, mp_fahrenheit, bp_fahrenheit
        for section in sections_list:
            heading = section.get('TOCHeading', '')

            if 'melting point' in heading.lower():
                for info in section.get('Information', []):
                    for markup in info.get('Value', {}).get('StringWithMarkup', []):
                        text = markup.get('String', '')
                        match = re.search(r'[-+]?\d+(?:\.\d+)?', text)
                        if match:
                            value = float(match.group(0))
                            if '°C' in text:
                                mp_celsius.append(value)
                            elif '°F' in text:
                                mp_fahrenheit.append(value)

            elif 'boiling point' in heading.lower():
                for info in section.get('Information', []):
                    for markup in info.get('Value', {}).get('StringWithMarkup', []):
                        text = markup.get('String', '')
                        match = re.search(r'[-+]?\d+(?:\.\d+)?', text)
                        if match:
                            value = float(match.group(0))
                            if '°C' in text:
                                bp_celsius.append(value)
                            elif '°F' in text:
                                bp_fahrenheit.append(value)
                                

            # Recursive descent into subsections if available
            if 'Section' in section:
                search_sections(section['Section'])

    if 'Record' in data and 'Section' in data['Record']:
        search_sections(data['Record']['Section'])

    # Calculate average temperatures in Celsius and Fahrenheit
    avg_mp_celsius = math.ceil(sum(mp_celsius) / len(mp_celsius)) if mp_celsius else None
    avg_bp_celsius = math.ceil(sum(bp_celsius) / len(bp_celsius)) if bp_celsius else None
    avg_mp_fahrenheit = math.ceil(sum(mp_fahrenheit) / len(mp_fahrenheit)) if mp_fahrenheit else None
    avg_bp_fahrenheit = math.ceil(sum(bp_fahrenheit) / len(bp_fahrenheit)) if bp_fahrenheit else None

    return avg_mp_celsius, avg_bp_celsius, avg_mp_fahrenheit, avg_bp_fahrenheit

def compound_state(mp_celsius: Optional[float], bp_celsius: Optional[float], mp_fahrenheit: Optional[float], bp_fahrenheit: Optional[float]) -> str:
    """
    Determine the physical state (solid, liquid, gas) of a compound at room temperature
    based on its melting and boiling points in both Celsius and Fahrenheit.

    The decision is first made using Celsius values, if available. If not, Fahrenheit values
    are used. If neither provides a conclusive result, 'state unknown' is returned.

    Args:
        mp_celsius (Optional[float]): Melting point in degrees Celsius.
        bp_celsius (Optional[float]): Boiling point in degrees Celsius.
        mp_fahrenheit (Optional[float]): Melting point in degrees Fahrenheit.
        bp_fahrenheit (Optional[float]): Boiling point in degrees Fahrenheit.

    Returns:
        str: The physical state at room temperature ("solid", "liquid", "gas", or "unknown").
    """
    room_temp_celsius = 20  # Celsius
    room_temp_fahrenheit = 68  # Fahrenheit

    # Determine state based on average melting and boiling points in Celsius
    if mp_celsius is not None and bp_celsius is not None:
        if room_temp_celsius < mp_celsius:
            state_celsius = 'solid'
        elif mp_celsius <= room_temp_celsius < bp_celsius:
            state_celsius = 'liquid'
        elif room_temp_celsius >= bp_celsius:
            state_celsius = 'gas'
    elif mp_celsius is not None:
        state_celsius = 'solid' if room_temp_celsius < mp_celsius else 'liquid'
    elif bp_celsius is not None:
        state_celsius = 'gas' if room_temp_celsius >= bp_celsius else 'liquid'
    else:
        state_celsius = 'unknown'

        # Determine state based on average melting and boiling points in Fahrenheit
    if mp_fahrenheit is not None and bp_fahrenheit is not None:
        if room_temp_fahrenheit < mp_fahrenheit:
            state_fahrenheit = 'solid'
        elif mp_fahrenheit <= room_temp_fahrenheit < bp_fahrenheit:
            state_fahrenheit = 'liquid'
        elif room_temp_fahrenheit >= bp_fahrenheit:
            state_fahrenheit = 'gas'
    elif mp_fahrenheit is not None:
        state_fahrenheit = 'solid' if room_temp_fahrenheit < mp_fahrenheit else 'liquid'
    elif bp_fahrenheit is not None:
        state_fahrenheit = 'gas' if room_temp_fahrenheit >= bp_fahrenheit else 'liquid'
    else:
        state_fahrenheit = 'unknown'

    # Decide on state based on both Celsius and Fahrenheit data
    if state_celsius != "unknown":
        return state_celsius
    elif state_fahrenheit != "unknown":
        return state_fahrenheit
    else:
        return "unknown"


def prioritize_pictograms(pictogram_list: List[str]) -> List[str]:
    """Sort a list of GHS pictograms based on their hazard priority.

    Priority is assigned based on hazard severity, with lower numbers indicating higher severity.
    Unknown pictograms are deprioritized to the end.

    Args:
        pictogram_list (List[str]): A list of GHS pictogram names.

    Returns:
        List[str]: The input list sorted by predefined hazard priority.
    """
    # Define priority order (lower number = higher priority)
    pictogram_priority = {
    "Explosive": 1,
    "Compressed Gas": 1,
    "Oxidizer": 2,
    "Flammable": 3,
    "Corrosive": 4,
    "Health Hazard": 5,
    "Acute Toxic": 5,  # Same level as Health Hazard
    "Irritant": 6,
    "Environmental Hazard": 6  # Same level as Irritant
}
    # Sort based on priority (fallback to high number if pictogram is unknown)
    return sorted(pictogram_list, key=lambda x: pictogram_priority.get(x, 99))


def is_chemically_compatible(
    existing_pictograms: List[str],
    new_pictograms: List[str],
    existing_acid_base_class: str,
    new_acid_base_class: str,
    group_name: str
) -> bool:
    """
    Determines whether two chemical compounds are chemically compatible for storage
    based on pictograms, acid/base classification and storage group.

    Parameters:
    - existing_pictograms: List of hazard pictograms for the existing compound.
    - new_pictograms: List of hazard pictograms for the compound to be added.
    - existing_acid_base_class: Acid/base classification of the existing compound ('acid', 'base', etc.).
    - new_acid_base_class: Acid/base classification of the new compound.
    - group_name: Name of the storage group (used for applying group-specific rules).

    Returns:
    - True if the compounds are compatible, False otherwise.
    """
    # Rule 1: Acid/base incompatibility
    if ("acid" in existing_acid_base_class and "base" in new_acid_base_class) or \
       ("base" in existing_acid_base_class and "acid" in new_acid_base_class):
        return False

    # Rule 2: Incompatible pictograms
    incompatible_pairs = [
        ("Flammable", "Oxidizer"),
        ("Flammable", "Corrosive"),
        ("Corrosive", "Oxidizer")
    ]
    for pic1 in existing_pictograms:
        for pic2 in new_pictograms:
            if (pic1, pic2) in incompatible_pairs or (pic2, pic1) in incompatible_pairs:
                return False

    # Rule 3: Acid + Corrosive + Acute Toxic or Health Hazard
    if "acid" in existing_acid_base_class and "Corrosive" in existing_pictograms:
        if "Acute Toxic" in new_pictograms or "Health Hazard" in new_pictograms:
            return False
    if "acid" in new_acid_base_class and "Corrosive" in new_pictograms:
        if "Acute Toxic" in existing_pictograms or "Health Hazard" in existing_pictograms:
            return False

    # Rule 5: Group-based restrictions (overrides compound-level checks)
    if group_name:
        group_name = group_name.lower()

        if group_name == "oxidizer":
            if "Flammable" in new_pictograms or "Corrosive" in new_pictograms:
                return False
        if group_name == "flammable" or group_name == "pyrophoric":
            if "Oxidizer" in new_pictograms or "Corrosive" in new_pictograms:
                return False
        if "corrosive" in group_name or "irritant" in group_name:
            if "Oxidizer" in new_pictograms or "Flammable" in new_pictograms:
                return False
        if "toxicity" in group_name or group_name == "cmr_stot":
            if "acid" in new_acid_base_class:
                return False
        if "acid" in group_name:
            if "Health Hazard" in new_pictograms or "Acute Toxic" in new_pictograms:
                return False

    return True

# Keep default_group(), default_group_gas(), and initialize_storage_groups() as before
def default_group():
    return {"solid": [], "liquid": [], "unknown": [], "gas": []}

def default_group_gas():
    return {"gas": [], "solid": [], "liquid": [], "unknown": []}

def initialize_storage_groups() -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
    return {
        "none": default_group(),
        "hazardous_environment": default_group(),
        "acute_toxicity": default_group(),
        "cmr_stot": default_group(),
        "toxicity_2_3": default_group(),
        "acid_corrosive_1": default_group(),
        "acid_irritant": default_group(),
        "base_corrosive_1": default_group(),
        "base_irritant": default_group(),
        "pyrophoric": default_group(),
        "flammable": default_group(),
        "oxidizer": default_group(),
        "explosive": default_group(),
        "compressed_gas": default_group_gas(),
        "nitric_acid": default_group()
    }

def chemsort_multiple_order(compounds: List[Dict[str, Any]], storage_groups: Dict[str, Dict[str, List[Dict[str, Any]]]]) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
    """
    Sorts a list of chemical compounds into appropriate storage groups based on 
    GHS pictograms, hazard statements, physical state, and acid/base classification.

    The function processes compounds by assigning them into predefined groups 
    if they are compatible with it and the compounds contained in the group.
    If a compound is not compatible with any known group, it is either added to a matching 
    custom group or a new custom group is created.

    Parameters:
        compounds (List[Dict[str, Any]]): 
            A list of dictionaries where each dictionary represents a compound with keys like 
            'name', 'sorted_pictograms', 'hazard_statements', 'acid_base_class', and 'state_room_temp'.
        
        storage_groups (Dict[str, Dict[str, List[Dict[str, Any]]]]): 
            A dictionary representing existing storage groups. Each key is a group name, 
            mapping to another dictionary with keys 'solid', 'liquid', and 'gas', each containing 
            a list of compatible compounds.
            Empty at first use of the function.

    Returns:
        Dict[str, Dict[str, List[Dict[str, Any]]]]: 
            Updated dictionary of storage groups with the compounds sorted into the appropriate categories.
    """
    phrases_hazard = [
        "may cause genetic defects", "cancer", "may damage fertility", "causes damage to organs"
    ]
    phrases_flam = [
        "catches fire spontaneously", "in contact with water emits", "may react explosively"
    ]

    pictogram_priority = {
        "Explosive": 1,
        "Oxidizer": 2,
        "Flammable": 3,
        "Corrosive": 4,
        "Acute Toxic": 5,
        "Health Hazard": 5,
        "Irritant": 6,
        "Environmental Hazard": 6,
        "Compressed Gas": 1
    }

    def compound_priority(compound):
        if compound["sorted_pictograms"]:
            return pictogram_priority.get(compound["sorted_pictograms"][0], 100)
        return 100

    compounds = sorted(compounds, key=compound_priority)

    custom_group_counter = 1
    custom_group_prefix = "custom_storage_"
    custom_groups = [key for key in storage_groups if key.startswith(custom_group_prefix)]
    
    def is_compatible_with_group(group_name, compound, group_dict=storage_groups):
        """
        Determines if a compound is compatible with a storage group.
        
        Args:
            group_name (str): Name of the storage group.
            compound (Dict): Dictionary containing compound information.
            group_dict (Dict): Dictionary containing all storage groups.
            
        Returns:
            bool: True if compatible, False otherwise.
        """
        # Get the available state keys for this group
        available_states = group_dict[group_name].keys()
        
        # Check if all available state lists are empty
        if all(not group_dict[group_name][state_key] for state_key in available_states):
            return is_chemically_compatible(
                [],                           # existing_pictograms
                compound["sorted_pictograms"], # new_pictograms
                "",                           # existing_acid_base_class
                compound["acid_base_class"],   # new_acid_base_class
                group_name                    # group_name
            )

        # Otherwise, check compatibility against all existing compounds in all states
        for state_key in available_states:
            compounds_in_group = group_dict[group_name][state_key]
            for existing in compounds_in_group:
                if not is_chemically_compatible(
                    existing["sorted_pictograms"],    # existing_pictograms
                    compound["sorted_pictograms"],    # new_pictograms
                    existing["acid_base_class"],      # existing_acid_base_class
                    compound["acid_base_class"],      # new_acid_base_class
                    group_name                        # group_name
                ):
                    return False
        return True

    for compound in compounds:
        chemical = compound["name"]
        sorted_pictograms = compound["sorted_pictograms"]
        hazard_statements = compound["hazard_statements"]
        acid_base_class = compound["acid_base_class"]
        state = compound["state_room_temp"]
        state_key = 'liquid' if 'liquid' in state else 'solid' if 'solid' in state else 'gas' if "gas" in state else "unknown"

        all_statements = " ".join(hazard_statements).lower()
        sorted_successfully = False

        if sorted_pictograms:
            first_picto = sorted_pictograms[0]

            if chemical.lower() == "nitric acid":
                storage_groups["nitric_acid"][state_key].append(compound)
                sorted_successfully = True

            elif first_picto == "Compressed Gas":
                storage_groups["compressed_gas"][state_key].append(compound)
                sorted_successfully = True

            elif first_picto == "Explosive":
                storage_groups["explosive"][state_key].append(compound)
                sorted_successfully = True

            elif first_picto == "Oxidizer":
                group = "oxidizer"
                if is_compatible_with_group(group, compound):
                    storage_groups["oxidizer"][state_key].append(compound)
                    sorted_successfully = True

            elif first_picto == "Flammable":
                group = "pyrophoric" if any(p in all_statements for p in phrases_flam) else "flammable"
                if is_compatible_with_group(group, compound):
                    storage_groups[group][state_key].append(compound)
                    sorted_successfully = True

            elif first_picto == "Corrosive":
                is_base = "base" in acid_base_class
                is_acid = "acid" in acid_base_class
                is_severe = "causes severe skin burns and eye damage" in all_statements
                group = None
                if is_base:
                    group = "base_corrosive_1" if is_severe else "base_irritant"
                elif is_acid:
                    group = "acid_corrosive_1" if is_severe else "acid_irritant"
                if group and is_compatible_with_group(group, compound):
                    storage_groups[group][state_key].append(compound)
                    sorted_successfully = True

            elif first_picto in ["Acute Toxic", "Health Hazard"]:
                if "fatal" in all_statements or "toxic" in all_statements:
                    group = "acute_toxicity"
                elif any(p in all_statements for p in phrases_hazard):
                    group = "cmr_stot"
                else:
                    group = "toxicity_2_3"
                if is_compatible_with_group(group, compound):
                    storage_groups[group][state_key].append(compound)
                    sorted_successfully = True

            elif first_picto in ["Irritant", "Environmental Hazard"]:
                group = "hazardous_environment" if "toxic to aquatic life" in all_statements else "none"
                if is_compatible_with_group(group, compound):
                    storage_groups[group][state_key].append(compound)
                    sorted_successfully = True

        if not sorted_successfully:
            if not sorted_pictograms:
                group = "none"
                if is_compatible_with_group(group, compound):
                    storage_groups[group][state_key].append(compound)
                    sorted_successfully = True
            else:
                for custom_group in custom_groups:
                    if is_compatible_with_group(custom_group, compound):
                        storage_groups[custom_group][state_key].append(compound)
                        sorted_successfully = True
                        break

        if not sorted_successfully:
            while True:
                new_group_name = f"{custom_group_prefix}{custom_group_counter}"
                if new_group_name not in storage_groups:
                    break
                custom_group_counter += 1
            custom_groups.append(new_group_name)
            storage_groups[new_group_name] = default_group()
            storage_groups[new_group_name][state_key].append(compound)

    return storage_groups