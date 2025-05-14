import re
import json
import requests
import urllib.parse
from rdkit import Chem
import pubchempy as pcp
from typing import List, Tuple
from rdkit import Chem


def get_compound_safety_data(compound_name: str, debug: bool = False) -> Tuple[str, List[str], List[str]]:
    """
    Fetches the PubChem CID, GHS pictograms, and hazard statements with percentages for a chemical compound.

    Args:
        compound_name (str): The common name of the compound (e.g., "acetone").
        debug (bool): If True, prints debug messages.

    Returns:
        Tuple[str, List[str], List[str]]: A tuple containing:
            - CID as a string (or empty string on failure)
            - List of pictogram names
            - List of hazard statements with percentages, deduplicated by hazard code
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
    


def get_name_and_smiles(cid: str) -> tuple[str, str, str]:
    """Return the Record Title (generic name), IUPAC name, and SMILES from a given CID using PubChemPy."""
    compound = pcp.Compound.from_cid(cid)

    record_title = compound.synonyms[0] if compound.synonyms else "Unknown"
    iupac_name = compound.iupac_name or "Unknown"
    smiles = compound.isomeric_smiles or compound.canonical_smiles or "Unknown"

    return record_title, iupac_name, smiles

def classify_acid_base(name: str, iupac_name: str, smiles: str, ghs_statements: list[str]) -> str:
    """Classify compound as acid or base based on name, IUPAC, SMILES structure, and GHS hazard statements."""
    name = name.lower()
    iupac_name = iupac_name.lower()

    result = []

    # Checks acid/base indicators in Name, IUPAC name, GHS hazard statements
    full_name = name + " " + iupac_name

    if any("H290" in stmt or "corrosive to metals" in stmt.lower() for stmt in ghs_statements):
        result.append("Unsure (from GHS H290)")

    if "acid" in full_name:
        result.append("acid")

    if any(base_word in full_name for base_word in ["hydroxide", "amine", "ammonia", "amide"]):
        result.append("base")


    # Substructure Matching with SMARTS (RDKit)
    acid_smarts = {
        "Carboxylic acid": "[CX3](=O)[OX2H1]",  # COOH
        "Sulfonic acid": "S(=O)(=O)[OH]",       # SO3H
        "Phenol": "c[OH]",                      # OH on aromatic ring
    }

    base_smarts = {
        "Ammonia":        "[NX3;H3]",            # NH3
        "Amide": "[NX3][CX3](=O)[#6]",
        "Urea-like": "[NX3][CX3](=O)[NX3]",
        "Primary amine": "[NX3;H2][CX4]",        
        "Secondary amine": "[NX3;H1][CX4][CX4]",
        "Tertiary amine": "[NX3]([CX4])([CX4])",
        "Imidazole-like": "n1cncc1",
        "Aniline": "c1ccc(cc1)[NH2]",
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"

    found_acid = [label for label, smarts in acid_smarts.items() if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))]
    found_base = [label for label, smarts in base_smarts.items() if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))]

    if found_acid:
        result.append("acid")
    if found_base:
        result.append("basic")

    if not result:
        return "unknown"


    acid_base_class = tuple(result)
    
    return acid_base_class


def prioritize_pictograms(pictogram_list: list):
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

def is_compatible_picto_ab(existing_pictograms, new_pictograms, existing_acid_base_class, new_acid_base_class):
    """
    Checks if two chemicals (represented by their pictograms and acid/base class)
    are compatible for storage.
    """

    # Rule 1: Acid/base incompatibility
    if ("acid" in existing_acid_base_class and "base" in new_acid_base_class) or \
       ("base" in existing_acid_base_class and "acid" in new_acid_base_class):
        return False

    # Rule 2: Pictogram incompatibilities
    incompatible_pairs = [
        ("Flammable", "Oxidizer"),
        ("Corrosive", "Flammable"),
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

    # All checks passed
    return True

def chemsort_multiple_order(compounds: list):
    """
    Sort multiple compounds into storage groups based on pictograms, hazard statements,
    and acid/base incompatibilities using a compatibility function.
    
    Compounds are first sorted by the priority of their highest priority pictogram
    before assigning to groups. Handles the incompatibility of certain pictograms (e.g., 
    Flammable and Corrosive) and ensures compounds are stored in the correct groups.
    """

    storage_groups = {
        "none": [],
        "hazardous_environment": [],
        "acute_toxicity": [],
        "cmr_stot": [],
        "toxicity_2_3": [],
        "acid_corrosive_1": [],
        "acid_irritant": [],
        "base_corrosive_1": [],
        "base_irritant": [],
        "pyrophoric": [],
        "flammable": [],
        "oxidizer": [],
        "explosive": [],
        "compressed_gas": [],
        "base_corrosive_flammable": [],
        "acid_corrosive_flammable": [],
        "corrosive_flammable": [],
        "oxidizer_flammable": [],
        "base_oxidizer_corrosive": [],
        "acid_oxidizer_corrosive": [],
        "no_category": [],
        "nitric_acid": []
    }

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

    # Check if two compounds' pictograms are incompatible
    def is_incompatible(pictograms_1, pictograms_2):
        incompatible_pairs = [
            ("Flammable", "Corrosive"),  # Flammable and Corrosive are incompatible
            ("Flammable", "Oxidizer"),   # Flammable and Oxidizer are incompatible
            ("Corrosive", "Oxidizer")    # Corrosive and Oxidizer are incompatible
        ]
        # Check if any pair of pictograms is incompatible
        for p1 in pictograms_1:
            for p2 in pictograms_2:
                if (p1, p2) in incompatible_pairs or (p2, p1) in incompatible_pairs:
                    return True
        return False

    def compound_priority(compound):
        if compound["sorted_pictograms"]:
            return pictogram_priority.get(compound["sorted_pictograms"][0], 100)
        return 100

    compounds = sorted(compounds, key=compound_priority)

    for compound in compounds:
        chemical = compound["name"]
        sorted_pictograms = compound["sorted_pictograms"]
        hazard_statements = compound["hazard_statements"]
        acid_base_class = compound["acid_base_class"]

        all_statements = " ".join(hazard_statements).lower()
        sorted_successfully = False

        i = 0
        while not sorted_successfully and i < len(sorted_pictograms):
            pictogram = sorted_pictograms[i]

            def is_compatible_with_group(group_name):
                # Check if the compound can be compatible with a group by comparing their pictograms
                for existing in storage_groups[group_name]:
                    if is_incompatible(existing["pictograms"], [pictogram]):
                        return False  # Incompatible pictograms found
                # Check for acid/base compatibility
                return all(
                    is_compatible_picto_ab(
                        existing["pictograms"],
                        [pictogram],
                        existing["acid_base_class"],
                        acid_base_class
                    )
                    for existing in storage_groups[group_name]
                )


            if chemical.lower() == "nitric acid":
                storage_groups["nitric_acid"].append({
                    "name": chemical,
                    "pictograms": sorted_pictograms,
                    "acid_base_class": acid_base_class
                })
                sorted_successfully = True
                
            elif pictogram == "Compressed Gas":
                storage_groups["compressed_gas"].append({
                    "name": chemical,
                    "pictograms": sorted_pictograms,
                    "acid_base_class": acid_base_class
                })
                sorted_successfully = True

            elif pictogram == "Explosive":
                storage_groups["explosive"].append({
                    "name": chemical,
                    "pictograms": sorted_pictograms,
                    "acid_base_class": acid_base_class
                })
                sorted_successfully = True

            elif pictogram == "Oxidizer":
                if "Flammable" in sorted_pictograms:
                    storage_groups["oxidizer_flammable"].append({
                    "name": chemical,
                    "pictograms": sorted_pictograms,
                    "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True
                    
                elif "Corrosive" in sorted_pictograms and "base" in acid_base_class:
                    storage_groups["base_oxidizer_corrosive"].append({
                    "name": chemical,
                    "pictograms": sorted_pictograms,
                    "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True
                    
                    
                elif "Corrosive" in sorted_pictograms and "acid" in acid_base_class:
                    storage_groups["acid_oxidizer_corrosive"].append({
                    "name": chemical,
                    "pictograms": sorted_pictograms,
                    "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True
                    
                
                elif is_compatible_with_group("oxidizer"):
                    storage_groups["oxidizer"].append({
                    "name": chemical,
                    "pictograms": sorted_pictograms,
                    "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True

############
            elif pictogram == "Flammable":
                # First, check for base_corrosive_flammable or acid_corrosive_flammable
                if "Corrosive" in sorted_pictograms and not "base" in acid_base_class and not "acid" in acid_base_class:
                    storage_groups["corrosive_flammable"].append({
                        "name": chemical,
                        "pictograms": sorted_pictograms,
                        "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True

                    
                elif "base" in acid_base_class and "Corrosive" in sorted_pictograms:
                    storage_groups["base_corrosive_flammable"].append({
                        "name": chemical,
                        "pictograms": sorted_pictograms,
                        "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True
                elif "acid" in acid_base_class:
                    storage_groups["acid_corrosive_flammable"].append({
                        "name": chemical,
                        "pictograms": sorted_pictograms,
                        "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True
                else:
                    # If neither base nor acid, just categorize as flammable
                    group = "pyrophoric" if any(p in all_statements for p in phrases_flam) else "flammable"
                    storage_groups[group].append({
                        "name": chemical,
                        "pictograms": sorted_pictograms,
                        "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True

            elif pictogram == "Corrosive":
                # Now check for corrosive base or acid after flammable categorization
                is_base = "base" in acid_base_class
                is_acid = any("acid" in item.lower() for item in acid_base_class)
                is_severe = "causes severe skin burns and eye damage" in all_statements

                if is_base:
                    # Categorize corrosive base into base_corrosive_1 or base_irritant
                    group = "base_corrosive_1" if is_severe else "base_irritant"
                    if is_compatible_with_group(group):
                        storage_groups[group].append({
                            "name": chemical,
                            "pictograms": sorted_pictograms,
                            "acid_base_class": acid_base_class
                        })
                        sorted_successfully = True

                elif is_acid:
                    # Categorize corrosive acid into acid_corrosive_1 or acid_irritant
                    group = "acid_corrosive_1" if is_severe else "acid_irritant"
                    if is_compatible_with_group(group):
                        storage_groups[group].append({
                            "name": chemical,
                            "pictograms": sorted_pictograms,
                            "acid_base_class": acid_base_class
                        })
                        sorted_successfully = True

            elif pictogram in ["Acute Toxic", "Health Hazard"]:
                if "fatal" in all_statements or "toxic" in all_statements:
                    group = "acute_toxicity"
                elif any(p in all_statements for p in phrases_hazard):
                    group = "cmr_stot"
                else:
                    group = "toxicity_2_3"
                if is_compatible_with_group(group):
                    storage_groups[group].append({
                        "name": chemical,
                        "pictograms": sorted_pictograms,
                        "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True

            elif pictogram in ["Irritant", "Environmental Hazard"]:
                if "toxic to aquatic life with long lasting effects" in all_statements:
                    group = "hazardous_environment"
                elif "harmful" in all_statements:
                    group = "none"
                else:
                    group = "none"
                if is_compatible_with_group(group):
                    storage_groups[group].append({
                        "name": chemical,
                        "pictograms": sorted_pictograms,
                        "acid_base_class": acid_base_class
                    })
                    sorted_successfully = True

            i += 1

        if not sorted_successfully:
            storage_groups["no_category"].append({
                "name": chemical,
                "pictograms": sorted_pictograms,
                "acid_base_class": acid_base_class
            })

    return storage_groups


chemicals = ["acetic acid", 
             "nitric acid", 
             "nitrous oxide", 
             "triethylamine", 
             "sodium hydroxide", 
             "hydrogen peroxide", 
             "acetone", 
             "chlorine", 
             "ethanolamine",
             "ammonia",
             "toluene",
             "xylene",
             "water",
             "glucose"]
chemicals_to_sort = []

for chemical in chemicals:
    cid, pictos, hazard_statements = get_compound_safety_data(chemical)
    name, iupac, smiles = get_name_and_smiles(cid)
    acid_base_class = classify_acid_base(name, iupac, smiles, hazard_statements)
    sorted_pictos = prioritize_pictograms(pictos)
    
    
    chemicals_to_sort.append({
    "name": name,
    "sorted_pictograms": sorted_pictos,
    "hazard_statements": hazard_statements,
    "acid_base_class": acid_base_class
    })
    
      
sorted_compounds = chemsort_multiple_order(chemicals_to_sort)


for group, chemicals in sorted_compounds.items():
    if chemicals:
        formatted_chemicals = "\n  ".join(str(d) for d in chemicals)
        print(f"{group}:\n  {formatted_chemicals}")
