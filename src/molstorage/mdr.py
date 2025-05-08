import re
import requests
import urllib.parse
from rdkit import Chem
import pubchempy as pcp
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager
from typing import List, Tuple


def get_cid(query: str) -> str:
    """Resolve a name or SMILES to a PubChem CID."""
    encoded = urllib.parse.quote(query, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/cids/TXT"
    r = requests.get(url)
    r.raise_for_status()
    cid = r.text.strip()
    if not cid:
        raise ValueError(f"No CID found for '{query}'")
    print(f"Found CID: {cid}")
    return cid


def get_name_and_smiles(cid: str) -> tuple[str, str, str]:
    """Return the Record Title (generic name), IUPAC name, and SMILES from a given CID using PubChemPy."""
    compound = pcp.Compound.from_cid(cid)
    record_title = compound.synonyms[0] if compound.synonyms else "Unknown"
    iupac_name = compound.iupac_name or "Unknown"
    smiles = compound.isomeric_smiles or compound.canonical_smiles or "Unknown"
    return record_title, iupac_name, smiles


def get_safety_info(cid_name: str) -> Tuple[List[str], List[str]]:
    """Scrape GHS hazard statements and pictograms from PubChem Safety and Hazards section."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid_name}"
    opts = Options()
    opts.add_argument("--headless=new")
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=opts)

    try:
        driver.get(url)
        print(f"[INFO] Loading {url}")
        WebDriverWait(driver, 10).until(lambda d: d.execute_script("return document.readyState") == "complete")
        if "404" in driver.title:
            raise ValueError(f"Invalid CID or compound not found: {cid_name}")

        section = driver.find_elements(By.CSS_SELECTOR, "section#Safety-and-Hazards")
        if not section:
            return ["No hazard statements"], ["No pictograms"]

        pictograms = []
        for el in section[0].find_elements(By.CSS_SELECTOR, "div.captioned.inline-block"):
            caption = el.get_attribute("data-caption") or ""
            caption = caption.strip()
            if re.fullmatch(r"[\d\s\-\/]+", caption) or any(char.isdigit() for char in caption):
                continue
            if caption:
                pictograms.append(caption)
        if not pictograms:
            pictograms.append("No pictograms")

        statements = []
        def extract_statements(tag_selector: str) -> List[str]:
            found = []
            for tag in section[0].find_elements(By.CSS_SELECTOR, tag_selector):
                text = tag.text.strip()
                if re.match(r"H\d{3}", text):
                    found.append(text)
            return found

        statements = extract_statements("p")
        if not statements:
            statements = extract_statements("div")

        if not statements:
            statements = ["No hazard statements"]

        unique_statements = []
        seen = set()
        for s in statements:
            if s not in seen:
                seen.add(s)
                unique_statements.append(s)

        return unique_statements, pictograms

    finally:
        driver.quit()


def classify_acid_base(name: str, iupac_name: str, smiles: str, ghs_statements: list[str]) -> str:
    """Classify compound as acid or base based on name, IUPAC, SMILES structure, and GHS hazard statements."""
    name = name.lower()
    iupac_name = iupac_name.lower()
    smiles_upper = smiles.upper()

    result = []

    # Checks for Name, IUPAC, GHS Statements
    full_name = name + " " + iupac_name
    if any("H290" in stmt or "corrosive to metals" in stmt.lower() for stmt in ghs_statements):
        result.append("Acid/Base (from GHS H290)")

    if "acid" in full_name:
        result.append("Acid (from name)")

    if any(base_word in full_name for base_word in ["hydroxide", "amine", "ammonia", "amide"]):
        result.append("Base (from name)")

    if name.endswith("ide") or iupac_name.endswith("ide"):
        result.append("Possibly base (from suffix 'ide')")

    # Substructure Matching with SMARTS (RDKit)
    acid_smarts = {
        "Carboxylic acid": "[CX3](=O)[OX2H1]",  # COOH
        "Sulfonic acid": "S(=O)(=O)[OH]",       # SO3H
        "Phenol": "c[OH]",                      # OH on aromatic ring
    }

    base_smarts = {
        "Ammonia": "[NX3;H3]",                   # NH3
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
        result.append("Acidic groups: " + ", ".join(found_acid))
    if found_base:
        result.append("Basic groups: " + ", ".join(found_base))

    if not result:
        return "Unknown (no clear acid/base features)"

    return " | ".join(result)


def extract_classification(description: str) -> str:
    """Simplify acid/base classification from descriptive string."""
    desc_lower = description.lower()
    if "acid" in desc_lower and "base" in desc_lower:
        return "Both"
    elif "acid" in desc_lower:
        return "Acid"
    elif "base" in desc_lower:
        return "Base"
    else:
        return "Unknown"
    

from bs4 import BeautifulSoup
import requests
from io import BytesIO
from pdfminer.high_level import extract_text
import re

class SDSLookupError(Exception):
    """Custom exception for SDS lookup failures."""
    pass

def get_fisher_sds_link(chemical_name):
    """
    Retrieve the Safety Data Sheet (SDS) URL for a given chemical from the Fisher Scientific website.

    This function takes the name of a chemical, formats it for use in a Fisher Scientific SDS search,
    sends a GET request to the site, parses the HTML response, and extracts the first SDS link found.

    Args:
        chemical_name (str): The name of the chemical to search for.

    Returns:
        str or None: The full URL to the SDS PDF or page if found; otherwise, None.
        
     Raises:
        ValueError: If the input chemical name is empty.
        requests.exceptions.RequestException: For network-related errors (connection issues, timeouts, etc.).
        SDSLookupError: If the page could not be retrieved or no SDS link was found.
    """
    
    if not chemical_name.strip():
        raise ValueError("Chemical name must not be empty.")

    # Format chemical name for URL
    search_query = chemical_name.replace(" ", "+")
    search_url = f"https://www.fishersci.com/us/en/catalog/search/sds?selectLang=EN&store=&msdsKeyword={search_query}"
    headers = {"User-Agent": "Mozilla/5.0"}

    try:
        response = requests.get(search_url, headers=headers, timeout=10)
        response.raise_for_status()  # Raise an HTTPError for bad status codes
    except requests.exceptions.RequestException as e:
        raise requests.exceptions.RequestException(
            f"Failed to retrieve SDS page: {e}"
        ) from e

    soup = BeautifulSoup(response.text, "html.parser")
    link_divs = soup.find_all("div", class_="catlog_items")

    for div in link_divs:
        a_tag = div.find("a", href=True)
        if a_tag and "store/msds" in a_tag["href"]:
            full_url = "https://www.fishersci.com" + a_tag["href"]
            return full_url

    raise SDSLookupError(f"No SDS link found for chemical: '{chemical_name}'")

class SDSExtractionError(Exception):
    """Custom exception for SDS section extraction failures."""
    pass

def extract_sds_sections(pdf_url):
    """
    Extracts specific sections (7, 10, and 2) from a Safety Data Sheet (SDS) in PDF format.

    This function downloads a PDF from a given URL, extracts text from it, and then processes
    the text to extract the following sections:
    - Section 7: Handling and Storage
    - Section 10: Stability and Reactivity
    - Section 2: Hazard(s) Identification (along with hazard-category pairs)

    Args:
    - pdf_url (str): URL of the PDF file to extract sections from.

    Returns:
    - tuple: A tuple containing the extracted sections:
        - section_7_text (list): List of lines from Section 7.
        - section_10_text (list): List of lines from Section 10.
        - hazard_dict (dict): Dictionary mapping hazards to their categories from Section 2.
        - incompatible_mat (dict): Dictionary with incompatible materials from Section 10.

    Raises:
    - ValueError: If the PDF URL is invalid or empty.
    - requests.exceptions.RequestException: If there is an issue downloading the PDF.
    - SDSExtractionError: If there is an issue extracting text from the PDF.
    """
    
    if not pdf_url or not pdf_url.startswith("http"):
        raise ValueError("Invalid or empty PDF URL provided.")

    headers = {"User-Agent": "Mozilla/5.0"}
    
    # Step 1: Download pdf
    try:
        print(f"[INFO] Downloading PDF from: {pdf_url}")
        response = requests.get(pdf_url, headers=headers, timeout=15)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise requests.exceptions.RequestException(f"Failed to download PDF: {e}") from e

    # Step 2: Extract text from pdf
    try:
        pdf_file = BytesIO(response.content)
        text = extract_text(pdf_file)
    except Exception as e:
        raise SDSExtractionError(f"Failed to extract text from PDF: {e}") from e

    # Step 3: Process lines
    full_text = extract_text(pdf_file)
    lines = full_text.splitlines()
    
    # Extract section 7: Handling and storage
    section_7_text = []
    section_7 = False
    
    for line in lines:
        line_clean = line.strip().lower()

        # Start capturing at Section 7
        if "7. handling and storage" in line_clean:
            section_7 = True
            section_7_text.append(line)
            continue

        # Stop capturing at Section 8
        if section_7 and ("8. exposure controls / personal protection" in line_clean):
            break

        if section_7:
            section_7_text.append(line) # Every line in section 7 gets appended


    # Extract section 10: Stability and reactivity
    section_10_text = []
    section_10 = False
    for line in lines:
        line_clean = line.strip().lower()

        # Start capturing at Section 10
        if "10. stability and reactivity" in line_clean:
            section_10 = True
            section_10_text.append(line)
            continue

        # Stop capturing at Section 11
        if section_10 and ("11. toxicological information" in line_clean):
            break

        if section_10:
            section_10_text.append(line) # Every line in section 10 gets appended
         
         
    # Extract section 2 and format dangers:categories
    # Step 3: Process lines
    section_2_text = []
    section_2 = False

    for line in lines:
        line_clean = line.strip().lower()

        if "2. hazard" in line_clean:  # allow partial matches
            section_2 = True
            section_2_text.append(line)
            continue

        if section_2 and ("label elements" in line_clean):
            break

        if section_2:
            section_2_text.append(line)

    if not section_2_text:
        raise SDSExtractionError("Section 2 not found in PDF.")
    
       # Find the start of the hazard section
    empty_line_idx = [i for i, line in enumerate(section_2_text) if not line]
    hazard_start = empty_line_idx[0] + 1

   # Dynamically find the index of the first "Category" line
    category_start = None
    for i in range(hazard_start, len(section_2_text)):
        if section_2_text[i].strip().startswith("Category"):
            category_start = i
            break

    if category_start is None:
        return("No categories found in section 2.")

    # Find hazard_end as the last non-empty line before the first "Category"
    hazard_end = category_start - 1
    while hazard_end > hazard_start and not section_2_text[hazard_end].strip():
        hazard_end -= 1

    # Extract hazards
    hazards = section_2_text[hazard_start:hazard_end + 1]

    # Extract categories
    categories = [line for line in section_2_text[category_start:] if line.strip().startswith("Category")] 
    

    # Map hazards to categories (assuming they're in the same order, which they should be)
    hazard_dict = {}

    i = 0
    j = 0
    while i < len(categories) and j < len(hazards):
        hazard = hazards[j]
    
        # Check for specific target organ toxicity (single exposure) or (repeated exposure)
        if "specific target organ toxicity" in hazard.lower():
            target_organs = hazards[j + 1]
            combined_hazard = f"{hazard.strip()} {target_organs.strip()}"
            hazard_dict[combined_hazard] = categories[i]
            i += 1
            j += 2  # Skip the next line as it's already paired with the current one

        else:
           # For regular hazards, just add the hazard with the category
           hazard_dict[hazard] = categories[i]
           i += 1
           j += 1
           
           
        # Get the incompatible materials from section 10 (usually always here)
        key = "Incompatible Materials"
        for i, line in enumerate(section_10_text):
            if line.strip().lower() == key.lower():
                next_line = section_10_text[i + 2].strip()
                
                # Split into a list
                if next_line:
                    materials = [m.strip().capitalize() for m in re.split(r',\s*', next_line)]
                    incompatible_mat = {key: materials}
                    
                
            
    return section_7_text, section_10_text, hazard_dict, incompatible_mat

def main():
    # Sample PDF URL (replace this with a valid PDF URL to test)
    pdf_url = "https://www.fishersci.com/store/msds?partNumber=S181100&productDescription=silver-nitrate-cert-acs-g&vendorId=VN00033897&keyword=true&countryCode=US&language=en"
    try:
        section_7_text, section_10_text, hazard_dict, incompatible_mat = extract_sds_sections(pdf_url)

        # Output extracted sections
        print("\n[INFO] Extracted Section 7")
        for line in section_7_text:
            print(line)

        print("\n[INFO] Extracted Section 10")
        for line in section_10_text:
            print(line)

        print("\n[INFO] Extracted Section 2")
        for hazard, category in hazard_dict.items():
            print(f"{hazard}: {category}")
            print()
            
        for key, materials in incompatible_mat.items():
            print(f"{key}:")
            for material in materials:
                print(f"- {material}")
    
    except ValueError as e:
        print(f"[ERROR] {e}")
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] {e}")
    except SDSExtractionError as e:
        print(f"[ERROR] {e}")
    except Exception as e:
        print(f"[ERROR] An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()



def can_be_stored_together(products: list[tuple[str, list[str], str, str]]) -> bool:
    """Checks if two products can be stored together based on pictograms, physical state, and acid/base classification."""
    if len(products) != 2:
        raise ValueError("Function supports comparison between exactly two products.")

    name1, pictos1, state1, smiles1 = products[0]
    name2, pictos2, state2, smiles2 = products[1]

    name1_lower = name1.lower()
    name2_lower = name2.lower()
    pictos1_set = set([p.title() for p in pictos1])
    pictos2_set = set([p.title() for p in pictos2])

    print(f"[INFO] {name1} ({state1}) Pictograms: {pictos1_set}")
    print(f"[INFO] {name2} ({state2}) Pictograms: {pictos2_set}")

    if state1 != state2:
        print(f"[ALERT] Incompatible physical states: {state1} and {state2}.")
        return False

    if name1_lower == name2_lower:
        if "Explosive" in pictos1_set or "Compressed Gas" in pictos1_set:
            print(f"[ALERT] Explosive or compressed gas must be stored alone.")
            return False
        return True

    for name, pictos in [(name1, pictos1_set), (name2, pictos2_set)]:
        if "Explosive" in pictos or "Compressed Gas" in pictos:
            print(f"[ALERT] {name} contains 'Explosive' or 'Compressed Gas' and must be stored alone.")
            return False

    incompatibles = [
        {"Oxidizer", "Flammable"},
        {"Corrosive", "Flammable"},
    ]

    for rule in incompatibles:
        if any(p in pictos1_set for p in rule) and any(p in pictos2_set for p in rule):
            if rule.issubset(pictos1_set) or rule.issubset(pictos2_set):
                print(f"[CONFLICT] Incompatibility detected between '{name1}' and '{name2}' due to: {rule}")
                return False

    # Acid/Base Classification using descriptive output
    desc1 = classify_acid_base(name1, name1, smiles1, pictos1)
    desc2 = classify_acid_base(name2, name2, smiles2, pictos2)

    print(f"[INFO] Acid/Base classification for {name1}: {desc1}")
    print(f"[INFO] Acid/Base classification for {name2}: {desc2}")

    class1 = extract_classification(desc1)
    class2 = extract_classification(desc2)

    acidic_corrosive_1 = class1 == "Acid" and "Corrosive" in pictos1_set
    acidic_corrosive_2 = class2 == "Acid" and "Corrosive" in pictos2_set
    toxic_or_health_1 = "Acute Toxic" in pictos1_set or "Health Hazard" in pictos1_set
    toxic_or_health_2 = "Acute Toxic" in pictos2_set or "Health Hazard" in pictos2_set

    # Explicit check for corrosive and flammable conflict
    if "Corrosive" in pictos1_set and "Flammable" in pictos2_set:
        print(f"[CONFLICT] '{name1}' (Corrosive) should not be stored with '{name2}' (Flammable).")
        return False
    if "Flammable" in pictos1_set and "Corrosive" in pictos2_set:
        print(f"[CONFLICT] '{name1}' (Flammable) should not be stored with '{name2}' (Corrosive).")
        return False

    # Handle the case where one or both compounds are classified as "Unknown"
    if class1 == "Unknown" or class2 == "Unknown":
        print(f"[INFO] One or both compounds have an 'Unknown' classification, but safety hazards are still checked.")
        return True

    if class1 == "Acid" and class2 == "Base":
        print(f"[CONFLICT] Acidic '{name1}' should not be stored with basic '{name2}'.")
        return False
    if class1 == "Base" and class2 == "Acid":
        print(f"[CONFLICT] Basic '{name1}' should not be stored with acidic '{name2}'.")
        return False

    return True


def main():
    print("🔍 Enter two compound names or SMILES strings to check compatibility.\n")
    
    # Get user input
    query1 = input("Enter first compound name or SMILES: ").strip()
    query2 = input("Enter second compound name or SMILES: ").strip()

    product_queries = [query1, query2]
    products = []

    # Process each product query
    for query in product_queries:
        try:
            cid = get_cid(query)  # Get the PubChem CID for the compound
            name, iupac, smiles = get_name_and_smiles(cid)  # Get the name, IUPAC, and SMILES string
            statements, pictos = get_safety_info(cid)  # Get safety info (hazard statements, pictograms)
            
            # Now ask for the physical state (outside of the try-except block)
            state = input(f"Enter the physical state of '{name}' (solid or liquid): ").strip().lower()
            
            # Append all four data points: name, pictograms, state, and SMILES
            products.append((name, pictos, state, smiles))

            # Classify the compound as an acid or base
            classification = classify_acid_base(name, name, "", pictos)
            print(f"[INFO] Classification for {name}: {classification}")
        
        except Exception as e:
            print(f"[ERROR] Problem with '{query}': {e}")
            continue  # Skip to the next product in case of error

    # If both products are collected, check compatibility
    if len(products) == 2:
        print(f"\n[INFO] All safety data collected.")
        if can_be_stored_together(products):
            print("✅ These products can be stored together.")
        else:
            print("❌ These products should NOT be stored together.")
    else:
        print("❌ Could not retrieve safety data for both products.")

if __name__ == "__main__":
    main()