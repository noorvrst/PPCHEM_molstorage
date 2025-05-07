from rdkit import Chem
from rdkit.Chem import Descriptors
import re  # Import the 're' module for regular expressions
import requests  # Import the 'requests' module to make HTTP requests
from selenium import webdriver  # Import webdriver from Selenium to control a web browser
from selenium.webdriver.chrome.service import Service  # Import Service to manage ChromeDriver
from selenium.webdriver.chrome.options import Options  # Import Options to configure Chrome browser settings
from selenium.webdriver.common.by import By  # Import By to specify HTML element locating strategy
from selenium.webdriver.support.ui import WebDriverWait  # Import WebDriverWait to wait for elements to load
from selenium.webdriver.support import expected_conditions as EC  # Import EC for expected conditions in wait
from webdriver_manager.chrome import ChromeDriverManager  # Import ChromeDriverManager to auto-install ChromeDriver
from bs4 import BeautifulSoup
from io import BytesIO
from pdfminer.high_level import extract_text
import urllib.parse  # Import urllib.parse to encode query strings for URLs
import pubchempy as pcp


# Function to get CID from a name or SMILES string
def get_cid(query: str) -> str:
    """Resolve a name or SMILES string to a PubChem CID."""
    is_smiles = any(c in query for c in "=#[]()123456789\\/")  # Check if query contains SMILES characters

    if is_smiles:
        encoded = urllib.parse.quote(query, safe="")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/cids/TXT"
    else:
        encoded = urllib.parse.quote(query, safe="")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/cids/TXT"

    r = requests.get(url)
    r.raise_for_status()
    cid = r.text.strip()
    if not cid:
        raise ValueError(f"No CID found for '{query}'")
    return cid


# Function to get compound name from CID
def get_compound_name(cid: str) -> str:
    """Fetch the compound name from PubChem given its CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
    r = requests.get(url)
    r.raise_for_status()
    data = r.json()
    return data["PropertyTable"]["Properties"][0]["IUPACName"]


# Function to get the SMILES string from PubChem
def get_smiles(cid: str) -> str:
    """Fetch the SMILES string for a given PubChem CID (Compound ID)."""
    try:
        compound = pcp.Compound.from_cid(cid)  # Fetch the compound data using the CID
        return compound.isomeric_smiles  # Get the isomeric SMILES string
    except Exception as e:
        print(f"[ERROR] Could not retrieve SMILES for CID {cid}: {e}")
        return ""


def classify_acid_base(name: str, iupac_name: str, smiles: str, ghs_statements: list[str]) -> str:
    """Classify compound as acid or base based on name, IUPAC, SMILES structure, and GHS hazard statements."""
    name = name.lower()
    iupac_name = iupac_name.lower()
    smiles_upper = smiles.upper()

    result = []

    # --- Checks for Name, IUPAC, GHS Statements ---
    full_name = name + " " + iupac_name

    if any("H290" in stmt or "corrosive to metals" in stmt.lower() for stmt in ghs_statements):
        result.append("Acid/Base (from GHS H290)")

    if "acid" in full_name:
        result.append("Acid (from name)")

    if any(base_word in full_name for base_word in ["hydroxide", "amine", "ammonia", "amide"]):
        result.append("Base (from name)")

    if name.endswith("ide") or iupac_name.endswith("ide"):
        result.append("Possibly base (from suffix 'ide')")

    if any(group in smiles_upper for group in ["COOH", "C(=O)OH", "SO3H"]):
        result.append("Acid (from SMILES text)")

    if any(group in smiles_upper for group in ["NH2", "NH3", "NH", "OH"]) and "COOH" not in smiles_upper:
        result.append("Base (from SMILES text)")

    # --- Substructure Matching with SMARTS (RDKit) ---
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

    from rdkit import Chem

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


# Function to get safety pictograms from PubChem
def get_safety_pictograms(cid: str) -> list[str]:
    """Scrape GHS pictogram names, filtering out number-like captions like '1-3-0'."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    opts = Options()
    opts.headless = True  # Run Chrome in headless mode, no GUI

    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=opts)

    try:
        driver.get(url)

        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.CSS_SELECTOR, "section#Safety-and-Hazards"))
        )
        section = driver.find_element(By.CSS_SELECTOR, "section#Safety-and-Hazards")

        raw = []
        for el in section.find_elements(By.CSS_SELECTOR, "div.captioned.inline-block"):
            caption = el.get_attribute("data-caption") or ""
            caption = caption.strip()

            if re.fullmatch(r"[\d\s\-\/]+", caption):
                continue
            if re.search(r"\d", caption):
                continue

            if caption:
                raw.append(caption)

        return raw

    finally:
        driver.quit()


# Function to get Fisher SDS link
def get_fisher_sds_link(chemical_name):
    search_query = chemical_name.replace(" ", "+")
    search_url = f"https://www.fishersci.com/us/en/catalog/search/sds?selectLang=EN&store=&msdsKeyword={search_query}"
    headers = {"User-Agent": "Mozilla/5.0"}

    response = requests.get(search_url, headers=headers)
    if response.status_code != 200:
        print("Failed to fetch the page.")
        return None

    soup = BeautifulSoup(response.text, "html.parser")
    link_divs = soup.find_all("div", class_="catlog_items")

    for div in link_divs:
        a_tag = div.find("a", href=True)
        if a_tag and "store/msds" in a_tag["href"]:
            full_url = "https://www.fishersci.com" + a_tag["href"]
            return full_url

    return None


# Function to extract specific sections from the SDS PDF
def extract_sds_sections(pdf_url):
    headers = {"User-Agent": "Mozilla/5.0"}
    
    response = requests.get(pdf_url, headers=headers)
    if response.status_code != 200:
        print("Failed to download the PDF.")
        return None, None
    
    pdf_file = BytesIO(response.content)
    text = extract_text(pdf_file)
    
    section7 = extract_section(text, "7. HANDLING AND STORAGE", "8.")
    section10 = extract_section(text, "10. STABILITY AND REACTIVITY", "11.")
    
    return section7, section10


# Helper function to extract a specific section of the SDS
def extract_section(text, start_marker, end_marker):
    start_idx = text.find(start_marker)
    if start_idx == -1:
        return None
    
    end_idx = text.find(end_marker, start_idx)
    if end_idx == -1:
        section_text = text[start_idx:]
    else:
        section_text = text[start_idx:end_idx]
    
    section_text = section_text.replace(start_marker, "").strip()
    section_text = re.sub(r'\n\s*\n', '\n\n', section_text)
    
    return section_text


# Main program
if __name__ == "__main__":
    print("🔍 Enter the name or SMILES of a compound to retrieve its safety and storage information.\n")
    query = input("Enter compound name or SMILES: ").strip()

    try:
        # Step 1: PubChem - Get safety pictograms
        cid = get_cid(query)
        compound_name = get_compound_name(cid)
        print(f"\n🧪 Molecule: {compound_name}")
        pictos = get_safety_pictograms(cid)

        # Fetch IUPAC name and SMILES for classification
        iupac_name = get_compound_name(cid)
        smiles = ''  # Optionally, you can get the SMILES from PubChem API here if needed
        
        # Step 2: Classify as Acid/Base
        classification = classify_acid_base(compound_name, iupac_name, smiles, pictos)
        print("\n⚖️ Acid/Base Classification:", classification)

        print("\n✅ GHS Safety Pictograms:")
        for p in pictos:
            print("  -", p)

        # Step 3: Fisher SDS - Get sections 7 and 10
        print("\n📄 Fetching SDS (Safety Data Sheet)...")
        sds_link = get_fisher_sds_link(compound_name)
        if sds_link:
            print(f"[INFO] SDS PDF found: {sds_link}")
            section7, section10 = extract_sds_sections(sds_link)

            if section7:
                print("\n📦 Section 7 - Handling and Storage:\n")
                print(section7)
            else:
                print("⚠️ Section 7 not found.")

            if section10:
                print("\n🧪 Section 10 - Stability and Reactivity:\n")
                print(section10)
            else:
                print("⚠️ Section 10 not found.")
        else:
            print("❌ SDS PDF link not found.")
    except Exception as e:
        print("❌ Error:", e)


def can_be_stored_together(products: list[tuple[str, list[str], str]]) -> bool:
    """Checks if two products can be stored together based on pictograms, physical state, and acid/base classification (if Oxidizer)."""
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

    # 🚫 Rule: Do not store a solid and a liquid together
    if state1 != state2:
        print(f"[ALERT] Incompatible physical states: {state1} and {state2}.")
        return False

    # Check if either product is an Oxidizer
    if "Oxidizer" in pictos1_set or "Oxidizer" in pictos2_set:
        # Classify both products as Acid or Base only if they are Oxidizers
        classification1 = classify_acid_base(name1, name1, smiles1, pictos1)
        classification2 = classify_acid_base(name2, name2, smiles2, pictos2)

        print(f"[INFO] Classification for {name1}: {classification1}")
        print(f"[INFO] Classification for {name2}: {classification2}")

    # Case: Same product
    if name1_lower == name2_lower:
        if "Explosive" in pictos1_set or "Compressed Gas" in pictos1_set:
            print(f"[ALERT] Explosive or compressed gas must be stored alone.")
            return False
        return True

    # Case: Either one must be stored alone
    for name, pictos in [(name1, pictos1_set), (name2, pictos2_set)]:
        if "Explosive" in pictos or "Compressed Gas" in pictos:
            print(f"[ALERT] {name} contains 'Explosive' or 'Compressed Gas' and must be stored alone.")
            return False

    # Known incompatibility rules
    incompatibles = [
        {"Oxidizer", "Flammable"},
        {"Corrosive", "Flammable"},
    ]

    for rule in incompatibles:
        if any(p in pictos1_set for p in rule) and any(p in pictos2_set for p in rule):
            if not (rule.issubset(pictos1_set) or rule.issubset(pictos2_set)):
                print(f"[CONFLICT] Incompatibility detected between '{name1}' and '{name2}' due to: {rule}")
                return False

    # Handle acid/base incompatibilities more carefully by checking classification
    classification1 = classify_acid_base(name1, name1, smiles1, pictos1)
    classification2 = classify_acid_base(name2, name2, smiles2, pictos2)

    # 🚫 New Rule: Acidic + Corrosive must not be stored with Acute Toxic or Health Hazard
    acidic_corrosive_1 = "Acid" in classification1 and "Corrosive" in pictos1_set
    acidic_corrosive_2 = "Acid" in classification2 and "Corrosive" in pictos2_set
    toxic_or_health_1 = "Acute Toxic" in pictos1_set or "Health Hazard" in pictos1_set
    toxic_or_health_2 = "Acute Toxic" in pictos2_set or "Health Hazard" in pictos2_set
    if (acidic_corrosive_1 and toxic_or_health_2) or (acidic_corrosive_2 and toxic_or_health_1):
        print(f"[CONFLICT] Acidic and corrosive product cannot be stored with acute toxic or health hazard product.")
        return False

    # Prevent unknown classifications from being classified as base or acid in conflict checks
    if "Unknown" in classification1 or "Unknown" in classification2:
        print(f"[INFO] Unable to determine acid/base classification for {name1} or {name2}. No acid/base incompatibility checked.")
        return True

    # If both are acid and base, do not store together
    if "Acid" in classification1 and "Base" in classification2:
        print(f"[CONFLICT] Acidic compound '{name1}' should not be stored with basic '{name2}'.")
        return False
    if "Base" in classification1 and "Acid" in classification2:
        print(f"[CONFLICT] Basic compound '{name1}' should not be stored with acidic '{name2}'.")
        return False
    return True


# MAIN PROGRAM
if __name__ == "__main__":
    print("🔍 Enter two compound names or SMILES strings to check storage compatibility.\n")
    
    products = []

    for i in range(2):
        query = input(f"Enter compound #{i+1} name or SMILES: ").strip()

        try:
            cid = get_cid(query)  # Assuming get_cid() can resolve the query to a CID
            compound_name = get_compound_name(cid)  # Get the compound name from the CID
            pictos = get_safety_pictograms(cid)  # Get the GHS pictograms for safety
            smiles = get_smiles(cid)  # Get the SMILES string for the compound
            
            state = input(f"Enter the physical state of '{compound_name}' (solid or liquid): ").strip().lower()

            # Append all four data points: name, pictograms, state, and SMILES
            products.append((compound_name, pictos, state, smiles))

            # Classify the compound as an acid or base
            classification = classify_acid_base(compound_name, compound_name, "", pictos)
            print(f"[INFO] Classification for {compound_name}: {classification}")
        except Exception as e:
            print(f"[ERROR] Problem with '{query}': {e}")
            continue

    if len(products) == 2:
        print(f"\n[INFO] All safety data collected.")
        if can_be_stored_together(products):
            print("✅ These products can be stored together.")
        else:
            print("❌ These products should NOT be stored together.")
    else:
        print("❌ Could not retrieve safety data for both products.")
