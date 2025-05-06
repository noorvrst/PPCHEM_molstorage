import re
import requests
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager
import urllib.parse
from rdkit import Chem
import pubchempy as pcp
from typing import List, Tuple
from rdkit import Chem


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
        # Build URL from CID or compound name
        driver.get(url)
        print(f"[INFO] Loading {url}")

        # Wait for page to fully load
        WebDriverWait(driver, 10).until(lambda d: d.execute_script("return document.readyState") == "complete")
        if "404" in driver.title:
            raise ValueError(f"Invalid CID or compound not found: {cid_name}")

        # Check for Safety and Hazards section
        section = None
        elements = driver.find_elements(By.CSS_SELECTOR, "section#Safety-and-Hazards")
        if elements:
            section = elements[0]
        else:
            return ["No hazard statements"], ["No pictograms"]

        # --- Extract pictograms ---
        pictograms = []
        for el in section.find_elements(By.CSS_SELECTOR, "div.captioned.inline-block"):
            caption = el.get_attribute("data-caption") or ""
            caption = caption.strip()
            if re.fullmatch(r"[\d\s\-\/]+", caption):
                continue
            if any(char.isdigit() for char in caption):
                continue
            if caption:
                pictograms.append(caption)
        if not pictograms:
            pictograms.append("No pictograms")

        # --- Extract hazard statements ---
        statements = []

        def extract_statements(tag_selector: str) -> List[str]:
            found = []
            for tag in section.find_elements(By.CSS_SELECTOR, tag_selector):
                text = tag.text.strip()
                if re.match(r"H\d{3}", text):
                    found.append(text)
            return found

        # Priority: try <p> first
        statements = extract_statements("p")

        # If <p> tags found nothing, try <div>
        if not statements:
            statements = extract_statements("div")

        # If still none found
        if not statements:
            statements = ["No hazard statements"]

        # Remove duplicates, preserve order
        seen = set()
        unique_statements = []
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


def main():
    # List of compound names to test
    compounds = ["acetic acid", "sodium hydroxide", "ammonia", "benzene", "caffeine", "water"]
    compounds_2 = ["ammonia"]

    for compound in compounds:
        print(f"\n{'=' * 40}")
        print(f"🔍 Analyzing compound: {compound}")
        try:
            cid = get_cid(compound)
            name, iupac, smiles = get_name_and_smiles(cid)
            hazards, pictos = get_safety_info(cid)
            acid_base = classify_acid_base(name, iupac, smiles, hazards)

            print(f"CID: {cid}")
            print(f"Name: {name}")
            print(f"IUPAC: {iupac}")
            print(f"SMILES: {smiles}")
            
            print("Pictograms:")
            for pic in pictos:
                print(f"  - {pic}")

            print("Hazard Statements:")
            for h in hazards:
                print(f"  - {h}")

            print(f"Acid/Base Classification: {acid_base}")
            
        except Exception as e:
            print(f"❌ Error processing '{compound}': {e}")

if __name__ == "__main__":
    main()