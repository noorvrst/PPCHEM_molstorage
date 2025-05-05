from rdkit import Chem
from rdkit.Chem import Descriptors

import re  # (Import the 're' module for regular expressions)
import requests  # (Import the 'requests' module to make HTTP requests)
from selenium import webdriver  # (Import webdriver from Selenium to control a web browser)
from selenium.webdriver.chrome.service import Service  # (Import Service to manage ChromeDriver)
from selenium.webdriver.chrome.options import Options  # (Import Options to configure Chrome browser settings)
from selenium.webdriver.common.by import By  # (Import By to specify HTML element locating strategy)
from selenium.webdriver.support.ui import WebDriverWait  # (Import WebDriverWait to wait for elements to load)
from selenium.webdriver.support import expected_conditions as EC  # (Import EC for expected conditions in wait)
from webdriver_manager.chrome import ChromeDriverManager  # (Import ChromeDriverManager to auto-install ChromeDriver)

from bs4 import BeautifulSoup
import requests
from io import BytesIO
from pdfminer.high_level import extract_text
import re
import urllib.parse  # (Import urllib.parse to encode query strings for URLs)

# Function to get CID from a name or SMILES string
def get_cid(query: str) -> str:
    """Resolve a name or SMILES string to a PubChem CID."""
    is_smiles = any(c in query for c in "=#[]()123456789\\/")  # (Check if query contains SMILES characters)

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
    print(f"[INFO] Found CID: {cid}")
    return cid

# Function to get compound name from CID
def get_compound_name(cid: str) -> str:
    """Fetch the compound name from PubChem given its CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
    r = requests.get(url)
    r.raise_for_status()
    data = r.json()
    return data["PropertyTable"]["Properties"][0]["IUPACName"]

# Function to get GHS pictograms
def get_safety_pictograms(cid: str) -> list[str]:
    """Scrape GHS pictogram names, filtering out number-like captions like '1-3-0'."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    opts = Options()
    opts.headless = True  # Run Chrome in headless mode, no GUI

    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=opts)

    try:
        driver.get(url)
        print(f"[INFO] Loading {url}")

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

        print(f"[INFO] Pictograms: {raw}")
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

    print("No SDS link found.")
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
    if not section7:
        section7 = extract_section(text, "7. Handling and storage", "8.")
    
    section10 = extract_section(text, "10. STABILITY AND REACTIVITY", "11.")
    if not section10:
        section10 = extract_section(text, "10. Stability and reactivity", "11.")
    
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

        print("\n✅ GHS Safety Pictograms:")
        for p in pictos:
            print("  -", p)

        # Step 2: Fisher SDS - Get sections 7 and 10
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
    
# function to test the compatibility of two molecules
def can_be_stored_together(products: list[tuple[str, list[str]]]) -> bool:
    """Checks if two products can be stored together based on GHS pictogram incompatibilities."""
    if len(products) != 2:
        raise ValueError("Function supports comparison between exactly two products.")

    name1, pictos1 = products[0]
    name2, pictos2 = products[1]

    name1_lower = name1.lower()
    name2_lower = name2.lower()
    pictos1_set = set([p.title() for p in pictos1])
    pictos2_set = set([p.title() for p in pictos2])

    print(f"[INFO] {name1} Pictograms: {pictos1_set}")
    print(f"[INFO] {name2} Pictograms: {pictos2_set}")

    # Case: Same product
    if name1_lower == name2_lower:
        if "Explosive" in pictos1_set or "Compressed Gas" in pictos1_set:
            print(f"[ALERT] Two identical products contain 'Explosive' or 'Compressed Gas'. Must be stored alone.")
            return False
        print("[INFO] Products are identical and safe to store together.")
        return True

    # Case: Either one must be stored alone
    for name, pictos in [(name1, pictos1_set), (name2, pictos2_set)]:
        if "Explosive" in pictos or "Compressed Gas" in pictos:
            print(f"[ALERT] Product {name} contains 'Explosive' or 'Compressed Gas' and must be stored alone.")
            return False

    # Known incompatibility rules
    incompatibles = [
        {"Oxidizer", "Flammable"},
        {"Corrosive", "Flammable"},
        {"Corrosive", "Health Hazard"},
        {"Corrosive", "Acute Toxic"},
    ]

    for rule in incompatibles:
        if any(p in pictos1_set for p in rule) and any(p in pictos2_set for p in rule):
            if not (rule.issubset(pictos1_set) or rule.issubset(pictos2_set)):
                print(f"[CONFLICT] Incompatibility detected between '{name1}' and '{name2}' due to: {rule}")
                return False

    return True
# MAIN PROGRAM
if __name__ == "__main__":
    print("🔍 Enter two compound names or SMILES strings to check compatibility.\n")
    query1 = input("Enter first compound name or SMILES: ").strip()
    query2 = input("Enter second compound name or SMILES: ").strip()

    product_queries = [query1, query2]

    products = []
    for query in product_queries:
        try:
            cid = get_cid(query)
            compound_name = get_compound_name(cid)
            pictos = get_safety_pictograms(cid)
            products.append((compound_name, pictos))
        except Exception as e:
            print(f"[ERROR] Problem with '{query}': {e}")
            continue

    if len(products) == 2:
        print(f"\n[INFO] All pictograms collected.")
        if can_be_stored_together(products):
            print("✅ The products can be stored together.")
        else:
            print("❌ The products should NOT be stored together.")
    else:
        print("❌ Could not retrieve safety data for both products.")
