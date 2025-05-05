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

import requests  # (Duplicate import of requests, can be removed)
import urllib.parse  # (Import urllib.parse to encode query strings for URLs)

# Make sure to install all the necessary packages (selenium, requests)
# Gets the cid ("id") of the name or smile of the molecule and then returns its chemical safety information. (2 functions work together)

def get_cid(query: str) -> str:  # (Define function to get CID from a name or SMILES string)
    """Resolve a name or SMILES string to a PubChem CID."""
    is_smiles = any(c in query for c in "=#[]()123456789\\/")  # (Check if query contains SMILES characters)

    if is_smiles:  # (If the query is a SMILES string)
        encoded = urllib.parse.quote(query, safe="")  # (URL-encode the SMILES string)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/cids/TXT"  # (Build URL for SMILES to CID)
    else:  # (If the query is a compound name)
        encoded = urllib.parse.quote(query, safe="")  # (URL-encode the name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/cids/TXT"  # (Build URL for name to CID)

    r = requests.get(url)  # (Send GET request to PubChem)
    r.raise_for_status()  # (Raise error if request failed)
    cid = r.text.strip()  # (Get response text and strip whitespace)
    if not cid:  # (If no CID found)
        raise ValueError(f"No CID found for '{query}'")  # (Raise error indicating no result)
    print(f"[INFO] Found CID: {cid}")  # (Print found CID)
    return cid  # (Return the CID)

def get_compound_name(cid: str) -> str:
    """Fetch the compound name from PubChem given its CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
    r = requests.get(url)
    r.raise_for_status()
    data = r.json()
    return data["PropertyTable"]["Properties"][0]["IUPACName"]


def get_safety_pictograms(cid: str) -> list[str]:  # (Define function to get safety pictograms from PubChem using CID)
    """Scrape GHS pictogram names, filtering out number-like captions like '1-3-0'."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"  # (Construct the compound URL)
    opts = Options()  # (Create Chrome options object)
    opts.headless = True  # (Run Chrome in headless mode, no GUI)

    service = Service(ChromeDriverManager().install())  # (Automatically install and configure ChromeDriver)
    driver = webdriver.Chrome(service=service, options=opts)  # (Launch Chrome browser with given options)

    try:
        driver.get(url)  # (Navigate to the compound page)
        print(f"[INFO] Loading {url}")  # (Print info message)

        WebDriverWait(driver, 20).until(  # (Wait up to 20 seconds for element to load)
            EC.presence_of_element_located((By.CSS_SELECTOR, "section#Safety-and-Hazards"))  # (Wait for Safety section)
        )
        section = driver.find_element(By.CSS_SELECTOR, "section#Safety-and-Hazards")  # (Find the Safety section)

        raw = []  # (Initialize empty list to store captions)
        for el in section.find_elements(By.CSS_SELECTOR, "div.captioned.inline-block"):  # (Iterate over pictogram elements)
            caption = el.get_attribute("data-caption") or ""  # (Get caption text or empty string)
            caption = caption.strip()  # (Remove surrounding whitespace)

            # Skip things that look like '1-3-0' or mostly digits/dashes
            if re.fullmatch(r"[\d\s\-\/]+", caption):  # (Skip if caption is only numbers/dashes)
                continue

            # Skip captions that contain any digits
            if re.search(r"\d", caption):  # (Skip if any number appears in the caption)
                continue

            if caption:  # (If caption is not empty)
                raw.append(caption)  # (Add caption to the list)

        print(f"[INFO] Pictograms: {raw}")  # (Print the extracted pictograms)
        return raw  # (Return the list of pictograms)

    finally:
        driver.quit()  # (Close the browser whether an error occurred or not)

if __name__ == "__main__":
    print("🔍 Enter the name or SMILES of a compound to retrieve its safety information.\n")
    query = input("Enter compound name or SMILES: ").strip()  # (Prompt user for input and strip whitespace)
    try:
        cid = get_cid(query)  # (Get CID from the input query)
        compound_name = get_compound_name(cid)  # (Get the compound's name from PubChem)
        print(f"\n🧪 Molecule: {compound_name}")  # (Print the compound name)

        pics = get_safety_pictograms(cid)  # (Get safety pictograms for the CID)
        print("\n✅ Final safety pictograms:")
        for p in pics:
            print("  -", p)
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