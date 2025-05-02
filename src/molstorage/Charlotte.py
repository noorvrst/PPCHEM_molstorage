import re
import requests
import urllib.parse
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

# Function to get the CID of the compound
def get_cid(query: str) -> str:
    """Resolve a name or SMILES string to a PubChem CID."""
    is_smiles = any(c in query for c in "=#[]()123456789\\/")
    encoded = urllib.parse.quote(query, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{'smiles' if is_smiles else 'name'}/{encoded}/cids/TXT"
    r = requests.get(url)
    r.raise_for_status()
    cid = r.text.strip()
    if not cid:
        raise ValueError(f"No CID found for '{query}'")
    print(f"[INFO] Found CID for {query}: {cid}")
    return cid

# Function to get the safety pictograms
def get_safety_pictograms(cid: str) -> list[str]:
    """Retrieve the GHS pictograms, filtering out numeric or invalid captions."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    opts = Options()
    opts.headless = True
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

            if re.fullmatch(r"[\d\s\-\/A-Z]+", caption) and not re.search(r"[a-z]", caption):
                continue

            if caption:
                raw.append(caption)

        print(f"[INFO] Pictograms extracted: {raw}")
        return raw

    finally:
        driver.quit()

# Function to check incompatibility between two products and include product names
def can_be_stored_together(products: list[tuple[str, list[str]]]) -> bool:
    """Checks if two products can be stored together based on GHS pictogram incompatibilities."""
    if len(products) != 2:
        raise ValueError("Function supports comparison between exactly two products.")

    name1, pictos1 = products[0]
    name2, pictos2 = products[1]

    # Rule: Explosive or Gas Cylinder must be stored alone
    for name, pictos in products:
        if "Explosive" in pictos or "Gas Cylinder" in pictos:
            print(f"[ALERT] Product {name} contains 'Explosive' or 'Gas Cylinder' and must be stored alone.")
            return False

    pictos1_set = set(pictos1)
    pictos2_set = set(pictos2)
    print(f"[INFO] {name1} Pictograms: {pictos1_set}")
    print(f"[INFO] {name2} Pictograms: {pictos2_set}")

    # Incompatibility rules (between products only)
    incompatibles = [
        {"Oxidizer", "Flammable"},
        {"Corrosive", "Flammable"},
        {"Corrosive", "Acute Toxic"},
        {"Corrosive", "Health Hazard"},
        {"Acute Toxic", "Health Hazard"},
        {"Corrosive", "Toxic"},
        {"Oxidizer", "Toxic"},
        {"Health Hazard", "Flame"},
        {"Toxic", "Health Hazard"},
    ]

    for rule in incompatibles:
        if any(p in pictos1_set for p in rule) and any(p in pictos2_set for p in rule):
            if not (rule.issubset(pictos1_set) or rule.issubset(pictos2_set)):
                print(f"[CONFLICT] Incompatibility detected between '{name1}' and '{name2}' due to: {rule}")
                return False

    return True

# MAIN PROGRAM
if __name__ == "__main__":
    product_names = ["trinitrotoluene", "glucose"]  # Change as needed

    products = []
    for name in product_names:
        try:
            cid = get_cid(name)
            pictos = get_safety_pictograms(cid)
            products.append((name, pictos))
        except Exception as e:
            print(f"[ERROR] Problem with '{name}': {e}")
            continue

    print(f"[INFO] All Pictograms collected: {[p[1] for p in products]}")
    if can_be_stored_together(products):
        print("✅ The products can be stored together.")
    else:
        print("❌ The products should NOT be stored together.")