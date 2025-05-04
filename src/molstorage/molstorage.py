import re
import requests
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

import requests
import urllib.parse

# Make sure to install all the necessary packages (selenium, requests)
# Gets the cid ("id") of the name or smile of the molecule and then returns its chemical safety information. (2 functions work together)

def get_cid(query: str) -> str:
    """Resolve a name or SMILES string to a PubChem CID."""
    # Heuristic: if the string contains structural characters, treat as SMILES
    is_smiles = any(c in query for c in "=#[]()123456789\\/")

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


def get_safety_pictograms(cid: str) -> list[str]:
    """Scrape GHS pictogram names, filtering out number-like captions like '1-3-0'."""
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

            # Skip things that look like '1-3-0' or mostly digits/dashes
            if re.fullmatch(r"[\d\s\-\/]+", caption):
                continue

            if caption:
                raw.append(caption)

        print(f"[INFO] Pictograms: {raw}")
        return raw

    finally:
        driver.quit()

if __name__ == "__main__":
    query = input("Enter compound name or SMILES: ").strip()
    try:
        cid = get_cid(query)
        pics = get_safety_pictograms(cid)
        print("\n✅ Final safety pictograms:")
        for p in pics:
            print("  -", p)
    except Exception as e:
        print("❌ Error:", e)