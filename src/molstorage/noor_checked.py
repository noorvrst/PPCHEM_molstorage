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
