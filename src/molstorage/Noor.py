from bs4 import BeautifulSoup
import requests
from io import BytesIO
from pdfminer.high_level import extract_text
import re

def get_fisher_sds_link(chemical_name):
    # Format chemical name for URL
    search_query = chemical_name.replace(" ", "+")
    search_url = f"https://www.fishersci.com/us/en/catalog/search/sds?selectLang=EN&store=&msdsKeyword={search_query}"
    headers = {"User-Agent": "Mozilla/5.0"}

    # Send GET request to Fisher SDS search page
    response = requests.get(search_url, headers=headers)
    if response.status_code != 200:
        print("Failed to fetch the page.")
        return None

    soup = BeautifulSoup(response.text, "html.parser")

    # Look for divs containing SDS links
    link_divs = soup.find_all("div", class_="catlog_items")

    for div in link_divs:
        a_tag = div.find("a", href=True)
        if a_tag and "store/msds" in a_tag["href"]:
            full_url = "https://www.fishersci.com" + a_tag["href"]
            return full_url

    print("No SDS link found.")
    return None

def extract_sds_sections(pdf_url):
    headers = {"User-Agent": "Mozilla/5.0"}
    
    # Download the PDF
    response = requests.get(pdf_url, headers=headers)
    if response.status_code != 200:
        print("Failed to download the PDF.")
        return None, None
    
    # Extract text from PDF
    pdf_file = BytesIO(response.content)
    text = extract_text(pdf_file)
    
    # Find Section 7: Handling and Storage
    section7 = extract_section(text, "7. HANDLING AND STORAGE", "8.")
    if not section7:
        section7 = extract_section(text, "7. Handling and storage", "8.")
    
    # Find Section 10: Stability and Reactivity
    section10 = extract_section(text, "10. STABILITY AND REACTIVITY", "11.")
    if not section10:
        section10 = extract_section(text, "10. Stability and reactivity", "11.")
    
    return section7, section10

def extract_section(text, start_marker, end_marker):
    # Find the start of the section
    start_idx = text.find(start_marker)
    if start_idx == -1:
        return None
    
    # Find the end of the section (next section marker)
    end_idx = text.find(end_marker, start_idx)
    if end_idx == -1:
        # If next section not found, take until end of text
        section_text = text[start_idx:]
    else:
        section_text = text[start_idx:end_idx]
    
    # Clean up the text
    section_text = section_text.replace(start_marker, "").strip()
    section_text = re.sub(r'\n\s*\n', '\n\n', section_text)  # Remove excessive newlines
    
    return section_text

def main():
    chemical = input("Enter chemical name: ")
    sds_link = get_fisher_sds_link(chemical)
    
    if not sds_link:
        print("Could not find SDS link.")
        return
    
    print(f"\nFound SDS at: {sds_link}")
    
    section7, section10 = extract_sds_sections(sds_link)
    
    print("\n=== Section 7: Handling and Storage ===")
    print(section7 if section7 else "Section 7 not found in SDS.")
    
    print("\n=== Section 10: Stability and Reactivity ===")
    print(section10 if section10 else "Section 10 not found in SDS.")

if __name__ == "__main__":
    main()