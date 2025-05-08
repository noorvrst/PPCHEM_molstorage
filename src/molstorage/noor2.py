from noor_checked import *

from noor_checked import extract_all_sections

section_7_text, section_10_text, hazard_dict = extract_all_sections("https://www.fishersci.com/store/msds?partNumber=D3720&productDescription=methylene-chlor-cert-acs-l&vendorId=VN00033897&keyword=true&countryCode=US&language=en")

def analyse_section(word,text)->None:
    # searches a text for a specific word
    total_text = " ".join(text)
    sentences = total_text.split(".")
    print("\n".join(sentences))
    for sentence in sentences:
        if word in sentence:
            print(sentence)
    return None

def extract_storage_specifications(text):
    dictio = handling_storage_section_7(text)
    for title, info in dictio.items():
        if title == "Storage":
            print(title)
            for sentence in info:
                print(f"- {sentence}")

extract_storage_specifications(section_7_text)

