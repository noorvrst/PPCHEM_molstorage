"""from noor_checked import *

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

"""

from Noor import *

titles7:list[str] = ["Handling","Storage"]
titles10:list[str] = ["Reactive Hazard","Stability","Conditions to Avoid","Incompatible Materials","Hazardous Decomposition Products","Hazardous Polymerization","Hazardous Reactions"]

def separate_section10_per_title(section,titles)->dict:
    dictio = {}
    for i,title in enumerate(titles):
        if title != titles[len(titles)-1]:
            if title in section:
                start_of_info_index = section.index(title)+len(title)
                end_of_info_index = section.index(titles10[i + 1])
                info = section[start_of_info_index:end_of_info_index]
                dictio[title] = info.replace("_","").replace("\n","")
        else:
            if title in section:
                start_of_info_index = section.index(title)+len(title)
                info = section[start_of_info_index:]
                dictio[title] = info.replace("_","").replace("\n","")
    return dictio

def main()->None:
    chemical1 = input("Enter first chemical name: ")
    chemical2 = input("Enter second chemical name: ")
    sds_link1= get_fisher_sds_link(chemical1)
    sds_link2= get_fisher_sds_link(chemical2)
    section7, section10_1 = extract_sds_sections(sds_link1)    
    section7, section10_2 = extract_sds_sections(sds_link2)  
    dictio_1=separate_section10_per_title(section10_1,titles10)
    dictio_2=separate_section10_per_title(section10_2,titles10)
    dictio_12 = {}
    for item in dictio_1:
        dictio_12[dictio_1[item]]=dictio_2[item]
    for key,value in dictio_12.items():
        print(f"\n{chemical1}: {key} \n{chemical2}: {value}\n")

if __name__=="__main__":
    main()
