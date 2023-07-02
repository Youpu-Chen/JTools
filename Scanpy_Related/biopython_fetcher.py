from Bio import Entrez


def readPubmedID(input_name):
    '''
    read and parse the pubmed ID
    '''
    pubmed_ids = []

    with open(input_name) as input:
        for item in input:
            pubmed_ids.append(item.strip())
    return pubmed_ids
    


def get_last_author_and_affiliations(pubmed_id):
    Entrez.email = "otis.hongpu@gmail.com"  # Replace with your email address

    # Fetch the PubMed record
    handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="medline", retmode="text")
    record = handle.read()
    handle.close()

    # Process the record to extract author affiliations and the last author
    affiliations = []
    last_author = None
    record_lines = record.splitlines()
    for line in record_lines:
        if line.startswith("AD  - "):
            affiliation = line[len("AD  - "):].strip()
            affiliations.append(affiliation)
        elif line.startswith("AU  - "):
            author = line[len("AU  - "):].strip()
            last_author = author

    return last_author, affiliations

if __name__ == "__main__":
    # # Test
    # pubmed_id = "19664746"  # Replace with the desired PubMed ID
    # last_author, affiliations = get_last_author_and_affiliations(pubmed_id)
    
    # print("Last Author:")
    # print(last_author)
    
    # print("\nAuthor Affiliations:")
    # for idx, affiliation in enumerate(affiliations, start=1):
    #     if idx > 1:
    #         break
    #     else:
    #         print(f"{idx}. {affiliation}")

    # Run
    # 1. parse the pubmed ID
    pubmed_ids = readPubmedID("./paper_pubmed_id.txt")
    # print(pubmed_ids)

    
    # 2. saving the pubmed info
    ## Note: the last author and affiliation were saved as dict
    output = open("./retrieve_withpubmed_id.txt", "w")

    for pubmed_id in pubmed_ids:
        # retrieve the raw last_author and affiliations
        last_author, affiliations = get_last_author_and_affiliations(pubmed_id)

        # retrieve the affiliation
        a_tag = ""
        for idx, affiliation in enumerate(affiliations, start=1):
            if idx > 1:
                break
            else:
                a_tag = affiliation
        
        print(f"{pubmed_id}: {last_author}, {a_tag}")
        output.write(f"{pubmed_id}\t{last_author}\t{a_tag}\n")
        a_tag = ""
    
    output.close()