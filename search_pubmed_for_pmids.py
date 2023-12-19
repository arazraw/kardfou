from Bio import Entrez

def search_pubmed(author, start_year, end_year):
    # Set your email here
    Entrez.email = "araz.rawshani@gu.se"

    # Define search query with date range
    query = f'"{author}"[Author] AND ("{start_year}/01/01"[Date - Publication] : "{end_year}/12/31"[Date - Publication])'

    # Perform the search
    handle = Entrez.esearch(db="pubmed", term=query, retmax=1000)
    record = Entrez.read(handle)
    handle.close()

    # Get the list of PMIDs
    pmids = record["IdList"]
    return pmids

# Example usage
author = "Redfors, Björn"
start_year = 1990
end_year = 2023
pmids = search_pubmed(author, start_year, end_year)
print(", ".join(pmids))