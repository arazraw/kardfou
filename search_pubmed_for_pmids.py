from Bio import Entrez
import time

def search_pubmed(author, start_year, end_year):
    # Set your email here
    Entrez.email = "araz.rawshani@gu.se"
    #Entrez.api_key = "abcdef012"

    # Define search query with date range
    query = f'"{author}"[Author] AND ("{start_year}/01/01"[Date - Publication] : "{end_year}/12/31"[Date - Publication])'

    # Perform the search
    handle = Entrez.esearch(db="pubmed", term=query, retmax=400)
    record = Entrez.read(handle)
    handle.close()

    # Get the list of PMIDs
    pmids = record["IdList"]

    # Introduce a delay of 1 second
    #time.sleep(1)

    return pmids

# Sven-Erik: 


author = "Ravn-Fischer, Annica"
start_year = 1990
end_year = 2023
pmids = search_pubmed(author, start_year, end_year)
print(author)
print(", ".join(pmids))