from Bio import Entrez
import time

def search_pubmed(author, start_year, end_year):
    # Set your email here
    Entrez.email = "araz.rawshani@gu.se"

    # Define search query with date range
    #query = f'"{author}"[Author] AND ("{start_year}/01/01"[Date - Publication] : "{end_year}/12/31"[Date - Publication])'
    query = f'"{author}"[Author]'

    # Perform the search
    handle = Entrez.esearch(db="pubmed", term=query, retmax=500)
    record = Entrez.read(handle)
    handle.close()

    # Get the list of PMIDs
    pmids = record["IdList"]

    # Introduce a delay of 1 second
    #time.sleep(1)

    return pmids

# Sven-Erik: 

author = "Bergh, Niklas"
start_year = 1990
end_year = 2023
pmids = search_pubmed(author, start_year, end_year)
print(author)
print(", ".join(pmids))