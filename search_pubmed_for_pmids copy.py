from Bio import Entrez
#import time

def search_pubmed(author):
    # Set your email here
    Entrez.email = "araz.rawshani@gu.se"

    # Define search query
    query = f'"{author}"[Author]'

    # Perform the search
    handle = Entrez.esearch(db="pubmed", term=query, retmax=500)
    record = Entrez.read(handle)
    handle.close()

    # Get the list of PMIDs
    pmids = record["IdList"]

    # Introduce a delay of 1 second if needed
    #time.sleep(1)

    return pmids

# List of authors to search for
# JQ Smith
authors = [
    "Wideqvist, Maria",
    "Abu Al Chay, Moner",
    "Frisk Torell, Matilda",
    "Novo, Mirza",
    "Jha, Sandeep",
    "Valeljung, Inger",
    "Mellberg, Tomas",
    "Björkenstam, Marie",
    "Bartfay, Sven Erik",
    "Backelin, Charlotte",
    "Martinsson, Andreas",
    "Taha, Amar",
    "Bene, Orsolya",
    "Vahedi, Farzad",
    "Rubulis, Aigars",
    "Herczku, Csaba",
    "Bollano, Entela",
    "Bergfeldt, Lennart",
    "Sigurjonsdottir, Runa",
    "Lundgren, Peter",
    "Smith, Gustav",
    "Hjalmarsson, Clara",
    "Andersson, Bert",
    "Bergh, Niklas",
    "Omerovic, Elmir",
    "Redfors, Björn",
    "Dworeck, Christian",
    "Odell, Annika",
    "Andréen, Sofie",
    "Ljungman, Charlotta",
    "Dahlberg, Pia",
    "Axelsson, Karl-Jonas",
    "Angerås, Oskar",
    "Råmunddal, Truls",
    "Odenstedt, Jacob",
    "Romeo, Stefano",
    "Pirazzi, Carlo",
    "Bobbio, Emanuele",
    "Karason, Kristjan",
    "Bentzel, Sara",
    "Ravn-Fischer, Annica",
    "Rawshani, Araz",
    "Skoglund, Kristofer",
    "Lachonius, Maria"
]


for author in authors:
    pmids = search_pubmed(author)
    print(f"Author: {author}")
    if pmids:
        print(", ".join(pmids))
    else:
        print("No publications found.")
    print("-" * 20)  # Separator for readability
