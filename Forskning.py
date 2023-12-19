# Libraries ========================================================================================

from requests import get
import re

import streamlit as st
import sqlite3
from pymed import PubMed
import pandas as pd
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, DataReturnMode
import base64

import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx


# Settings ========================================================================================
st.set_page_config(layout="wide",
                   page_title="VO Kardiologi - FOU & ALF Rapport")
st.title ("FOU & ALF Rapport från VO Kardiologi")
st.subheader("Kardiologen SU/Sahlgrenska")
st.write('\n')
st.write('Välkommen till Kardiologens FOU & ALF rapport. Här kan du registrera artiklar och se statistik över dessa. Du kan även se en graf över hur författarna samarbetar med varandra internet och externt. Du kan även registrera utbildnings- och handledningsaktiviteter.')


# Colors
purple_color = 'rgb(132, 112, 255)' #8470FF
yellow_color = 'rgb(255, 215, 0)' #FFD700
dark_yellow_color = 'rgb(255, 215, 0)'
light_blue_color = '#90D2DC'


# Custom CSS to hide the GitHub icon
hide_menu_style = """
        <style>
        #MainMenu {visibility: hidden;}
        #GithubIcon {visibility: hidden;}
        footer {visibility: hidden;}
        footer {visibility: hidden;}
        header {visibility: hidden;}
        </style>
        """
st.markdown(hide_menu_style, unsafe_allow_html=True)


# List of substrings to include in some plots
included_substrings = ["Araz Rawshani", "Dan Ioanes", "Redfors", "Hjalmarsson", "Skoglund", "Truls", "Bergh", "Ljungman",
                       "Bergfeldt", "Bert Andersson", "Bollano", "Omerovic", "Petur Petursson", "Oskar Angerås", "Oscar Angerås", "Oskar Angeras", "Oscar Angeras",
                       "Sebastian", "Dworeck", "Sigurjonsdottir", "Lachonius", "Sara Bentzel", "Kristjan Karason", "Emanuele Bobbio",
                       "Stefano Romeo", "Jakob Odensted", "Karl-Jonas Axelsson", "Pia Dahlberg", "Antros Louca", "Mohammed Munir",
                       "Charlotta Ljungman", "Sofie Andréen", "Annika Odell", "Gustav Smith"]


# Functions ============================================================================================================

# Fetch citation data from OpenCitations
api_token = "e7dce73f-7ffd-48a6-95c9-7b1403e07d0a"

#def fetch_citation_data(doi):
#    api_call = f"https://opencitations.net/index/coci/api/v1/citations/{doi}"
#    headers = {"Authorization": f"Bearer {api_token}"}
#    response = get(api_call)
#
#    if response.status_code == 200:
#        data = response.json()
#        # Filter out items where the DOI is the same as the queried DOI
#        citation_dois = [item.get('citing', 'N/A') for item in data if item.get('citing') != doi]
#        return len(citation_dois), ", ".join(citation_dois)
#    else:
#        return 0, ""
def fetch_citation_data(doi):
    api_call = f"https://opencitations.net/index/coci/api/v1/citations/{doi}"
    headers = {"Authorization": f"Bearer {api_token}"}
    response = get(api_call)

    if response.status_code == 200:
        data = response.json()
        citation_dois = [item_citing for item_citing in (item.get('citing') for item in data) if item_citing != doi]
        return len(citation_dois), ", ".join(citation_dois)
    else:
        return 0, ""


# Extract and clean author names from the JSON-like format
#def extract_author_names(authors_str):
#    authors_list = re.findall(r"'lastname': '([^']+)', 'firstname': '([^']+)'", authors_str)
#    unique_authors = set([(lastname, firstname) for lastname, firstname in authors_list])
#    cleaned_authors = [f"{firstname} {lastname}" for lastname, firstname in unique_authors]
#    return cleaned_authors
def extract_author_names(authors_str):
    cleaned_authors = {f"{firstname} {lastname}" for lastname, firstname in re.findall(r"'lastname': '([^']+)', 'firstname': '([^']+)'", authors_str)}
    return list(cleaned_authors)


#def create_author_network(authors_data):
#    G = nx.Graph()
#    for authors_str in authors_data:
#        authors = [author.strip() for author in authors_str[0].split(",")]
#        for i in range(len(authors)):
#            for j in range(i+1, len(authors)):
#                G.add_edge(authors[i], authors[j])
#    return G
def create_author_network(authors_data):
    G = nx.Graph()
    for authors_str in authors_data:
        authors = (author.strip() for author in authors_str[0].split(","))
        for i, author in enumerate(authors):
            for other_author in islice(authors, i+1, None):
                G.add_edge(author, other_author)
    return G


# Function to setup AgGrid with specified column order and widths
def setup_aggrid(dataframe):
    gb = GridOptionsBuilder.from_dataframe(dataframe)
    gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)  # Enable pagination with 20 rows per page
    gb.configure_selection('single', use_checkbox=False)  # Single row selection

    # Set the flex for all columns
    for col in dataframe.columns:
        if col == "Title":
            gb.configure_column(col, flex=2)  # Set a higher flex value for the "Title" column
        else:
            gb.configure_column(col, flex=1)

    gridOptions = gb.build()
    return AgGrid(dataframe, gridOptions=gridOptions, fit_columns_on_grid_load=True, key='grid')

# CONNECT AND QUERY DATABASE ============================================================================================
# Connect to the SQLite database
conn = sqlite3.connect("studies.db")
c = conn.cursor()

# Query for the total number of studies
c.execute("SELECT COUNT(*) FROM studies")
total_studies = c.fetchone()[0]
# Query for the total number of citations
c.execute("SELECT SUM(Citation_Count) FROM studies")
total_citations = c.fetchone()[0] or 0  # Handle None result

# Streamlit columns for cards
col1, col2 = st.columns(2)

# Card for Total Number of Studies
col1.markdown(f"""
<div style="background-color:#FFD700;padding:20px;border-radius:10px;">
    <h5 style="color:#595125;">Antal publicerade studier: {total_studies}</h5>
</div>
""", unsafe_allow_html=True)

# Card for Total Number of Citations
col2.markdown(f"""
<div style="background-color:#8470FF;padding:20px;border-radius:10px;">
    <h5 style="color:#2a2355;">Antal genererade citat: {total_citations}</h5>
</div>
""", unsafe_allow_html=True)




# Sidebar for input =====================================================================================================
st.sidebar.header("Registrera artiklar")
pmids_input = st.sidebar.text_input("Ange ett eller flera (separerade med kommatecken) PubMed ID (PMID).")
pmid_list = [pmid.strip() for pmid in pmids_input.split(",") if pmid.strip()]

# Add password field
password = st.sidebar.text_input("Ange lösenord:", type='password')

if st.sidebar.button("Importera"):
    # Check if the password is correct before submitting
    if password == "correct":  # Replace with the correct password
        for pmid in pmid_list:
            try:
                pubmed = PubMed(tool="PubMedSearcher", email="your_email@example.com")  # Replace with your email
                query = pubmed.query(pmid, max_results=1)
                article = next(query, None)

                if article:
                    # Retrieve DOI
                    doi = article.doi or "N/A"
                    pub_year = article.publication_date.year if article.publication_date else "N/A"

                    # Fetch citation data using the DOI
                    if doi != "N/A":
                        citation_count, citation_dois = fetch_citation_data(doi)
                        if citation_count is not None:
                            title = article.title or "N/A"
                            authors_raw = str(article.authors) if hasattr(article, "authors") else "N/A"
                            authors = ", ".join(extract_author_names(authors_raw))
                            journal = article.journal or "N/A"

                            # Save to database
                            c.execute('''
                                INSERT OR REPLACE INTO studies
                                (DOI, PMID, Title, Authors, Journal, Year, Citation_Count, Citation_DOIs)
                                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                            ''', (doi, pmid, title, authors, journal, pub_year, citation_count, citation_dois))
                            conn.commit()

                            st.sidebar.success(f"Din artikel med PMID {pmid} har lagts till och sparats!")
                        else:
                            st.sidebar.error(f"Hittade inte DOI data {doi}")
                    else:
                        st.sidebar.warning(f"DOI saknas för PMID {pmid}")
                else:
                    st.sidebar.warning(f"Inget PMID hittades: {pmid}")
            except Exception as e:
                st.sidebar.error(f"Ett fel inträffade: {e}")
    else:
        st.sidebar.error("Fel lösenord. Försök igen.")

st.write('\n')

# Studies Table using AgGrid ======================================================
st.subheader("Samtliga studier")
st.write("Klicka på en rad för att se detaljar. Du kan filtrera med hjälp av varje kolumns filterknappar.")
c.execute("SELECT * FROM studies")
data = c.fetchall()
if data:
    df = pd.DataFrame(data, columns=["DOI", "Pubmed ID", "Title", "Authors", "Journal", "Year", "Citation_Count", "Citation_DOIs"])
    df = df.rename(columns={
        "Pubmed ID": "Pubmed ID",
        "Title": "Titel",
        "Authors": "Författare",
        "Journal": "Tidskrift",
        "Year": "År",
        "Citation_Count": "Antal citat",
        "Citation_DOIs": "Citerande_Artiklar"
    })
    df = df[["Titel", "År", "Författare", "Tidskrift", "DOI", "Pubmed ID", "Antal citat", "Citerande_Artiklar"]]
    grid_response = setup_aggrid(df)

    # Check if a row is selected
    if grid_response['selected_rows']:
        selected_row = grid_response['selected_rows'][0]  # Get the first selected row

        # Display the contents of the selected row in a nicer format
        st.write("Detaljer för vald studie:")
        st.markdown(f"**Titel:** {selected_row['Titel']}")
        st.markdown(f"**Författare:** {selected_row['Författare']}")
        st.markdown(f"**Tidskrift:** {selected_row['Tidskrift']}")
        st.markdown(f"**DOI:** {selected_row['DOI']}")
        st.markdown(f"**Pubmed ID:** {selected_row['Pubmed ID']}")
        st.markdown(f"**Antal citat:** {selected_row['Antal citat']}")

        # Convert DOI links in "Citerande Artiklar" to clickable links
        citing_articles = selected_row['Citerande_Artiklar'].split(', ')
        citing_links = [f"[{doi}](https://doi.org/{doi})" for doi in citing_articles if doi]
        st.markdown("**Citerande Artiklar:** " + ', '.join(citing_links))

    # Add a download button for the DataFrame
    csv = df.to_csv(index=False)
    st.download_button(
        label="Ladda ner CSV fil",
        data=csv,
        file_name="Artiklar_Kardiologen_Sahlgrenska.csv",
        mime="text/csv",
    )

else:
    st.info("Inga studier importerade ännu.")


# CHARTS ================================================================================================================

# Bar Chart of Published Journals
st.subheader("Tidskrifter")

c.execute("SELECT Journal, COUNT(*) FROM studies GROUP BY Journal")
journal_data = c.fetchall()
if journal_data:
    df_journal = pd.DataFrame(journal_data, columns=["Journal", "Count"])

    # Sort by count and keep the top 35
    df_journal = df_journal.sort_values(by="Count", ascending=False).head(35)

    # Truncate journal names to 30 characters
    df_journal['Journal'] = df_journal['Journal'].apply(lambda x: x[:40] + '...' if len(x) > 40 else x)

    fig_journal = px.bar(df_journal, y="Journal", x="Count", labels={"Count": "Antal artiklar"}, color_discrete_sequence=[purple_color])
    st.plotly_chart(fig_journal)
else:
    st.info("Inga data tillgängliga.")




# Number of publications per author =============================================
st.subheader("Antal publikationer per författare")
st.write("""
Antal publikationer är inte ett mått på forskningskvalitet. Antal publikationer korrelerar med forskningstid och forskningsområde.
""")

# Checkbox to apply filter
apply_filter = st.checkbox("Visa endast Kardiologens anställda", key="author_publications")

# Fetch author data from the database
c.execute("SELECT Authors FROM studies")
author_data = c.fetchall()

# Flatten the list of authors from all rows and count occurrences
authors = [author for authors_str in author_data for author in authors_str[0].split(", ")]
author_occurrences = pd.Series(authors).value_counts().reset_index()
author_occurrences.columns = ["Author", "Publications"]

# Apply filter if checked
if apply_filter:
    author_occurrences = author_occurrences[author_occurrences['Author'].apply(lambda author: any(substring in author for substring in included_substrings))]

# Create the bar plot
fig_author_publications = px.bar(author_occurrences, y="Author", x="Publications", color_discrete_sequence=[purple_color], labels={"Publications": "Antal publikationer"})
fig_author_publications.update_layout(yaxis={'categoryorder':'total ascending'}, height=800)

# Display the plot
st.plotly_chart(fig_author_publications)





# Number of citations per author =============================================
st.subheader("Antal citeringar per författare")
st.write("""
Antal citat är inte ett mått på forskningskvalitet. Antal citat korrelerar starkare med forskningsområde, snarare än kvalitet.
""")

# Checkbox to apply filter
apply_filter = st.checkbox("Visa endast Kardiologens anställda", key="author_citations")

# Fetch author data along with citation counts from the database
c.execute("SELECT Authors, Citation_Count FROM studies")
author_citation_data = c.fetchall()

# Process data to accumulate total citations per author
author_citations = {}
for authors_str, citation_count in author_citation_data:
    authors = authors_str.split(", ")
    for author in authors:
        author_citations[author] = author_citations.get(author, 0) + citation_count

# Convert to DataFrame
df_author_citations = pd.DataFrame(list(author_citations.items()), columns=["Author", "Citations"])

# Apply filter if checked
if apply_filter:
    df_author_citations = df_author_citations[df_author_citations['Author'].apply(lambda author: any(substring in author for substring in included_substrings))]

# Remove authors with zero citations
df_author_citations = df_author_citations[df_author_citations['Citations'] != 0]

df_author_citations.sort_values(by="Citations", ascending=False, inplace=True)

# Create the bar plot
fig_author_citations = px.bar(df_author_citations, y="Author", x="Citations", color_discrete_sequence=[purple_color], labels={"Citations": "Antal citeringar"})
fig_author_citations.update_layout(yaxis={'categoryorder':'total ascending'}, height=800)

# Display the plot
st.plotly_chart(fig_author_citations)



# TRENDER ========================================
st.subheader("Trender i publikationer och citeringar")
df_articles = pd.read_sql_query("SELECT * FROM studies", conn)

# Aggregate data
yearly_data = df_articles.groupby('Year').agg(
    Publications=pd.NamedAgg(column="Title", aggfunc="count"),  # Counting the number of publications per year
    Citations=pd.NamedAgg(column="Citation_Count", aggfunc="sum")    # Summing the citations per year
).reset_index()

# Create subplots: use 'secondary_y' for the second y-axis (citations)
fig = make_subplots(specs=[[{"secondary_y": True}]])

# Add bar chart for publications
fig.add_trace(
    go.Bar(x=yearly_data['Year'], y=yearly_data['Publications'], name="Publikationer", marker=dict(color=purple_color)),
    secondary_y=False,
)

# Add line chart for citations
fig.add_trace(
    go.Scatter(
        x=yearly_data['Year'], 
        y=yearly_data['Citations'], 
        name="Citeringar", 
        mode='lines+markers', 
        line=dict(color=yellow_color),  # Line color
        marker=dict(
            color=dark_yellow_color,  # Darker yellow for dots
            size=10  # Twice as large
        )
    ),
    secondary_y=True,
)

# Set x-axis title
fig.update_xaxes(title_text="Year")

# Set y-axes titles
fig.update_yaxes(title_text="Publikationer", secondary_y=False)
fig.update_yaxes(title_text="Citeringar", secondary_y=True)

# Show the figure
st.plotly_chart(fig)



# Author Network Graph =============================================
# Define colors for matching and non-matching authors
st.subheader("Forskarnätverk")
match_color = purple_color
non_match_color = yellow_color

# Function to create the author network graph
# Function to create the author network graph with an optional filter for a selected author
def create_author_network_graph(author_data, selected_author=None):
    G = nx.Graph()
    for authors_str in author_data:
        authors = [author.strip() for author in authors_str[0].split(",")]
        for i in range(len(authors)):
            for j in range(i+1, len(authors)):
                if selected_author is None or selected_author in [authors[i], authors[j]]:
                    G.add_edge(authors[i], authors[j])
    return G

# Fetch all authors for the dropdown
c.execute("SELECT Authors FROM studies")
author_data = c.fetchall()
all_authors = set([author for authors_str in author_data for author in authors_str[0].split(", ")])

# Streamlit dropdown for selecting an author
#selected_author = st.selectbox("Välj en forskare för att se dennes samarbeten.", [""] + sorted(all_authors), key="author_network_graph")
selected_author = st.selectbox("Välj en forskare (lämna tomt för att se alla)", [""] + sorted(all_authors), index=sorted(all_authors).index("Sebastian Völz"), key="author_network_graph")
st.write("Här kan du se hur forskarna samarbetar med varandra internet och externt. Lila cirklar är Kardiologens anställda. Gula cirklar är externa forskare.")
# Create the network graph based on the selected author
G = create_author_network_graph(author_data, selected_author if selected_author else None)

pos = nx.spring_layout(G)  # Positioning of the nodes

# Nodes
node_x, node_y = [], []
node_text = []  # For labels
node_color = []  # For colors

for node in G.nodes():
    x, y = pos[node]
    node_x.append(x)
    node_y.append(y)
    node_text.append(node)  # Add author name

    # Check if the author matches any of the substrings in 'included_substrings'
    if any(substring in node for substring in included_substrings):
        node_color.append(match_color)
    else:
        node_color.append(non_match_color)

node_trace = go.Scatter(
    x=node_x, 
    y=node_y, 
    mode='markers+text', 
    text=node_text, 
    hoverinfo='text',
    marker=dict(
        showscale=False, 
        colorscale='YlGnBu', 
        size=15,
        color=node_color  # Use the color list
    ),
    textposition='top center', 
    textfont=dict(family='Arial', size=12, color='black')
)

# Edges
edge_x, edge_y = [], []
for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_x.extend([x0, x1, None])
    edge_y.extend([y0, y1, None])
edge_trace = go.Scatter(
    x=edge_x, 
    y=edge_y, 
    line=dict(width=0.5, color='#888'), 
    hoverinfo='none', 
    mode='lines'
)

# Create figure and add traces
fig = go.Figure(
    data=[edge_trace, node_trace], 
    layout=go.Layout(
        showlegend=False, 
        hovermode='closest', 
        margin=dict(b=0, l=0, r=0, t=0)
    )
)

# Display the graph
st.plotly_chart(fig)


# Remove rows from the database by PMID
st.write("\n")
st.write("\n")
st.subheader("Ta bort artikel från databasen")
pmid_to_remove = st.text_input("Ange ett PMID för att ta bort artikeln från databasen:")
remove_password = st.text_input("Bekräfta med lösenord", type='password')

if st.button("Ta bort artikel"):
    # Check if the password is correct before removing the row
    if remove_password == "correct":  # Replace with the correct password
        try:
            # Check if the PMID exists in the database
            c.execute("SELECT COUNT(*) FROM studies WHERE PMID = ?", (pmid_to_remove,))
            row_count = c.fetchone()[0]

            if row_count > 0:
                # Remove the row with the matching PMID from the database
                c.execute("DELETE FROM studies WHERE PMID = ?", (pmid_to_remove,))
                conn.commit()

                st.success(f"Rad med PMID {pmid_to_remove} har tagits bort från databasen.")
            else:
                st.warning(f"Ingen rad hittades med PMID {pmid_to_remove} i databasen.")
        except Exception as e:
            st.error(f"Ett fel inträffade: {e}")
    else:
        st.error("Fel lösenord. Försök igen.")


# Close the database connection
conn.close()