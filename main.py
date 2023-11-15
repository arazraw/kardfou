import streamlit as st
from requests import get

def fetch_citation_data(identifier):
    api_call = f"https://opencitations.net/api/v1/metadata/{identifier}"
    response = get(api_call)

    if response.status_code == 200:
        data = response.json()
        if data:
            citations = data[0].get('citation', '').split('; ')
            citation_count = len(citations) if citations[0] != '' else 0
            return citation_count, "; ".join(citations)
        else:
            return 0, ""
    else:
        return None, None


# Function to display the data
def display_data(data):
    for item in data:
        st.write(f"**Title:** {item.get('title', 'N/A')}")
        st.write(f"**Authors:** {item.get('author', 'N/A')}")
        st.write(f"**Publication Date:** {item.get('pub_date', 'N/A')}")
        st.write(f"**Venue:** {item.get('venue', 'N/A')}")
        st.write(f"**Volume:** {item.get('volume', 'N/A')}")
        st.write(f"**Issue:** {item.get('issue', 'N/A')}")
        st.write(f"**Pages:** {item.get('page', 'N/A')}")
        st.write(f"**DOI:** {item.get('doi', 'N/A')}")

        citations = item.get('citation', '').split('; ')
        st.write(f"**Number of Citations:** {len(citations)}")
        st.write("**Citations:**")
        for citation in citations:
            st.write(citation)

        st.write("----")

# Streamlit app
def main():
    st.title("OpenCitations Citation Fetcher")

    # Input for DOI
    identifier = st.text_input("Enter DOI (format: doi:10.xxxx/xxxx):")

    # Button to fetch citations
    if st.button("Fetch Citations"):
        if identifier:
            data = fetch_citation_data(identifier)
            if data:
                display_data(data)
            else:
                st.error("No data found or an error occurred.")
        else:
            st.error("Please enter a DOI.")

if __name__ == "__main__":
    main()
