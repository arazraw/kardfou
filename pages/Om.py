import streamlit as st

st.set_page_config(page_title="VO Kardiologi", page_icon="favicon.ico")
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

st.image("logo.png", width=180)
st.title ("FOU Rapport - VO Kardiologi")
st.header ("Om rapporten")
st.write("Rapporten är en sammanställning av forskning, handledning och undervisning på VO Kardiologi. Den innehåller information från alla år, med möjlighet att filtrera på år och forskare. Rapportens källkod finns att tillgå.")
st.write("---")
st.write("Verksamhetschef: Kristofer Skoglund (kristofer.skoglund@vgregion.se)")
st.write("Information om applikationen: Araz Rawshani (araz.rawshani@vgregion.se)")
st.subheader("Kardiologens FOU-Råd")
st.write("""
         Kristofer Skoglund, Araz Rawshani, Gustav Smith, Annika Odell, Sofie Andréen,\n
         Oskar Angerås, Entela Bollano, Annica Ravn-Fischer, Elmir Omerovic, Björn Redfors,\n
         Clara Hjalmarsson, Charlotta Ljungman, Stefano Romeo, Pia Dahlberg,\n
         Bert Andersson, Lennart Bergfeldt, Runa Sigurjonsdottir,
         """)

# Link to website of ALF Västra Götaland (https://www.alfvastragotaland.se/)
st.subheader("Länkar")
# Button with link to https://www.alfvastragotaland.se/
if st.button("ALF Västra Götaland"):
    st.write("https://www.alfvastragotaland.se/")
#Button with link to https://www.sahlgrenska.se/omraden/omrade-6/kardiologi/
if st.button("VO Kardiologi"):
    st.write("https://www.sahlgrenska.se/omraden/omrade-6/kardiologi/")