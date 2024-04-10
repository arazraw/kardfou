import streamlit as st
import pandas as pd
import os
import matplotlib.pyplot as plt
import sqlite3

st.set_page_config(page_title="VO Kardiologi - Handledning & Undervisning", page_icon="favicon.ico",)

st.image("logo.png", width=180)
st.title ("Handledning & Undervisning på VO Kardiologi")

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


# FUNCTIONS AND CONSTANTS ======================================================
secret_value = st.secrets["admin_pw"]

# Colors
purple_color = '#8470FF'
yellow_color = 'rgb(255, 215, 0)' #FFD700
dark_yellow_color = 'rgb(255, 215, 0)'
light_blue_color = '#90D2DC'

# Insert data in sqlite database
def insert_data(data):
    conn = sqlite3.connect('teaching.db')
    cursor = conn.cursor()
    cursor.execute('''
        INSERT OR REPLACE INTO teaching (
            email, 
            associate_professor, 
            combination_profession, 
            percentage_research, 
            Year_of_dissertation, 
            active_researchers, 
            Year_registered_phd_students, 
            current_main_supervisor, 
            current_co_supervisor, 
            previous_main_supervisor, 
            previous_co_supervisor, 
            participate_nationally, 
            participate_internationally, 
            participate_regionally, 
            pil_hpe_101, 
            pil_hpe_102, 
            pil_hpe_103, 
            pil_hpe_201, 
            pil_hpe_301, 
            Handledning_av_forskarstuderande, 
            Handledning_i_forskarutbildning, 
            use_prom, 
            use_prom_examples
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
        (
            data['email'], 
            data['associate_professor'],
            data['combination_profession'],
            data['percentage_research'],
            data['Year_of_dissertation'],
            data['active_researchers'],
            data['Year_registered_phd_students'],
            data['current_main_supervisor'],
            data['current_co_supervisor'],
            data['previous_main_supervisor'],
            data['previous_co_supervisor'],
            data['participate_nationally'],
            data['participate_internationally'],
            data['participate_regionally'],
            data['pil_hpe_101'],
            data['pil_hpe_102'],
            data['pil_hpe_103'],
            data['pil_hpe_201'],
            data['pil_hpe_301'],
            data['handledning_forskar'],
            data['handledning_utbildning'],
            data['use_prom'],
            data['use_prom_examples']
        )
    )
    conn.commit()
    conn.close()



def fetch_data():
    conn = sqlite3.connect('teaching.db')
    df = pd.read_sql_query("SELECT * FROM teaching", conn)
    conn.close()
    return df


# Sidebar for data input
st.sidebar.title("Lägg till utbildning & handledning")
with st.sidebar.form(key='researcher_form'):
    email = st.text_input('E-post')
    associate_professor = 1 if st.checkbox("Docent") else 0
    combination_profession = 1 if st.checkbox("Kombinationstjänst") else 0
    active_researchers = 1 if st.checkbox("Är aktiv forskare") else 0
    percentage_research = st.number_input("Andel forskning i tjänst", 0, 100)
    Year_of_dissertation = st.number_input("År för disputation", 1900, 2023)
    current_main_supervisor = st.number_input("Antal huvudhandledarskap nu", 0)
    current_co_supervisor = st.number_input("Antal bihandledarskap nu", 0)
    previous_main_supervisor = st.number_input("Antal huvudhandledarskap tidigare", 0)
    previous_co_supervisor = st.number_input("Antal bihandledarskap tidigare", 0)
    Year_registered_phd_students = st.text_input("År då doktorander registrerats")
    participate_nationally = st.number_input("Genomfört nationell presentation under året", 0)
    participate_internationally = st.number_input("Genomfört internationell presentation under året", 0)
    participate_regionally = st.number_input("Genomfört regional presentation under året", 0)
    pil_hpe_101 = 1 if st.checkbox("PIL101 eller HPE101") else 0
    pil_hpe_102 = 1 if st.checkbox("PIL102 eller HPE102") else 0
    pil_hpe_103 = 1 if st.checkbox("PIL103 eller HPE103") else 0
    pil_hpe_201 = 1 if st.checkbox("PIL201 eller HPE201") else 0
    pil_hpe_301 = 1 if st.checkbox("PIL301 eller HPE301") else 0
    handledning_forskar = 1 if st.checkbox("Handledning av forskarstuderande") else 0
    handledning_utbildning = 1 if st.checkbox("Handledning i forskarutbildning") else 0
    use_prom = 1 if st.checkbox("Använder PROM") else 0
    use_prom_examples = st.text_area("Exempel på PROM")
    password = st.text_input('Lösenord', type='password')
    submit_button = st.form_submit_button("Spara")


if submit_button:
    if password == secret_value:
        # Collecting form data
        new_data = {
            "email": email, 
            "associate_professor": associate_professor,
            "combination_profession": combination_profession,
            "percentage_research": percentage_research,
            "Year_of_dissertation": Year_of_dissertation,
            "active_researchers": active_researchers,
            "Year_registered_phd_students": Year_registered_phd_students,
            "current_main_supervisor": current_main_supervisor,
            "current_co_supervisor": current_co_supervisor,
            "previous_main_supervisor": previous_main_supervisor,
            "previous_co_supervisor": previous_co_supervisor,
            "participate_nationally": participate_nationally,
            "participate_internationally": participate_internationally,
            "participate_regionally": participate_regionally,
            "pil_hpe_101": pil_hpe_101,
            "pil_hpe_102": pil_hpe_102,
            "pil_hpe_103": pil_hpe_103,
            "pil_hpe_201": pil_hpe_201,
            "pil_hpe_301": pil_hpe_301,
            "handledning_forskar": handledning_forskar,
            "handledning_utbildning": handledning_utbildning,
            "use_prom": use_prom,
            "use_prom_examples": use_prom_examples
        }

        # Inserting data into the database
        insert_data(new_data)
        st.success("Data har sparats!")
    else:
        st.error("Fel lösenord. Prova igen.")

# Dashboard
df = fetch_data()
if not df.empty:    
    # Using f-strings for string concatenation
    st.markdown(f"**Docenter:** {df['associate_professor'].sum()}")
    st.markdown(f"**Kombinationstjänster:** {df['combination_profession'].sum()}")
    st.markdown(f"**Aktiva forskare:** {df['active_researchers'].sum()}")
    st.markdown(f"**Antal huvudhandledarskap nu:** {df['current_main_supervisor'].sum()}")
    st.markdown(f"**Antal bihandledarskap nu:** {df['current_co_supervisor'].sum()}")
    st.markdown(f"**Antal huvudhandledarskap tidigare:** {df['previous_main_supervisor'].sum()}")
    st.markdown(f"**Antal bihandledarskap tidigare:** {df['previous_co_supervisor'].sum()}")
    
    st.markdown("---")
    # Aggregate Data
    course_counts = {
        'PIL101 / HPE101': df['pil_hpe_101'].sum(),
        'PIL102 l HPE102': df['pil_hpe_102'].sum(),
        'PIL103 l HPE103': df['pil_hpe_103'].sum(),
        'PIL201 l HPE201': df['pil_hpe_201'].sum(),
        'PIL301 l HPE301': df['pil_hpe_301'].sum(),
        'Handl av forskarstud': df['Handledning_av_forskarstuderande'].sum(),
        'Handl i forskarutb': df['Handledning_i_forskarutbildning'].sum()
    }

    # Create Bar Charts
    fig, ax = plt.subplots()
    ax.bar(course_counts.keys(), course_counts.values(), color=purple_color)
    ax.set_xlabel('')
    ax.set_ylabel('Antal')
    ax.set_title('Genomförda kurser')
    plt.xticks(rotation=45)
    st.pyplot(fig)

    #st.markdown("---")
    # Create the histogram
    #fig, ax = plt.subplots()
    #ax.hist(df['percentage_research'], bins=range(0, 101, 10), color=purple_color)
    #ax.set_xlabel('Andel forskning')
    #ax.set_ylabel('Antal anställda')
    #ax.set_title('Forskning (%) som del av anställning')

    # Display the Histogram in Streamlit
    #st.pyplot(fig)

    #st.markdown("---")
    #st.subheader("Rådata")
    #st.dataframe(df)
else:
    st.warning("Ingen data hittades.")


st.markdown("---")

# Add input fields for email and password
delete_email = st.text_input('Radera rad med e-post')
delete_password = st.text_input('Lösenord', type='password')

# Add a button to submit the deletion request
delete_button = st.button("Ta bort rad")

if delete_button:
    if delete_password == secret_value:
        # Check if the email exists in the DataFrame
        if delete_email in df['email'].values:
            # Delete the row from the database
            conn = sqlite3.connect('teaching.db')
            cursor = conn.cursor()
            cursor.execute("DELETE FROM teaching WHERE email = ?", (delete_email,))
            conn.commit()
            conn.close()
            
            # Update the DataFrame to reflect the changes
            df = df[df['email'] != delete_email]
            
            st.success(f"Rad med e-post '{delete_email}' har raderats.")
        else:
            st.error(f"E-post '{delete_email}' hittades inte.")
    else:
        st.error("Fel lösenord. Försök igen.")
