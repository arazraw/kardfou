import sqlite3

# Connect tos the SQLite database
conn = sqlite3.connect("studies.db")
c = conn.cursor()

# Create the "studies" table with the desired schema
c.execute('''
    CREATE TABLE IF NOT EXISTS studies (
        DOI TEXT, 
        PMID TEXT PRIMARY KEY, 
        Title TEXT, 
        Authors TEXT, 
        Journal TEXT,
        Year INTEGER,
        Citation_Count INTEGER,
        Citation_DOIs TEXT
    )
''')

# Commit the changes and close the database connection
conn.commit()
conn.close()
