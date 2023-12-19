import sqlite3

# Connect to SQLite database (it will be created if it doesn't exist)
conn = sqlite3.connect('teaching.db')

# Create a cursor object using the cursor() method
cursor = conn.cursor()

# SQL query to create a table with the correct data types
create_table_query = '''
CREATE TABLE IF NOT EXISTS teaching (
    email TEXT PRIMARY KEY,
    associate_professor INTEGER,  -- BOOLEAN represented as INTEGER
    combination_profession INTEGER,  -- BOOLEAN
    percentage_research REAL,
    Year_of_dissertation INTEGER,
    active_researchers INTEGER,  -- BOOLEAN
    Year_registered_phd_students INTEGER,
    current_main_supervisor INTEGER,
    current_co_supervisor INTEGER,
    previous_main_supervisor INTEGER,
    previous_co_supervisor INTEGER,
    participate_nationally INTEGER,  -- BOOLEAN
    participate_internationally INTEGER,  -- BOOLEAN
    participate_regionally INTEGER,  -- BOOLEAN
    pil_hpe_101 INTEGER,  -- BOOLEAN
    pil_hpe_102 INTEGER,  -- BOOLEAN
    pil_hpe_103 INTEGER,  -- BOOLEAN
    pil_hpe_201 INTEGER,  -- BOOLEAN
    pil_hpe_301 INTEGER,  -- BOOLEAN
    Handledning_av_forskarstuderande INTEGER,  -- BOOLEAN
    Handledning_i_forskarutbildning INTEGER,  -- BOOLEAN
    use_prom INTEGER,  -- BOOLEAN
    use_prom_examples TEXT
);
'''

# Execute the SQL command
cursor.execute(create_table_query)

# Commit the changes and close the connection
conn.commit()
conn.close()
