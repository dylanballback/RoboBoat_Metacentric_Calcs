import sqlite3
import adsk.core, adsk.fusion, traceback

# Function to create the SQLite database and the hull_data table.
# This function needs to be run only once.
def create_database():
    # Connect to SQLite database. If it doesn't exist, it will be created.
    conn = sqlite3.connect('boat_hull_data.db')
    cursor = conn.cursor()
    # Create a new table named 'hull_data' if it doesn't already exist.
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS hull_data (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            hull_type TEXT,
            hull_size_category TEXT,
            offset INTEGER,
            volume REAL,
            length REAL,
            width REAL,
            cg_x REAL,
            cg_y REAL,
            cg_z REAL
        )
    ''')
    conn.commit()
    conn.close()

# Function to insert data into the hull_data table.
def insert_data(hull_type, hull_size_category, offset, volume, length, width, cg_x, cg_y, cg_z):
    # Connect to the SQLite database.
    conn = sqlite3.connect('boat_hull_data.db')
    cursor = conn.cursor()
    # Insert data into the table.
    cursor.execute('''
        INSERT INTO hull_data (hull_type, hull_size_category, offset, volume, length, width, cg_x, cg_y, cg_z)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (hull_type, hull_size_category, offset, volume, length, width, cg_x, cg_y, cg_z))
    conn.commit()
    conn.close()


# Function to retrieve data for a specific hull based on its id.
def get_data_for_hull(hull_id):
    # Connect to the SQLite database.
    conn = sqlite3.connect('boat_hull_data.db')
    cursor = conn.cursor()
    # Select and retrieve data for the specified hull_id.
    cursor.execute('SELECT * FROM hull_data WHERE id = ?', (hull_id,))
    rows = cursor.fetchall()
    conn.close()
    return rows

# Function to get aggregate statistics for hulls.
def get_aggregate_statistics():
    # Connect to the SQLite database.
    conn = sqlite3.connect('boat_hull_data.db')
    cursor = conn.cursor()
    # Select and retrieve aggregate data grouped by hull_type.
    cursor.execute('SELECT hull_type, AVG(volume), AVG(length), AVG(width) FROM hull_data GROUP BY hull_type')
    rows = cursor.fetchall()
    conn.close()
    return rows


# Call create_database only once manually to set up the database.
# create_database()


#Example on how to use in main script.
"""
# Main function where the Fusion 360 script runs.
def run(context):
    try:
        # [Your existing run function code goes here]
        # After each data point is calculated:
        insert_data(hull_type, hull_size_category, offset, volume, length, width, cg_x, cg_y, cg_z)

        # [Rest of your run function code goes here]

    except:
        # Error handling
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

"""