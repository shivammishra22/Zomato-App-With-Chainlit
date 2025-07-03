import pandas as pd
from docx import Document
import numpy as np

# === USER INPUTS ===
file_path = r"C:\Users\shivam.mishra2\Downloads\New_Psur_File\marketing_exposure_tables.xlsx"
ddd_value = 10  # Replace with actual DDD value
country_name = "South Africa"
medicine = "Esomeprazole"
place = "South Africa"
date = "2020-01-01"

# === DATA LOADING ===
df = pd.read_excel(file_path, engine='openpyxl')

# === DATA CLEANING ===
df["Strength in mg"] = df["Strength in mg"].astype(str).str.replace("mg", "").str.strip()
df["Strength in mg"] = pd.to_numeric(df["Strength in mg"], errors='coerce')

pack_column = next((col for col in ["Pack", "Packs"] if col in df.columns), None)
if pack_column:
    df[pack_column] = pd.to_numeric(df[pack_column].astype(str).str.replace(",", "").str.extract(r'(\d+)')[0], errors='coerce').fillna(0).astype(int)

if "Pack size" in df.columns:
    pack_size_extracted = df["Pack size"].astype(str).str.replace("'", "").str.extract(r'(\d+)\s*[xX]\s*(\d+)')
    df["Pack size"] = pd.to_numeric(pack_size_extracted[0], errors='coerce').fillna(1).astype(int) * pd.to_numeric(pack_size_extracted[1], errors='coerce').fillna(1).astype(int)

col = "Number of tablets / Capsules/Injections"
df[col] = df[col].astype(str).str.replace(",", "").str.split(":").str[-1].str.strip()
df[col] = pd.to_numeric(df[col], errors='coerce')

df["Delivered quantity (mg)"] = pd.to_numeric(df["Delivered quantity (mg)"].astype(str).str.replace(",", ""), errors='coerce')

df.rename(columns={"Product": "Molecule"}, inplace=True)

# === CALCULATIONS ===
df["Sales Figure (mg) or period/Volume of sales (in mg)"] = df["Number of tablets / Capsules/Injections"] * df["Strength in mg"]
df["Patients Exposure (PTY) for period"] = df["Sales Figure (mg) or period/Volume of sales (in mg)"] / (ddd_value * 365)
df["Patients Exposure (PTY) for period"] = df["Patients Exposure (PTY) for period"].round(0)

# === FILTERING ===

# Filter data
df_country = df[df["Country"] == "South Africa"].reset_index(drop=True)
df_non_country = df[df["Country"] != "South Africa"].reset_index(drop=True)

# Calculate totals
country_total_exposure = df_country["Patients Exposure (PTY) for period"].sum()
non_country_total_exposure = df_non_country["Patients Exposure (PTY) for period"].sum()

# Create total rows
country_total_row = pd.DataFrame([{
    "Country": "Total",
    "Patients Exposure (PTY) for period": country_total_exposure
}])

non_country_total_row = pd.DataFrame([{
    "Country": "Total",
    "Patients Exposure (PTY) for period": non_country_total_exposure
}])

# Append total rows
df_country = pd.concat([df_country, country_total_row], ignore_index=True)
df_non_country = pd.concat([df_non_country, non_country_total_row], ignore_index=True)
# Replace NaN values with blank space only in rows where Country is 'Total'
df.loc[df['Country'] == 'Total'] = df.loc[df['Country'] == 'Total'].fillna("")

# Replace all NaN values with blank space
df.replace(np.nan, '', inplace=True)


# === WORD DOCUMENT GENERATION ===
doc = Document()
doc.add_heading("5.3 Cumulative and Interval Patient Exposure from Marketing Experience", level=1)

summary_text = (
    f"The MAH obtained initial MA for their generic formulation of {medicine} in {place} on {date}.\n"
    f"The post-authorization exposure to {medicine} during the cumulative period was {int(total_exposure)} patients\n"
    f"({place}: {int(country_exposure)} and Non {place}: {int(non_country_exposure)}) treatment days approximately and presented in Table 3."
)
doc.add_paragraph(summary_text)

def add_table_to_doc(document, dataframe, title):
    document.add_heading(title, level=2)
    table = document.add_table(rows=1, cols=len(dataframe.columns))
    table.style = 'Table Grid'

    hdr_cells = table.rows[0].cells
    for i, col_name in enumerate(dataframe.columns):
        hdr_cells[i].text = str(col_name)

    for _, row in dataframe.iterrows():
        row_cells = table.add_row().cells
        for i, item in enumerate(row):
            row_cells[i].text = str(item)

add_table_to_doc(doc, df_country, f"{country_name} Data")
add_table_to_doc(doc, df_non_country, f"Non-{country_name} Data")

doc.save("Esomeprazole_Exposure.docx")
print("✅ Word document saved as 'Esomeprazole_Exposure.docx'")
