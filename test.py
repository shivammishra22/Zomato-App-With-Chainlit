

import os
import re
import pandas as pd
from docx import Document

def extract_table_after_text(doc, search_text):
    pattern = re.compile(re.escape(search_text), re.IGNORECASE)
    found_index = None

    # Find the paragraph index that matches the search text
    for i, para in enumerate(doc.paragraphs):
        if pattern.search(para.text):
            found_index = i
            break

    if found_index is None:
        print(f"❌ Text not found: {search_text}")
        return []

    # Locate the table that comes after the matched paragraph
    para_counter = 0
    table_counter = 0
    for block in doc.element.body:
        if block.tag.endswith('p'):
            para_counter += 1
        elif block.tag.endswith('tbl'):
            if para_counter > found_index:
                table = doc.tables[table_counter]
                break
            table_counter += 1
    else:
        print(f"❌ No table found after: {search_text}")
        return []

    # Extract table data
    table_data = []
    for row in table.rows:
        table_data.append([cell.text.strip() for cell in row.cells])

    # Remove duplicate header rows
    if len(table_data) > 1 and table_data[0] == table_data[1]:
        table_data.pop(1)

    return table_data
def save_multiple_tables_to_excel(tables_dict, excel_path):
    os.makedirs(os.path.dirname(excel_path), exist_ok=True)
    with pd.ExcelWriter(excel_path) as writer:
        for sheet_name, table_data in tables_dict.items():
            if table_data:
                df = pd.DataFrame(table_data[1:], columns=table_data[0])
                df.to_excel(writer, sheet_name=sheet_name[:31], index=False)  # Excel sheet name limit
    print(f"✅ All tables saved to {excel_path}")

# === USAGE ===
docx_path = r"C:\Users\shivam.mishra2\Downloads\ALL_PSUR_File\PSUR_all _Data\Esomeprazole ADCO_19-Apr-2017 to 31-Dec-2024\Draft\DRA\Data request form_Esomeprazole ADCO.docx"
excel_path = r"C:\Users\shivam.mishra2\Downloads\New_Psur_File\marketing_exposure_tables.xlsx"

search_texts = [
    "Sales data of Esomeprazole 20-40 mg gastro-resistant tablets for the period: 19 April 2019 to 31 December 2024",
    "Cumulative sales data sale required",
    "Distribution details for Esomeprazole 40 mg tablets",
    
]

doc = Document(docx_path)
tables_dict = {}

for text in search_texts:
    table_data = extract_table_after_text(doc, text)
    if table_data:
        tables_dict[text[:30]] = table_data  # Sheet name trimmed to 30 chars

if tables_dict:
    save_multiple_tables_to_excel(tables_dict, excel_path)

# Second File

import pandas as pd
from docx import Document

# === USER INPUTS ===
file_path = r"C:\Users\shivam.mishra2\Downloads\New_Psur_File\marketing_exposure_tables.xlsx"
ddd_value = 10  # Replace with actual DDD value
country_name = "South Africa"
medicine = "Esomeprazole"
place = "South Africa"
date = "2020-01-01"

# === LOAD DATA ===
df = pd.read_excel(file_path, engine='openpyxl')

# === CLEANING AND CONVERSION ===

# Handle 'Pack' or 'Packs'
possible_pack_columns = ["Pack", "Packs"]
pack_column = next((col for col in possible_pack_columns if col in df.columns), None)
if pack_column:
    df[pack_column] = pd.to_numeric(df[pack_column].astype(str).str.extract(r'(\d+)')[0], errors='coerce').fillna(0).astype(int)
else:
    print("❌ Neither 'Pack' nor 'Packs' column found.")

# Clean 'Pack size'
if "Pack size" in df.columns:
    pack_size_extracted = df["Pack size"].astype(str).str.replace("'", "").str.extract(r'(\d+)\s*[xX]\s*(\d+)')
    df["Pack size"] = pd.to_numeric(pack_size_extracted[0], errors='coerce').fillna(1).astype(int) * pd.to_numeric(pack_size_extracted[1], errors='coerce').fillna(1).astype(int)

# Clean numeric columns
numeric_columns = [
    "Number of tablets / Capsules/Injections",
    "Delivered quantity (mg)",
    "Strength in mg"
]
for col in numeric_columns:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col].astype(str).str.extract(r'(\d+)')[0], errors='coerce').fillna(0).astype(int)

# === CALCULATIONS ===
df["Sales Figure (mg) or period/Volume of sales (in mg)"] = df["Number of tablets / Capsules/Injections"] * df["Strength in mg"]
df["Patients Exposure (PTY) for period"] = df["Sales Figure (mg) or period/Volume of sales (in mg)"] / (ddd_value * 365)
df["Patients Exposure (PTY) for period"] = df["Patients Exposure (PTY) for period"].round(0)

# === FILTERING ===
df_country = df[df['Country'] == country_name].reset_index(drop=True)
df_non_country = df[df['Country'] != country_name].reset_index(drop=True)

total_exposure = df["Patients Exposure (PTY) for period"].sum()
country_exposure = df_country["Patients Exposure (PTY) for period"].sum()
non_country_exposure = df_non_country["Patients Exposure (PTY) for period"].sum()

# === WORD DOCUMENT ===
doc = Document()
doc.add_heading('5.3 Cumulative and Interval Patient Exposure from Marketing Experience', level=1)

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
