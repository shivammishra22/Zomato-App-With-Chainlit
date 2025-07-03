import os
import pandas as pd
from docx import Document

def extract_specific_table(docx_path, keywords):
    doc = Document(docx_path)
    matched_table = None

    for table in doc.tables:
        for row in table.rows:
            row_text = [cell.text.strip() for cell in row.cells]
            if any(keyword in cell for cell in row_text for keyword in keywords):
                matched_table = table
                break
        if matched_table:
            break

    if matched_table:
        extracted_data = []
        for row in matched_table.rows:
            extracted_data.append([cell.text.strip() for cell in row.cells])

        # Remove duplicate header rows if the first two rows are identical
        if len(extracted_data) > 1 and extracted_data[0] == extracted_data[1]:
            extracted_data.pop(1)

        return extracted_data
    else:
        print("❌ No matching table found.")
        return []

def save_table_to_excel(table_data, excel_path):
    os.makedirs(os.path.dirname(excel_path), exist_ok=True)
    df = pd.DataFrame(table_data[1:], columns=table_data[0])  # Use first row as header
    df.to_excel(excel_path, sheet_name='Extracted_Table', index=False)
    print(f"✅ Specific table saved to {excel_path}")

# ✅ Replace with actual paths
docx_path = r"C:\Users\shivam.mishra2\Downloads\embedding\01PSUR\Data request form.docx"
excel_path = r"C:\Users\shivam.mishra2\Downloads\New_Psur_File\intial_table1.xlsx"

# Keywords to identify the specific table
keywords = ["Molecular Product", "Study Number", "Test Product Name","Active comparator name","TestProduct","Active Comparator","Placebo","Total"]
# Extract and save
table_data = extract_specific_table(docx_path, keywords)
if table_data:
    save_table_to_excel(table_data, excel_path)
