import os
import re
import pandas as pd
from docx import Document

# === STEP 1: Extract table from DOCX and save to Excel ===
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
        extracted_data = [[cell.text.strip() for cell in row.cells] for row in matched_table.rows]
        if len(extracted_data) > 1 and extracted_data[0] == extracted_data[1]:
            extracted_data.pop(1)
        return extracted_data
    else:
        print("❌ No matching table found.")
        return []

def save_table_to_excel(table_data, excel_path):
    os.makedirs(os.path.dirname(excel_path), exist_ok=True)
    df = pd.DataFrame(table_data[1:], columns=table_data[0])
    df.to_excel(excel_path, sheet_name='Extracted_Table', index=False)
    print(f"✅ Table saved to Excel: {excel_path}")
    return excel_path  # Return path for next step

# === STEP 2: Generate Summary DOCX from Excel ===
def generate_summary_doc(excel_path, output_doc_path):
    df = pd.read_excel(excel_path, engine='openpyxl')
    columns = pd.MultiIndex.from_tuples([
        ('Molecule/Product', ''), ('Study Number', ''), ('Study Title', ''), ('Test product name', ''),
        ('Active comparator name', ''), ('No. of subjects enrolled', 'Test Product'),
        ('No. of subjects enrolled', 'Active Comparator'), ('No. of subjects enrolled', 'Placebo'),
        ('No. of subjects enrolled', 'Total'), ('Gender distribution of subjects enrolled', 'Male'),
        ('Gender distribution of subjects enrolled', 'Female'), ('Age distribution of subjects enrolled', '<18 years'),
        ('Age distribution of subjects enrolled', '18-65 years'), ('Age distribution of subjects enrolled', '>65 years'),
        ('Racial distribution of subjects enrolled', 'Asian'), ('Racial distribution of subjects enrolled', 'Black'),
        ('Racial distribution of subjects enrolled', 'Caucasian'), ('Racial distribution of subjects enrolled', 'Other'),
        ('Racial distribution of subjects enrolled', 'Unknown')
    ])
    df.columns = columns
    df_new = df.iloc[2:].copy()
    df_new.reset_index(drop=True, inplace=True)

    numeric_cols = [col for col in columns if col[1]]
    df_new[numeric_cols] = df_new[numeric_cols].apply(pd.to_numeric, errors='coerce')

    target_columns = [
        ('Age distribution of subjects enrolled', '<18 years'),
        ('Age distribution of subjects enrolled', '18-65 years'),
        ('Age distribution of subjects enrolled', '>65 years'),
        ('Racial distribution of subjects enrolled', 'Asian'),
        ('Racial distribution of subjects enrolled', 'Black'),
        ('Racial distribution of subjects enrolled', 'Caucasian')
    ]
    non_zero_columns = [col for col in target_columns if df_new[col].sum() != 0]

    doc = Document()
    doc.add_heading('1 Jubilant Generics Limited', level=1)
    doc.add_heading('2 Levetiracetam Periodic Safety Update Report', level=2)
    doc.add_paragraph('Reporting period: 30-Nov-2024 to 30-Nov-2024')
    doc.add_heading('5 ESTIMATED EXPOSURE AND USE PATTERNS', level=2)
    doc.add_heading('5.1 General considerations', level=3)
    doc.add_paragraph(
        'For clinical trials patient exposure can be accurately calculated because dosage and duration of treatment are clearly known. '
        'In terms of post marketing use patient exposure cannot be accurately calculated for certain reasons such as varying dosage and '
        'duration of treatment as well as changing or unknown patient compliance.'
    )
    doc.add_heading('5.2 Cumulative Subject Exposure in Clinical Trials', level=3)

    total_subjects = int(df_new[('No. of subjects enrolled', 'Total')].sum())
    num_studies = len(df_new)
    age_group = non_zero_columns[0][1] if non_zero_columns else "Unknown"
    gender_group = 'Male' if df_new[('Gender distribution of subjects enrolled', 'Male')].sum() > 0 else 'Female'
    medicine_name = "levetiracetam"
    summary_text = (
        f"Jubilant as MAH has not conducted any Clinical Trials. However, Jubilant has conducted {num_studies} BA/BE study with {medicine_name} till the DLP of the PSUR 30-Nov-2024 "
        f"and cumulative subject exposure in the completed clinical trial were {total_subjects} subjects.\n\n"
        f"Of these {total_subjects}, all were {gender_group} subjects of age distribution between {age_group} years. "
        f"Cumulative subject exposure to {medicine_name} in BA/BE studies is given in the table below:"
    )
    doc.add_paragraph(summary_text)
    os.makedirs(os.path.dirname(output_doc_path), exist_ok=True)
    doc.save(output_doc_path)
    print(f"✅ Summary document saved to: {output_doc_path}")
    return output_doc_path  # Return for merging

# === STEP 3: Clean original DOCX ===
def delete_tables_with_keywords(doc, keywords):
    tables_to_delete = []
    for table in doc.tables:
        table_text = " ".join(cell.text for row in table.rows for cell in row.cells)
        if any(re.search(keyword, table_text, re.IGNORECASE) for keyword in keywords):
            tables_to_delete.append(table)
    for table in tables_to_delete:
        tbl = table._element
        tbl.getparent().remove(tbl)

def extract_text_before_stop_line(doc, stop_line):
    full_text = "\n".join([para.text for para in doc.paragraphs])
    match = re.search(rf"^(.*?)\b{re.escape(stop_line)}\b", full_text, re.DOTALL | re.IGNORECASE)
    return match.group(1) if match else ""

def delete_multiple_paragraphs(doc, lines_to_delete_text):
    lines_to_delete = [line.strip() for line in lines_to_delete_text.splitlines() if line.strip()]
    for para in doc.paragraphs:
        if para.text.strip() in lines_to_delete:
            p = para._element
            p.getparent().remove(p)
            p._p = p._element = None

def remove_blank_paragraphs(doc):
    for para in doc.paragraphs:
        if not para.text.strip():
            p = para._element
            p.getparent().remove(p)

def clear_text(doc):
    for para in doc.paragraphs:
        para.clear()

def process_document(doc_path, keywords, stop_line, output_cleaned_path):
    doc = Document(doc_path)
    delete_tables_with_keywords(doc, keywords)
    extracted_text = extract_text_before_stop_line(doc, stop_line)
    if extracted_text:
        delete_multiple_paragraphs(doc, extracted_text)
    remove_blank_paragraphs(doc)
    clear_text(doc)
    os.makedirs(os.path.dirname(output_cleaned_path), exist_ok=True)
    doc.save(output_cleaned_path)
    print(f"✅ Cleaned document saved to: {output_cleaned_path}")
    return output_cleaned_path

# === STEP 4: Merge two documents ===
def merge_documents(doc1_path, doc2_path, merged_path):
    doc1 = Document(doc1_path)
    doc2 = Document(doc2_path)
    merged = Document()
    for element in doc1.element.body:
        merged.element.body.append(element)
    for element in doc2.element.body:
        merged.element.body.append(element)
    merged.save(merged_path)
    print(f"✅ Documents merged to: {merged_path}")

# === MAIN EXECUTION ===
if __name__ == "__main__":
    # Base docx path
    docx_path = r"C:\Users\shivam.mishra2\Downloads\embedding\01PSUR\Data request form.docx"
    base_dir = os.path.dirname(docx_path).replace("embedding\\01PSUR", "New_Psur_File")
    os.makedirs(base_dir, exist_ok=True)

    excel_path = os.path.join(base_dir, "intial_table1.xlsx")
    section5_path = os.path.join(base_dir, "Section5.docx")
    final_table_path = os.path.join(base_dir, "final_Table.docx")
    merged_path = os.path.join(base_dir, "mergedFINAL_document.docx")

    # Keywords
    table_keywords = ["Molecular Product", "Study Number", "Test Product Name", "Active comparator name", "TestProduct", "Active Comparator", "Placebo", "Total"]
    delete_keywords = ["Country", "Signal term"]
    stop_line = "Pooled literature Data: PVG Department"

    # Run all steps
    table_data = extract_specific_table(docx_path, table_keywords)
    if table_data:
        excel_generated = save_table_to_excel(table_data, excel_path)
        section5_generated = generate_summary_doc(excel_generated, section5_path)
        final_cleaned = process_document(docx_path, delete_keywords, stop_line, final_table_path)
        merge_documents(section5_generated, final_cleaned, merged_path)
