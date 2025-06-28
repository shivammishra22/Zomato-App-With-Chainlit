import os
import pandas as pd
from docx import Document
from config.styling_config import DocumentStyling

# ------------------ Step 1: Extract Specific Table from Word ------------------ #
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

        if len(extracted_data) > 1 and extracted_data[0] == extracted_data[1]:
            extracted_data.pop(1)

        return extracted_data
    else:
        print("❌ No matching table found.")
        return []

# ------------------ Step 2: Save Table to Excel ------------------ #
def save_table_to_excel(table_data, excel_path):
    os.makedirs(os.path.dirname(excel_path), exist_ok=True)
    df = pd.DataFrame(table_data[1:], columns=table_data[0])
    df.to_excel(excel_path, sheet_name='Extracted_Table', index=False)
    print(f"✅ Specific table saved to {excel_path}")

# ------------------ Step 3: Process Excel Data ------------------ #
def process_excel_data(excel_path):
    df = pd.read_excel(excel_path, engine='openpyxl')

    df.columns = [
        "Molecule/ Product", "Study Number", "Study Title", "Test product name",
        "Active comparator name", "Test\nProduct", "Active Comparator", "Placebo", "Total",
        "Male", "Female", "<18 years", "18-65 years", ">65 years",
        "Asian", "Black", "Caucasian", "Other", "Unknown"
    ]

    df = df[2:].reset_index(drop=True)

    columns_to_convert = [
        "Test\nProduct", "Active Comparator", "Placebo", "Total",
        "Male", "Female", "<18 years", "18-65 years", ">65 years",
        "Asian", "Black", "Caucasian", "Other", "Unknown"
    ]
    df[columns_to_convert] = df[columns_to_convert].apply(pd.to_numeric, errors='coerce')

    target_columns = ["<18 years", "18-65 years", ">65 years", "Asian", "Black", "Caucasian"]
    non_zero_columns = [col for col in target_columns if df[col].sum() != 0]

    return df, non_zero_columns

# ------------------ Step 4: Generate Word Document ------------------ #
def generate_word_document(df, non_zero_columns, word_path):
    doc = Document()

    # Heading
    DocumentStyling.create_split_heading(doc.add_paragraph(), '1 Jubilant Generics Limited')
    DocumentStyling.create_split_heading(doc.add_paragraph(), '2 Levetiracetam Periodic Safety Update Report')

    # Reporting Period
    para = doc.add_paragraph()
    run = para.add_run('Reporting period: 30-Nov-2024 to 30-Nov-2024')
    DocumentStyling.apply_bold_content_style(run)

    # Section 5 Title
    DocumentStyling.create_split_heading(doc.add_paragraph(), '5 ESTIMATED EXPOSURE AND USE PATTERNS')
    DocumentStyling.create_split_subheading(doc.add_paragraph(), '5.1 General considerations')

    para = doc.add_paragraph()
    run = para.add_run(
        'For clinical trials patient exposure can be accurately calculated because dosage and duration of treatment are clearly known. '
        'In terms of post marketing use patient exposure cannot be accurately calculated for certain reasons such as varying dosage and '
        'duration of treatment as well as changing or unknown patient compliance.'
    )
    DocumentStyling.apply_content_style(run)

    DocumentStyling.create_split_subheading(doc.add_paragraph(), '5.2 Cumulative Subject Exposure in Clinical Trials')

    total_subjects = int(df["Total"].sum())
    num_studies = len(df)
    age_group = non_zero_columns[0] if len(non_zero_columns) > 0 else "Unknown"
    gender_group = non_zero_columns[1] if len(non_zero_columns) > 1 else "Unknown"

    summary_text = (
        f"Jubilant as MAH has not conducted any Clinical Trials. However, Jubilant has conducted {num_studies} BA/BE study with levetiracetam till the DLP of the PSUR 30-Nov-2024 "
        f"and cumulative subject exposure in the completed clinical trial were {total_subjects} subjects.\n\n"
        f"Of these {total_subjects}, all were {gender_group} male subjects of age distribution between {age_group} years. "
        "Cumulative subject exposure to levetiracetam in BA/BE studies is given in the table below:"
    )

    para = doc.add_paragraph()
    run = para.add_run(summary_text)
    DocumentStyling.apply_content_style(run)

    doc.add_paragraph("\n")

    # Table
    table = doc.add_table(rows=1, cols=len(df.columns))
    table.style = DocumentStyling.TABLE_STYLE if hasattr(DocumentStyling, 'TABLE_STYLE') else 'Table Grid'

    hdr_cells = table.rows[0].cells
    for i, column_name in enumerate(df.columns):
        run = hdr_cells[i].paragraphs[0].add_run(str(column_name))
        DocumentStyling.apply_table_header_style(hdr_cells[i])

    for _, row in df.iterrows():
        row_cells = table.add_row().cells
        for i, cell_value in enumerate(row):
            run = row_cells[i].paragraphs[0].add_run(str(cell_value))
            DocumentStyling.apply_content_style(run)

    doc.save(word_path)
    print(f"✅ Word document saved to {word_path}")

# ------------------ Main Execution ------------------ #
docx_path = r"C:\Users\shivam.mishra2\Downloads\Data request form.docx"
excel_path = r"C:\Users\shivam.mishra2\Downloads\embedding\01PSUR\Levetiracetam_Exposure.xlsx"
word_path = r"C:\Users\shivam.mishra2\Downloads\embedding\01PSUR\Levetirace_final.docx"

keywords = ["Molecular Product", "Study Number", "Test Product Name", "Active comparator name",
            "TestProduct", "Active Comparator", "Placebo", "Total"]

table_data = extract_specific_table(docx_path, keywords)
if table_data:
    save_table_to_excel(table_data, excel_path)
    df, non_zero_columns = process_excel_data(excel_path)
    generate_word_document(df, non_zero_columns, word_path)
