import re
from docx import Document
from docx.shared import RGBColor, Pt

# Keywords to identify the correct table (case-insensitive, flexible)
HEADER_KEYWORDS = [
    r"(?i)molecular/?\s*product",
    r"(?i)study number",
    r"(?i)study title",
    r"(?i)test\s*product name",
    r"(?i)active comparator name"
]

# Function to check if a table contains all header keywords
def table_matches_keywords(table):
    # Combine all cell texts in the first row to simulate the header
    header_row_text = " ".join(cell.text.strip() for cell in table.rows[0].cells)
    # Check if all keywords are found in the combined header text
    return all(re.search(pattern, header_row_text) for pattern in HEADER_KEYWORDS)

# Function to copy table content to another doc
def copy_table_content(source_table, target_doc):
    # Create a new table with the same number of columns
    new_table = target_doc.add_table(rows=0, cols=len(source_table.columns))
    new_table.style = source_table.style

    for row in source_table.rows:
        new_row = new_table.add_row()
        for i, cell in enumerate(row.cells):
            new_cell = new_row.cells[i]
            new_cell.text = cell.text
            # Optional: Copy basic formatting from the source run
            if cell.paragraphs and cell.paragraphs[0].runs:
                src_run = cell.paragraphs[0].runs[0]
                trg_para = new_cell.paragraphs[0]
                trg_run = trg_para.runs[0]
                trg_run.bold = src_run.bold
                trg_run.italic = src_run.italic
                trg_run.underline = src_run.underline
                trg_run.font.size = src_run.font.size
                trg_run.font.color.rgb = src_run.font.color.rgb

# Main function to extract and append the table
def append_matching_table_with_data(source_path, target_path, output_path):
    source_doc = Document(source_path)
    target_doc = Document(target_path)

    matched = False
    for table in source_doc.tables:
        if table_matches_keywords(table):
            matched = True
            # Add a heading before the table (optional)
            heading = target_doc.add_paragraph("Extracted Table with Header and Data")
            run = heading.runs[0]
            run.font.bold = True
            run.font.size = Pt(14)
            run.font.color.rgb = RGBColor(0, 0, 0)

            # Copy the full table (headers + data)
            copy_table_content(table, target_doc)
            break  # Only the first matching table is copied

    if matched:
        target_doc.save(output_path)
        print(f"✅ Matching table extracted and saved to: {output_path}")
    else:
        print("❌ No matching table found with the specified header keywords.")

# Example usage
append_matching_table_with_data(
    source_path="source.docx",
    target_path="target.docx",
    output_path="output_with_matched_table.docx"
)
