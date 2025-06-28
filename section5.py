import re
import os
import tempfile
from docx import Document
from docx.shared import Inches
from docx2pdf import convert
from pdf2image import convert_from_path

# Regex keywords to match table headers
HEADER_KEYWORDS = [
    r"(?i)molecular/?\s*product",
    r"(?i)study number",
    r"(?i)study title",
    r"(?i)test\s*product name",
    r"(?i)active comparator name"
]

# Check if a table matches header keywords
def table_matches_keywords(table):
    header_text = " ".join(cell.text.strip() for cell in table.rows[0].cells)
    return all(re.search(pattern, header_text) for pattern in HEADER_KEYWORDS)

# Extract matching tables into a new temporary docx file
def extract_matching_tables(source_docx, temp_filtered_docx):
    doc = Document(source_docx)
    new_doc = Document()

    matched = False
    for table in doc.tables:
        if table_matches_keywords(table):
            matched = True
            new_doc.add_paragraph("Extracted Table with Matching Header:")
            new_table = new_doc.add_table(rows=0, cols=len(table.columns))
            new_table.style = table.style

            for row in table.rows:
                new_row = new_table.add_row()
                for i, cell in enumerate(row.cells):
                    new_row.cells[i].text = cell.text

            new_doc.add_page_break()

    if matched:
        new_doc.save(temp_filtered_docx)
        return True
    return False

# Convert filtered docx to PDF, then to images, then insert into target docx
def insert_matching_table_images(source_docx, target_docx, output_docx):
    with tempfile.TemporaryDirectory() as tmpdir:
        filtered_docx = os.path.join(tmpdir, "filtered.docx")

        # Step 1: Extract matching tables
        has_match = extract_matching_tables(source_docx, filtered_docx)
        if not has_match:
            print("❌ No matching tables found with given headers.")
            return

        # Step 2: Convert filtered docx to PDF
        filtered_pdf = os.path.join(tmpdir, "filtered.pdf")
        convert(filtered_docx, filtered_pdf)

        # Step 3: Convert PDF to images (screenshots)
        images = convert_from_path(filtered_pdf)

        # Step 4: Insert images into target docx
        target = Document(target_docx)
        for i, img in enumerate(images):
            img_path = os.path.join(tmpdir, f"table_{i}.png")
            img.save(img_path)
            target.add_picture(img_path, width=Inches(6))
            target.add_page_break()

        target.save(output_docx)
        print(f"✅ Matching tables added as images to: {output_docx}")

# Example usage
insert_matching_table_images(
    source_docx="source.docx",
    target_docx="target.docx",
    output_docx="output_with_table_screenshots.docx"
)
