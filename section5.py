import re
import os
from docx import Document
from docxcompose.composer import Composer

# Regex header patterns to match
HEADER_KEYWORDS = [
    r"(?i)molecular/?\s*product",
    r"(?i)study number",
    r"(?i)study title",
    r"(?i)test\s*product name",
    r"(?i)active comparator name"
]

# Check if table matches the regex patterns
def table_matches_keywords(table):
    header_text = " ".join(cell.text.strip() for cell in table.rows[0].cells)
    return all(re.search(pattern, header_text) for pattern in HEADER_KEYWORDS)

# Extract content (paragraphs and tables) from source if any table matches
def extract_matching_content(source_path, temp_output_path):
    source_doc = Document(source_path)
    new_doc = Document()
    matched = False

    for block in iter_block_items(source_doc):
        if isinstance(block, TableWrapper):
            if table_matches_keywords(block.table):
                matched = True
                new_doc.add_paragraph("Matching Table Section:")
                copy_table(block.table, new_doc)
        elif isinstance(block, ParagraphWrapper):
            if matched:
                new_doc.add_paragraph(block.paragraph.text)

    if matched:
        new_doc.save(temp_output_path)
    return matched

# Util: Wrap document block items
from docx.table import Table
from docx.text.paragraph import Paragraph

class ParagraphWrapper:
    def __init__(self, paragraph):
        self.paragraph = paragraph

class TableWrapper:
    def __init__(self, table):
        self.table = table

def iter_block_items(parent):
    for child in parent.element.body.iterchildren():
        if child.tag.endswith('}p'):
            yield ParagraphWrapper(Paragraph(child, parent))
        elif child.tag.endswith('}tbl'):
            yield TableWrapper(Table(child, parent))

# Copy table to target doc (without changing formatting)
def copy_table(source_table, target_doc):
    new_table = target_doc.add_table(rows=0, cols=len(source_table.columns))
    new_table.style = source_table.style

    for row in source_table.rows:
        new_row = new_table.add_row()
        for i, cell in enumerate(row.cells):
            new_row.cells[i].text = cell.text

# Merge result into base document using Composer (preserves formatting)
def append_docx_with_matching_page(base_path, source_path, output_path):
    temp_docx = "temp_matching.docx"
    matched = extract_matching_content(source_path, temp_docx)

    if not matched:
        print("❌ No matching content found with regex headers.")
        return

    base_doc = Document(base_path)
    composer = Composer(base_doc)
    matching_doc = Document(temp_docx)

    composer.append(matching_doc)
    composer.save(output_path)

    os.remove(temp_docx)
    print(f"✅ Matching content appended and saved to: {output_path}")

# Example usage
append_docx_with_matching_page(
    base_path="main_report.docx",
    source_path="append_page.docx",  # Full source document to scan
    output_path="merged_output.docx"
)
