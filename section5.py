from docx import Document
from docx.shared import RGBColor, Pt
from copy import deepcopy

# Function to copy table from source to target document
def append_table_with_format(source_path, target_path, output_path, table_index=0):
    source_doc = Document(source_path)
    target_doc = Document(target_path)

    # Select the specific table from source
    table = source_doc.tables[table_index]

    # Append a heading before the table (optional)
    heading = target_doc.add_paragraph("Appended Table with Sub-Headings")
    run = heading.runs[0]
    run.font.bold = True
    run.font.size = Pt(14)
    run.font.color.rgb = RGBColor(0, 0, 0)

    # Append table
    new_table = target_doc.add_table(rows=0, cols=len(table.columns))
    new_table.style = table.style

    for row in table.rows:
        new_row = new_table.add_row()
        for i, cell in enumerate(row.cells):
            new_cell = new_row.cells[i]
            new_cell.text = cell.text
            # Optional: Copy font styles (bold, size, etc.) if needed
            if cell.paragraphs and cell.paragraphs[0].runs:
                src_run = cell.paragraphs[0].runs[0]
                trg_para = new_cell.paragraphs[0]
                trg_run = trg_para.runs[0]
                trg_run.bold = src_run.bold
                trg_run.italic = src_run.italic
                trg_run.underline = src_run.underline
                trg_run.font.size = src_run.font.size
                trg_run.font.color.rgb = src_run.font.color.rgb

    # Save output
    target_doc.save(output_path)
    print(f"✅ Table successfully appended to: {output_path}")

# Example usage
append_table_with_format(
    source_path="source.docx",
    target_path="target.docx",
    output_path="merged_output.docx",
    table_index=0  # Index of the table in source file
)
