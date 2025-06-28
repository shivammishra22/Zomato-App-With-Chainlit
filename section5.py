import os
import tempfile
import subprocess
import platform
from docx import Document
from docxcompose.composer import Composer
from docx2pdf import convert
from PyPDF2 import PdfReader, PdfWriter

def convert_pdf_to_docx(pdf_path, output_docx_path):
    # Path to LibreOffice
    if platform.system() == "Windows":
        soffice = r"C:\Program Files\LibreOffice\program\soffice.exe"
    else:
        soffice = "libreoffice"

    subprocess.run([
        soffice,
        "--headless",
        "--convert-to", "docx",
        "--outdir", os.path.dirname(output_docx_path),
        pdf_path
    ], check=True)

def extract_last_page_and_append(base_path, append_path, output_path):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Step 1: Convert base docx to PDF
        base_pdf = os.path.join(tmpdir, "base.pdf")
        convert(base_path, base_pdf)

        # Step 2: Extract last page of PDF
        reader = PdfReader(base_pdf)
        writer = PdfWriter()
        writer.add_page(reader.pages[-1])  # Last page only
        last_page_pdf = os.path.join(tmpdir, "last_page.pdf")
        with open(last_page_pdf, "wb") as f:
            writer.write(f)

        # Step 3: Convert last page PDF back to DOCX
        convert_pdf_to_docx(last_page_pdf, tmpdir)
        last_page_docx = os.path.join(tmpdir, "last_page.docx")

        # Step 4: Use docxcompose to append last page to append_path.docx
        append_doc = Document(append_path)
        composer = Composer(append_doc)
        composer.append(Document(last_page_docx))
        composer.save(output_path)

        print(f"✅ Last page from '{base_path}' appended to '{append_path}' → saved as '{output_path}'")

# Example usage
extract_last_page_and_append(
    base_path="main_report.docx",       # Extract last page from this
    append_path="append_this.docx",     # Append to this
    output_path="merged_result.docx"    # Final merged output
)
