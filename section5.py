import os
import tempfile
from docx import Document
from docxcompose.composer import Composer
from PyPDF2 import PdfReader, PdfWriter
from docx2pdf import convert
import subprocess
import platform

def extract_last_pages_from_pdf(pdf_path, output_path, num_pages=2):
    reader = PdfReader(pdf_path)
    writer = PdfWriter()

    total = len(reader.pages)
    for i in range(total - num_pages, total):
        writer.add_page(reader.pages[i])

    with open(output_path, "wb") as f:
        writer.write(f)

def convert_pdf_to_docx(pdf_path, output_path):
    # LibreOffice must be installed
    if platform.system() == "Windows":
        soffice_cmd = r'C:\Program Files\LibreOffice\program\soffice.exe'
    else:
        soffice_cmd = 'libreoffice'

    subprocess.run([
        soffice_cmd,
        '--headless',
        '--convert-to', 'docx',
        '--outdir', os.path.dirname(output_path),
        pdf_path
    ], check=True)

def append_last_pages(base_docx, source_docx, output_docx, pages_to_append=2):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Step 1: Convert source DOCX to PDF
        temp_pdf = os.path.join(tmpdir, "source.pdf")
        convert(source_docx, temp_pdf)

        # Step 2: Extract last N pages from PDF
        last_pages_pdf = os.path.join(tmpdir, "last_pages.pdf")
        extract_last_pages_from_pdf(temp_pdf, last_pages_pdf, num_pages=pages_to_append)

        # Step 3: Convert last 2 pages back to DOCX
        convert_pdf_to_docx(last_pages_pdf, tmpdir)
        extracted_docx = os.path.join(tmpdir, "last_pages.docx")

        # Step 4: Append to base DOCX
        base = Document(base_docx)
        composer = Composer(base)
        composer.append(Document(extracted_docx))
        composer.save(output_docx)
        print(f"✅ Last {pages_to_append} pages appended to: {output_docx}")

# Example usage
append_last_pages(
    base_docx="main_report.docx",
    source_docx="source_pages.docx",
    output_docx="merged_output.docx",
    pages_to_append=2
)
