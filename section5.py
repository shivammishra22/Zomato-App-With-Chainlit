from docx import Document
from config.styling_config import DocumentStyling

# Create a new Document
doc = Document()

# Utility functions using DocumentStyling
def add_styled_section_heading(text):
    para = doc.add_paragraph()
    DocumentStyling.create_split_heading(para, text)

def add_styled_subheading(text):
    para = doc.add_paragraph()
    DocumentStyling.create_split_subheading(para, text)

def add_styled_paragraph(text):
    para = doc.add_paragraph()
    run = para.add_run(text)
    DocumentStyling.apply_content_style(run)

# ----------------- Section 7 ----------------- #
add_styled_section_heading("7 SUMMARIES OF SIGNIFICANT FINDINGS FROM CLINICAL TRIALS DURING THE REPORTING INTERVAL")

# 7.1
add_styled_subheading("7.1 Completed clinical trials")
add_styled_paragraph("MAH did not conduct or completed any clinical trials during the reporting period covered by this report.")

# 7.2
add_styled_subheading("7.2 Ongoing clinical trials")
add_styled_paragraph("No clinical trials were ongoing during the period covered by this report.")

# Metadata
add_styled_paragraph(
    "Version status: Final v 1.0 Confidential LEVE-UA-002\n"
    "Version date: 18-Feb-2025 Page 39 of 65\n\n"
    "Jubilant Generics Limited Levetiracetam Periodic Safety Update Report\n"
    "Reporting period: 30-Nov-2021 to 30-Nov-2024"
)

# 7.3
add_styled_subheading("7.3 Long term follow up")
add_styled_paragraph("No long-term studies were conducted during the period covered by this report.")

# 7.4
add_styled_subheading("7.4 Other therapeutic use of medicinal product")
add_styled_paragraph(
    "No such information is available as the MAH did not conduct any programmes that follow a specific protocol, "
    "with solicited reporting as per ICH-E2D (e.g., expanded access programmes, compassionate use programmes, "
    "particular patient use, and other organised data collection)."
)

# 7.5
add_styled_subheading("7.5 New safety data related to fixed combination therapies")
add_styled_paragraph("No clinical study with levetiracetam in other combination was conducted during the reporting interval.")

# ----------------- Section 8 ----------------- #
add_styled_section_heading("8 FINDINGS FROM NON-INTERVENTIONAL STUDIES")
add_styled_paragraph(
    "The MAH has not conducted any non-interventional study related to the current product under consideration "
    "since obtaining initial granting of MA. Hence, no such information is available."
)

# Save the document
doc.save("sec7_8.docx")
print("✅ Document 'sec7_8.docx' created using styling_config.")
