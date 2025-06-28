from docx import Document
from docx.shared import RGBColor

# Create a new Document
doc = Document()

# Function to add heading with black color
def add_heading(doc, text, level):
    heading = doc.add_heading(text, level=level)
    run = heading.runs[0]
    run.font.color.rgb = RGBColor(0, 0, 0)

# Section 7
add_heading(doc, "7 SUMMARIES OF SIGNIFICANT FINDINGS FROM CLINICAL TRIALS DURING THE REPORTING INTERVAL", level=2)

# Subsection 7.1
add_heading(doc, "7.1 Completed clinical trials", level=3)
doc.add_paragraph("MAH did not conduct or completed any clinical trials during the reporting period covered by this report.")

# Subsection 7.2
add_heading(doc, "7.2 Ongoing clinical trials", level=3)
doc.add_paragraph("No clinical trials were ongoing during the period covered by this report.")

# Metadata
doc.add_paragraph(
    "Version status: Final v 1.0 Confidential LEVE-UA-002\n"
    "Version date: 18-Feb-2025 Page 39 of 65\n\n"
    "Jubilant Generics Limited Levetiracetam Periodic Safety Update Report\n"
    "Reporting period: 30-Nov-2021 to 30-Nov-2024"
)

# Subsection 7.3
add_heading(doc, "7.3 Long term follow up", level=3)
doc.add_paragraph("No long-term studies were conducted during the period covered by this report.")

# Subsection 7.4
add_heading(doc, "7.4 Other therapeutic use of medicinal product", level=3)
doc.add_paragraph(
    "No such information is available as the MAH did not conduct any programmes that follow a specific protocol, "
    "with solicited reporting as per ICH-E2D (e.g., expanded access programmes, compassionate use programmes, "
    "particular patient use, and other organised data collection)."
)

# Subsection 7.5
add_heading(doc, "7.5 New safety data related to fixed combination therapies", level=3)
doc.add_paragraph("No clinical study with levetiracetam in other combination was conducted during the reporting interval.")

# Section 8
add_heading(doc, "8 FINDINGS FROM NON-INTERVENTIONAL STUDIES", level=2)
doc.add_paragraph(
    "The MAH has not conducted any non-interventional study related to the current product under consideration "
    "since obtaining initial granting of MA. Hence, no such information is available."
)

# Save the document
doc.save("sec7_8.docx")
