import pandas as pd
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_PARAGRAPH_ALIGNMENT
from docx.oxml.ns import qn

# Update this path to your actual Excel file
excel_path = r"C:\Users\shivam.mishra2\Downloads\New_Psur_File\intial_table.xlsx"
df = pd.read_excel(excel_path, engine='openpyxl')

# Define multi-index columns
columns = pd.MultiIndex.from_tuples([
    ('Molecule/Product', ''),
    ('Study Number', ''),
    ('Study Title', ''),
    ('Test product name', ''),
    ('Active comparator name', ''),
    ('No. of subjects enrolled', 'Test Product'),
    ('No. of subjects enrolled', 'Active Comparator'),
    ('No. of subjects enrolled', 'Placebo'),
    ('No. of subjects enrolled', 'Total'),
    ('Gender distribution of subjects enrolled', 'Male'),
    ('Gender distribution of subjects enrolled', 'Female'),
    ('Age distribution of subjects enrolled', '<18 years'),
    ('Age distribution of subjects enrolled', '18-65 years'),
    ('Age distribution of subjects enrolled', '>65 years'),
    ('Racial distribution of subjects enrolled', 'Asian'),
    ('Racial distribution of subjects enrolled', 'Black'),
    ('Racial distribution of subjects enrolled', 'Caucasian'),
    ('Racial distribution of subjects enrolled', 'Other'),
    ('Racial distribution of subjects enrolled', 'Unknown')
])

df.columns = columns
df_new = df.iloc[2:].copy()
df_new.reset_index(drop=True, inplace=True)

# Convert columns to numeric
columns_to_convert = [
    ('No. of subjects enrolled', 'Test Product'),
    ('No. of subjects enrolled', 'Active Comparator'),
    ('No. of subjects enrolled', 'Placebo'),
    ('No. of subjects enrolled', 'Total'),
    ('Gender distribution of subjects enrolled', 'Male'),
    ('Gender distribution of subjects enrolled', 'Female'),
    ('Age distribution of subjects enrolled', '<18 years'),
    ('Age distribution of subjects enrolled', '18-65 years'),
    ('Age distribution of subjects enrolled', '>65 years'),
    ('Racial distribution of subjects enrolled', 'Asian'),
    ('Racial distribution of subjects enrolled', 'Black'),
    ('Racial distribution of subjects enrolled', 'Caucasian'),
    ('Racial distribution of subjects enrolled', 'Other'),
    ('Racial distribution of subjects enrolled', 'Unknown')
]
df_new[columns_to_convert] = df_new[columns_to_convert].apply(pd.to_numeric, errors='coerce')

# Identify non-zero demographic columns
target_columns = [
    ('Age distribution of subjects enrolled', '<18 years'),
    ('Age distribution of subjects enrolled', '18-65 years'),
    ('Age distribution of subjects enrolled', '>65 years'),
    ('Racial distribution of subjects enrolled', 'Asian'),
    ('Racial distribution of subjects enrolled', 'Black'),
    ('Racial distribution of subjects enrolled', 'Caucasian')
]
non_zero_columns = [col for col in target_columns if df_new[col].sum() != 0]

# Create Word document
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
age_group = non_zero_columns[0][1] if len(non_zero_columns) > 0 else "Unknown"
gender_group = 'Male' if df_new[('Gender distribution of subjects enrolled', 'Male')].sum() > 0 else 'Female'

summary_text = (
    f"Jubilant as MAH has not conducted any Clinical Trials. However, Jubilant has conducted {num_studies} BA/BE study with levetiracetam till the DLP of the PSUR 30-Nov-2024 "
    f"and cumulative subject exposure in the completed clinical trial were {total_subjects} subjects.\n\n"
    f"Of these {total_subjects}, all were {gender_group} subjects of age distribution between {age_group} years. "
    "Cumulative subject exposure to levetiracetam in BA/BE studies is given in the table below:"
)
doc.add_paragraph(summary_text)

df.iloc[2:]
