import os
import re
import pandas as pd
import numpy as np
from docx import Document
import requests
from bs4 import BeautifulSoup
import urllib3

# Suppress InsecureRequestWarning
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def fetch_ddd_fallback(medicine, code):
    # Use your earlier fallback BeautifulSoup code
    if pd.isna(code):
        return np.nan
    url = f"https://atcddd.fhi.no/atc_ddd_index/?code={code}"
    try:
        response = requests.get(url, verify=False, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, "html.parser")
        ddd_values = soup.find_all("td", align="right")
        for td in ddd_values:
            value = td.get_text(strip=True)
            if value.replace('.', '', 1).isdigit():
                print(f"Fetched fallback DDD for {medicine} ({code}): {value}")
                return float(value)
        print(f"No DDD value found on fallback site for {medicine} ({code})")
        return np.nan
    except Exception as e:
        print(f"Error during fallback DDD fetch for {medicine} ({code}): {e}")
        return np.nan

# Your existing extract_table_after_text, save_table_to_excel, generate_fallback_doc, and calculate_exposure_and_generate_doc remain unchanged.

# === MAIN EXECUTION ===

if __name__ == "__main__":
    docx_path = r"C:\Users\shivam.mishra2\Downloads\ALL_PSUR_File\PSUR_all _Data\Olanzapine PSUR_South Africa_29-Sep-17 to 31-Mar-25\Draft\DRA\Data request form_olanzapine.docx"
    search_text = "Cumulative sales data sale required"
    excel_output_path = r"C:\Users\shivam.mishra2\Downloads\New_Psur_File\marketing_exposure_tables.xlsx"

    ddd_excel_path = "drug_code_map_with_ddd.xlsx"  # Path to your DDD Excel (from earlier code)
    country = "South Africa"
    medicine = "Esomeprazole"
    place = "South Africa"
    date = "2020-01-01"

    # STEP 1: Read the DDD Excel
    try:
        ddd_df = pd.read_excel(ddd_excel_path)
    except Exception as e:
        print(f"⚠️ Could not read DDD Excel: {e}")
        ddd_df = pd.DataFrame()

    try:
        if not os.path.exists(docx_path):
            raise FileNotFoundError(f"File not found: {docx_path}")
        doc = Document(docx_path)
        table_data = extract_table_after_text(doc, search_text)

        # --- DDD VALUE LOGIC START ---
        # Try Excel first
        ddd_row = ddd_df[ddd_df["Drug Name"].str.lower() == medicine.lower()] if not ddd_df.empty else pd.DataFrame()
        if not ddd_row.empty and not pd.isna(ddd_row.iloc[0]["DDD Value"]):
            ddd_value = ddd_row.iloc[0]["DDD Value"]
            print(f"✅ DDD value found in Excel: {ddd_value}")
        else:
            # Try fallback using code if available
            if not ddd_row.empty:
                code = ddd_row.iloc[0]["Drug Code"]
            else:
                # Try to get Drug Code from the Excel if name not found
                code = np.nan
                if "Drug Code" in ddd_df.columns:
                    possible = ddd_df[ddd_df["Drug Name"].str.lower().str.contains(medicine.lower())]
                    if not possible.empty:
                        code = possible.iloc[0]["Drug Code"]
            ddd_value = fetch_ddd_fallback(medicine, code)
            if pd.isna(ddd_value):
                print("❌ DDD value not found anywhere. Will use fallback document.")
        # --- DDD VALUE LOGIC END ---

        if table_data and not pd.isna(ddd_value):
            save_table_to_excel(table_data, excel_output_path)
            calculate_exposure_and_generate_doc(excel_output_path, ddd_value, country, medicine, place, date)
        else:
            print("⚠️ Table not found after search text or DDD is missing. Generating fallback document.")
            generate_fallback_doc(medicine)
    except Exception as e:
        print(f"⚠️ Error: {e}")
        generate_fallback_doc(medicine)
