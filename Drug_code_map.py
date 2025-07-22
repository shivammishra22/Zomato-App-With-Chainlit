import pandas as pd

# === Step 1: Create Hash Map for Product to Dosage Form ===
product_dosage_map = {
    "Esomeprazole": "Gastro-resistant",
    "Zipola 5": "Film coated Tablet",
    "Zipola 10": "Film coated Tablet",
    "Jubilonz OD10": "Oro dispersible tablet",
    "Jubilonz OD5": "Oro dispersible tablet",
    "SCHIZOLANZ": "Oro dispersible tablet",
    "Olanzapine film coated tablets": "Film coated Tablet"
}

# === Step 2: Read Excel File ===
excel_path = r"C:\Users\shivam.mishra2\Downloads\New_Psur_File\marketing_exposure_tables.xlsx"

try:
    df = pd.read_excel(excel_path, engine='openpyxl')
except Exception as e:
    print(f"❌ Error reading Excel: {e}")
    exit()

# === Step 3: Clean and Match Product Names ===
def map_dosage(product_name):
    for key in product_dosage_map:
        if key.lower() in str(product_name).lower():
            return product_dosage_map[key]
    return ""

df["Dosage Form (Units)"] = df["Product"].apply(map_dosage)

# === Step 4: Save Updated Excel File ===
try:
    df.to_excel(excel_path, index=False)
    print(f"✅ Updated Excel saved with 'Dosage Form (Units)' column at: {excel_path}")
except Exception as e:
    print(f"❌ Error saving Excel: {e}")
    
