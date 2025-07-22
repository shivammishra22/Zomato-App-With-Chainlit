import pandas as pd
import requests
from bs4 import BeautifulSoup
import urllib3
import numpy as np

# Suppress only the InsecureRequestWarning
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Read the Excel file (make sure the path is correct)
df = pd.read_excel("drug_code_map.xlsx")  # The input Excel file

def fetch_ddd(drug_name, drug_code):
    url = f"https://atcddd.fhi.no/atc_ddd_index/?code={drug_code}"
    try:
        response = requests.get(url, verify=False, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, "html.parser")
        ddd_values = soup.find_all("td", align="right")
        for td in ddd_values:
            value = td.get_text(strip=True)
            # Only consider if the value is numeric
            if value.replace('.', '', 1).isdigit():
                print(f"DDD Value for {drug_name} ({drug_code}): {value}")
                return float(value)
        print(f"No numeric DDD value found for {drug_name} ({drug_code})")
        return np.nan
    except Exception as e:
        print(f"Error for {drug_name} ({drug_code}): {e}")
        return np.nan

# Add a new column to the DataFrame with the fetched DDD values
df['DDD Value'] = df.apply(lambda row: fetch_ddd(row['Drug Name'], row['Drug Code']), axis=1)

# Save to a new Excel file
df.to_excel("drug_code_map_with_ddd.xlsx", index=False)
print("Results saved in 'drug_code_map_with_ddd.xlsx'")
