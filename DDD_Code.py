import requests
from bs4 import BeautifulSoup
import urllib3
from drug_code_map import drug_code_map

# Suppress only the InsecureRequestWarning
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def fetch_ddd(drug_name):
    # Get the ATC code from the drug name
    code = drug_code_map.get(drug_name.lower())
    if not code:
        print(f"Drug '{drug_name}' not found in the map.")
        return float('nan')

    # Construct the URL
    url = f"https://atcddd.fhi.no/atc_ddd_index/?code={code}"

    # Send request
    response = requests.get(url, verify=False)
    soup = BeautifulSoup(response.content, "html.parser")

    # Find the DDD value
    ddd_values = soup.find_all("td", align="right")
    for td in ddd_values:
        value = td.get_text(strip=True)
        # Check if the value is numeric
        if value.replace('.', '', 1).isdigit():
            print(f"DDD Value for {drug_name} ({code}): {value}")
            return float(value)

    print(f"No numeric DDD value found for {drug_name} ({code})")
    return float('nan')

# Example usage
ddd_results = {}
for drug in drug_code_map:
    ddd_results[drug] = fetch_ddd(drug)

# Optional: print the final dictionary
print(ddd_results)
