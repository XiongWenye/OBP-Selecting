import pandas as pd

csv_file_path = 'Compound_OBP_binding.csv'
data = pd.read_csv(csv_file_path)

compound_names = [
    'butanoic acid, ethyl ester / ethyl butyrate',
    'dicyclohexyl ketone',
    'undecanal',
    'acetic acid',
    '2-heptanone',
    'hexanal',
    'cyclopentanone',
    '3-hydroxy-2-butanone',
    '2-pentanone'
]

def find_binding_data(data, compound_names):
    binding_data = {}
    for compound in compound_names:
        compound_data = data[data['Compound name'].str.contains(compound, case=False, na=False)]
        if not compound_data.empty:
            binding_data[compound] = compound_data.iloc[0].to_dict()
    return binding_data

def process_binding_value(value):
    if isinstance(value, str) and '>' in value:
        return float(value.replace('>', '')) + 0.1
    try:
        return float(value)
    except ValueError:
        return None

binding_data = find_binding_data(data, compound_names)

def select_best_obp(binding_data):
    obp_scores = {}
    for compound, data in binding_data.items():
        for obp, value in data.items():
            if obp not in ['CAS-number', 'Compound name'] and pd.notna(value):
                processed_value = process_binding_value(value)
                if processed_value is not None:
                    obp_scores[obp] = obp_scores.get(obp, 0) + processed_value
    sorted_obps = sorted(obp_scores.items(), key=lambda item: item[1], reverse=True)
    return sorted_obps

best_obps = select_best_obp(binding_data)

best_obp_df = pd.DataFrame(best_obps, columns=['OBP', 'Score'])
output_file_path = 'Best_OBP_Combination.csv'
best_obp_df.to_csv(output_file_path, index=False)

print(f"Best OBP combination has been saved to {output_file_path}")