# import json
# import os

# def process_geojson(input_file, output_file):
#     with open(input_file, 'r', encoding='utf-8') as f:
#         data = json.load(f)
    
#     new_features = []
    
#     for feature in data['features']:
#         new_feature = {
#             'type': 'Feature',
#             'geometry': feature['geometry'],
#             'properties': {
#                 'name': feature['properties']['name']
#             }
#         }
#         new_features.append(new_feature)
    
#     new_geojson = {
#         'type': 'FeatureCollection',
#         'features': new_features
#     }
    
#     with open(output_file, 'w', encoding='utf-8') as f:
#         json.dump(new_geojson, f, ensure_ascii=False, indent=2)

# # Directory containing the script and the input files
# input_dir = os.path.dirname(os.path.abspath(__file__))

# # Create a new directory to save the output files
# output_dir = os.path.join(input_dir, 'test')
# os.makedirs(output_dir, exist_ok=True)

# # Process each file in the input directory
# for filename in os.listdir(input_dir):
#     if filename.endswith('.json'):
#         input_file = os.path.join(input_dir, filename)
#         output_file = os.path.join(output_dir, filename)
#         process_geojson(input_file, output_file)
#         print(f'Processed {filename}')

import json
import os

def process_geojson(input_files, output_file):
    all_features = []
    
    for input_file in input_files:
        with open(input_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        for feature in data['features']:
            new_feature = {
                'type': 'Feature',
                'geometry': feature['geometry'],
                'properties': {
                    'name': feature['properties']['name']
                }
            }
            all_features.append(new_feature)
    
    combined_geojson = {
        'type': 'FeatureCollection',
        'features': all_features
    }
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(combined_geojson, f, ensure_ascii=False, indent=2)

# Directory containing the script and the input files
input_dir = os.path.dirname(os.path.abspath(__file__))

# Collect all input files
input_files = [os.path.join(input_dir, filename) for filename in os.listdir(input_dir) if filename.endswith('.json')]

# Define the output file path
output_file = os.path.join(input_dir, 'test.json')

# Process the input files and generate the combined output file
process_geojson(input_files, output_file)
print(f'Combined GeoJSON file created: {output_file}')

