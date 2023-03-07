import argparse
import csv
import os

parser = argparse.ArgumentParser(description='Aggregate CSV files in multiple directories.')
parser.add_argument('input_directory', type=str, help='Path to the root directory containing input CSV files.')
parser.add_argument('output_file', type=str, help='Path to the output CSV file.')
args = parser.parse_args()

all_csv_files = []

for directory, subdirectories, files in os.walk(args.input_directory):
    for file in files:
        if file.endswith('_confindr_report.csv'):
            file_path = os.path.join(directory, file)
            all_csv_files.append(file_path)

with open(args.output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    first_file = True
    for file_path in all_csv_files:
        with open(file_path, 'r') as infile:
            reader = csv.reader(infile)
            skip_header = False
            for row in reader:
                if skip_header and not first_file:
                    writer.writerow(row)
                elif first_file:
                    writer.writerow(row)
                    first_file = False
                skip_header = True
