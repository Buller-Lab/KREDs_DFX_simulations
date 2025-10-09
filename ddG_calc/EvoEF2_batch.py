import os
import subprocess
import csv
import argparse

def repair_pdb(input_pdb, output_folder):
    repair_command = f"EvoEF2/EvoEF2 --command=RepairStructure --pdb={input_pdb}"
    subprocess.run(repair_command, shell=True)

    # Extract filename without extension
    file_name = os.path.splitext(os.path.basename(input_pdb))[0]
    repair_file = f"{file_name}_Repair.pdb"

    # Move repair file to output folder
    repair_file_path = os.path.join(os.getcwd(), repair_file)
    output_path = os.path.join(output_folder, repair_file)
    os.rename(repair_file_path, output_path)
    return output_path

def compute_stability(input_pdb):
    stability_command = f"EvoEF2/EvoEF2 --command=ComputeStability --pdb={input_pdb}"
    result = subprocess.check_output(stability_command, shell=True, text=True)
    
    # Extract stability value from output
    for line in result.split('\n'):
        if line.startswith("Total"):
            stability = line.split()[2]
            return stability

def process_pdb_files(input_folder, output_folder, output_csv):
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Process each pdb file in input folder
    with open(output_csv, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Filename", "Stability"])

        for file in os.listdir(input_folder):
            if file.endswith(".pdb"):
                pdb_file = os.path.join(input_folder, file)
                repair_file = repair_pdb(pdb_file, output_folder)
                stability = compute_stability(repair_file)
                csv_writer.writerow([file, stability])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files with EvoEF2")
    parser.add_argument("--input_folder", default="EvoEF2_input", 
                       help="Input folder containing PDB files (default: EvoEF2_input)")
    parser.add_argument("--output_folder", default="EvoEF2_output", 
                       help="Output folder for repaired PDB files and CSV (default: EvoEF2_output)")
    args = parser.parse_args()

    output_csv = os.path.join(args.output_folder, "EvoEF2.csv")
    process_pdb_files(args.input_folder, args.output_folder, output_csv)
