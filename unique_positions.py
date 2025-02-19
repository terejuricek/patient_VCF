from pathlib import Path

my_path = Path("/Parth/to/vcf_positions/testing/")  # Convert to Path object
folder_pattern = ["EP_*", "LP_*"]

# all folder names with the pattern without the need to use glob
patients = [folder.name for folder in my_path.iterdir() if any(folder.match(pattern) for pattern in folder_pattern)]
print(patients)

# create a new folder in my_path to store the output files uding Path module
output_dir = my_path / "output" 
output_dir.mkdir(parents=True, exist_ok=True)

for folder in patients:
    patient_dir = my_path / folder
    # create a new txt file with the patient name where the unique positions will be stored
    output_file = output_dir / f"{folder}.txt"
    unique_positions = set()
    for vcf_file in patient_dir.glob("*.vcf"):
        # open and process the vcf file - itereate over the lines and extract the chrom and pos
        # add the unique positions to the txt file
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                columns = line.strip().split('\t')
                if len(columns) < 2:
                    continue
                chrom, pos = columns[0], columns[1]
                # add the pair LIST [chrom, pos] to the set
                unique_positions.add((chrom, pos))
    # write the unique positions to the output file
    with open(output_file, 'w') as f:
        for chrom, pos in sorted(unique_positions, key=lambda x: (x[0], (x[1]))):
            f.write(f"{chrom}\t{pos}\n")
