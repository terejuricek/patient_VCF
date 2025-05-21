# additional script to asses what posistions need to be found in CRAM

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import fnmatch

my_path = Path("/Users/terezajurickova/Desktop/vcf_positions/testing/")
folder_pattern = ["EP_*", "LP_*"]

# filter folders by pattern
patients = [folder.name for folder in my_path.iterdir() if folder.is_dir() and any(fnmatch.fnmatch(folder.name, pattern) for pattern in folder_pattern)]
print(patients)

# create output directory
output_dir = my_path / "output"
output_dir.mkdir(parents=True, exist_ok=True)

for folder in patients:
    patient_dir = my_path / folder
    output_file = output_dir / f"{folder}.txt"
    unique_positions = set()

    for vcf_file in patient_dir.glob("*.vcf"):
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                columns = line.strip().split('\t')
                if len(columns) < 2:
                    continue
                chrom, pos = columns[0], columns[1]
                unique_positions.add((chrom, pos))

    with open(output_file, 'w') as f:
        for chrom, pos in sorted(unique_positions, key=lambda x: (x[0], int(x[1]))):
            f.write(f"{chrom}\t{pos}\n")