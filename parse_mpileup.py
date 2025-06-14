import re
from collections import Counter

def clean_bases(bases):
    # Remove start-of-read (^.) and end-of-read ($)
    bases = re.sub(r'\^.', '', bases)
    bases = bases.replace('$', '')
    
    # Remove insertions and deletions using regex
    bases = re.sub(r'[\+\-](\d+)[ACGTNacgtn]+', '', bases)
    
    return bases

def parse_mpileup_line(line):
    fields = re.split(r'\s+', line.strip())
    chrom, pos, ref_base, depth, bases, _ = fields[:6]
    
    bases = clean_bases(bases)

    # Convert . and , to reference base
    base_counts = Counter()
    for base in bases:
        if base in '.,':
            base_counts[ref_base.upper()] += 1
        elif base.upper() in 'ACGTN':
            base_counts[base.upper()] += 1

    return chrom, pos, ref_base.upper(), int(depth), base_counts

def determine_genotype(base_counts, ref_base):
    total = sum(base_counts.values())
    if total == 0:
        return './.'

    sorted_bases = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
    alt_base = next((b for b, _ in sorted_bases if b != ref_base), None)
    if not alt_base:
        return '0/0'

    alt_freq = base_counts[alt_base] / total
    if alt_freq >= 0.8:
        return '1/1'
    elif alt_freq >= 0.2:
        return '0/1'
    else:
        return '0/0'

# Test single line
line = "chr1 1000 A 10 gg.-1cc**C+3aaaG**^]T <<<<<<<<<"
chrom, pos, ref, depth, counts = parse_mpileup_line(line)
gt = determine_genotype(counts, ref)
print(f"{chrom}\t{pos}\t{ref}\t{gt}")


def process_mpileup(input_file):
    """Read and process SAMtools mpileup output file."""
    try:
        with open(input_file, 'r') as f:
            for line in f:
                chrom, pos, ref_base, depth, base_counts = parse_mpileup_line(line)
                genotype = determine_genotype(base_counts, ref_base)

                alt_base = max((b for b in base_counts if b != ref_base.upper()), 
                               key=lambda b: base_counts[b], default=".")
                
                print(f"{chrom}\t{pos}\t.\t{ref_base.upper()}\t{alt_base}\t.\t.\t.\tGT\t{genotype}")
    except FileNotFoundError:
        chrom, pos, ref_base, depth, base_counts = parse_mpileup_line(input_file)
        genotype = determine_genotype(base_counts, ref_base)

        alt_base = max((b for b in base_counts if b != ref_base.upper()), 
                       key=lambda b: base_counts[b], default=".")

        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")
        print(f"{chrom}\t{pos}\t.\t{ref_base.upper()}\t{alt_base}\t.\t.\t.\tGT\t{genotype}")

print("Processing single line input:")
print(process_mpileup("chr1 1000 A 10 gg.-1cc**C+3aaaG**^]T <<<<<<<<<"))