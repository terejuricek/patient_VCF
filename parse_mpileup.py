import re

input = "chr1  1000    A       10      .,Aa^].$aG      <<<<<<<<<"
input = "chr1  1000    A       10      gg.cc**C+3aaaG**^]T      <<<<<<<<<"
print("i+nput".find("+"))

base_string = "cc+3aaaTT-2aaG-1aCC+2a^]T"



def parse_bases(base_string, ref_base):

    base_counts = {'A': 0, 
                   'C': 0, 
                   'G': 0, 
                   'T': 0, 
                   'N': 0}
    
    # RULES for parsing:
        # \+[0-9]+[ACGTNacgtn*#]+
        # \-[0-9]+[ACGTNacgtn*#]+
        # \^.
        # $
    
    indel = 0
    while indel != -1:
        insert = base_string.find("+")
        # print(insert)
        if insert != -1:
            skip = int(base_string[insert+1])
            base_string = str(base_string[:insert]) + str(base_string[insert + 2 + skip:])
            # print(base_string)

        delet = base_string.find("-")
        if delet != -1:
            skip = int(base_string[delet+1])    
            base_string = str(base_string[:delet]) + str(base_string[delet + 2 + skip:])
            # print(base_string)
        
        if insert == -1 and delet == -1:
            indel = -1
        
        
    base_string = re.sub(r'\^.', '', base_string)                       # removes the ^ character and one imidiately following base - quality
    base_string = base_string.replace('$', '')                          # removes the $ character
    
    base_string = re.sub(r'[\+\-](\d+)[ACGTNacgtn]+', '', base_string)  # removes <> and anything in between - insertion and deletion characters

    # print(base_string)
    
    for base in base_string:
        if base in '.,':
            base_counts[ref_base.upper()] += 1
        elif base.upper() in base_counts:
            base_counts[base.upper()] += 1
    
    return base_counts
    
def compute_frequencies(counts):
    total = sum(counts.values())
    return {allele: count / total for allele, count in counts.items()} if total else {}


print(parse_bases(input.split()[4], input.split()[2]))
print(compute_frequencies(parse_bases(input.split()[4], input.split()[2])))


##################################################################################################################################




import sys
import re
from collections import Counter

def parse_mpileup_line(line):
    """Extract relevant information from a SAMtools mpileup line."""
    fields = line.strip().split(" ")
    chrom, pos, ref_base, depth, bases, qual_scores = fields[:6]
    print(bases)
    indel = 0
    while indel != -1:
        insert = bases.find("+")
        # print(insert)
        if insert != -1:
            skip = int(bases[insert+1])
            bases = str(bases[:insert]) + str(bases[insert + 2 + skip:])
            # print(bases)

        delet = bases.find("-")
        if delet != -1:
            skip = int(bases[delet+1])    
            bases = str(bases[:delet]) + str(bases[delet + 2 + skip:])
            # print(bases)
        
        if insert == -1 and delet == -1:
            indel = -1
    
    # Count occurrences of each base
    bases = re.sub(r'\^.', '', bases)  # Remove start-of-read markers
    bases = re.sub(r'[\$,]', '', bases)  # Remove end-of-read markers
    
    base_counts = Counter(bases.upper())  # Normalize to uppercase
    print(chrom, pos, ref_base, depth, base_counts)
    return chrom, pos, ref_base, depth, base_counts

def determine_genotype(base_counts, ref_base):
    """Calculate allele frequencies and determine genotype."""
    total_bases = sum(base_counts.values())
    
    if total_bases == 0:
        return "./."  # No calls
    
    sorted_bases = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Identify most frequent base (alternative allele)
    alt_base = sorted_bases[0][0] if sorted_bases[0][0] != ref_base else sorted_bases[1][0]
    alt_freq = base_counts.get(alt_base, 0) / total_bases

    # Assign genotype based on allele frequency thresholds
    if alt_freq >= 0.8:
        genotype = "1/1"  # Homozygous alternate
    elif alt_freq >= 0.2:
        genotype = "0/1"  # Heterozygous
    else:
        genotype = "0/0"  # Homozygous reference
    
    return genotype

def process_mpileup(input_file):
    """Read and process SAMtools mpileup output file."""
    print("#CHROM\tPOS\tREF\tGT")
    try:
        with open(input_file, 'r') as f:
            for line in f:
                chrom, pos, ref_base, depth, base_counts = parse_mpileup_line(line)
                genotype = determine_genotype(base_counts, ref_base)
                print(f"{chrom}\t{pos}\t{ref_base}\t{genotype}")
    except FileNotFoundError:
        chrom, pos, ref_base, depth, base_counts = parse_mpileup_line(input_file)
        genotype = determine_genotype(base_counts, ref_base)
        print(f"{chrom}\t{pos}\t{ref_base}\t{genotype}")

# if __name__ == "__main__":
#     if len(sys.argv) != 2:
#         print("Usage: python variant_call.py <mpileup_file>")
#         sys.exit(1)

#     process_mpileup(sys.argv[1])

input = "chr1 1000 A 10 gg.-1cc**C+3aaaG**^]T <<<<<<<<<"

print("#########################")
print(parse_mpileup_line(input))
print("#########################")
print(process_mpileup(input))
print("#########################")
