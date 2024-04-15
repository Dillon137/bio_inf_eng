filename = "rosalind_revp.txt"

DNA_comp = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

RNA_to_proteins = {
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L', 'UAC': 'Y',
    'UAU': 'Y', 'UAA': '*', 'UAG': '*', 'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W', 'CUA': 'L',
    'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P', 'CAC': 'H', 'CAU': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'AUA': 'I', 'AUC': 'I', 'AUU': 'I',
    'AUG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T', 'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R', 'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'GCA': 'A',
    'GCC': 'A', 'GCG': 'A', 'GCU': 'A', 'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G',
    'GGG': 'G', 'GGU': 'G'
}

mono_iso_mass = {
    'A':   71.03711,
    'C':   103.00919,
    'D':   115.02694,
    'E':   129.04259,
    'F':   147.06841,
    'G':   57.02146,
    'H':   137.05891,
    'I':   113.08406,
    'K':   128.09496,
    'L':   113.08406,
    'M':   131.04049,
    'N':   114.04293,
    'P':   97.05276,
    'Q':   128.05858,
    'R':   156.10111,
    'S':   87.03203,
    'T':   101.04768,
    'V':   99.06841,
    'W':   186.07931,
    'Y':   163.06333,
}


def read_file(filename):
    with open(f"./data/{filename}", 'r') as f:
        text = f.readlines()
        f.close()
    text = [i.strip() for i in text]
    return text


def DNA_to_RNA(dna):
    dna = dna.replace('T', 'U')
    return dna


def RNA_to_DNA(rna):
    rna = rna.replace('U', 'T')
    return rna


def DNA_compliment(dna):
    trans = str.maketrans(DNA_comp)
    dna_c = dna.translate(trans)[::-1]
    return dna_c


def find_motif(string, motif):
    motifs = [i+1 for i in range(len(string)) if string.startswith(motif, i)]
    return motifs


def find_shared_motifs(strings):
    c = 1
    start = 0
    match = False
    longest = ""
    for i, m in enumerate(strings[0]):
        test = strings[0][start:c]

        for string in strings[1:]:
            if string.find(test) == -1:
                start += 1
                c += 1
                match = False
                break

            match = True

        if match:
            longest = strings[0][start:c]
            c += 1

    return longest


def locate_restriction_sites(dna):
    sites = []
    r_comp = DNA_compliment(dna)

    for n in range(4, 13):
        for i in range(len(dna)-(n-1)):
            if r_comp[-(i+n):-(i)] != dna[i:i+n]:
                continue

            sites.append([i+1, n])

    return sites


def remove_exons(dna, exons):
    for s in exons:
        dna = dna.replace(s, "")

    return dna


def protein_translation(rna):
    groups = [rna[i:i+3] for i in range(0, len(rna), 3)]
    proteins = [RNA_to_proteins[group] for group in groups]
    return proteins


def calc_protein_mass(protein_str):
    mass = 0
    for p in protein_str:
        mass += mono_iso_mass[p]
    return mass.__round__(3)


def h_dist(dna_1, dna_2):
    c = 0
    for idx, base in enumerate(dna_1):
        if base != dna_2[idx]:
            c += 1
    return c


text = read_file(filename)
text = ''.join(text)
text = text.split(">Rosalind_")[1:]
text = [i[4:] for i in text]

dna = text[0]

for i in locate_restriction_sites(dna):
    print(*i)