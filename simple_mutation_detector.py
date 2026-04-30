from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner  # FIX 7: was missing


def inputFile():
    return input("\nEnter the filepath of the DNA file: ")


def readTheFile(filePath):
    sequence = ""  # FIX 1: initialize sequence before using it
    try:
        with open(filePath, "r") as f:
            for line in f:
                if not line.startswith(">"):  # FIX 2: was checking '<' instead of '>'
                    sequence += line.strip()

        sequence = sequence.upper()

        if set(sequence) - set("ATCG"):  # FIX 3: ATCG must be a string in quotes
            print("\nWarning: Non-DNA characters found")

        return sequence

    except FileNotFoundError:
        print("\nFilepath not found")
        return None


def runBlast(sequence):
    print("\nRunning BLAST...")

    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nt",
        sequence=sequence   # FIX 4: was using undefined 'dna_sequence', use parameter name
    )

    with open("blast_result.xml", "w") as out_file:
        out_file.write(result_handle.read())

    print("\nBLAST completed and saved to file")


def get_top_hit():
    try:
        with open("blast_result.xml", "r") as f:
            blast_records = NCBIXML.parse(f)

            for record in blast_records:
                if record.alignments:
                    top = record.alignments[0]

                    print("\nTop Match Found:")
                    print("\nTitle:", top.title)
                    print("\nAccession:", top.accession)

                    return top.accession

    except:
        print("File empty or not found")

    return None


# FIX 5: fetch_reference() was indented inside get_top_hit() — moved to top level
def fetch_reference(accession):
    Entrez.email = "mallap894@lynchburg.edu"

    print("\nFetching reference DNA from database...")

    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text"
    )

    record = SeqIO.read(handle, "fasta")

    reference_seq = str(record.seq).upper()  # FIX 6: was stored as 'seq_name' but returned 'reference_seq'

    print("\nReference DNA retrieved")
    print("\nReference length:", len(reference_seq))

    return reference_seq


def find_mutations_aligned(reference, patient):

    print("\nPerforming sequence alignment...")

    aligner = PairwiseAligner()
    aligner.mode = 'local'

    alignments = aligner.align(reference, patient)

    best = alignments[0]

    aligned_ref_seq = str(best.target)
    aligned_pat_seq = str(best.query)

    print("Alignment complete")

    mutations = []

    min_len = min(len(aligned_ref_seq), len(aligned_pat_seq))

    for i in range(min_len):
        r = aligned_ref_seq[i]
        p = aligned_pat_seq[i]

        if r == '-' or p == '-':
            mutations.append((i, r, p))
        elif r != p:
            mutations.append((i, r, p))

    return mutations, aligned_ref_seq, aligned_pat_seq


def main():

    # FIX 8: main() was calling wrong function names — now matches actual function names
    file_path = inputFile()
    dna = readTheFile(file_path)

    if not dna:
        return

    print("\nDNA Loaded")
    print("Length:", len(dna))
    print("Preview:", dna[:50] + "...")

    runBlast(dna)

    accession = get_top_hit()

    if not accession:
        return

    reference_dna = fetch_reference(accession)

    mutations, aligned_ref, aligned_pat = find_mutations_aligned(reference_dna, dna)

    print("\nMutation Analysis")
    print("Total mutations found:", len(mutations))

    if mutations:
        print("\nFirst 10 mutations:")
        for m in mutations[:10]:
            ref_base = m[1] if m[1] != '-' else '(gap)'
            pat_base = m[2] if m[2] != '-' else '(gap)'
            print(f"  Position {m[0]}: {ref_base} -> {pat_base}")
    else:
        print("No mutations detected — sequences match perfectly.")


if __name__ == "__main__":
    main()  # FIX 9: stray period removed from main().