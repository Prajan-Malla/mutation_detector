python
# --------------------------------------------
# IMPORTS
# --------------------------------------------

# BioPython tool to convert DNA -> protein
from Bio.Seq import Seq  

# Import YOUR existing mutation detector functions
# This keeps your project modular (good design!)
from simple_mutation_detector import (
    readTheFile,
    runBlast,
    get_top_hit,
    fetch_reference,
    find_mutations_aligned
)


# --------------------------------------------
# DNA → PROTEIN TRANSLATION
# --------------------------------------------
def translate_dna(seq):
    """
    Converts a DNA sequence into a protein sequence.

    Why this matters:
    DNA mutations only become biologically meaningful
    when they change the resulting protein.

    Example:
    DNA: ATG → Protein: M (Methionine)

    Returns:
        Protein sequence as a string
    """
    return str(Seq(seq).translate())


# --------------------------------------------
# MUTATION CLASSIFICATION (PROTEIN LEVEL)
# --------------------------------------------
def classify_protein_changes(ref_dna, pat_dna):
    """
    Compares reference DNA and patient DNA at the protein level.

    Steps:
    1. Convert both DNA sequences into proteins
    2. Compare amino acids one by one
    3. Detect type of mutation

    Types:
    - Silent: No protein change (usually harmless)
    - Missense: Amino acid changes (may affect function)
    - Nonsense: Creates STOP codon (very serious)
    - STOP lost: Protein becomes longer (can be harmful)
    """

    # Convert both sequences to proteins
    ref_protein = translate_dna(ref_dna)
    pat_protein = translate_dna(pat_dna)

    effects = []

    # Only compare up to shortest length to avoid index errors
    min_len = min(len(ref_protein), len(pat_protein))

    for i in range(min_len):
        r = ref_protein[i]   # reference amino acid
        p = pat_protein[i]   # patient amino acid

        # If amino acids differ → mutation affects protein
        if r != p:

            # STOP codon gained → protein stops early
            if p == "*":
                effects.append((i, "Nonsense mutation (STOP gained)"))

            # STOP codon lost → protein keeps going
            elif r == "*":
                effects.append((i, "STOP lost"))

            # Otherwise → amino acid changed
            else:
                effects.append((i, "Missense mutation"))

    # If proteins are identical → mutation is silent
    if ref_protein == pat_protein:
        effects.append(("All", "Silent mutation"))

    return effects


# --------------------------------------------
# FRAMESHIFT DETECTION
# --------------------------------------------
def detect_frameshift(dna_seq):
    """
    Detects frameshift mutations.

    Key idea:
    DNA is read in groups of 3 (codons).
    If length is not divisible by 3 → reading frame shifts.

    Example:
    ATG-AAA-GGG (valid)
    ATG-AAG-GG... (shifted → everything changes)

    Frameshift = VERY HIGH IMPACT
    """
    return len(dna_seq) % 3 != 0


# --------------------------------------------
# SIGNIFICANCE SCORING SYSTEM
# --------------------------------------------
def mutation_significance_score(ref_dna, pat_dna):
    """
    Assigns a numerical score based on how dangerous mutations are.

    Logic:
    - Frameshift → +5 (very severe)
    - Nonsense → +5 (protein stops early)
    - Missense → +2 (moderate effect)

    Higher score = more biologically significant mutation
    """

    score = 0

    # Check frameshift first (biggest impact)
    if detect_frameshift(pat_dna):
        score += 5

    # Analyze protein-level effects
    effects = classify_protein_changes(ref_dna, pat_dna)

    for e in effects:
        if "Nonsense" in e[1]:
            score += 5
        elif "Missense" in e[1]:
            score += 2

    return score, effects


# --------------------------------------------
# MAIN PIPELINE (CONNECTS EVERYTHING)
# --------------------------------------------
def main():
    """
    This is the full workflow:

    1. Load patient DNA
    2. Run BLAST → find closest known sequence
    3. Fetch reference DNA
    4. Align sequences → find mutations
    5. Analyze biological significance
    """

    # Ask user for DNA file
    file_path = input("\nEnter DNA file path: ")

    # Read DNA sequence from file
    dna = readTheFile(file_path)

    if not dna:
        return

    print("\nRunning base mutation detector...")

    # Step 1: Run BLAST (find similar known sequence)
    runBlast(dna)

    # Step 2: Get best match accession ID
    accession = get_top_hit()
    if not accession:
        return

    # Step 3: Fetch reference DNA from database
    reference_dna = fetch_reference(accession)

    # Step 4: Align sequences + detect mutations
    mutations, aligned_ref, aligned_pat = find_mutations_aligned(reference_dna, dna)

    print("\n--- BASIC MUTATION RESULTS ---")
    print("Total mutations:", len(mutations))


    # --------------------------------------------
    # SIGNIFICANCE ANALYSIS
    # --------------------------------------------
    print("\n--- SIGNIFICANCE ANALYSIS ---")

    # Calculate score + mutation effects
    score, effects = mutation_significance_score(reference_dna, dna)

    # Check frameshift
    if detect_frameshift(dna):
        print("⚠️ Frameshift mutation detected (HIGH IMPACT)")

    print("\nProtein-level effects:")

    # Show first few mutations to avoid clutter
    for e in effects[:10]:
        print(f"  Position {e[0]}: {e[1]}")

    print("\nFinal Score:", score)

    # Interpret score for user
    if score >= 5:
        print("🔥 High-impact mutation")
    elif score >= 2:
        print("⚠️ Moderate-impact mutation")
    else:
        print("✅ Likely benign mutation")


# --------------------------------------------
# PROGRAM ENTRY POINT
# --------------------------------------------
# This ensures the script runs only when executed directly
if __name__ == "__main__":
    main()
