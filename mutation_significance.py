from Bio.Seq import Seq
from simple_mutation_detector import (
    readTheFile,
    runBlast,
    get_top_hit,
    fetch_reference,
    find_mutations_aligned
)

# ----------------------------
# TRANSLATION
# ----------------------------
def translate_dna(seq):
    return str(Seq(seq).translate())

# ----------------------------
# MUTATION CLASSIFICATION
# ----------------------------
def classify_protein_changes(ref_dna, pat_dna):
    ref_protein = translate_dna(ref_dna)
    pat_protein = translate_dna(pat_dna)

    effects = []

    min_len = min(len(ref_protein), len(pat_protein))

    for i in range(min_len):
        r = ref_protein[i]
        p = pat_protein[i]

        if r != p:
            if p == "*":
                effects.append((i, "Nonsense mutation (STOP gained)"))
            elif r == "*":
                effects.append((i, "STOP lost"))
            else:
                effects.append((i, "Missense mutation"))

    if ref_protein == pat_protein:
        effects.append(("All", "Silent mutation"))

    return effects

# ----------------------------
# FRAMESHIFT CHECK
# ----------------------------
def detect_frameshift(dna_seq):
    return len(dna_seq) % 3 != 0

# ----------------------------
# SIGNIFICANCE SCORING
# ----------------------------
def mutation_significance_score(ref_dna, pat_dna):
    score = 0

    # Frameshift = high impact
    if detect_frameshift(pat_dna):
        score += 5

    effects = classify_protein_changes(ref_dna, pat_dna)

    for e in effects:
        if "Nonsense" in e[1]:
            score += 5
        elif "Missense" in e[1]:
            score += 2

    return score, effects

# ----------------------------
# MAIN PIPELINE (WRAPPER)
# ----------------------------
def main():
    file_path = input("\nEnter DNA file path: ")

    dna = readTheFile(file_path)
    if not dna:
        return

    print("\nRunning base mutation detector...")
    runBlast(dna)

    accession = get_top_hit()
    if not accession:
        return

    reference_dna = fetch_reference(accession)

    mutations, aligned_ref, aligned_pat = find_mutations_aligned(reference_dna, dna)

    print("\n--- BASIC MUTATION RESULTS ---")
    print("Total mutations:", len(mutations))

    # ----------------------------
    # SIGNIFICANCE ANALYSIS
    # ----------------------------
    print("\n--- SIGNIFICANCE ANALYSIS ---")

    score, effects = mutation_significance_score(reference_dna, dna)

    if detect_frameshift(dna):
        print("⚠️ Frameshift mutation detected (HIGH IMPACT)")

    print("\nProtein-level effects:")
    for e in effects[:10]:
        print(f"  Position {e[0]}: {e[1]}")

    print("\nFinal Score:", score)

    if score >= 5:
        print("🔥 High-impact mutation")
    elif score >= 2:
        print("⚠️ Moderate-impact mutation")
    else:
        print("✅ Likely benign mutation")


if __name__ == "__main__":
    main()