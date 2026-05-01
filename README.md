
 Simple DNA Mutation Detector (BLAST + Alignment)

This project is a Python-based bioinformatics tool that detects DNA mutations by:

1. Reading a patient DNA sequence from a file  
2. Running NCBI BLAST to find the closest reference sequence  
3. Fetching the reference DNA sequence from the NCBI database  
4. Aligning patient and reference sequences  
5. Detecting and reporting mutations (substitutions and gaps)

It uses **Biopython** and real-time access to the NCBI nucleotide database.

---

## ⚙️ Features

-  Reads FASTA or raw DNA sequence files  
-  Runs online BLAST search (NCBI BLASTn)  
-  Fetches reference sequence using accession ID  
-  Performs pairwise sequence alignment  
-  Detects:
- Base substitutions (A → G, C → T, etc.)
- Insertions/deletions (gaps)
- 📊 Outputs mutation positions and differences  

---

Example Output

DNA Loaded
Length: 1200
Preview: ATGCGTACGTTAGC...

Running BLAST...
BLAST completed and saved to file

Top Match Found:
Title: Homo sapiens hemoglobin beta gene
Accession: NM_000518

Reference DNA retrieved
Reference length: 1600

Performing sequence alignment...
Alignment complete

Mutation Analysis
Total mutations found: 12

First 10 mutations:
Position 266: G -> A
Position 293: C -> T
Position 310: T -> (gap)


---



##  Requirements

Install dependencies:

```bash
pip install biopython
````

## How to Run

1. Save your DNA file (FASTA or raw sequence)

2. Run the script:

```bash
python mutation_detector.py
```

3. Enter file path when prompted:

```
Enter the filepath of the DNA file:
```

---

## 🔬 Workflow Diagram

```
DNA File
   ↓
Read Sequence
   ↓
NCBI BLAST (find best match)
   ↓
Fetch Reference Genome
   ↓
Pairwise Alignment
   ↓
Mutation Detection
   ↓
Results Output
```

---

## Key Libraries Used

* Bio.Blast (NCBIWWW, NCBIXML) → BLAST search
* Bio.Entrez → Fetch sequences from NCBI
* Bio.SeqIO → Parse FASTA files
* Bio.Align.PairwiseAligner → Sequence alignment

---

## ⚠️ Limitations

* if the e score is less, the mutation is inaccurate.

---

## Author : Prajan Malla

