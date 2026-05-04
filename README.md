🧬 Simple Mutation Detector – Bioinformatics Analysis Tool
📌 Project Overview

The Simple Mutation Detector is a Python-based bioinformatics pipeline designed to identify and analyze genetic mutations by comparing a patient DNA sequence against a reference genome. The system integrates NCBI BLAST, sequence alignment algorithms, and protein-level translation to assess the biological significance of detected mutations.

This project demonstrates applied skills in computational biology, sequence analysis, and data-driven biological interpretation, bridging raw DNA data to meaningful mutation insights.

⚙️ Key Features
📂 Reads DNA sequences from FASTA/text files
🌐 Integrates NCBI BLAST (blastn) for sequence matching
🧬 Retrieves reference genomes using NCBI Entrez API
🔬 Performs pairwise sequence alignment to detect variations
⚡ Identifies mutations including:
Substitutions
Insertions
Deletions
🧪 Converts DNA sequences into protein sequences using BioPython
🧠 Classifies mutation effects:
Silent mutations
Missense mutations
Nonsense mutations (STOP gained)
STOP loss mutations
⚠️ Detects frameshift mutations
📊 Computes a mutation severity score to estimate biological impact
🧠 Technical Highlights

This project applies key concepts in computational biology:

Sequence alignment using dynamic comparison of biological sequences
Codon-based translation from DNA to protein
Mutation impact classification based on protein-level changes
Integration of external biological databases (NCBI)
Handling of real-world biological data formats (FASTA/GenBank)
🛠️ Tools & Technologies
Python 3
BioPython
NCBI BLAST API
NCBI Entrez Database
Pairwise Sequence Alignment
🔬 Workflow Pipeline
Input DNA sequence from file
Run BLAST search against NCBI database
Retrieve closest reference genome
Align patient vs reference sequence
Detect base-level mutations
Translate DNA into protein sequence
Classify mutation type and biological impact
Generate mutation severity score
📊 Outcome

The system outputs:

Total number of mutations
Detailed mutation positions
Protein-level impact analysis
Severity classification (benign → high-impact)

This allows rapid assessment of whether a mutation is likely neutral or biologically significant.

🎯 Skills Demonstrated
Bioinformatics pipeline design
API integration (NCBI Entrez & BLAST)
Sequence alignment and analysis
Object-based biological data processing
Python software engineering
Computational genetics interpretation
🚀 Future Improvements
Web-based interface for easier usage
Graphical mutation visualization
Support for multiple gene comparison
Machine learning model for mutation impact prediction
Exportable clinical-style reports (PDF/CSV)
