INTNet
=====
Multi-task multilabel deep neural networks for identification and classification of integrons.

INTNet is designed to identify and classify integron integrases, predict bacterial hosts, and associate ARGs (multi-label) with integrons. This versatile tool supports a range of input types including:
* Long Amino Acid Sequences (Full Length/Contigs)
* Long Nucleotide Sequences
* Short Amino Acid Reads (30-50 aa)
* Short Nucleotide Reads (100-150 nt)

All inputs should be in FASTA format.

**INTNet Components**\
INTNet comprises two specialized models to accommodate different read lengths:
* **INTNet-s**: Optimized for short reads, enhancing prediction accuracy for sequences ranging from 30 to 50 amino acids or 100 to 150 nucleotides.
* **INTNet-l**: Tailored for long sequences, ensuring robust predictions for full-length contigs or long nucleotide sequences.




![alt text](https://github.com/patience111/INTNet/blob/master/pics/INTNet_workflow.jpg)

Installation

