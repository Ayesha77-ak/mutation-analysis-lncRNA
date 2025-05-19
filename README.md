# 🧬 lncRNA Mutation Analysis Using 1000 Genomes Data

This project performs a bioinformatics analysis to identify and visualize mutations in long non-coding RNA regions using genomic variant data, simulating a workflow inspired by 1000 Genomes Project data. It utilizes variant calls in VCF format and lncRNA annotations from GENCODE to detect overlapping mutation events.

---

## 📌 Objective

To screen and quantify genomic variants located in annotated lncRNA regions, and visualize their distribution across chromosomes. This kind of analysis is relevant for understanding mutation clustering and functional impacts on non-coding regulatory regions—highly relevant to projects in comparative genomics and regulatory RNA biology.

---

## 🧰 Tools & Libraries Used

- `R` (version ≥ 4.0.0)
- `GenomicRanges` – for genomic interval operations
- `rtracklayer` – for importing GTF annotations
- `VariantAnnotation` – for VCF parsing and annotation
- `ggplot2` – for plotting mutation distribution

---



