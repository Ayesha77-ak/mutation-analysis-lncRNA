# lncRNA_mutation_analysis_1000G.R
# Author: Ayesha Ikram
# Description: Analyzes overlaps between lncRNA regions and variant calls (VCF) from genomic data

# Load libraries
library(GenomicRanges)       # For handling genomic coordinates
library(rtracklayer)         # To read GTF annotation files
library(VariantAnnotation)   # To read VCF files
library(ggplot2)             # For plotting

# ---- Step 1: Define file paths ----
vcf_file <- "example_data/NA12878_chr22.vcf"   # Example variant file (simulate 1000G)
gtf_file <- "example_data/gencode.v39.annotation.gtf"   # GTF file with lncRNA annotations

# ---- Step 2: Load VCF file ----
vcf <- readVcf(vcf_file, "hg38")  # Replace with your genome build if different
gr_variants <- rowRanges(vcf)    # Convert variants to GRanges object

# ---- Step 3: Load GTF file and extract lncRNAs ----
gtf <- import(gtf_file)
lncRNAs <- gtf[gtf$type == "transcript" & gtf$transcript_biotype == "lncRNA"]

# ---- Step 4: Identify overlapping mutations in lncRNAs ----
overlaps <- findOverlaps(gr_variants, lncRNAs)
variant_hits <- gr_variants[queryHits(overlaps)]
lncRNA_hits <- lncRNAs[subjectHits(overlaps)]

# ---- Step 5: Create results dataframe ----
results_df <- data.frame(
  Chromosome = as.character(seqnames(variant_hits)),
  Position = start(variant_hits),
  lncRNA_ID = mcols(lncRNA_hits)$transcript_id,
  Gene = mcols(lncRNA_hits)$gene_name
)

# ---- Step 6: Save results ----
dir.create("results", showWarnings = FALSE)
write.csv(results_df, "results/lncRNA_mutation_hits.csv", row.names = FALSE)

# ---- Step 7: Plot mutation counts by chromosome ----
dir.create("plots", showWarnings = FALSE)
results_df$Chromosome <- factor(results_df$Chromosome)
p <- ggplot(results_df, aes(x = Chromosome)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Mutation Counts in lncRNA Regions by Chromosome",
       x = "Chromosome", y = "Mutation Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/lncRNA_mutation_distribution.png", plot = p)

# ---- Step 8: Console output ----
cat("✅ Total overlapping mutations in lncRNAs:", nrow(results_df), "\n")
cat("📄 Results saved to: results/lncRNA_mutation_hits.csv\n")
cat("📊 Plot saved to: plots/lncRNA_mutation_distribution.png\n")
