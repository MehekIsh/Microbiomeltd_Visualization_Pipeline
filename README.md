## Plot_Microbiome
plot_microbiome is an R helper function that automates a full alpha-diversity workflow for microbiome studies. It takes a feature (ASV/OTU) table, taxonomy table, sample metadata, and an optional phylogenetic tree, then generates publication-ready alpha-diversity figures (TIFF) and an Excel workbook of all calculated indices. The function handles grouping by any metadata column, runs pairwise Wilcoxon tests with Benjamini–Hochberg correction, annotates only significant comparisons, and trims extreme outliers with a 3×IQR rule. Outputs include boxplots for observed richness, Shannon, Simpson, and Chao1, as well as wide/long alpha tables and tidy summaries of significant group differences, making it easy to drop the results directly into downstream analysis or manuscripts.
## Pipeline Documentation and Parameters
plot_microbiome <- function(
feature_table,
taxonomy_file,
metadata_file,
tree_file = NULL,
output_dir = "microbiome_outputs",
plot_types = c("abundance","rarefaction","alpha","beta"),
tax_ranks = c("Phylum","Class","Order","Family","Genus","Species"),
dpi = 300,
width = 10,
height = 6,
group_by = NULL,
alpha_show_p = TRUE
)
### Parameters (what each file/option does)
● feature_table: Tab-delimited counts per ASV/OTU by sample. Sample IDs must match the metadata.
● taxonomy_file: Taxonomic annotations for features. Required for abundance plots; not strictly needed for alpha, but safe to include.
● metadata_file: Tab-delimited sample metadata. Must contain the column named in group_by (e.g., “Group”).
● tree_file: Newick tree. Used for phylogenetic metrics if enabled; optional for the current alpha run.
● output_dir: Folder where figures, logs, and Excel workbooks are written.
● plot_types: Which analyses to produce. You’re using "alpha" here.
● group_by: Metadata column that defines groups for plotting and statistics. 

### Expected outputs (for plot_types = "alpha")
● Figures (TIFF)
● alpha_shannon_observed_p001.tiff, alpha_simpson_chao1_p001.tiff,
alpha_chao1_overall_box_p001.tiff
● The same three files with _p05.tiff suffix
● Boxplots use translucent fills, only statistically significant pairwise differences are
annotated (Wilcoxon, BH-adjusted), brackets are evenly spaced, and extreme outliers
are removed using a 3×IQR rule. A right-hand note on each figure states these criteria.
● Excel workbook
● alpha_diversity_combined.xlsx containing:
○ Alpha_All_Indices: one row per sample with the following columns, in order:
○ observed, chao1, diversity_inverse_simpson, diversity_gini_simpson,
diversity_shannon, diversity_fisher, diversity_coverage, evenness_camargo,
evenness_pielou, evenness_simpson, evenness_evar, evenness_bulla,
dominance_dbp, dominance_dmn, dominance_absolute, dominance_relative,

dominance_simpson, dominance_core_abundance, dominance_gini,
rarity_log_modulo_skewness, rarity_low_abundance, rarity_rare_abundance.
○ When the microbiome R package is available these are computed directly;
otherwise, everything feasible is filled, and the rest are marked NA.
○ Alpha_Long and Alpha_Wide data tables.
○ Observed_Wilcoxon_BH, Shannon_Wilcoxon_BH, Simpson_Wilcoxon_BH,
Chao1_Wilcoxon_BH: BH-adjusted pairwise Wilcoxon p-value matrices,
computed on the trimmed data to match the plots.
○ Significant_Pairs_p001 and Significant_Pairs_p05: tidy tables of the significant
comparisons retained in the figures for each alpha threshold.

### If something needs a fix
If anything errors, looks cramped, or a sheet is missing, please contact me:
●  The exact console message,
● A screenshot or listing of the files written to output_dir,
● The first ~10 lines of your metadata1.txt (including headers).

### Quick checks that solve most issues:
● The group_by column exists in the metadata and is spelled exactly the same.
● Sample IDs in the metadata match those in the feature table.
● Files are tab-delimited without stray quotes.
● If the microbiome package is not installed, some indices in Alpha_All_Indices will be NA,
which is expected.
