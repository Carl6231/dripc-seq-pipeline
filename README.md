# DRIPc-seq representative workflow

This repository provides a standardized representative workflow for documenting the main analytical steps used in DRIPc-seq processing and differential R-loop peak analysis. It is intended to improve transparency and reproducibility, but it should not be interpreted as a complete parameter-level reconstruction of every project-specific analysis detail underlying the manuscript-reported primary results.

For differential R-loop peak analysis, the DiffBind workflow exports all reported peaks and filters significant differential R-loop-enriched regions using an FDR threshold by default, together with the stated fold-change threshold. Legacy command-line aliases may be retained only for backward compatibility with older wrapper scripts; the implemented significance criterion in the revised workflow is FDR-based.
