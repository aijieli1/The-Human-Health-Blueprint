<img src="blueprint.jpg" alt="blueprint" width="800"/>

# Paper Figures Code

This folder contains the code used to generate the figures for the paper:

**Paper Title**  
Author(s): Aijie Li❤

## Link to paper
TBA (will be updated once the paper is public)  

## Contents

Figure numbering follows the main-figure order in the Cardiovascular Diabetology submission manuscript.

| Main figure | Figure title | Available scripts |
| --- | --- | --- |
| Fig. 1 | Genetic architecture of cardiometabolic diseases and metabolic traits | Panel b: [Figure1b_code1.r](Figure1b_code1.r) (Genomic SEM factor-path component) and [Figure1b_code2.r](Figure1b_code2.r) (genetic-correlation heatmap). No script for panel a is included. |
| Fig. 2 | Molecular features of CMD subtypes | Panel a: [Figure2a.r](Figure2a.r) (circular chromosome/SNV/gene plot). Panel b: [Figure2b.r](Figure2b.r) (GO/KEGG pathway and gene plot). |
| Fig. 3 | Tissue and cell type-specific network of cardiometabolic subtypes | No script currently included. |
| Fig. 4 | Multilayer genetic relationships among cardiometabolic subtypes | No script currently included. |
| Fig. 5 | MR-supported directional network involving CAD, HTN and T2D | No script currently included. |
| Fig. 6 | Candidate therapeutic associations with CMD subtypes | [Figure6.r](Figure6.r) (drug-gene-subtype plot). |

`Figure1b_code1.r` generates full Genomic SEM factor-path diagrams; it is retained as a component script for Fig. 1b. The manuscript shows simplified factor loadings in Fig. 1b and full factor-path diagrams in Supplementary Figs. 5b and 6b. The scripts generate figure components rather than all final assembled panels.

Renamed scripts: `Figure2b_code1.r` → `Figure1b_code1.r`; `Figure2b_code2.r` → `Figure1b_code2.r`; `Figure4a.r` → `Figure2a.r`; `Figure4b.r` → `Figure2b.r`. `Figure6.r` retains its name. Script contents are unchanged.

## Requirements
- Python 3.8+
- R version 4.3.2+

## Usage
1. Clone this repository:
   ```bash
   git clone https://github.com/aijieli1/The-Human-Health-Blueprint.git

If you have any questions about this repo, please open an issue or contact **aijie li**.  
- [Email](mailto:your_email@example.com)
