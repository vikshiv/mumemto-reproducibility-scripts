# Reproducibility scripts for **mumemto**

Mumemto is a tool for analyzing pangenome sequence collections. It identifies **maximal unique/exact matches (multi-MUMs and multi-MEMs)** present across a collection of sequences. Mumemto can **visualize** pangenome synteny, **identify misassemblies**, and provide a unifiying structure to a pangenome.

Preprint available at [https://doi.org/10.1101/2025.01.05.631388](https://doi.org/10.1101/2025.01.05.631388).

This repo contains scripts to recreate the results reported in the preprint. The main software repo is available at [https://github.com/vikshiv/mumemto](https://github.com/vikshiv/mumemto).

## How to use this repo
Each figure has a collection of scripts to reproduce the corresponding figure in the paper. The assembly data for each experiment can be found here: [onedrive link]. For each bash script, change the path at the top to match the data path.

See references below for the original assembly sources. HPRC assemblies were scaffolded with [RagTag](https://github.com/malonge/RagTag).

### Figure 2
First run the two bash scripts in `graph_construction` to build the Minigraph-Cactus comparison graph and the Mumemto-based graphs. Then to generate Fig 2H, use the jupyter notebook in `fig2`.

### Figure 3
First generate the multi-MUM information using mumemto with the bash script in `fig3`. Then use the two jupyter notebooks to recreate Fig3.

### Figure 4
First generate the multi-MUM information using mumemto with the bash script in `fig4`. Then use the included jupyter notebooks to recreate Fig4.

## References
- W.-W. Liao, M. Asri, J. Ebler, D. Doerr, M. Haukness, G. Hickey, S. Lu, J. K. Lucas, J. Monlong, H. J. Abel, et al., “A draft human pangenome reference,” Nature, vol. 617, no. 7960, pp. 312–324, 2023.
- D. Tang, Y. Jia, J. Zhang, H. Li, L. Cheng, P. Wang, Z. Bao, Z. Liu, S. Feng, X. Zhu, et al., “Genome evolution and diversity of wild and cultivated potatoes,” Nature, vol. 606, no. 7914, pp. 535–541, 2022.
- Potato Database, http://solomics.agis.org.cn/potato/, Accessed: 2024-12-10.
- G. Hickey, J. Monlong, J. Ebler, A. M. Novak, J. M. Eizenga, Y. Gao, T. Marschall, H. Li, and B. Paten, “Pangenome graph construction from genome alignments with Minigraph-Cactus,” Nature Biotechnology, pp. 1–11, 2023.
