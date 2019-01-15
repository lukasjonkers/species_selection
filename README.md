# species_selection
Data and R code to assess planktonic foraminifera species importance and transfer function performance
described in "Sensitivity to species selection indicates the effect of nuisance variables on marine microfossil transfer functions"
Lukas Jonkers and Michal Kucera, Climate of the Past Discussions: https://doi.org/10.5194/cp-2018-107

Contents:

DATA:

dat.RDS: ForCenS planktonic foraminifera assemblages (Siccha and Kucera, 2017, Sci. Data, https://doi.org/10.1038/sdata.2017.109) and SST data (World Ocean Atlas 2001, Stephens et al. 2002).

fossil_data.RDS: fossil planktonic foraminifera assemblages.
Core MD95-2040: https://doi.pangaea.de/10.1594/PANGAEA.66811;
core M35003-4: https://doi.pangaea.de/10.1594/PANGAEA.55756;
MARGO LGM: https://doi.pangaea.de/10.1594/PANGAEA.227329.

CODE:

taxon_ranking.R: determine species ranking (calls random_Mat.R)

many_TFs.R: build transfer functions (calls many_MAT.R and many_WA.R)

assess_ranking.R: explore species ranking

reconstruct_SST.R: predict SST based on fossil assemblages. Use increasing number of transfer functions based on species ranking

Figs.R: figures 1 to 8
