# Simulation studies for studying the orientation of SNPs in MR-Egger
This repository contains simulation and analysis codes for the paper 
"Potential problems with Egger regression in Mendelian randomization".

### Files description
* `main_simulation.R`: Summary-level simulations described in the main text, including Simulation (a), (b) and (c) under both balanced and directional pleiotropy.
* `indiv_TSLS_PMR.R`: Individual-level simulations described in Supplementary Section S5 and S6. 
* `weakIV_simulation.R`: Summary-level simulations with irrelevant IVs described in Supplementary Section S7.
* `intercept_simulation.R`: Simulations to study the intercept test in MR-Egger described in Supplementary Section S8
* `mr_egger_func.R`: Functions to perform MR-Egger with a specific coding scheme, and with many random coding schemes.
* `/real_data_res/48pair_egger_intercept.RData`: Real data results of MR-Egger intercept using random coding schemes.
* `/real_data_res/48pair_default_egger_intercept.RData`: Real data results of MR-Egger intercept using the default coding.
