# PPbenchmark
Notebooks for Palaeoproteomics Benchmarking paper:

Rodriguez Palomo, I.; Nair, B.; Chiang, Y.; Dekker, J.; Dartigues, B.; Mackie, M.; Evans, M.; Macleod, R.; Olsen, J. V.; Collins, M. J. Benchmarking the identification of a single degraded protein to explore optimal search strategies for ancient proteins. Peer Community Journal, Volume 4 (2024), article no. e107. https://doi.org/10.24072/pcjournal.491

Notebooks:
1. `process_results.ipynb` This notebook collects and aggregates PSMs data from the tested proteomics tools. It also calculates q_values, adds retention times when not provided, peptide position and normalises protein IDs
2. `map_denovo.ipynb` Maps *de novo* peptides to the BLG sequence and calculates position accuracy
3. `make_plots.ipynb` For the database search tools, takes the data from `process_results.ipynb` and produces figures related to these tools. It also takes the physicochemical properties files and produces the position specific plot.
4. `plot_positions.Rmd` (and associated `.md`) Plots the BLG PSM position coverage and accuracy (for *de novo*)
5. `amyloid.pred.R` small R script to produce the amyloid propensity for `AmyloGram` and `appnn`

- For Python, the notebooks use `pyteomics` and `lxml` to read pepxml and mzML files. `seaborn`, `upsetplot` and `matplotlib` are used` for plotting.
- In the the Rmd file, plots are produced with `ggplot`
