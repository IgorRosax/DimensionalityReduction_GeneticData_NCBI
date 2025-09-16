
# HSLocalMDS: Experiments and Results for Dimensionality Reduction in Genetic Data

This repository contains the scripts, results, and visualizations supporting the article:

**"HSLocalMDS: An R library for dimensionality reduction based on Hyperbolic Smoothing Local Multidimensional Scaling"**

**Authors:** Igor da Silva, Fabiano Oliveira, Vinicius Xavier

In this work, we introduce the `HSLocalMDS` R package, which implements the Hyperbolic Smoothing Local Multidimensional Scaling (HSLocalMDS) method for dimensionality reduction. The experiments in this repository compare HSLocalMDS with other dimensionality reduction techniques using real genetic datasets from the NCBI GEO database.

## Repository Structure

- **01.GEOquery.R**: Download and preprocess GEO data.
- **02.GEOquerySumarise.R**: Summarize and overview the obtained data.
- **03.Script_RD.R**: Run dimensionality reduction methods, including HSLocalMDS.
- **04.script_evaluate_DRs.R**: Evaluation of the methods.
- **05.plot_accuracy.R**: Generate plots for accuracy and other metrics.
- **output/**: Analysis results, comparative tables, and plots.
  - PNG and CSV files with metrics such as Agreement Rate, Accuracy, Sensitivity, Specificity, Precision, Recall, and F1.
  - Subfolders for each dataset (e.g., GSE14020, GSE18842, etc.) containing specific results.
- **QVisVis/**: Auxiliary scripts for qualitative visualization and analysis.
- **RD_Results/**: Intermediate results and `.RData` files generated during execution.

## Datasets Used

The experiments use several GEO datasets: `GSE14020`, `GSE28735`, `GSE18842`, `GSE35988`, `GSE21034`, `GSE44076`, and `GSE29272`. Each dataset results is stored in specific subfolders within `output/` and `RD_Results/`.

## How to Reproduce the Results

1. Install the required R packages, including `HSLocalMDS` (see the article or package documentation for installation instructions), `GEOquery`, and other dependencies for dimensionality reduction and visualization.
2. Run the scripts in the suggested order above using R.
3. Results will be automatically saved in the `output/` and `RD_Results/` folders.

## Results and Visualizations

- Comparative tables of dimensionality reduction method performance, highlighting the results of HSLocalMDS versus other methods.
- Plots for accuracy, sensitivity, specificity, precision, recall, F1, and agreement rate.
- Qualitative visualizations of clusters and projections of reduced data.

## Citation and Credits

If you use these scripts or results, please cite the article:

> Igor da Silva, Fabiano Oliveira, Vinicius Xavier. "HSLocalMDS: An R library for dimensionality reduction based on Hyperbolic Smoothing Local Multidimensional Scaling". 2025.

This repository was developed to support academic experiments and the evaluation of dimensionality reduction methods on public genetic datasets.
