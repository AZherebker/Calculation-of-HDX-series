# Calculation of HDX Series

The code and the data for the paper `Decreasing the uncertainty in the comparison of molecular fingerprints of organic aerosols with H/D exchange mass spectrometry`.

### Folder and File Descriptions

- **`data/`**: This folder contains the raw data files used for the analysis.
  - **`HDX_1.csv`**: Mass spectrometry HDX data for sample 1.
  - **`HDX_2.csv`**: Mass spectrometry HDX data for sample 2.
  - **`Test-samples.csv`**: A test sample file with formulae lists used for the HDX calculations.

- **`notebooks/`**: Contains the Jupyter notebook used for processing the data and performing the HDX series calculations.
  - **`HDX_project_GIT.ipynb`**: The Jupyter notebook where the HDX calculations are performed, including data processing, analysis, and results generation.

- **`results/`**: Stores the output and results generated from the notebook.
  - **`Results-HDX.csv`**: The processed output after calculating HDX series from the input data.

- **`README.md`**: This file, which explains the structure and purpose of each folder and file.

## How to Run the Code

1. **Set Up the Environment**:
   - Ensure that you have all necessary dependencies installed. You can use PDM or `pip` to install them.
   - Required libraries: `numpy`, `pandas`, `matplotlib`, `tqdm`.

2. **Load the Notebook**:
   - Open the `HDX_project_GIT.ipynb` notebook located in the `notebooks/` folder.
   - Follow the instructions in the notebook to process the data.

3. **Input Data**:
   - The input data files are located in the `data/` folder. These include the HDX samples and test sample data.

4. **Output**:
   - After running the notebook, the results will be generated and saved in the `results` folder.
