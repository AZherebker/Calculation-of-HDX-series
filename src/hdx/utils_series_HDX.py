# Imports
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_overlaps(err_in_ppm: float, parent: pd.DataFrame, hdx: pd.DataFrame) -> pd.DataFrame:
    """
    Optimized version of checking overlaps in m/z values between two datasets within a given error margin (in ppm).

    Parameters:
        err_in_ppm (float): Error margin in parts per million (ppm).
        parent (pd.DataFrame): DataFrame containing 'mz' column for parent ion data.
        hdx (pd.DataFrame): DataFrame containing 'mz' column for HDX ion data.

    Returns:
        pd.DataFrame: A copy of the 'hdx' DataFrame with overlapping peaks removed.
    """
    # Calculate lower and upper bounds of m/z ranges for each parent
    parent['mz_lower'] = parent['mz'] - (err_in_ppm / 10**6) * parent['mz']
    parent['mz_upper'] = parent['mz'] + (err_in_ppm / 10**6) * parent['mz']
    
    # Perform a cross-join between parent and hdx datasets
    # Use pd.merge to create a cartesian product of both datasets
    merged_df = pd.merge(parent[['mz', 'mz_lower', 'mz_upper']], hdx[['mz']], how='cross')

    # Filter rows where hdx.mz falls within the parent's m/z bounds
    overlapping_peaks = merged_df[(merged_df['mz_lower'] < merged_df['mz_y']) & (merged_df['mz_y'] < merged_df['mz_upper'])]['mz_y'].unique()

    logging.info(f'Within {err_in_ppm} ppm found {len(overlapping_peaks)} peaks overlapping between parent ({len(parent["mz"])} peaks) and hdx ({len(hdx["mz"])} peaks)')

    # Filter out overlapping peaks from hdx dataset
    hdx_without_subtr = hdx[~hdx["mz"].isin(overlapping_peaks)]

    return hdx_without_subtr

# Function to find isotopic series in the data
def find_series(
    data: pd.DataFrame, 
    hdx_without_subtr: pd.DataFrame, 
    Err: float, 
    mz_accrcy: float, 
    d: float
) -> tuple[pd.DataFrame, list]:
    """
    Finds isotopic series in the 'data' based on comparison with the 'hdx_without_subtr' dataset.
    The function computes possible isotopic shifts and matches them with allowed errors.

    Parameters:
        data (pd.DataFrame): The parent ion data, must contain 'mz', 'o', 'h', 'c', 'n', and 's' columns.
        hdx_without_subtr (pd.DataFrame): HDX data with overlapping peaks removed.
        Err (float): Error tolerance for mass differences.
        mz_accrcy (float): Allowed mass accuracy for matching m/z values.
        d (float): Step size for calculating isotopic shifts.

    Returns:
        tuple: 
            - pd.DataFrame: Updated parent dataset with HDX information added.
            - list: A list containing DataFrames for each isotopic series found.
    """
    ErMax = Err
    ErMin = -Err

    hdx_list = []
    search = []
    hdx_fin = []

    series_list = []

    # Iterate through each m/z value in the parent dataset
    for m in tqdm(range(len(data["mz"])), desc="Finding series"):
        # Define m/z search range
        x1 = data["mz"][m] - mz_accrcy
        x2 = data["mz"][m] + d * data["o"][m] + mz_accrcy

        # Filter hdx data based on m/z search range
        s = hdx_without_subtr[(x1 < hdx_without_subtr['mz']) & (hdx_without_subtr['mz'] < x2)].copy()

        # Calculate isotopic shift (n) for each m/z value within the region
        s["n"] = round((s["mz"] - data["mz"][m]) / d)
        
        # Calculate errors based on the expected isotopic shifts
        s["err"] = s["mz"] - data["mz"][m] - s["n"] * d
        
        # Filter based on error tolerance
        s_f = s[(ErMin < s["err"]) & (s["err"] < ErMax)]

        if not s_f.empty:
            # Check if the number of shifts is within allowed range for O and S atoms, and H
            if data["o"][m] + data["s"][m] >= max(s_f["n"]) + 1 and max(s_f["n"]) + 1 <= data["h"][m]:
                hdx_list.append(int(max(s_f["n"])))
            elif data["o"][m] + data["s"][m] < max(s_f["n"]) + 1 or max(s_f["n"]) + 1 > data["h"][m]:
                search = s_f["n"].tolist()
                for i in range(len(search)):
                    if data["o"][m] + data["s"][m] >= search[-(i + 1)] + 1 and search[-(i + 1)] + 1 <= data["h"][m]:
                        hdx_list.append(int(search[-(i + 1)]))
                        break
                if data["o"][m] + data["s"][m] < search[-(i + 1)] + 1 or search[-(i + 1)] > data["h"][m]:
                    hdx_list.append('NA')
                    hdx_fin.append('NA')
            
            # Further filtering and series determination
            if hdx_list[-1] != 'NA':
                if s_f["n"].tolist()[0] >= 2:
                    hdx_list.append('NA')
                    hdx_fin.append('NA')
                if s_f["n"].tolist()[0] < 2:
                    s_fi = s_f[(s_f["n"] <= hdx_list[-1])].copy()
                    filter_n = s_fi["n"].tolist()
                    answer = []
                    answer.append(filter_n[0])
                    for i in range(len(filter_n) - 1):
                        if filter_n[i + 1] - filter_n[i] <= 1:
                            answer.append(filter_n[i + 1])
                        else:
                            break
                    hdx_fin.append(int(answer[-1]))
                    s_fin = s_f[(s_f["n"] <= hdx_fin[-1])].copy()
        else:
            hdx_list.append('NA')
            hdx_fin.append('NA')

        # Adding metadata to the series list
        if hdx_fin[-1] != 'NA':
            s_fin["Mp"] = data["mz"][m]
            s_fin["c"] = data["c"][m]
            s_fin["h"] = data["h"][m]
            s_fin["o"] = data["o"][m]
            s_fin["nit"] = data["n"][m]
            s_fin["s"] = data["s"][m]
            series_list.append(s_fin)

    # Calculate H-labile data based on final isotopic shifts
    Hl = []
    for item in hdx_fin: 
        try: 
            Hl.append(item + 1) 
        except: 
            Hl.append('NA')

    logging.info(len(hdx_fin))

    # Using assign to add columns in a cleaner way
    data_w_hdx = data.assign(HDX=hdx_fin, H_Labile=Hl)

    # Return the modified DataFrame
    return data_w_hdx, series_list

