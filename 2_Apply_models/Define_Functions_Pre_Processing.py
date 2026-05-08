import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter



# Apply SNV to the spectra
def apply_snv(df):

    wl_cols = [col for col in df.columns if "wl" in col]
    
    snv_array = df[wl_cols].apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=1), axis=1)
    
    snv_array.columns = wl_cols
    
    df_non_wl = df.drop(columns=wl_cols)
    
    df_snv = pd.concat([df_non_wl, snv_array], axis=1)
    
    return df_snv








def apply_sg(df, window_length, polyorder, deriv):
    wl_cols = [c for c in df.columns if "wl" in c]
    X = df[wl_cols].values
    
    X_sg = savgol_filter(
        X,
        window_length=window_length,
        polyorder=polyorder,
        deriv=deriv,
        axis=1
    )
    
    df_out = df.copy()
    df_out[wl_cols] = X_sg
    
    return df_out



def apply_binning(df, bin_size):

    wl_cols = [c for c in df.columns if c.startswith("wl_")]
    wl_vals = np.array([int(c.split("_")[1]) for c in wl_cols])

    order = np.argsort(wl_vals)
    wl_cols = np.array(wl_cols)[order]
    wl_vals = wl_vals[order]

    X = df[wl_cols].values
    n_samples, n_wl = X.shape

    n_bins = n_wl // bin_size  

    X = X[:, :n_bins * bin_size]
    wl_vals = wl_vals[:n_bins * bin_size]

    X_binned = X.reshape(n_samples, n_bins, bin_size).mean(axis=2)

    wl_centers = wl_vals.reshape(n_bins, bin_size)[:, bin_size // 2]

    new_cols = [f"wl_{wl}" for wl in wl_centers]

    df_out = df.drop(columns=wl_cols).copy()
    
    df_out = pd.concat(
    [df_out, pd.DataFrame(X_binned, index=df_out.index, columns=new_cols)],
    axis=1
)

    return df_out



def apply_preprocessing_pipeline(
    df,
    snv,
    sg,
    window_length,
    polyorder,
    deriv,
    bin_size
):

    df_out = df.copy()

    # Binning
    if pd.notna(bin_size):
        df_out = apply_binning(df_out, bin_size)

    # SNV
    if snv:
        df_out = apply_snv(df_out)

    # Savitzky-Golay
    if sg:
        df_out = apply_sg(df_out, window_length, polyorder, deriv)

    return df_out




def split_by_fold(full_df, fold_df, fold_col, test_fold):


    df = full_df.merge(
        fold_df[['ID', fold_col]],
        on='ID',
        how='inner'
    )

    test_df  = df[df[fold_col] == test_fold].drop(columns=[fold_col])
    train_df = df[df[fold_col] != test_fold].drop(columns=[fold_col])

    return train_df, test_df







def s_factor(alphas: np.ndarray, logit_factor: float) -> np.ndarray:
    if not (0.1 <= logit_factor <= 5.0):
        raise ValueError("logit_factor must lie between 0.1 and 5.0.")

    alphas = np.asarray(alphas)


    alpha_trans = (alphas ** logit_factor) / (
        (alphas ** logit_factor) + ((1 - alphas) ** logit_factor)
    )

    return alpha_trans