# In this script we prepare the hyperparameter search grid for the different models

from pathlib import Path
import pandas as pd
import numpy as np
import itertools

# Set directory
ROOT = Path(__file__).resolve().parents[1]

#--------------------------------------------------------------
# Define search space for preprocessing, for all models
snv_options = [True, False]
sg_options = [True, False]
window_lengths = [3,7,13]
polyorders = [1, 2]
derivs = [1, 2]
bin_sizes = [5, 10, 15, 20]


grid = []


for snv, sg, bin_size in itertools.product(
    snv_options, sg_options, bin_sizes
):
    if sg:
        for wl, po, d in itertools.product(window_lengths, polyorders, derivs):
            if po >= d:
                grid.append({
                    "SNV": snv,
                    "SG": True,
                    "window_length": wl,
                    "polyorder": po,
                    "deriv": d,
                    "bin_size": bin_size
                })
    else:
        grid.append({
            "SNV": snv,
            "SG": False,
            "window_length": None,
            "polyorder": None,
            "deriv": None,
            "bin_size": bin_size
        })

print(f"Total valid preprocessing pipelines: {len(grid)}")

#--------------------------------------------------------------
# Define Grid for TabPFN

ST = [0.9, 1, 1.1, 1.2, 1.3, 1.4]
ABS = [True, False]

grid_TabPFN = []

for base_cfg, st, abs_flag in itertools.product(grid, ST, ABS):
    cfg = base_cfg.copy()
    cfg.update({
        "ST": st,
        "ABS": abs_flag
    })
    grid_TabPFN.append(cfg)

print(f"Total valid preprocessing + TabPFN pipelines: {len(grid_TabPFN)}")

df_grid_TabPFN = pd.DataFrame(grid_TabPFN)

df_grid_TabPFN.to_csv(
    ROOT / "Data" / "Preprocessing_Grid" / "preprocessing_grid_TabPFN.csv",
    index=False
)

#--------------------------------------------------------------
# Define Grid for TabICL

Factor = np.linspace(0.8, 1.35, 12)

grid_TabICL = []

for base_cfg, factor in itertools.product(grid, Factor):
    cfg = base_cfg.copy()
    cfg.update({
        "Factor": factor
    })
    grid_TabICL.append(cfg)

print(f"Total valid preprocessing + TabICL pipelines: {len(grid_TabICL)}")

df_grid_TabICL = pd.DataFrame(grid_TabICL)

df_grid_TabICL.to_csv(
    ROOT / "Data" / "Preprocessing_Grid" / "preprocessing_grid_TabICL.csv",
    index=False
)


#--------------------------------------------------------------
# Define Grid for PLSQR

n_components = [5, 10, 25]
Factor = [0.9, 1, 1.1, 1.3]

grid_PLSQR = []

for base_cfg, n_components, factor in itertools.product(grid, n_components, Factor):
    cfg = base_cfg.copy()
    cfg.update({
        "n_components": n_components,
        "Factor": factor
    })
    grid_PLSQR.append(cfg)

print(f"Total preprocessing + PLSQR pipelines: {len(grid_PLSQR)}")

df_grid_PLSQR = pd.DataFrame(grid_PLSQR)

df_grid_PLSQR.to_csv(
    ROOT / "Data" / "Preprocessing_Grid" / "preprocessing_grid_PLSQR.csv",
    index=False
)


#--------------------------------------------------------------
# Define Grid for QRF

max_features = [0.35, 0.5,0.65]
Factor = [0.9, 1, 1.1, 1.3]

grid_QRF = []

for base_cfg, max_features, factor in itertools.product(grid, max_features, Factor):
    cfg = base_cfg.copy()
    cfg.update({
        "max_features": max_features,
        "Factor": factor
    })
    grid_QRF.append(cfg)

print(f"Total preprocessing + QRF pipelines: {len(grid_QRF)}")

df_grid_QRF = pd.DataFrame(grid_QRF)

df_grid_QRF.to_csv(
    ROOT / "Data" / "Preprocessing_Grid" / "preprocessing_grid_QRF.csv",
    index=False
)