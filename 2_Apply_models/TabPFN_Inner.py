# This is a script run on a hpc cluster for the tuning. Ask viacheslav.barkov@uni-osnabrueck.de for cluster implementation
# We renamed a few files, if there are any linking bugs, contact us, it may be due to old foler mismatches
from pathlib import Path
import argparse

import properscoring as ps
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from tabpfn import TabPFNRegressor

from Define_Functions_Pre_Processing import (
    apply_binning,
    apply_preprocessing_pipeline,
    apply_sg,
    apply_snv,
    split_by_fold,
)


def drop_first_columns(df: pd.DataFrame, n: int = 3) -> pd.DataFrame:
    """Always drop the first n columns (default = 3), if present."""
    if df.shape[1] > n:
        return df.iloc[:, n:].copy()
    return df.iloc[:, 0:0].copy()


def interval_score_95_from_quantiles(y, q, quantiles_list, alpha=0.05):
    y = np.asarray(y)
    q = np.asarray(q)
    taus = np.asarray(quantiles_list)

    idx_lower = np.where(np.isclose(taus, alpha / 2))[0][0]
    idx_upper = np.where(np.isclose(taus, 1 - alpha / 2))[0][0]

    q_left = q[:, idx_lower]
    q_right = q[:, idx_upper]

    sharpness = q_right - q_left
    calibration = (np.maximum(0.0, q_left - y) + np.maximum(0.0, y - q_right)) * (
        2.0 / alpha
    )

    return float(np.mean(sharpness + calibration))


def crps_from_quantiles_step_cdf(y, q, quantiles_list):
    y = np.asarray(y)
    q = np.asarray(q)
    taus = np.asarray(quantiles_list)

    # Probability masses from quantile gaps
    base_weights = np.diff(np.concatenate(([0.0], taus, [1.0])))  # (K+1,)

    # Ensemble support (duplicate last quantile for right tail)
    y_pred_ensemble = np.concatenate([q, q[:, -1:]], axis=1)  # (n_samples, K+1)

    # Repeat weights for all samples
    weights_matrix = np.tile(base_weights, (y_pred_ensemble.shape[0], 1))

    # CRPS
    crps = ps.crps_ensemble(y, y_pred_ensemble, weights=weights_matrix)

    return float(np.mean(crps))


def run_iteration(
    i,
    hyper_row,
    Params,
    target,
    Full_Dataset,
    Outer_Fold,
    soil_properties,
    quantiles_list,
    L,
):
    # Results caching
    # Create a unique path for this specific iteration
    cache_dir = Path("Results/TabPFN/Inner_CV/.cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"preds_{target}_fold_{i}_hyper_{hyper_row}.npz"

    # Check cache
    if cache_file.exists():
        # FAST PATH: Load existing predictions
        data = np.load(cache_file)
        y_true_all = data["y_true"]
        y_pred_mean_all = data["y_pred_mean"]
        y_pred_quantiles_all = data["y_pred_quantiles"]

    else:
        training_dataset, test_dataset = split_by_fold(
            Full_Dataset, Outer_Fold, "Outer_Fold", i
        )

        window_length = int(Params["window_length"]) if Params["SG"] else None
        polyorder = int(Params["polyorder"]) if Params["SG"] else None
        deriv = int(Params["deriv"]) if Params["SG"] else None

        y_true_list = []
        y_pred_mean_list = []
        y_pred_quantiles_list = []

        for l in range(1, L + 1):
            inner_fold_path = rf"Data/Folds/Inner_CV/Nested_Outer_{i}.csv" # COULD BE A BUG HERE; I CHANGED IT LATER DUE TO CHANGE IN STRUCTURE OF FILE NAMING 
            Inner_Fold = pd.read_csv(inner_fold_path)

            training_dataset_tune, calibration_dataset_tune = split_by_fold(
                training_dataset, Inner_Fold, "Inner_Fold", l
            )

            training_dataset_tune_processed = apply_preprocessing_pipeline(
                training_dataset_tune,
                snv=Params["SNV"],
                sg=Params["SG"],
                window_length=window_length,
                polyorder=polyorder,
                deriv=deriv,
                bin_size=Params["bin_size"],
            )

            calibration_dataset_tune_processed = apply_preprocessing_pipeline(
                calibration_dataset_tune,
                snv=Params["SNV"],
                sg=Params["SG"],
                window_length=window_length,
                polyorder=polyorder,
                deriv=deriv,
                bin_size=Params["bin_size"],
            )

            training_dataset_tune_ready = drop_first_columns(
                training_dataset_tune_processed
            )
            calibration_dataset_tune_ready = drop_first_columns(
                calibration_dataset_tune_processed
            )

            X_train_tune = training_dataset_tune_ready.drop(
                columns=soil_properties, errors="ignore"
            )
            X_calibration_tune = calibration_dataset_tune_ready.drop(
                columns=soil_properties, errors="ignore"
            )

            y_train_tune = training_dataset_tune_ready[target].values
            y_calibration_tune = calibration_dataset_tune_ready[target].values

            model_tune = TabPFNRegressor(
                softmax_temperature=Params["ST"],
                average_before_softmax=Params["ABS"],
                n_estimators=4,
                random_state=42,
            )

            model_tune.fit(X_train_tune, y_train_tune)

            y_pred_tune = model_tune.predict(
                X_calibration_tune,
                output_type="main",
                quantiles=quantiles_list,
            )

            y_true_list.append(y_calibration_tune)
            y_pred_mean_list.append(np.asarray(y_pred_tune["mean"]).ravel())
            y_pred_quantiles_list.append(np.stack(y_pred_tune["quantiles"], axis=1))

        # Concatenate results from lists
        y_true_all = np.concatenate(y_true_list)
        y_pred_mean_all = np.concatenate(y_pred_mean_list)
        y_pred_quantiles_all = np.vstack(y_pred_quantiles_list)

        # Save to cache
        np.savez_compressed(
            cache_file,
            y_true=y_true_all,
            y_pred_mean=y_pred_mean_all,
            y_pred_quantiles=y_pred_quantiles_all,
        )

    # Metrics
    rmse = np.sqrt(mean_squared_error(y_true_all, y_pred_mean_all))
    mae = mean_absolute_error(y_true_all, y_pred_mean_all)
    r2 = r2_score(y_true_all, y_pred_mean_all)

    interval_score_95 = interval_score_95_from_quantiles(
        y_true_all,
        y_pred_quantiles_all,
        quantiles_list,
    )

    interval_score_90 = interval_score_95_from_quantiles(
        y_true_all,
        y_pred_quantiles_all,
        quantiles_list,
        alpha=0.1,
    )

    crps_output = crps_from_quantiles_step_cdf(
        y_true_all,
        y_pred_quantiles_all,
        quantiles_list,
    )

    return pd.DataFrame(
        [
            {
                "Interval_Score_95": interval_score_95,
                "Interval_Score_90": interval_score_90,
                "CRPS": crps_output,
                "RMSE": rmse,
                "MAE": mae,
                "R2": r2,
                "Outer_Fold": i,
                "Hyper_Row": hyper_row,
                **Params.to_dict(),
            }
        ]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run TabPFN optimization for a specific soil target."
    )
    parser.add_argument(
        "--target",
        type=str,
        required=True,
        help="The name of the soil property target (e.g., SOM, Clay, Sand).",
    )

    args = parser.parse_args()
    target = args.target

    Full_Dataset = pd.read_csv(r"Data/Full_Dataset.csv")
    Outer_Fold = pd.read_csv(r"Data/Folds/Outer/Outer_Folds.csv")
    TabPFN_Grid = pd.read_csv(r"Data/Preprocessing_Grid/preprocessing_grid_TabPFN.csv")

    soil_properties = [
        "Sand",
        "Total_C",
        "SOM",
        "Clay",
        "pH_SMP",
        "P_Meh3",
        "K_Meh3",
        "Al_Meh3",
    ]

    core = np.arange(0.005, 1.00, 0.005)
    tails = [0.001, 0.999]

    quantiles_list = np.sort(np.r_[core, tails])

    k = 10
    L = 10
    # soil_properties_target = ["SOM", "Clay", "P_Meh3", "K_Meh3"]

    # --- 3. EXECUTION LOOP ---
    # We keep the target loop serial (outermost), but parallelize the Fold/Grid search.

    if target not in soil_properties:
        print(f"WARNING: '{target}' is not in the standard list of soil properties.")

    print(f"Processing Target: {target}")

    tasks = []

    for i in range(1, k + 1):
        for hyper_row in range(len(TabPFN_Grid)):
            Params = TabPFN_Grid.iloc[hyper_row]

            tasks.append(
                delayed(run_iteration)(
                    i,
                    hyper_row,
                    Params,
                    target,
                    Full_Dataset,
                    Outer_Fold,
                    soil_properties,
                    quantiles_list,
                    L,
                )
            )

    print(f"  -> Running {len(tasks)} tasks in parallel...")
    results = Parallel(n_jobs=16, verbose=5)(tasks)

    tuning_overview_df = pd.concat(results, ignore_index=True)

    out_folder = Path("Results/TabPFN/Inner_CV/")
    out_folder.mkdir(parents=True, exist_ok=True)
    out_path = out_folder / f"Tuning_Overview_{target}_IS.csv"

    tuning_overview_df.to_csv(out_path, index=False)
    print(f"  -> Saved to {out_path}")
