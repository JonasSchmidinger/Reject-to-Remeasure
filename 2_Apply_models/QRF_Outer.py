from pathlib import Path
import re
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from quantile_forest import RandomForestQuantileRegressor
from Define_Functions_Pre_Processing import (
    apply_snv,
    apply_sg,
    apply_binning,
    apply_preprocessing_pipeline,
    split_by_fold,
    s_factor
)

# Set directory
ROOT = Path(__file__).resolve().parents[1]

# Define soil properties
soil_properties = ["Total_C","SOM",	"Sand","Clay","pH_SMP",	"P_Meh3","K_Meh3","Al_Meh3"]


# Define quantile taus
quantiles_list = [0.001] \
    + [round(x, 3) for x in np.arange(0.005, 1.00, 0.005)] \
    + [0.999]

quantiles_list

# Drop ID columns
def drop_first_columns(df: pd.DataFrame, n: int = 3) -> pd.DataFrame:
    """Always drop the first n columns (default = 3), if present."""
    if df.shape[1] > n:
        return df.iloc[:, n:].copy()
    return df.iloc[:, 0:0].copy()

all_results = []  # metrics across all datasets/targets


# Store results
all_results = []  

Full_Dataset = pd.read_csv(ROOT / "Data" /  "Full_Dataset.csv")
Full_Dataset

Outer_Fold = pd.read_csv(ROOT / "Data" / "Folds" / "Outer" / "Outer_Folds.csv")
Outer_Fold


TEST_RUN = False   

if TEST_RUN:
    k = 2
    soil_properties_target = ['Sand']  
else:    
    k = 10
    soil_properties_target = ["Total_C","SOM",	"Sand","Clay","pH_SMP",	"P_Meh3","K_Meh3","Al_Meh3"]
    #soil_properties_target = ["Al_Meh3"]


for target in soil_properties_target:
    print(target)
    
    id_all = []
    y_true_all = []
    y_pred_all = []
    y_quantiles_all_folds = []
    
    
    tuning_overview_df = pd.read_csv(ROOT / "Results" / "QRF" / "Inner_CV" / f"Tuning_Overview_{target}_IS.csv")
    
    for i in range(1, k + 1):
        print(i)
        
        training_dataset, test_dataset = split_by_fold(Full_Dataset,Outer_Fold,"Outer_Fold",i)
        
        best_row = (
          tuning_overview_df
              .query("Outer_Fold == @i")
              .sort_values("CRPS", ascending=True)
              .iloc[0]
        )
        
        window_length = int(best_row["window_length"]) if best_row["SG"] else None
        polyorder  = int(best_row["polyorder"]) if best_row["SG"] else None
        deriv  = int(best_row["deriv"]) if best_row["SG"] else None

        
        training_dataset_processed = apply_preprocessing_pipeline(
            training_dataset,
            snv=best_row["SNV"],
            sg=best_row["SG"],
            window_length=window_length,
            polyorder=polyorder,
            deriv=deriv,
            bin_size=best_row["bin_size"]
        )
        
        test_dataset_processed = apply_preprocessing_pipeline(
            test_dataset,
            snv=best_row["SNV"],
            sg=best_row["SG"],
            window_length=window_length,
            polyorder=polyorder,
            deriv=deriv,
            bin_size=best_row["bin_size"]
        )
        
        training_dataset_ready = drop_first_columns(training_dataset_processed)        
        test_dataset_ready = drop_first_columns(test_dataset_processed)
        # training_dataset_ready = drop_first_columns(training_dataset_processed,n = 2)        
        # test_dataset_ready = drop_first_columns(test_dataset_processed,n = 2)
        
        
        X_train = training_dataset_ready.drop(columns=soil_properties, errors='ignore')
        X_test  = test_dataset_ready.drop(columns=soil_properties, errors='ignore')

        y_train = training_dataset_ready[target].values
        y_test  = test_dataset_ready[target].values
        
 
        model = RandomForestQuantileRegressor(
            n_estimators=500,
            max_features=float(best_row["max_features"]),
            random_state=42
        )
        model.fit(X_train, y_train)
        
        y_pred_mean = model.predict(X_test)
        
        quantiles_list_scaled = s_factor(quantiles_list, best_row["Factor"])
        
        fold_quantiles = np.column_stack([
            model.predict(X_test, quantiles=q_scaled) for q_scaled in quantiles_list_scaled
        ])

        y_pred_all.extend(y_pred_mean.tolist())
        y_true_all.extend(y_test.tolist())
        y_quantiles_all_folds.append(fold_quantiles)
        
    
    out_dir = ROOT / "Results" / "QRF" / "Outer"

    # True values
    pd.DataFrame({f"{target}_target": y_true_all}).to_csv(
        rf"{out_dir}\{target}_true.csv",
        index=False
    )
    
    # Mean predictions
    pd.DataFrame({f"{target}_pred": y_pred_all}).to_csv(
        rf"{out_dir}\{target}_pred.csv",
        index=False
    )
    
    # Quantiles
    y_quantiles_all = np.concatenate(y_quantiles_all_folds, axis=0)
    qcols = [f"q{q}" for q in quantiles_list]
    
    pd.DataFrame(y_quantiles_all, columns=qcols).to_csv(
        rf"{out_dir}\{target}_pred_quantiles.csv",
        index=False
    )
    
    # Metrics on concatenated out-of-fold predictions
    r2  = r2_score(y_true_all, y_pred_all)
    mse = mean_squared_error(y_true_all, y_pred_all)
    mae = mean_absolute_error(y_true_all, y_pred_all)
    
    all_results.append({
        "Target": target,
        "R2": r2,
        "MSE": mse,
        "MAE": mae
    })
    
    print(f"{target:>10s}  R2={r2: .3f}  RMSE={mse**0.5: .3f}  MAE={mae: .3f}")
    
    
        
# Summary table
id_all