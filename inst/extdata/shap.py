import xgboost
import pandas as pd
import shap

def shap_py(xgb_model, train_data, out_file):
    print("Reading model...")
    model_xgb_test = xgboost.Booster()
    model_xgb_test.load_model(xgb_model)

    X = pd.read_csv(train_data, index_col = 0, sep = '\t')
    del X['binary_celltype']
    
    # Explain
    print("Calculating SHAP values...")
    explainer = shap.Explainer(model_xgb_test)
    shap_values = explainer(X)

    shap_values_pd = pd.DataFrame.from_records(shap_values.values)
    shap_values_pd.index = X.index
    shap_values_pd.columns = X.columns

    print("Saving SHAP values...")
    shap_values_pd.to_csv(out_file, sep='\t', index=True)

    print("Done")

