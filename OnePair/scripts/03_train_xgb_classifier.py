#!/usr/bin/env python3
"""
Script: 03_train_xgb_classifier.py
Purpose: Train XGBoost binary classifier for variant quality prediction.
         Performs stratified train/val/test split, handles class imbalance,
         and saves model with metadata.

Dependencies: pandas, numpy, xgboost, scikit-learn, matplotlib
"""

import json
import sys
from collections import Counter
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    roc_auc_score,
)


# =============================================================================
# METRICS FUNCTIONS
# =============================================================================
def summarize_metrics(y_true: np.ndarray, y_pred: np.ndarray, label: str = "") -> Dict:
    """
    Compute comprehensive classification metrics.
    
    Parameters
    ----------
    y_true : np.ndarray
        True labels
    y_pred : np.ndarray
        Predicted labels
    label : str
        Dataset label for logging
    
    Returns
    -------
    Dict
        Dictionary of metric names and values
    """
    cm = confusion_matrix(y_true, y_pred)
    tn, fp, fn, tp = cm.ravel()
    
    dict_metrics = {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "mcc": matthews_corrcoef(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred),
        "recall": recall_score(y_true, y_pred),
        "specificity": (tn / (tn + fp)) if (tn + fp) else float('nan'),
        "FOR": (fn / (fn + tn)) if (fn + tn) else float('nan'),
        "f1": f1_score(y_true, y_pred)
    }
    
    print(f"\n--- {label} Set Metrics ---")
    print(f"Confusion matrix:\n{cm}")
    for k, v in dict_metrics.items():
        print(f"  {k:20s}: {v:.4f}")
    
    return dict_metrics


# =============================================================================
# TRAINING FUNCTION
# =============================================================================
def train_xgboost_model(
    fp_input: str,
    dir_output: str,
    model_suffix: str,
    test_size: float = 0.05,
    val_size: float = 0.05,
    random_state: int = 123) -> Tuple[XGBClassifier, Dict, Dict, Dict]:

    """
    Train XGBoost classifier with stratified train/val/test split.
    
    Parameters
    ----------
    fp_input : str
        Path to preprocessed training CSV
    dir_output : str
        Output directory for model and artifacts
    model_suffix : str
        Suffix for output file naming
    test_size : float
        Fraction for test set
    val_size : float
        Fraction for validation set
    random_state : int
        Random seed
    
    Returns
    -------
    Tuple[XGBClassifier, Dict, Dict, Dict]
        Trained model and metrics dictionaries
    """
    
    print("\nSUBSTEP 1: LOADING PREPROCESSED DATA")
    
    
    df = pd.read_csv(fp_input)
    df['GT'] = df['GT'].astype('category')
    df['GINeighbors'] = df['GINeighbors'].astype('category')
    print(f"Loaded data shape: {df.dtypes}")
    print(f"Columns: {list(df.columns)}")
    
    # Separate features and labels
    y = df['label']
    X = df.drop(['label'], axis=1)
    
    print(f"Feature shape: {X.shape}")
    print(f"Label distribution:\n{y.value_counts()}")
    
    
    print("\nSUBSTEP 2: TRAIN/VAL/TEST SPLIT")
    
    
    # First split: separate test set
    X_train_val, X_test, y_train_val, y_test = train_test_split(
        X, y,
        test_size=test_size,
        stratify=y,
        random_state=random_state
    )
    
    # Second split: separate validation from training
    val_size_rel = val_size / (1 - test_size)
    X_train, X_val, y_train, y_val = train_test_split(
        X_train_val, y_train_val,
        test_size=val_size_rel,
        stratify=y_train_val,
        random_state=random_state
    )
    
    print(f"Train shape: {X_train.shape}")
    print(f"Val shape:   {X_val.shape}")
    print(f"Test shape:  {X_test.shape}")
    
    # Class distribution
    print(f"\nTrain label counts: {Counter(y_train)}")
    print(f"Val label counts:   {Counter(y_val)}")
    print(f"Test label counts:  {Counter(y_test)}")
    
    # Calculate class imbalance
    pos_count = (y_train == 1).sum() if y_train.dtype != bool else (y_train == True).sum()
    neg_count = (y_train == 0).sum() if y_train.dtype != bool else (y_train == False).sum()
    
    
    print(f"\nClass imbalance handling:")
    print(f"  Positive samples: {pos_count:,}")
    print(f"  Negative samples: {neg_count:,}")    
    
    print("\nSUBSTEP 3: TRAINING XGBOOST CLASSIFIER")
    print(df.head)
    # Configure and train model
    model_xgb = XGBClassifier(
        objective='binary:logistic',
        tree_method='hist',
        use_categorical=True,
        device='cuda',
        enable_categorical=True,
        eval_metric=['logloss'],
        random_state=random_state,
        n_estimators=5000,
        early_stopping_rounds=50,
        learning_rate=0.1,
        max_depth=12,
        min_child_weight=10,
        gamma=1,
        subsample=0.7,
        colsample_bytree=1,
        reg_lambda=10,
        reg_alpha=10,
    )
    
    print("Training XGBoost with early stopping...")
    model_xgb.fit(
        X_train, y_train,
        eval_set=[(X_train, y_train), (X_val, y_val)],
        verbose=False
    )
    
    print(f"Training complete. Best iteration: {model_xgb.best_iteration}")
    
    
    print("\nSUBSTEP 4: EVALUATING MODEL PERFORMANCE")
    
    
    # Make predictions
    y_train_pred = model_xgb.predict(X_train)
    y_val_pred = model_xgb.predict(X_val)
    y_test_pred = model_xgb.predict(X_test)
    
    # Compute metrics
    dict_train_metrics = summarize_metrics(y_train, y_train_pred, "Training")
    dict_val_metrics = summarize_metrics(y_val, y_val_pred, "Validation")
    dict_test_metrics = summarize_metrics(y_test, y_test_pred, "Test")
    
    
    print("\nSUBSTEP 5: SAVING MODEL AND ARTIFACTS")
    
    
    # Create output directory
    dir_output = Path(dir_output)
    dir_output.mkdir(parents=True, exist_ok=True)
    
    # Timestamp for filenames
    tz_jst = timezone(timedelta(hours=9))
    dt_local = datetime.now(tz_jst)
    timestamp = dt_local.strftime("%Y%m%d_%H%M")
    
    # Save model as JSON booster
    fp_model = dir_output / f"{model_suffix}_{timestamp}_booster.json"
    model_xgb.save_model(fp_model)
    print(f"Saved model to: {fp_model}")
    
    # Save metadata
    dict_metadata = {
        'timestamp': timestamp,
        'xgboost_version': xgb.__version__,
        'input_file': str(fp_input),
        'model_suffix': model_suffix,
        'random_state': random_state,
        'dataset_sizes': {
            'train': len(X_train),
            'val': len(X_val),
            'test': len(X_test)
        },
        'hyperparameters': model_xgb.get_params(),
        'best_iteration': int(model_xgb.best_iteration),
        'metrics': {
            'train': dict_train_metrics,
            'val': dict_val_metrics,
            'test': dict_test_metrics
        },
        'feature_names': list(X.columns)
    }
    
    fp_metadata = dir_output / f"{model_suffix}_{timestamp}_metadata.json"
    with open(fp_metadata, 'w') as f:
        json.dump(dict_metadata, f, indent=2, default=str)
    print(f"Saved metadata to: {fp_metadata}")
    
    # Save feature importance plots
    dir_plots = dir_output.parent / 'outputs/plots'
    dir_plots.mkdir(exist_ok=True)
    
    for importance_type in ['gain', 'weight', 'cover']:
        plt.figure(figsize=(10, 6))
        xgb.plot_importance(
            model_xgb, 
            importance_type=importance_type,
            values_format='{v:.0f}',
            title=f'Feature Importance ({importance_type.capitalize()})',
            xlabel=importance_type.capitalize(),
            ylabel='Features'
        )
        fp_plot = dir_plots / f"{model_suffix}_feature_importance_{importance_type}.png"
        plt.tight_layout()
        plt.savefig(fp_plot, dpi=300)
        plt.close()

    print(f"Saved feature importance plots to: {fp_plot.parent}")    
    return model_xgb, dict_train_metrics, dict_val_metrics, dict_test_metrics


# =============================================================================
# MAIN
# =============================================================================
def main() -> int:
    """Main execution function."""
    WD = Path(__file__).parent.parent
    INPUT_PATH = WD / 'data' / 'pg1kg30x6.bisnv.neighbor2all.annot2.samplewise.V.10M_labeled_mlready.csv'
    OUTPUT_DIR = WD / 'models'
    MODEL_SUFFIX = 'SB-OnePair-VDVC'
    TEST_SIZE = 0.05
    VAL_SIZE = 0.05
    RANDOM_STATE = 123

    print(f"Input: {INPUT_PATH}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Model suffix: {MODEL_SUFFIX}")

    # Validate input
    if not INPUT_PATH.exists():
        print(f"[e] Input file not found: {INPUT_PATH}")
        return 1
    
    try:
        print(f"XGBoost version: {xgb.__version__}")
        
        train_xgboost_model(
            fp_input=str(INPUT_PATH),
            dir_output=str(OUTPUT_DIR),
            model_suffix=MODEL_SUFFIX,
            test_size=TEST_SIZE,
            val_size=VAL_SIZE,
            random_state=RANDOM_STATE
        )
        
        print("Model training completed successfully")
        return 0
        
    except Exception as e:
        print(f"Error during model training: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
