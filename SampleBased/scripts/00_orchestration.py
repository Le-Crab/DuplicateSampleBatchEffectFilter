#!/home/vlad/venvs/vlads-py312/bin/python
"""
Batch effect filtering pipeleine based on paired samples and ML.
Pipeline Orchestration script that executes workflow steps sequentially.
Intended for cloning and execution with some prior configuration.
"""

import subprocess
import sys
from pathlib import Path

def check_dependencies():
    """Verify Python dependencies are installed."""
    required = ["numpy", "pandas", "xgboost", "matplotlib" ]  # adjust as needed
    missing = []
    
    for package in required:
        try:
            __import__(package)
        except ImportError:
            missing.append(package)
    
    if missing:
        print("Missing dependencies. Install with:")
        print(f"pip install {' '.join(missing)}")
        sys.exit(1)

def run_step(script_name, step_number):
    """Execute a single pipeline step with basic error handling."""
    print(f"\n{'='*50}")
    print(f"STEP {step_number}: Executing {script_name}")
    print(f"{'='*50}")
    
    script_path = Path(__file__).parent / script_name
    
    if not script_path.exists():
        print(f"ERROR: {script_name} not found at {script_path}")
        sys.exit(1)
    
    try:
        # Run python scripts with the current interpreter; run shell scripts with bash;
        # run R scripts with Rscript.
        if script_path.suffix == '.sh':
            cmd = ['bash', str(script_path)]
        elif script_path.suffix == '.R':
            cmd = ['Rscript', str(script_path)]
        else:
            cmd = [sys.executable, str(script_path)]

        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False  # Stream output to console in real-time
        )
        print(f"✓ {script_name} completed successfully")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ {script_name} failed with exit code {e.returncode}")
        sys.exit(1)


def main():
    """Execute the pipeline: label, preprocess, train, predict, filter, arrange, plot."""
    
    print("Starting Pipeline")
    print(f"Working directory: {Path(__file__).parent.absolute()}")
    
    # check for tools
    check_dependencies()
    # Define pipeline steps in execution order
    pipeline_steps = [
        ("01_weak_sup_labeling.py", 1),   # Step 1: Labeling of training data
        ("02_prepare_training_data.py", 2),       # Step 2: Data cleaning/preparation
        ("03_train_xgb_classifier.py", 3),       # Step 3: Train XGBoost classifier
        ("04_predict_variants.py", 4),       # Step 4: Predict erroneous variants
        ("05_filter_bcf_files.py", 5),       # Step 5: Filter VCF/BCF files
        ("06_merge_bcf_files.py", 6),       # Step 6: Merge filtered BCF files
        ("07_bcf2pca.sh", 7),       # Step 7: BCF -> PCA/GWAS pipeline
        ("08_vec2pca.R", 8),       # Step 8: PCA plot from .vec files
        ("09_mhqq.R", 9),       # Step 9: Manhattan + QQ plots from GWAS .assoc
    ]
    
    # Execute each step sequentially
    for script_name, step_num in pipeline_steps:
        run_step(script_name, step_num)
    
    print(f"\n{'='*50}")
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()