#!/home/vlad/venvs/vlads-py312/bin/python
"""scripts/00_orchestration.py

LociBased pipeline orchestrator.

Runs the numbered workflow steps sequentially:
    01) prepare ML-ready training table
    02) train XGBoost model
    03) predict discordant/erroneous variant IDs
    04) filter ONE multi-sample VCF/BCF using predicted IDs
    05) convert filtered BCF -> PLINK/PCA/GWAS
    06) PCA plotting
    07) Manhattan + QQ plotting
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
        ("01_prepare_training_data.py", 1),  # ML-ready training table
        ("02_train_xgb_classifier.py", 2),   # train model
        ("03_predict_variants.py", 3),       # predict erroneous IDs
        ("04_filter_bcf_file.py", 4),        # filter multi-sample BCF/VCF
        ("05_bcf2pca.sh", 5),                # PCA/GWAS
        ("06_vec2pca.R", 6),                 # PCA plots
        ("07_mhqq.R", 7),                    # Manhattan/QQ plots
    ]
    
    # Execute each step sequentially
    for script_name, step_num in pipeline_steps:
        run_step(script_name, step_num)
    
    print(f"\n{'='*50}")
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()