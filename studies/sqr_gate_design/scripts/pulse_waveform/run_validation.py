"""Wrapper to run validate_results.py and capture output to a file."""
import sys
import os
import traceback
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
LOG_FILE = STUDY_DIR / "data" / "validation_output.txt"

logf = open(LOG_FILE, "w", encoding="utf-8")

try:
    sys.stdout = logf
    sys.stderr = logf

    os.chdir(STUDY_DIR)
    sys.path.insert(0, str(SCRIPT_DIR))

    from validate_results import main
    main()
except Exception:
    traceback.print_exc(file=logf)
finally:
    logf.flush()
    logf.close()
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    print("Wrapper complete. Check data/validation_output.txt")
