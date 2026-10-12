"""Portable roots and opt-in guard; no scientific computation at import time."""
from pathlib import Path
import os
import sys

INPUT_ROOT = Path(os.environ.get('KSS_INPUT_ROOT', '/input/ur/data')).resolve()
OUTPUT_ROOT = Path(os.environ.get('KSS_OUTPUT_ROOT', '/output/ur/data')).resolve()

def require_run():
    if '--help' in sys.argv or '--run' not in sys.argv:
        print('Use ../run.py --help. Direct execution requires --run and KSS_INPUT_ROOT/KSS_OUTPUT_ROOT.')
        raise SystemExit(0 if '--help' in sys.argv else 2)
    if INPUT_ROOT == OUTPUT_ROOT or INPUT_ROOT in OUTPUT_ROOT.parents or OUTPUT_ROOT in INPUT_ROOT.parents:
        raise SystemExit('Input and output roots must be separate, non-nested directories.')
    stage = {'rf_ecological_context.py':'rf_lopo', 'composition_numeric.py':'composition',
             'breadth_rf_depth.py':'breadth_rf_depth'}[Path(sys.argv[0]).name]
    if (OUTPUT_ROOT/stage).exists():
        raise SystemExit('Refusing an existing stage output: '+str(OUTPUT_ROOT/stage))
    sys.argv.remove('--run')
