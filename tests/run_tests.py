"""Minimal runner for environments without pytest:  python3 tests/run_tests.py   (pytest users: python -m pytest tests/ -v)"""
import sys, os, glob, importlib, unittest, traceback
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); ok = skip = bad = 0
for f in sorted(glob.glob(os.path.join(here, 'test_*.py'))):
    mod = importlib.import_module(os.path.basename(f)[:-3])
    for n in sorted(x for x in dir(mod) if x.startswith('test_')):
        try: getattr(mod, n)(); ok += 1; print('PASS', mod.__name__, n)
        except unittest.SkipTest as e: skip += 1; print('SKIP', mod.__name__, n, '-', e)
        except Exception: bad += 1; print('FAIL', mod.__name__, n); traceback.print_exc()
print(f'{ok} passed, {skip} skipped, {bad} failed'); sys.exit(1 if bad else 0)
