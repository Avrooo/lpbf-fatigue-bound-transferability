"""Write results/run_manifest.json: Python and package versions, input checksums, git commit (if any), options and completion status.
Usage: python3 scripts/write_manifest.py STATUS [OPTIONS...]   (STATUS is 'complete' or 'failed')"""
import sys, os, json, hashlib, platform, subprocess, datetime, importlib.metadata as md
m = {'status': sys.argv[1], 'options': sys.argv[2:], 'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
     'python': platform.python_version(), 'platform': platform.platform(),
     'packages': {p: md.version(p) for p in ('numpy', 'scipy', 'pandas', 'openpyxl', 'matplotlib', 'python-docx')},
     'inputs_sha256': {f: hashlib.sha256(open('data/' + f, 'rb').read()).hexdigest() for f in ('FatigueData-AM2022.xlsx', 'FatigueData-AM2022.json') if os.path.exists('data/' + f)},
     'random_seeds': 'fixed on the command line of every step in run_all.sh (reproduce 11, learning curves 51/52, within-paper 100 splits seed 5, coverage_conservatism 11, digitization seeds 1-5)'}
try: m['git_commit'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], stderr=subprocess.DEVNULL).decode().strip()
except Exception: m['git_commit'] = None
os.makedirs('results', exist_ok=True); json.dump(m, open('results/run_manifest.json', 'w'), indent=2); print('wrote results/run_manifest.json', m['status'])
