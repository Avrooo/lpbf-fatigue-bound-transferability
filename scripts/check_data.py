"""Verify the two input files against data/SHA256SUMS. Usage: python3 scripts/check_data.py [DATA_DIR] (default: data).
Exit status 1 on a missing file or a checksum mismatch, unless the environment variable ALLOW_DATA_MISMATCH=1 is set
(then a mismatch is reported but the run continues; use only to inspect a different upstream version, never for the reported results)."""
import sys, os, hashlib
D = sys.argv[1] if len(sys.argv) > 1 else 'data'; bad = 0
for line in open(os.path.join(D, 'SHA256SUMS')):
    if not line.strip() or line.startswith('#'): continue
    want, name = line.split(); p = os.path.join(D, name.lstrip('*'))
    if not os.path.exists(p): print('MISSING', p); bad += 1; continue
    got = hashlib.sha256(open(p, 'rb').read()).hexdigest(); ok = got == want; print('OK      ' if ok else 'MISMATCH', name, got); bad += not ok
if bad and os.environ.get('ALLOW_DATA_MISMATCH') != '1': sys.exit('Input check failed. Run: python3 scripts/fetch_data.py data  (or place the original files in data/).')
