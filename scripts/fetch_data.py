"""Retrieve the public FatigueData-AM2022 files (Scientific Data 2023; figshare DOI 10.6084/m9.figshare.22337629) and check them.
Usage: python3 fetch_data.py OUTDIR
The figshare API lists the article's files; we download the Excel and JSON files and print SHA-256 sums. The analysis in this package used
files with the sums below; a different figshare version can give different sums, so a mismatch is a warning, not an error (compare the
printed counts: 1,610 S-N datasets, 15,146 S-N points, 2,460 parameter rows)."""
import sys, os, json, hashlib, urllib.request
EXPECTED = {'FatigueData-AM2022.xlsx': '95c0c41266f32a275bcef834bfd7674079861a2d00d3974e5e6fd0ad7100528b', 'FatigueData-AM2022.json': '36a5ed57538739c37056bc956880f7c89bed009090999ea214d42612cf2c40de'}
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
MANUAL = ('Automatic download failed. Download FatigueData-AM2022.xlsx and FatigueData-AM2022.json by hand from https://doi.org/10.6084/m9.figshare.22337629 '
          'and place them in ' + out + '/ (then run sh run_all.sh).')
try: meta = json.load(urllib.request.urlopen('https://api.figshare.com/v2/articles/22337629', timeout=60))
except Exception as e: sys.exit(f'{e}\n{MANUAL}')
print('title:', meta.get('title'), '| version:', meta.get('version'), '| license:', (meta.get('license') or {}).get('name'))
for f in meta['files']:
    name = f['name']
    if not name.lower().endswith(('.xlsx', '.json')): continue
    dest = os.path.join(out, name); print('downloading', name, f.get('size'), 'bytes'); urllib.request.urlretrieve(f['download_url'], dest)
    h = hashlib.sha256(open(dest, 'rb').read()).hexdigest(); exp = EXPECTED.get(name)
    print('  sha256', h, '| matches the analysed file' if h == exp else '| DIFFERS from the analysed file (check version)' if exp else '')
