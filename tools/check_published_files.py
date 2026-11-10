"""Verify curated native/GDS/results/reference hashes without Sonnet."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
def main():
 manifest=json.loads((ROOT/'final_designs/manifest.json').read_text(encoding='utf-8-sig'));issues=[]
 for row in manifest['files']:
  p=ROOT/row['path']
  if not p.is_file():issues.append('Missing: '+row['path']);continue
  if p.stat().st_size!=row['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:issues.append('Changed: '+row['path'])
 if issues:print('\n'.join(issues));return 1
 print(f"Verified {len(manifest['files'])} packaged artifact hashes. Full-filter evidence status is recorded in final_full_filter/metrics.json.");return 0
if __name__=='__main__':sys.exit(main())
