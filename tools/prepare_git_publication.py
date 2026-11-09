"""Refresh integrity hashes for the curated final-design folder; no Git filtering."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def main():
 base=ROOT/'final_designs';files=[]
 for folder in ['final_designs','final_old','reference_inputs']:
  for p in sorted((ROOT/folder).rglob('*')):
   if not p.is_file() or p==base/'manifest.json':continue
   if p.suffix in ['.log','.sid','.npz','.tmp'] or 'sondata' in p.parts:continue
   files.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 for name in ['reference_taper_geometry.png','reference_passband_plot.png']:
  p=ROOT/name;files.append(dict(path=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 result=dict(schema_version=3,purpose='Integrity manifest only; does not control Git tracking. Includes final_old history and supplied reference_inputs.',full_filter_status='4.92 GHz point complete; final full-band results pending. final_old is historical.',files=files,total_bytes=sum(f['bytes'] for f in files))
 (base/'manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(f"Recorded {len(files)} artifact hashes; {result['total_bytes']} bytes.")
if __name__=='__main__':main()
