"""Rebuild frozen selected Sonnet projects from repository-contained GDS inputs.

No EM calculation is launched. --verify-only needs no Sonnet installation.
The full-filter layout is regenerated from the current taper assembly's metal.
"""
from pathlib import Path
import argparse,copy,json,re,shutil,subprocess,xml.etree.ElementTree as ET
import numpy as np
import gdstk
ROOT=Path(__file__).resolve().parents[1]
FINAL=ROOT/'final_designs'
NS='https://www.sonnetsoftware.com/schema/project'
ET.register_namespace('',NS);ET.register_namespace('xsi','http://www.w3.org/2001/XMLSchema-instance')
def tag(n):return '{'+NS+'}'+n
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def points(p):return np.array([(float(x),float(y)) for x,y in re.findall(r'\(([^,]+),([^\)]+)\)',' '.join(q.text or '' for q in p.findall(tag('Points'))))])
def signature(e):
 if e is None:return None
 attrs=dict(e.attrib)
 # Native GDS import writes integral layer thicknesses as e.g. 12000.0.
 # Compare the physical value, retaining strict checks on all other fields.
 if 'Thickness' in attrs:attrs['Thickness']=float(attrs['Thickness'])
 return (e.tag,attrs,(e.text or '').strip(),[signature(q) for q in e])
def native_polys(path):return [gdstk.Polygon(points(p),layer=1) for p in ET.parse(path).iter(tag('PlanarPolygon'))]
def gds_polys(path):
 lib=gdstk.read_gds(str(path));assert np.isclose(lib.unit,1e-6);tops=lib.top_level();assert len(tops)==1
 return tops[0].get_polygons()
def flipped(polys,height):return [gdstk.Polygon(np.column_stack((p.points[:,0],height-p.points[:,1])),layer=1) for p in polys]
def same(a,b):return not gdstk.boolean(a,b,'xor',precision=.001)
def full_from_blocks():
 """Assemble ten copies of the saved regular body and mirrored taper metal."""
 meta=read(FINAL/'full_filter/parameters.json');source=gds_polys(FINAL/'taper_3/layout.gds');first=meta['first_regular_left_um'];uw=meta['regular_width_um'];gap=meta['regular_gap_um'];n=meta['n_regular'];height=1500
 width=2*first+n*uw+(n-1)*gap;assert width==meta['w_sim']
 left=gdstk.boolean(source,[gdstk.rectangle((0,0),(first,height))],'and',precision=.001,layer=1)
 right=[gdstk.Polygon(np.column_stack((width-p.points[:,0],p.points[:,1])),layer=1) for p in left]
 # Box-coordinate regular body: native y=456..1066, GDS y=434..1044.
 body=gdstk.boolean(source,[gdstk.rectangle((first,434),(first+uw,1044))],'and',precision=.001,layer=1)
 assert body
 parts=left+right+[gdstk.rectangle((first,1052),(width-first,height),layer=1),gdstk.rectangle((first,0),(width-first,434),layer=1)]
 for i in range(n):
  parts.extend(p.copy().translate(i*(uw+gap),0) for p in body)
  if i<n-1:
   mid=first+i*(uw+gap)+uw+gap/2;parts.append(gdstk.rectangle((mid-2,434),(mid+2,1052),layer=1))
 assert same(parts,gds_polys(FINAL/'full_filter/layout.gds')),'Full layout differs from checked-in input'
 return parts

def verify(folder):
 desc=read(folder/'design.json');project=folder/'project.sonx';tree=ET.parse(project);geo=tree.find(tag('Geometry'));size=geo.find(tag('Box')+'/'+tag('Size'));width,height=desc['box_um']
 assert [float(size.get('X')),float(size.get('Y'))]==[width,height]
 polys=gds_polys(folder/'layout.gds');assert same(flipped(polys,height),native_polys(project)),'GDS/native metal mismatch'
 assert all(np.allclose(p.points/2,np.round(p.points/2)) for p in polys),'Off-grid native metal'
 mesh=geo.find(tag('Box')+'/'+tag('NumCells'));assert [float(mesh.get('X')),float(mesh.get('Y'))]==[width/2,height/2]
 ports=[p for p in geo.iter(tag('Port'))];assert len(ports)==2
 assert all(float(p.find(tag('Impedance')).get('Resistance'))==50 for p in ports)
 if desc['name']=='full_filter':full_from_blocks()
 return dict(design=desc['name'],native_GDS_XOR_empty=True,grid_um=2,ports=2,full_filter_regenerated_from_blocks=desc['name']=='full_filter')

def rebuild(folder,dest,sonnet_bin):
 check=verify(folder);desc=read(folder/'design.json');width,height=desc['box_um'];dest.mkdir(parents=True,exist_ok=True)
 polys=full_from_blocks() if desc['name']=='full_filter' else gds_polys(folder/'layout.gds')
 lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys);lib.write_gds(str(dest/'layout.gds'))
 tree=ET.parse(folder/'project.sonx');original=copy.deepcopy(tree.getroot());geometry=tree.find(tag('Geometry'));old_ports=[]
 for poly in geometry.iter(tag('PlanarPolygon')):
  xy=points(poly)
  for port in poly.iter(tag('Port')):
   i=int(port.get('Vertex'));old_ports.append((copy.deepcopy(port),xy[i],xy[(i+1)%len(xy)]))
 layer=geometry.find('.//'+tag('PlanarTechLayer'));assert layer is not None
 if layer.find(tag('ImportExportMapping')) is None:ET.SubElement(layer,tag('ImportExportMapping'),GDSStreamDataTypes='1:0')
 tree.write(dest/'template.sonx',encoding='utf-8',xml_declaration=True)
 (dest/'metal.lay').write_text('Metal-1 metal_boundary 1 0 0 "SUPERCONDUCTOR" on 0\n',encoding='ascii')
 def run(args,name):
  q=subprocess.run([str(x) for x in args],cwd=dest,capture_output=True,text=True,timeout=90)
  (dest/name).write_text(q.stdout+q.stderr,encoding='utf-8');q.check_returncode()
 run([sonnet_bin/'gds.exe','-v','-i',dest/'template.sonx','-g','-l',dest/'metal.lay','-p','TOP','-u','um','-s',f'{width},{height}','-t','min','-C','2,2','-o',dest/'imported.sonx',dest/'layout.gds'],'import.log')
 result=ET.parse(dest/'imported.sonx');root=result.getroot();geo=root.find(tag('Geometry'));imported=list(geo.iter(tag('PlanarPolygon')))
 assert same([gdstk.Polygon(points(p)) for p in imported],native_polys(folder/'project.sonx')),'Imported metal mismatch'
 # Reattach exact original port edges; vertex order can change during import.
 for p in imported:
  objs=p.find(tag('Objects'))
  if objs is not None:
   for q in list(objs):
    if q.tag==tag('Port'):objs.remove(q)
 for port,a,b in old_ports:
  hits=[]
  for p in imported:
   xy=points(p)
   for i,(c,d) in enumerate(zip(xy,np.roll(xy,-1,axis=0))):
    if (np.allclose(a,c,atol=1e-6) and np.allclose(b,d,atol=1e-6)) or (np.allclose(a,d,atol=1e-6) and np.allclose(b,c,atol=1e-6)):hits.append((p,i))
  assert len(hits)==1,('Port edge ambiguous',port.attrib,len(hits));p,i=hits[0];port.set('Vertex',str(i));objs=p.find(tag('Objects'))
  if objs is None:objs=ET.SubElement(p,tag('Objects'))
  objs.append(port)
 # Preserve nongeometric solver settings and reference planes exactly.
 oldgeo=original.find(tag('Geometry'))
 for name in ['Units','Control','Sweeps','OutputFiles','Optimize']:
  old=original.find(tag(name));current=root.find(tag(name))
  if current is not None:root.remove(current)
  if old is not None:root.append(copy.deepcopy(old))
 oldbox=oldgeo.find(tag('Box'));newbox=geo.find(tag('Box'))
 for q in list(newbox):
  if q.tag==tag('RefPlane'):newbox.remove(q)
 for q in oldbox.findall(tag('RefPlane')):newbox.append(copy.deepcopy(q))
 for name in ['Dielectric','Conductor','DielectricMaterialModel','BoxCover']:
  assert [signature(q) for q in geo.iter(tag(name))]==[signature(q) for q in oldgeo.iter(tag(name))],('Changed material',name)
 for name in ['Units','Control','Sweeps']:assert signature(root.find(tag(name)))==signature(original.find(tag(name)))
 root.find(tag('OutputFiles')).set('Folder','.')
 result.write(dest/'project.sonx',encoding='utf-8',xml_declaration=True)
 run([sonnet_bin/'soncmd.exe','-Read',dest/'project.sonx'],'read.log')
 shutil.copy2(folder/'design.json',dest/'design.json')
 report=dict(**check,solver_started=False,native_reimport_and_read_passed=True,settings_preserved=True)
 (dest/'rebuild_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');return report

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify-only',action='store_true');parser.add_argument('--output-dir',type=Path,default=ROOT/'build/final_designs');parser.add_argument('--sonnet-bin',type=Path,default=Path(r'C:\Program Files\Sonnet Software\19.52.2025.261\bin'));args=parser.parse_args()
 reports=[]
 for d in read(FINAL/'designs.json')['designs']:
  folder=FINAL/d['name'];report=verify(folder) if args.verify_only else rebuild(folder,(args.output_dir/d['name']).resolve(),args.sonnet_bin);reports.append(report);print(json.dumps(report),flush=True)
 print('Verified selected inputs; no EM solver was launched.' if args.verify_only else 'Rebuilt and read every native project; no EM solver was launched.')
if __name__=='__main__':main()
