"""Read-only EM diagnostics; no new solver jobs and no reference updates."""
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper import taper3 as work
from geometry.taper.progressive_fitting import fit
from geometry.taper.staged_fitting import circuit


def profile_direct(case):
    source=work.NATIVE/case/'direct/generated_ncell_4.snp'
    out=work.OUT/case/'direct/C3_profile_older1percent'
    measured=json.loads((work.OUT/case/'inductance.json').read_text())
    ls=[measured['L0_nH']]+work.REFERENCE_FIT['L_nH']
    cs=[251.]+work.REFERENCE_FIT['C_fF']
    gs=[None]+work.REFERENCE_FIT['G_fF']
    free=[('G',1)]+[('C',i) for i in range(1,5)]+[('G',i) for i in range(2,5)]
    bounds=[(.1,3000.)]+[(.99*(cs if k=='C' else gs)[i],1.01*(cs if k=='C' else gs)[i]) for k,i in free[1:]]
    rows=[]
    for value in (0.,100.,176.,251.,400.,600.,1000.):
        cp=cs.copy();cp[0]=value
        result=fit(source,out/f'C3_{value:g}ff',ls,cp,gs,free,
            capacitance_bounds=bounds,title=f'Direct T3: C0 fixed at {value:g} fF')
        rows.append(dict(C3_fixed_fF=value,RMS=result['complex_RMS'],
            G_T3_T2_fF=result['G_fF'][1],delay_ps=result['delay_ps'],
            C_fF=result['C_fF'],G_fF=result['G_fF'],active_bounds=result['active_bounds']))
    (out/'profile.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2),flush=True)
    return rows


def fixture_check(case):
    """Check exported port-line Z0, without promoting a new acceptance model."""
    out=work.OUT/case/'direct/fixture_diagnostic';out.mkdir(exist_ok=True)
    net=read_snp(work.NATIVE/case/'direct/generated_ncell_4.snp')
    ix=np.linspace(0,len(net.f)-1,1501).astype(int);f=net.f[ix];obs=net.s[ix]
    ls=[json.loads((work.OUT/case/'inductance.json').read_text())['L0_nH']]+work.REFERENCE_FIT['L_nH']
    cs=[251.]+work.REFERENCE_FIT['C_fF'];gs=[None]+work.REFERENCE_FIT['G_fF']
    rows=[]
    for z in (50.,56.887):
        for objective in ('S21','all_S'):
            def prediction(x,frequency=f):
                cp=cs.copy();gp=gs.copy();cp[0],gp[1]=x[:2]
                return circuit(frequency,ls,cp,gp,left_ps=x[2],right_ps=x[3],line_z=z)
            def residual(x):
                d=prediction(x)-obs
                q=d[:,1,0] if objective=='S21' else d.reshape(-1)
                return np.r_[q.real,q.imag]
            candidates=[]
            for c in (10,176,251,600,990):
                for g in (100,153,250,500):
                    candidates.append(least_squares(residual,[c,g,8,5],
                        bounds=([0,.1,-40,-40],[1000,3000,40,40]),
                        x_scale=[251,153,10,10],max_nfev=180))
            best=min(candidates,key=lambda q:q.fun@q.fun)
            p=prediction(best.x,net.f);delta=p-net.s
            row=dict(line_z_ohm=z,objective=objective,x=best.x.tolist(),
                active_bounds=best.active_mask.tolist(),
                S21_RMS=float(np.sqrt(np.mean(abs(delta[:,1,0])**2))),
                S11_RMS=float(np.sqrt(np.mean(abs(delta[:,0,0])**2))),
                S22_RMS=float(np.sqrt(np.mean(abs(delta[:,1,1])**2))))
            rows.append(row);print('FIXTURE DIAGNOSTIC',row,flush=True)
            np.savez(out/f'{objective}_z{z:g}.npz',f=net.f,actual=net.s,model=p)
    (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


def internal_line_check(case):
    """Fixed geometry-derived top-rail length, diagnostic model only."""
    root=work.OUT/case/'response';out=root/'internal_line_diagnostic';out.mkdir(exist_ok=True)
    net=read_snp(work.NATIVE/case/'response/generated_ncell_4.snp')
    meta=json.loads((work.NATIVE/case/'response/geometry.json').read_text())
    s=meta['t3_settings']['body'];length=s['clearance']['d_right']+s['cap']['t_right']+s['meander']['w']/2
    tau=length*np.sqrt(6.5538)/299.792458  # um divided by speed in um/ps.
    ls=[json.loads((work.OUT/case/'inductance.json').read_text())['L0_nH']]+work.REFERENCE_FIT['L_nH']
    old=np.array(work.REFERENCE_FIT['C_fF']+work.REFERENCE_FIT['G_fF'][1:])
    ix=np.linspace(0,len(net.f)-1,1501).astype(int);f=net.f[ix];obs=net.s[ix,1,0]
    def calculate(freq,x,z,t):
        cp=[x[0]]+x[3:7].tolist();gp=[x[1],x[2]]+x[7:10].tolist()
        w=2*np.pi*freq;a=np.ones(len(w),complex);b=np.zeros(len(w),complex);c=b.copy();d=a.copy()
        for i,(l,cap,g) in enumerate(zip(ls,cp,gp)):
            zs=1/(1j*w*g*1e-15);b=b+a*zs;d=d+c*zs
            y=1j*w*cap*1e-15+1/(1j*w*l*1e-9);a=a+b*y;c=c+d*y
            if i==0:
                co=np.cos(w*t*1e-12);si=1j*np.sin(w*t*1e-12)
                a,b,c,d=a*co+b*si/z,a*si*z+b*co,c*co+d*si/z,c*si*z+d*co
        return 2/(a+b/50+c*50+d)*np.exp(-1j*w*x[-1]*1e-12)
    lower=np.r_[0,.1,.1,.99*old,-100];upper=np.r_[1000,3000,3000,1.01*old,100]
    scale=np.r_[251,261,153,old,20]
    rows=[]
    for z,t in ((50.,0.),(50.,tau),(56.887,tau)):
        def residual(x):
            q=calculate(f,x,z,t)-obs;return np.r_[q.real,q.imag]
        candidates=[]
        for cap in (100,251,600,900):
            for gj in (153,500):
                candidates.append(least_squares(residual,np.r_[cap,261,gj,old,5],bounds=(lower,upper),
                    x_scale=scale,max_nfev=400))
        best=min(candidates,key=lambda q:q.fun@q.fun);m=calculate(net.f,best.x,z,t)
        row=dict(line_z=z,line_delay_ps=t,line_length_um=length,x=best.x.tolist(),
            active_bounds=best.active_mask.tolist(),complex_RMS=float(np.sqrt(np.mean(abs(m-net.s[:,1,0])**2))))
        rows.append(row);print('INTERNAL LINE',row,flush=True)
    (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n');return rows


def check_node_separation(case):
    """Sever each coil's top neck numerically; capacitor nodes must separate.

    This is analysis on copies of polygons, never a geometry modification.
    DC connectivity alone cannot detect a short between grounded inductors.
    """
    import gdstk
    import xml.etree.ElementTree as ET
    from geometry.prepare_sonnet import tag,points
    native=work.NATIVE/case/'response'
    top=gdstk.read_gds(str(native/'geometry.gds')).top_level()[0]
    cuts=[];coils=[]
    def walk(cell,offset=(0,0)):
        for ref in cell.references:
            if ref.rotation or ref.x_reflection or ref.magnification!=1:raise ValueError('Unsupported transform')
            xy=np.array(offset)+ref.origin
            if 'COIL' in ref.cell.name or 'MEANDER' in ref.cell.name:
                bb=np.array(ref.cell.bounding_box())+xy
                x=float(bb[:,0].mean());y=float(bb[1,1])
                cuts.append(gdstk.rectangle((x-6,y-7),(x+6,y-3)))
                coils.append(dict(name=ref.cell.name,x=x,y_top=y))
            else:walk(ref.cell,xy)
    walk(top)
    polys=[]
    for q in ET.parse(native/'generated_ncell_4.sonx').iter(tag('PlanarPolygon')):
        xy=points(q);xy[:,1]=448-xy[:,1];polys.append(gdstk.Polygon(xy))
    separated=gdstk.boolean(polys,cuts,'not',precision=.001)
    result=dict(coils=coils,native_nodes_after_cut=len(separated),
        expected_nodes=len(coils)+2,pass_no_hidden_DC_bridge=len(separated)==len(coils)+2,
        scope='Native metal analyzed after severing only shunt coil top leads; production metal unchanged')
    (work.OUT/case/'node_separation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('NODE SEPARATION',result,flush=True);return result


if __name__=='__main__':
    profile_direct('w84_dr256_rightpad148_i36_j22')
