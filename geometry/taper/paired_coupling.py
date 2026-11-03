"""Cancel shared downstream response between restored and direct-feed tests.

This is a diagnostic of a local coupling block, not an automatic acceptance
of a frequency-independent capacitor. Lead assumptions and residuals are saved.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from circuit_common import read_snp
from geometry.prepare_sonnet import tag,points
from geometry.taper.bridged_fixtures import NATIVE,CASES,VERIFIED,ROOT
from geometry.taper.audit_response_settings import signature
from geometry.taper.run_bridged_stage import OUT,digest


def abcd(s):
    a,b,c,d=s[:,0,0],s[:,0,1],s[:,1,0],s[:,1,1]
    m=np.empty_like(s)
    m[:,0,0]=((1+a)*(1-d)+b*c)/(2*c)
    m[:,0,1]=50*((1+a)*(1+d)-b*c)/(2*c)
    m[:,1,0]=((1-a)*(1-d)-b*c)/(100*c)
    m[:,1,1]=((1-a)*(1+d)+b*c)/(2*c)
    return m


def line(f,length_um,z=56.8872):
    theta=2*np.pi*f*length_um*1e-6*np.sqrt(6.55377)/299792458
    m=np.empty((len(f),2,2),complex)
    m[:,0,0]=m[:,1,1]=np.cos(theta)
    m[:,0,1]=1j*z*np.sin(theta);m[:,1,0]=1j*np.sin(theta)/z
    return m


def reuse_t1_control():
    original=Path(VERIFIED['artifacts']['sonnet_project'])
    new=NATIVE/CASES[1]/'generated_ncell_2.sonx'
    a,b=signature(original),signature(new)
    checks={k:a[k]==b[k] for k in a}
    geometry=lambda path:[gdstk.Polygon(points(q)) for q in ET.parse(path).iter(tag('PlanarPolygon'))]
    checks['native_geometry_XOR_empty']=not bool(gdstk.boolean(geometry(original),geometry(new),'xor',precision=.001))
    status=json.loads((ROOT/'output_data/taper_exploration/finger536_t20_right310_input6_j56/run_status.json').read_text())
    snp=Path(VERIFIED['artifacts']['touchstone'])
    checks['archived_project_hash']=status['project_sha256']==digest(original)
    checks['archived_snp_hash']=status['snp_sha256']==digest(snp)
    out=OUT/CASES[1];out.mkdir(parents=True,exist_ok=True)
    report=dict(source_project=str(original),matching_prepared_project=str(new),source_snp=str(snp),
                checks=checks,reused_without_new_EM=True)
    (out/'control_reuse_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    if not all(checks.values()):raise ValueError(report)
    return snp,out


def analyze(full_snp,direct_snp,out,full_lead_um,direct_lead_um):
    full,direct=read_snp(full_snp),read_snp(direct_snp)
    mapping={round(f):i for i,f in enumerate(direct.f)}
    pairs=[(i,mapping[round(f)]) for i,f in enumerate(full.f) if round(f) in mapping]
    if len(pairs)<20:raise ValueError('Insufficient exact shared frequency samples')
    ia,ib=np.array(pairs).T;f=full.f[ia]
    # M_full * inverse(M_direct) cancels their common downstream block.
    relative=abcd(full.s[ia])@np.linalg.inv(abcd(direct.s[ib]))
    pi=np.linalg.inv(line(f,full_lead_um))@relative@line(f,direct_lead_um)
    b=pi[:,0,1];w=2*np.pi*f
    apparent=-1/(w*b.imag)*1e15
    r=f/5e9;basis=np.c_[-1/r,r]
    coefficient=np.linalg.lstsq(basis,b.imag,rcond=None)[0]
    prediction=basis@coefficient
    result=dict(method='Paired ABCD downstream cancellation; local pi block',samples=len(f),
                interval_GHz=[float(f[0]/1e9),float(f[-1]/1e9)],
                full_lead_um=full_lead_um,direct_lead_um=direct_lead_um,
                lead_assumptions=dict(Z_ohm=56.8872,Eeff=6.55377),
                median_apparent_Cg_fF=float(np.median(apparent)),
                apparent_Cg_range_fF=[float(apparent.min()),float(apparent.max())],
                Cg_after_series_L_separation_fF=float(1/(2*np.pi*5e9*coefficient[0])*1e15),
                residual_series_L_nH=float(coefficient[1]/(2*np.pi*5e9)*1e9),
                series_B_imaginary_fit_RMSE_ohm=float(np.sqrt(np.mean((b.imag-prediction)**2))),
                series_B_real_range_ohm=[float(b.real.min()),float(b.real.max())],
                acceptance_pending=True,
                limitation='Locality and effective lead models must be validated; this diagnostic alone is not a certified coupling extraction')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    np.savez(out/'paired_pi_data.npz',f=f,pi=pi,apparent_Cg_fF=apparent)
    (out/'paired_coupling.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    snp,out=reuse_t1_control()
    analyze(snp,NATIVE/CASES[0]/'generated_ncell_2.snp',out/'paired_input_diagnostic',344,372)
