"""No-EM diagnosis: infer the T2 shunt admittance demanded by frozen downstream metal.

This is an algebraic consistency check, not a capacitance fit or accepted design.
It never releases T1 or regular parameters. Known input coupling is a hypothesis.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.locked_taper2 import ROOT,FIXED,response


def downstream_abcd(f):
    w=2*np.pi*f
    m=np.tile(np.eye(2,dtype=complex),(len(f),1,1))
    for i,cg in enumerate([FIXED['Cg_input_fF'],FIXED['Cg_junction_fF'],FIXED['Cg_output_fF']]):
        a=np.tile(np.eye(2,dtype=complex),(len(f),1,1));a[:,0,1]=1/(1j*w*cg*1e-15);m=m@a
        if i<2:
            l=[FIXED['L_taper_nH'],FIXED['L_regular_nH']][i]*1e-9
            c=[FIXED['C_taper_fF'],FIXED['C_regular_fF']][i]*1e-15
            a=np.tile(np.eye(2,dtype=complex),(len(f),1,1));a[:,1,0]=1j*w*c+1/(1j*w*l);m=m@a
    return m


def required_shunt(f,s21,cg):
    w=2*np.pi*f;z=1/(1j*w*cg*1e-15);m=downstream_abcd(f)
    a,b,c,d=m[:,0,0],m[:,0,1],m[:,1,0],m[:,1,1]
    base=a+b/50+c*50+d+z*(c+d/50)
    coefficient=(50+z)*(a+b/50)
    delay=(FIXED['fixture_delay_left_ps']+FIXED['fixture_delay_right_ps'])*1e-12
    undelayed=s21*np.exp(1j*w*delay)
    return (2/undelayed-base)/coefficient


def required_junction(f,s21,c0,cg,l2):
    """Algebraic junction consistency at assumed T2 C values, never a free fit."""
    w=2*np.pi*f
    prefix=np.tile(np.eye(2,dtype=complex),(len(f),1,1));prefix[:,0,1]=1/(1j*w*cg*1e-15)
    shunt=np.tile(np.eye(2,dtype=complex),(len(f),1,1));shunt[:,1,0]=1j*w*c0*1e-15+1/(1j*w*l2*1e-9)
    prefix=prefix@shunt
    # Strip the fixed first series coupling from the downstream ABCD block.
    inverse=np.tile(np.eye(2,dtype=complex),(len(f),1,1));inverse[:,0,1]=-1/(1j*w*FIXED['Cg_input_fF']*1e-15)
    suffix=inverse@downstream_abcd(f)
    derivative=np.zeros_like(prefix);derivative[:,0,1]=1
    def denominator(m):return m[:,0,0]+m[:,0,1]/50+m[:,1,0]*50+m[:,1,1]
    delay=(FIXED['fixture_delay_left_ps']+FIXED['fixture_delay_right_ps'])*1e-12
    den=2/(s21*np.exp(1j*w*delay))
    return (den-denominator(prefix@suffix))/denominator(prefix@derivative@suffix)


def diagnose(name):
    out=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((out/'trial_locked.json').read_text());net=read_snp(trial['snp']);l=trial['fit']['L_t2_nH']
    mask=(net.f>=4.6e9)&(net.f<=5.5e9)
    f=net.f[mask];observed=net.s[mask,1,0]
    y=required_shunt(f,observed,459)
    required_c=(y.imag+1/(2*np.pi*f*l*1e-9))/(2*np.pi*f)*1e15
    # The identity must recover a synthetic ideal T2 exactly.
    sy=response(f,[100,459],l)[:,1,0]
    exact=required_shunt(f,sy,459)
    expected=1j*2*np.pi*f*100e-15+1/(1j*2*np.pi*f*l*1e-9)
    if not np.allclose(exact,expected,rtol=1e-8,atol=1e-10):raise ValueError('Admittance inversion failed synthetic check')
    peaks,_=find_peaks(abs(observed),prominence=.02)
    samples=[]
    junction=required_junction(f,observed,100,459,l)
    for frequency in [4.7,4.85,4.9,5.,5.126,5.3]:
        i=int(np.argmin(abs(f-frequency*1e9)))
        samples.append({'frequency_GHz':float(f[i]/1e9),'required_C0_fF':float(required_c[i]),'required_conductance_mS':float(y.real[i]*1e3),
                        'junction_C_fF_if_T2_at_targets':float(-1/(2*np.pi*f[i]*junction.imag[i])*1e15),
                        'junction_R_ohm_if_T2_at_targets':float(junction.real[i])})
    fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True)
    axes[0].plot(f/1e9,abs(observed),'r',label='Sonnet');axes[0].plot(f/1e9,abs(sy),'k--',label='Target, fixed downstream');axes[0].legend();axes[0].set_ylabel('|S21|')
    axes[1].plot(f/1e9,required_c);axes[1].axhline(100,color='k',ls='--');axes[1].set_ylabel('Required C0 (fF)');axes[1].set_ylim(-500,500)
    axes[2].plot(f/1e9,y.real*1e3);axes[2].axhline(0,color='k',ls='--');axes[2].set(ylabel='Required G (mS)',xlabel='Frequency (GHz)',ylim=(-10,10))
    for ax in axes:ax.grid(alpha=.3)
    fig.suptitle('Frozen downstream consistency check: input Cg = 459 fF');fig.tight_layout();fig.savefig(out/'required_shunt_diagnostic.png',dpi=160);plt.close(fig)
    result={'purpose':'Algebraic diagnostic, not physical extraction','fixed_input_Cg_fF':459,'synthetic_identity_pass':True,'actual_peaks_GHz':(f[peaks]/1e9).tolist(),'samples':samples}
    (out/'required_shunt_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':diagnose('locked_t1_idc30_span180_trace4_ground160_rise320_right560')
