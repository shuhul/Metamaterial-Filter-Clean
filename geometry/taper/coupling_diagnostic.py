"""Offline diagnostics for interactions omitted by the nearest-neighbor chain.

Never used to declare a target pass. Known downstream values are held fixed
to check whether input-to-t1 electric feedthrough or mutual L explains bias.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.second_experiment import L1,LR


def network(f,l2,x,kind):
    # x: C02, input Cg, t2-t1 Cg, extra interaction, left ps, right ps.
    w=2*np.pi*f
    y=np.zeros((len(f),5,5),dtype=complex)
    for node,cap in zip((1,2,3),(x[0],101.492048434,188.450588735)):
        y[:,node,node]+=1j*w*cap*1e-15
    inductance=np.diag([l2,L1,LR])*1e-9
    if kind=='mutual':inductance[0,1]=inductance[1,0]=x[3]*1e-9
    y[:,1:4,1:4]+=np.linalg.inv(inductance)[None,:,:]/(1j*w[:,None,None])
    pairs=[(0,1,x[1]),(1,2,x[2]),(2,3,8.229325259),(3,4,8.312405937)]
    if kind=='cross_cap':pairs.append((0,2,x[3]))
    for a,b,c in pairs:
        v=1j*w*c*1e-15
        y[:,a,a]+=v;y[:,b,b]+=v;y[:,a,b]-=v;y[:,b,a]-=v
    ports=[0,4];internal=[1,2,3]
    yp=y[:,ports][:,:,ports]-y[:,ports][:,:,internal]@np.linalg.solve(y[:,internal][:,:,internal],y[:,internal][:,:,ports])
    identity=np.eye(2)[None,:,:]
    s=(identity-50*yp)@np.linalg.inv(identity+50*yp)
    s[:,0,0]*=np.exp(-2j*w*x[4]*1e-12);s[:,1,1]*=np.exp(-2j*w*x[5]*1e-12)
    s[:,0,1]*=np.exp(-1j*w*(x[4]+x[5])*1e-12);s[:,1,0]=s[:,0,1]
    return s


def diagnose(snp,l2,output):
    net=read_snp(snp);idx=np.flatnonzero(net.f<=9e9)[::2]
    f=net.f[idx];obs=net.s[idx];records={}
    for kind in ('cross_cap','mutual'):
        lower=[20,100,10,0 if kind=='cross_cap' else -.5,0,0]
        upper=[250,1200,140,50 if kind=='cross_cap' else .5,20,20]
        def residual(x):
            s=network(f,l2,x,kind);a=3*(s[:,1,0]-obs[:,1,0]);b=s[:,0,0]-obs[:,0,0]
            return np.r_[a.real,a.imag,b.real,b.imag]
        starts=[[100,459,56.7,v,4,.6] for v in ([1,10,25] if kind=='cross_cap' else [-.1,0,.1])]
        results=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=500) for x in starts]
        r=min(results,key=lambda q:q.fun@q.fun);s=network(f,l2,r.x,kind)
        records[kind]={'x':r.x.tolist(),'success':bool(r.success),'cost':float(r.fun@r.fun),
                       's21_rmse':float(np.sqrt(np.mean(abs(s[:,1,0]-obs[:,1,0])**2))),
                       's11_rmse':float(np.sqrt(np.mean(abs(s[:,0,0]-obs[:,0,0])**2)))}
    Path(output).write_text(json.dumps(records,indent=2)+'\n')
    return records


def input_parasitic_response(f,x,l2,l1):
    """Diagnostic chain plus IDC series L and input-to-ground C (x[9:11])."""
    w=2*np.pi*f
    eye=np.tile(np.eye(2,dtype=complex),(len(f),1,1))
    m=eye.copy();m[:,1,0]=1j*w*x[10]*1e-15
    for i,cg in enumerate(x[3:7]):
        a=eye.copy();a[:,0,1]=1/(1j*w*cg*1e-15)
        if i==0:a[:,0,1]+=1j*w*x[9]*1e-9
        m=m@a
        if i<3:
            a=eye.copy();a[:,1,0]=1j*w*x[i]*1e-15+1/(1j*w*[l2,l1,LR][i]*1e-9)
            m=m@a
    if len(x)>11:
        # The native output end still contains its grounded regular meander.
        # It loads the output port even though it is not an interior mode.
        a=eye.copy();a[:,1,0]=1j*w*x[11]*1e-15+1/(1j*w*LR*1e-9)
        m=m@a
    a,b,c,d=m[:,0,0],m[:,0,1],m[:,1,0],m[:,1,1]
    den=a+b/50+c*50+d
    s=np.empty_like(m)
    s[:,0,0]=(a+b/50-c*50-d)/den
    s[:,1,1]=(-a+b/50-c*50+d)/den
    s[:,1,0]=s[:,0,1]=2/den
    s[:,0,0]*=np.exp(-2j*w*x[7]*1e-12)
    s[:,1,1]*=np.exp(-2j*w*x[8]*1e-12)
    s[:,1,0]*=np.exp(-1j*w*(x[7]+x[8])*1e-12)
    s[:,0,1]=s[:,1,0]
    return s


def diagnose_input_parasitics(output):
    """Offline model check only; does not replace target extraction or run EM."""
    output=Path(output);primary=json.loads((output/'fit_network.json').read_text())
    net=read_snp(output/'response.snp');records=[]
    lower=[20,60,150,100,10,2,2,0,0,0,0]
    upper=[300,160,230,1600,140,25,25,20,20,2,200]
    for stop in (9,6):
        ix=np.flatnonzero(net.f<=stop*1e9)[::2];f=net.f[ix];obs=net.s[ix]
        def residual(x):
            s=input_parasitic_response(f,x,primary['L_t2_nH'],primary['L_t1_nH'])
            a=3*(s[:,1,0]-obs[:,1,0]);b=s[:,0,0]-obs[:,0,0]
            return np.r_[a.real,a.imag,b.real,b.imag]
        starts=[np.r_[primary['x'],esl,cap] for esl,cap in ((0,0),(.4,60))]
        results=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=400,ftol=1e-10) for x in starts]
        r=min(results,key=lambda q:q.fun@q.fun)
        s=input_parasitic_response(f,r.x,primary['L_t2_nH'],primary['L_t1_nH'])
        record={'stop_GHz':stop,'x':r.x.tolist(),'success':bool(r.success),'active_bounds':r.active_mask.tolist(),
                's21_rmse':float(np.sqrt(np.mean(abs(s[:,1,0]-obs[:,1,0])**2))),
                's11_rmse':float(np.sqrt(np.mean(abs(s[:,0,0]-obs[:,0,0])**2)))}
        print(json.dumps(record),flush=True);records.append(record)
    (output/'input_parasitic_diagnostic.json').write_text(json.dumps({
        'purpose':'Model inadequacy diagnostic; not a replacement target fit',
        'parameter_order':['C02','C01','Cr','Cin','Cg21','Cg','Cout','delay_left_ps','delay_right_ps','IDC_ESL_nH','input_ground_fF'],
        'fits':records},indent=2)+'\n')
    return records


def diagnose_output_end(output):
    """Check the physically present output-end LC against all measured channels."""
    output=Path(output);primary=json.loads((output/'fit_network.json').read_text())
    net=read_snp(output/'response.snp');records=[]
    lower=[20,60,150,100,10,2,2,0,0,100]
    upper=[300,160,230,1600,140,25,25,20,20,300]
    for stop in (9,6):
        ix=np.flatnonzero(net.f<=stop*1e9)[::2];f=net.f[ix];obs=net.s[ix]
        def response(x):return input_parasitic_response(f,np.r_[x[:9],0,0,x[9]],primary['L_t2_nH'],primary['L_t1_nH'])
        def residual(x):
            s=response(x);a=3*(s[:,1,0]-obs[:,1,0]);b=s[:,0,0]-obs[:,0,0];c=s[:,1,1]-obs[:,1,1]
            return np.r_[a.real,a.imag,b.real,b.imag,c.real,c.imag]
        starts=[np.r_[primary['x'],cap] for cap in (180,220)]
        results=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=400,ftol=1e-10) for x in starts]
        r=min(results,key=lambda q:q.fun@q.fun);s=response(r.x)
        record={'stop_GHz':stop,'x':r.x.tolist(),'success':bool(r.success),'active_bounds':r.active_mask.tolist(),
                's21_rmse':float(np.sqrt(np.mean(abs(s[:,1,0]-obs[:,1,0])**2))),
                's11_rmse':float(np.sqrt(np.mean(abs(s[:,0,0]-obs[:,0,0])**2))),
                's22_rmse':float(np.sqrt(np.mean(abs(s[:,1,1]-obs[:,1,1])**2)))}
        print(json.dumps(record),flush=True);records.append(record)
    (output/'output_end_diagnostic.json').write_text(json.dumps({
        'purpose':'Physical output fixture model check; not a replacement target fit',
        'parameter_order':['C02','C01','Cr','Cin','Cg21','Cg','Cout','delay_left_ps','delay_right_ps','Cend_fF'],
        'fits':records},indent=2)+'\n')
    return records


def isolation_response(f,x,l_nh,cadd_ff):
    """Single LC with unequal feed couplings, port stray C and matched delays."""
    w=2*np.pi*f;eye=np.tile(np.eye(2,dtype=complex),(len(f),1,1))
    m=eye.copy();m[:,1,0]=1j*w*x[3]*1e-15
    a=eye.copy();a[:,0,1]=1/(1j*w*x[1]*1e-15);m=m@a
    a=eye.copy();a[:,1,0]=1j*w*(x[0]+cadd_ff)*1e-15+1/(1j*w*l_nh*1e-9);m=m@a
    a=eye.copy();a[:,0,1]=1/(1j*w*x[2]*1e-15);m=m@a
    a=eye.copy();a[:,1,0]=1j*w*x[4]*1e-15;m=m@a
    a,b,c,d=m[:,0,0],m[:,0,1],m[:,1,0],m[:,1,1];den=a+b/50+c*50+d
    s=np.empty_like(m);s[:,0,0]=(a+b/50-c*50-d)/den;s[:,1,1]=(-a+b/50-c*50+d)/den
    s[:,1,0]=s[:,0,1]=2/den
    s[:,0,0]*=np.exp(-2j*w*x[5]*1e-12);s[:,1,1]*=np.exp(-2j*w*x[6]*1e-12)
    s[:,1,0]*=np.exp(-1j*w*(x[5]+x[6])*1e-12);s[:,0,1]=s[:,1,0]
    return s


def diagnose_isolation(output):
    """Constrain C0 from all four existing C_add spectra; no EM or target refit."""
    output=Path(output);record=json.loads((output/'results.json').read_text());data=[]
    for point in record['sweep']:
        cap=point['C_add_fF'];net=read_snp(output/f'cadd_{cap:g}ff.snp');peak=point['frequency_Hz']
        ix=np.flatnonzero(abs(net.f-peak)<1e9)[::2]
        data.append((cap,net.f[ix],net.s[ix]))
    records=[]
    for count in (4,3):
        subset=data[-count:]
        def residual(x):
            errors=[]
            for cap,f,obs in subset:
                delta=isolation_response(f,x,record['L0_nH'],cap)-obs
                for i,j in ((1,0),(0,0),(1,1)):
                    errors.extend([delta[:,i,j].real,delta[:,i,j].imag])
            return np.concatenate(errors)
        lower=[1,.1,1,0,0,0,0];upper=[300,200,300,200,200,20,20]
        starts=[[c0,20,100,40,50,1,1] for c0 in (80,140)]
        results=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=350,ftol=1e-10) for x in starts]
        r=min(results,key=lambda q:q.fun@q.fun)
        record_fit={'number_of_loads':count,'x':r.x.tolist(),'success':bool(r.success),
                    'active_bounds':r.active_mask.tolist(),'real_residual_rmse':float(np.sqrt(np.mean(r.fun**2)))}
        print(json.dumps(record_fit),flush=True);records.append(record_fit)
    (output/'isolation_capacitance_diagnostic.json').write_text(json.dumps({
        'purpose':'Independent self-capacitance model diagnostic from archived sweeps; not an accepted extraction',
        'L_nH':record['L0_nH'],
        'parameter_order':['C0_fF','Cleft_fF','Cright_fF','C_port_left_fF','C_port_right_fF','delay_left_ps','delay_right_ps'],
        'fits':records},indent=2)+'\n')
    return records
