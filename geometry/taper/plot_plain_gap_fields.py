"""Validate and plot saved baseline fields; never starts an EM simulation."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import LineCollection
from geometry.prepare_sonnet import tag, points
from geometry.taper.plain_gap_study import NATIVE, OUT
from geometry.taper.current_csv import read


def native_polygons():
    root=ET.parse(NATIVE/'baseline_currents/generated_ncell_3.sonx').getroot()
    return [points(p) for p in root.iter(tag('PlanarPolygon'))]


def metal_mask(d,polygons):
    # Whole/Line exports use the project's displayed coordinates. Native XML
    # has LocalOrigin Y=1500 and YDirection=Negative: native Y=1500-export Y.
    xx,yy=np.meshgrid(d['x'],1500-d['y_export'])
    coords=np.column_stack([xx.ravel(),yy.ravel()])
    mask=np.array(gdstk.inside(coords,[gdstk.Polygon(p) for p in polygons]))
    return mask.reshape(xx.shape)


def main():
    folder=OUT/'baseline_currents';export=folder/'exports';polys=native_polygons()
    labels=['lower','middle','upper'];currents=[];charges=[];checks=[];leads=[]
    for label in labels:
        j=read(export/f'jxy_whole_{label}.csv');q=read(export/f'q_whole_{label}.csv')
        if 'Coulomb' not in ','.join(q['measurement']):raise ValueError('Not charge data')
        mask=metal_mask(j,polys);amp=np.sqrt(abs(j['jx'])**2+abs(j['jy'])**2)
        active=amp>amp.max()*1e-10
        invalid=int(np.count_nonzero(active&~mask))
        if invalid:raise ValueError(f'{invalid} nonzero current samples outside actual native metal')
        if not np.array_equal(j['x'],q['x']) or not np.array_equal(j['y_export'],q['y_export']):
            raise ValueError('Charge/current axes differ')
        qa=abs(q['scalar']);qactive=qa>qa.max()*1e-10
        qinvalid=int(np.count_nonzero(qactive&~mask))
        if qinvalid:raise ValueError(f'{qinvalid} nonzero charge samples outside actual native metal')
        currents.append((j,np.ma.masked_where(~mask|(amp==0),amp)))
        charges.append((q,np.ma.masked_where(~mask|(qa==0),qa)))
        checks.append(dict(label=label,frequency_GHz=j['frequency_Hz']/1e9,
            active_current_samples=int(active.sum()),current_outside_metal=invalid,
            active_charge_samples=int(qactive.sum()),charge_outside_metal=qinvalid))
        # Whole-export pixel centers land inside both 4 um vertical necks.
        # Compare the same Jy component and direction; this avoids assigning
        # a meaningless combined phase to the 2-D current vector.
        iy=np.flatnonzero(j['y_export']==1026)[0];samples={}
        for name,x in [('T1',1334),('unit',1706)]:
            ix=np.flatnonzero(j['x']==x)[0]
            if not mask[iy,ix]:raise ValueError('Phase sample is not on metal')
            value=j['jy'][iy,ix]
            samples[name]=dict(x_um=x,native_y_um=474,Jy_magnitude_Am=float(abs(value)),
                               phase_deg=float(np.angle(value,deg=True)))
        a=j['jy'][iy,np.flatnonzero(j['x']==1334)[0]]
        b=j['jy'][iy,np.flatnonzero(j['x']==1706)[0]]
        leads.append(dict(label=label,samples=samples,T1_relative_unit_phase_deg=float(np.angle(a/b,deg=True))))
    edges=[np.vstack([p,p[0]])*np.array([1,-1])+np.array([0,448]) for p in polys]
    fig,axs=plt.subplots(2,3,figsize=(20,7.2),sharex=True,sharey=True)
    fig.subplots_adjust(left=.07,right=.90,top=.87,bottom=.10,wspace=.13,hspace=.35)
    for r,(data,title,units) in enumerate([(currents,'Surface Current Magnitude','A/m'),
                                          (charges,'Surface Charge Magnitude','C/m²')]):
        maximum=max(float(a.max()) for _,a in data)
        norm=LogNorm(vmin=maximum*1e-4,vmax=maximum)
        for col,(d,a) in enumerate(data):
            ax=axs[r,col];y=d['y_export']-1052 # design Y = 448 - native Y
            mesh=ax.pcolormesh(d['x'],y,a,shading='nearest',cmap='inferno',norm=norm,rasterized=True)
            ax.add_collection(LineCollection(edges,colors='#8b9a9f',linewidths=.25,alpha=.6))
            ax.set_xlim(480,2250);ax.set_ylim(-640,36);ax.set_aspect('equal')
            ax.set_title(f"{['Lower Peak','Between Peaks','Upper Peak'][col]} · {d['frequency_Hz']/1e9:.4f} GHz",fontsize=15)
            if r==1:ax.set_xlabel('Horizontal Position (µm)',fontsize=12)
            if col==0:ax.set_ylabel(title+'\nVertical Position (µm)',fontsize=12)
            for x,name in [(966,'T2'),(1334,'T1'),(1706,'Unit'),(2064,'End')]:
                ax.text(x,-630,name,ha='center',fontsize=11,color='#32434b',
                        bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.5))
        cax=fig.add_axes([.92,.55 if r==0 else .13,.009,.27])
        fig.colorbar(mesh,cax=cax,label=f'{title} ({units})')
    fig.suptitle('Approved Layout: Current and Charge at the Two Resonances',fontsize=22)
    fig.savefig(folder/'current_charge_maps.png',dpi=170);fig.savefig(folder/'current_charge_maps.pdf');plt.close(fig)
    # Jy phase is a component phase, not a phase invented by combining Jx/Jy.
    fig,axs=plt.subplots(1,3,figsize=(20,5),sharex=True,sharey=True)
    for ax,(d,a),label in zip(axs,currents,labels):
        jy=d['jy'];phase=np.ma.masked_where(a.mask|(abs(jy)<abs(jy).max()*1e-3),np.angle(jy,deg=True))
        mesh=ax.pcolormesh(d['x'],d['y_export']-1052,phase,shading='nearest',cmap='twilight',vmin=-180,vmax=180)
        ax.add_collection(LineCollection(edges,colors='gray',linewidths=.25,alpha=.6))
        ax.set_xlim(480,2250);ax.set_ylim(-640,36);ax.set_aspect('equal');ax.set_title(label.capitalize())
        ax.set_xlabel('Horizontal Position (µm)')
    axs[0].set_ylabel('Vertical Position (µm)')
    fig.colorbar(mesh,ax=axs.tolist(),label='Jy phase (degrees)',shrink=.75)
    fig.suptitle('Vertical Surface-Current Phase · Same Port-1 Excitation',fontsize=20)
    fig.savefig(folder/'vertical_current_phase.png',dpi=170);plt.close(fig)
    report=dict(coordinate_mapping='native Y = 1500 - exported Y; design Y = exported Y - 1052',
        validation=checks,neck_current_phase_samples=leads,
        notes=['1 V / 50 ohm source at port 1, 50 ohm termination at port 2.',
               'Fields are not deembedded. No quantitative branch-current extraction: line exports failed metal-support validation.',
               'Common magnitude color scales across frequencies; four decades shown.',
               'Current/charge localization indicates mode participation, not unique lumped capacitances or loss.'])
    (folder/'field_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    for path in export.glob('export_*.json'):
        manifest=json.loads(path.read_text())
        valid=path.stem in ('export_jxy_whole','export_q_whole')
        manifest['accepted_for_field_analysis']=valid
        manifest['coordinate_validation']='field_validation.json (parent folder)' if valid else 'Rejected; see REJECTED_EXPORTS.md'
        path.write_text(json.dumps(manifest,indent=2)+'\n')
    (export/'REJECTED_EXPORTS.md').write_text('# Export validation\n\n'
        'Use `jxy_whole_*.csv` and `q_whole_*.csv` only.\n\n'
        'The older `charge_*.csv` are actually current: unknown Type="charge" silently falls back to JXY. '
        'Surface charge requires Type="q" and a Coulomb/sq m header.\n\n'
        'The older rectangular `jxy_*.csv` and `q_*.csv` have inconsistent Y coordinates and must not be interpreted. '
        'Line exports `jxy_leads_*.csv` also failed metal-support validation and are rejected. '
        'Whole-layout exports are validated against actual native metal support in field_validation.json.\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
