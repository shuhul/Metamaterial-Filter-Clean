"""Recover complete four-port direct records from a mixed legacy Sonnet export.

Changing component terminals can leave old two-port cache records in a native
export. Never flatten or reinterpret those rows as part of a four-port matrix.
Only explicitly requested complete 4x4 MA records are accepted.
"""
import numpy as np
import skrf


def four_port_direct(path,frequencies_GHz):
    lines=path.read_text().splitlines()
    headers=[q.upper().split() for q in lines if q.startswith('#')]
    if headers!=[['#','GHZ','S','MA','R','50.00000']]:
        raise ValueError('Unexpected native Touchstone units/normalization')
    blocks=[];current=[]
    for line in lines:
        if not line.strip() or line.startswith(('!','#')):continue
        values=list(map(float,line.split()))
        if not line[0].isspace():
            if current:blocks.append(current)
            current=values
        else:current.extend(values)
    if current:blocks.append(current)
    rows={}
    for values in blocks:
        if len(values)!=33:continue
        frequency=values[0]
        if frequency not in frequencies_GHz:continue
        if frequency in rows:raise ValueError('Duplicate direct frequency')
        ma=np.asarray(values[1:]).reshape(4,4,2)
        rows[frequency]=ma[:,:,0]*np.exp(1j*np.pi/180*ma[:,:,1])
    if set(rows)!=set(frequencies_GHz):raise ValueError('Missing complete requested four-port records')
    freq=sorted(rows)
    return skrf.Network(frequency=skrf.Frequency.from_f(np.array(freq)*1e9,unit='hz'),s=np.array([rows[f] for f in freq]),z0=50)
