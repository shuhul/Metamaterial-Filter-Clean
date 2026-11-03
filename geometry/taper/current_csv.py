"""Read Sonnet current exports without assuming coordinate labels are native."""
import csv
import numpy as np


def read(path):
    with open(path) as handle:rows=list(csv.reader(handle))
    result=dict(frequency_Hz=float(rows[2][1]),measurement=rows[7])
    section='scalar';data={};x=None
    for row in rows[8:]:
        if not row:continue
        first=row[0].strip()
        if first in ('X Directed','Y Directed'):
            # Sonnet labels scalar charge as "Y Directed" too. The quantity
            # header, not that legacy section name, determines its meaning.
            section='scalar' if 'Surface Charge' in result['measurement'] else ('jx' if first=='X Directed' else 'jy');continue
        if first=='X Position ->':
            x=np.array([float(v) for v in row[1:] if v.strip()]);continue
        if x is None:continue
        try:y=float(first)
        except ValueError:continue
        values=[complex(v.replace(' ','')) for v in row[1:] if v.strip()]
        if len(values)!=len(x):raise ValueError('CSV axis/data size mismatch')
        data.setdefault(section,[]).append((y,values))
    for kind,items in data.items():
        result[kind]=np.array([v for _,v in items]);result['y_export']=np.array([y for y,_ in items])
    result['x']=x
    return result
