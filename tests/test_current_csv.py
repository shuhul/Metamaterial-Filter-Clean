"""Regression for Sonnet's scalar-charge export with a misleading direction label."""
import csv
from pathlib import Path
import tempfile
import unittest
from geometry.taper.current_csv import read


class CurrentCsvTest(unittest.TestCase):
    def test_charge_is_scalar_even_when_section_says_y_directed(self):
        rows=[['VER : 2'],['soncmd'],['Frequency:','4.8441e9'],['Drive:'],
              ['theLevel:','0'],['Export Positions in:','UM'],['Export Steps:','4'],
              ['Surface Charge','Complex Form','Coulomb/sq m'],['Y Directed'],
              ['Y Position'],['X Position ->','2','6'],
              ['1498',' 1.2e-9 +3.4e-9j',' -5e-9 -6e-9j']]
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'charge.csv'
            with path.open('w',newline='') as handle:csv.writer(handle).writerows(rows)
            d=read(path)
        self.assertNotIn('jy',d)
        self.assertEqual(d['scalar'].shape,(1,2))
        self.assertEqual(d['scalar'][0,0],1.2e-9+3.4e-9j)
        self.assertEqual(d['scalar'][0,1],-5e-9-6e-9j)
        self.assertEqual(d['y_export'].tolist(),[1498])


if __name__=='__main__':unittest.main()
