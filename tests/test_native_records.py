import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from geometry.taper.native_records import four_port_direct


class NativeRecordTests(unittest.TestCase):
    def test_mixed_two_port_rows_cannot_pollute_four_port_matrix(self):
        with TemporaryDirectory() as folder:
            path=Path(folder)/'mixed.snp'
            path.write_text('# GHZ S MA R 50.00000\n4.65 1 0 0 0 0 0 1 0\n'+
                '5 '+ ' '.join(['1 0']*4)+'\n'+('  '+' '.join(['1 0']*4)+'\n')*3+
                '5.4 1 0 0 0 0 0 1 0\n')
            net=four_port_direct(path,[5])
            np.testing.assert_allclose(net.s,np.ones((1,4,4)))
            with self.assertRaises(ValueError):four_port_direct(path,[4.65,5])


if __name__=='__main__':unittest.main()
