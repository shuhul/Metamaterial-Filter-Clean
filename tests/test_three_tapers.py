"""Physical invariants for the new regular-cell-based T1 fixture."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import gdstk
from geometry.taper.three_tapers import BASE,build,validate,from_dict


class NewTaperTests(unittest.TestCase):
    def test_shared_ground_and_flat_top(self):
        for isolation in (False,True):
            cell,meta=build(BASE,isolation=isolation)
            polys=cell.get_polygons()
            self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),1)
            flat=gdstk.rectangle((0,0),(meta['w_sim'],448))
            self.assertFalse(gdstk.boolean([flat],polys,'not',precision=.001))
            self.assertEqual(BASE.d_top+BASE.cap.t_top+BASE.meander.l,618)

    def test_metadata_rebuild_validates_written_gds(self):
        cell,meta=build(BASE)
        self.assertEqual(from_dict(json.loads(json.dumps(meta['unit']))),BASE)
        with TemporaryDirectory() as temp:
            path=Path(temp)/'test.gds'
            lib=gdstk.Library(unit=1e-6,precision=1e-9)
            lib.add(cell,*cell.dependencies(True));lib.write_gds(str(path))
            result,_=validate(path,meta)
            self.assertEqual(result.n_cells,2)


if __name__=='__main__':unittest.main()
