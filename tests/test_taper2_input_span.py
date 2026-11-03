"""An independent feed span must preserve the entire IDC and downstream metal."""
import json
import unittest
from dataclasses import replace
import gdstk
from geometry.taper import taper2_plain as w

BASE='idc_rail8_dl256_pad148_w148_p200_i42n70_j42_r74tg2'


class InputSpanTests(unittest.TestCase):
    def settings(self):
        path=w.ROOT/'sonnet_files/taper_exploration'/w.CAMPAIGN/BASE/'response/geometry.json'
        saved=json.loads(path.read_text())
        return w.from_dict(saved['trial_settings']),saved,path

    def test_default_still_reproduces_saved_geometry_and_metadata(self):
        s,saved,path=self.settings();top,meta=w.build(s)
        self.assertEqual(meta,saved)
        original=gdstk.read_gds(str(path.with_suffix('.gds'))).top_level()[0]
        self.assertFalse(gdstk.boolean(top.get_polygons(),original.get_polygons(),'xor',precision=.001))

    def test_shorter_feed_only_translates_idc_and_downstream(self):
        s,_,_=self.settings();old,a=w.build(s)
        new,b=w.build(replace(s,coupling=replace(s.coupling,w_input_stub=344)))
        shift=a['cell_lefts'][1]-b['cell_lefts'][1]
        self.assertEqual(shift,196)
        self.assertEqual(a['unit'],b['unit'])
        self.assertEqual(a['IDC_clearance_um'],b['IDC_clearance_um'])
        start=a['cell_lefts'][1]-s.coupling.d-s.coupling.t_input
        clip=gdstk.rectangle((start,-1052),(a['w_sim'],448))
        p=gdstk.boolean(old.get_polygons(),[clip],'and',precision=.001)
        q=gdstk.boolean([p.copy().translate(shift,0) for p in new.get_polygons()],[clip],'and',precision=.001)
        self.assertFalse(gdstk.boolean(p,q,'xor',precision=.001))

    def test_invalid_span_rejected(self):
        s,_,_=self.settings()
        for width in (3,6):
            with self.subTest(width=width),self.assertRaises(ValueError):
                w.build(replace(s,coupling=replace(s.coupling,w_input_stub=width)))


if __name__=='__main__':unittest.main()
