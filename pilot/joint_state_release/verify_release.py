"""Check release provenance and exact scalar display values."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parent
manifest = json.loads((root / 'provenance.json').read_text())
for path, expected in manifest.items():
    actual = hashlib.sha256((root / 'source/titania-dielectric-redox' / path).read_bytes()).hexdigest()
    assert actual == expected, path
page = json.loads((root / 'page_data.json').read_text())
for key, name in [('reduction', 'reduction-o1-6-i2-2-z30-n32-119902.json'),
                  ('reoxidation', 'reox-o1.6-i2.2-z30-n320-r2.5-x0.05.json')]:
    raw = json.loads((root / 'results' / name).read_text())
    assert len(page[key]) == len(raw['checkpoints'])
    for shown, full in zip(page[key], raw['checkpoints']):
        for field, value in shown.items():
            if not isinstance(value, list):
                assert value == full[field], (key, field)
        if 'depth_nm' in shown:
            assert len(shown['depth_nm']) <= 120
            assert shown['depth_nm'][:12] == full['depth_nm'][:12]
            lookup = {depth: index for index, depth in enumerate(full['depth_nm'])}
            for field, value in shown.items():
                if isinstance(value, list):
                    assert value == [full[field][lookup[z]] for z in shown['depth_nm']], (key, field)
print(f'PASS: {len(manifest)} source hashes; exact scalar checkpoints; exact sampled profiles.')
