"""Verify release file bytes, frozen sources and complete license coverage.

Run with normal Python from an extracted versioned release. No third-party
dependency or BraTS data is required; no search is executed.
"""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent
def main():
    entries=[]
    for line in (ROOT/'MANIFEST.sha256').read_text(encoding='utf-8').splitlines():
        digest,relative=line.split('  ',1)
        p=(ROOT/relative).resolve()
        if not p.is_relative_to(ROOT):raise ValueError('Manifest escapes source root')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=digest:raise ValueError('Digest mismatch: '+relative)
        entries.append(relative)
    pins=json.loads((ROOT/'full_search/FROZEN_SOURCE_SHA256.json').read_text())
    for relative,digest in pins.items():
        p=ROOT/'full_search/frozen'/relative
        if hashlib.sha256(p.read_bytes()).hexdigest()!=digest:raise ValueError('Frozen mismatch: '+relative)
    scope=json.loads((ROOT/'LICENSE_SCOPE.json').read_text())
    for relative in entries+['MANIFEST.sha256']:
        if relative not in scope:raise ValueError('Missing license scope: '+relative)
    print(f'RELEASE_INTEGRITY_PASS: {len(entries)} files; {len(pins)} frozen entries; license scope complete')
if __name__=='__main__':
    try:main()
    except Exception as exc:
        print(f'RELEASE_INTEGRITY_FAIL: {exc}',file=sys.stderr)
        raise SystemExit(1) from exc
