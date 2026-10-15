"""Check shipped files against FILE_SHA256.json; no scientific computation."""
from pathlib import Path
import json,hashlib,sys
root=Path(__file__).resolve().parent
try:
    hashes=json.loads((root/'FILE_SHA256.json').read_text(encoding='utf-8'))
    for name,expected in hashes.items():
        path=(root/name).resolve()
        if not path.is_relative_to(root):raise ValueError('Manifest path escapes package')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('Digest mismatch: '+name)
    print('PACKAGE_SHA256_PASS: '+str(len(hashes))+' files')
except (OSError,ValueError,TypeError) as error:
    print('PACKAGE_SHA256_FAIL: '+str(error),file=sys.stderr)
    sys.exit(1)
