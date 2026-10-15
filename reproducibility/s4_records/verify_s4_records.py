"""Verify the frozen S4 records without requiring any crop/state NPZ.

python verify_s4_records.py
python verify_s4_records.py --project-root <source-project>

The optional mode checks original context, witness and stored-state hashes.
It reads inputs only. A failed code hash, boundary, Betti or SAFE check is fatal.
"""
from pathlib import Path
import argparse,json,hashlib
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def rank(cols):
    piv={}
    for c in cols:
        while c:
            k=c.bit_length()-1
            if k not in piv:piv[k]=c;break
            c^=piv[k]
    return len(piv)
parser=argparse.ArgumentParser();parser.add_argument('--project-root',type=Path);args=parser.parse_args()
e=json.loads((ROOT/'ATTACHMENT_RECALCULATION.json').read_text(encoding='utf-8'))
ledger=json.loads((ROOT/'ATTACHMENT_CELL_LEDGER.json').read_text(encoding='utf-8'))
assert sha(ROOT/'close_evidence_FROZEN.py')==e['matrix_code_sha256']
for cells,state in zip(ledger,e['states']):
    n=[len(s) for s in cells['cells_by_dimension']];d1,d2=cells['boundary_columns_f2']
    for c in d2:
        composed=0
        for i,d in enumerate(d1):
            if c>>i&1:composed^=d
        assert composed==0
    r1,r2=rank(d1),rank(d2)
    beta=[n[0]-r1,n[1]-r1-r2,n[2]-r2]
    assert beta==state['beta']
    # Frozen S4 rule: SAFE iff the attachment Betti vector is (1,0,0).
    # Use the recomputed vector, rather than trusting either stored field.
    expected_safe=beta==[1,0,0]
    if type(state['safe']) is not bool or state['safe']!=expected_safe:
        raise ValueError(f"SAFE mismatch for {state['state']}: "
                         f"stored={state['safe']!r}, recomputed={expected_safe}, beta={beta}")
    assert [r1,r2]==state['boundary_ranks_f2']
    print(state['state'],state['beta'],'SAFE='+str(state['safe']))
if args.project_root:
    import numpy as np
    p=args.project_root
    assert sha(p/e['ctx_file'])==e['ctx_sha256']
    assert sha(p/e['bundle_file'])==e['bundle_sha256']
    with np.load(p/e['bundle_file'],allow_pickle=False) as b:
        for i,state in enumerate(e['states']):
            h=hashlib.sha256(np.ascontiguousarray(b['wit_states'][i],dtype=np.uint8).tobytes()).hexdigest()
            assert h==state['witness_state_sha256']
    print('Original source-state hashes verified.')
print('All three frozen S4 results verified.')
