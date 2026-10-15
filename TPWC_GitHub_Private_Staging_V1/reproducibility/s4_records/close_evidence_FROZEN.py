"""Close the approved Narrative candidate without editing frozen sources.

Run: bundled-python -X utf8 code/close_evidence.py
Reads frozen context, stored witness states and production path rows. Replays
the path, computes cubical boundary ranks over F2, cross-checks the existing
probe, and appends Table S4 to a separate OR3 candidate. Assertions and explicit
input errors stop the build. No new search or theoretical result is introduced.
"""
from pathlib import Path
import csv, difflib, hashlib, itertools, json, re, shutil, subprocess, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if (ROOT/'FREEZE_MANIFEST.json').exists():
    raise RuntimeError('This version is frozen. Reproduce or revise in a new version directory.')
FC=ROOT.parent;PHASE=FC.parent;PROJECT=PHASE.parents[1]
NARR=FC/'R-Fig3_4_R4I_NARRATIVE'
SRC=NARR/'MANUSCRIPT_CANDIDATE'
CLOSURE=FC/'R-Fig3_4_R4H_closure/MANUSCRIPT_AND_OR_R4H/OR_candidates'
FROZEN=PROJECT/'MANUSCRIPT/PHASE_Q/FINAL'
MAN=ROOT/'MANUSCRIPT';OR=ROOT/'ONLINE_RESOURCES';E=ROOT/'evidence';D=ROOT/'diffs'
for p in (MAN,OR,E,D):p.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

# Record all approved candidate files, frozen Q assets and production files.
# All writes below are inside ROOT. Input manifests also cover the witness data.
protected=[p for folder in (NARR,FROZEN,PROJECT/'PRODUCTION',CLOSURE) for p in folder.rglob('*') if p.is_file()]
fd=PROJECT/'FIGURE_DATA/shard_1/FIG/BraTS2021_00110__k2xy_106_111_162'
protected += [fd/'ctx.npz',fd/'C00-F1__G26.npz',PHASE/'code/attach_along_witness.py',PHASE/'logs/recomputed_table_S2_effective.csv']
before={str(p.relative_to(PROJECT)):sha(p) for p in protected}
(E/'INPUT_SHA256.json').write_text(json.dumps(before,indent=2),encoding='utf-8')

try:
    ctx=np.load(fd/'ctx.npz',allow_pickle=False)
    bundle=np.load(fd/'C00-F1__G26.npz',allow_pickle=False)
    lo=ctx['crop_lo'];base=ctx['xt'].astype(np.uint8)
    coords=bundle['coords'];stored=bundle['wit_states']
    assert stored.shape==(5,180)
    assert bundle['status'].item()=='PURE_FOUND'
    idx=tuple((coords-lo).T)
    assert np.array_equal(base[idx],stored[0])
except (OSError,KeyError,ValueError) as exc:raise RuntimeError(f'Invalid stored witness bundle: {exc}') from exc
pathrows=[]
for f in sorted((PROJECT/'PRODUCTION').glob('shard_*/PATHS.csv')):
    with f.open(encoding='utf-8-sig',newline='') as stream:
        for line,r in enumerate(csv.DictReader(stream),start=2):
            if (r['system'],r['grid'],r['case_id'],r['block'])==('C00-F1','G26','BraTS2021_00110','k2xy@106,111,162'):
                pathrows.append(dict(r,source=str(f.relative_to(PROJECT)),csv_line=line))
pathrows.sort(key=lambda r:int(r['step']))
assert len(pathrows)==4 and [int(r['step']) for r in pathrows]==list(range(4))
crop=base.copy();crops=[crop.copy()]
for k,r in enumerate(pathrows):
    v=tuple(int(r[a]) for a in ('z','y','x'));isrem=r['type']=='REM'
    assert np.array_equal(coords[bundle['wit_vidx'][k]],v)
    assert bool(bundle['wit_is_rem'][k])==isrem
    assert bundle['wit_path'][k].tolist()==[int(isrem),*v]
    local=tuple(np.array(v)-lo)
    assert int(crop[local])==int(isrem),(k,'invalid flip occupancy')
    crop[local]=1-crop[local]
    assert np.array_equal(crop[idx],stored[k+1]),(k,'stored witness and replay disagree')
    assert r['safe']=='1' and bool(bundle['wit_flags'][k,0])
    crops.append(crop.copy())
(E/'ORIGINAL_WITNESS_PATH_ROWS.json').write_text(json.dumps(pathrows,indent=2),encoding='utf-8')

OFF=[o for o in itertools.product((-1,0,1),repeat=3) if o!=(0,0,0)]
TARGET=np.array((106,111,162))

def cells_from_neighbours(crop,v):
    """Build A as the union of closed contact cells, including their faces.

    Cell coordinates are local to Cv=[0,1]^3 in z,y,x order. Each factor is
    (0,0), (1,1) or (0,1). The crop contains the complete 26-neighbourhood.
    """
    local=v-lo
    assert np.all(local>=1) and np.all(local<np.array(crop.shape)-1)
    occupancy=[int(crop[tuple(local+o)]) for o in OFF]
    closed=set()
    for o,on in zip(OFF,occupancy):
        if not on:continue
        f=tuple((0,0) if a==-1 else (1,1) if a==1 else (0,1) for a in o)
        closed.update(itertools.product(*[[(a,b)] if a==b else [(a,a),(b,b),(a,b)] for a,b in f]))
    layers=[sorted(c for c in closed if sum(a!=b for a,b in c)==d) for d in range(3)]
    return occupancy,layers

def boundary(c):
    """Unoriented cellular boundary, appropriate over F2."""
    out=[]
    for axis,(a,b) in enumerate(c):
        if a!=b:
            for end in (a,b):
                face=list(c);face[axis]=(end,end);out.append(tuple(face))
    return out

def columns(lower,upper):
    lookup={c:i for i,c in enumerate(lower)}
    return [sum(1<<lookup[f] for f in boundary(c)) for c in upper]

def rank_f2(cols):
    """Sparse Gaussian elimination with Python integer bit vectors."""
    piv={}
    for col in cols:
        while col:
            pivot=col.bit_length()-1
            if pivot not in piv:piv[pivot]=col;break
            col^=piv[pivot]
    return len(piv)

records=[];cell_ledgers=[]
for k in range(3):
    occ,layers=cells_from_neighbours(crops[k],TARGET)
    b1,b2=columns(layers[0],layers[1]),columns(layers[1],layers[2])
    for col in b2:
        composed=0
        for j,b in enumerate(b1):
            if col>>j&1:composed^=b
        assert composed==0,'Boundary composition is nonzero'
    n=[len(c) for c in layers];r1,r2=rank_f2(b1),rank_f2(b2)
    beta=[n[0]-r1,n[1]-r1-r2,n[2]-r2]
    assert beta==([1,0,1] if k==0 else [1,0,0]),(k,beta)
    statehash=hashlib.sha256(np.ascontiguousarray(stored[k],dtype=np.uint8).tobytes()).hexdigest()
    crophash=hashlib.sha256(np.ascontiguousarray(crops[k],dtype=np.uint8).tobytes()).hexdigest()
    records.append(dict(state=f'Y{k}',bundle_array=f'wit_states[{k},:]',completed_moves=k,last_move='initial' if k==0 else 'REM '+str(k),
        target_zyx=TARGET.tolist(),beta=beta,safe=beta==[1,0,0],cells=n,boundary_ranks_f2=[r1,r2],
        witness_state_sha256=statehash,crop_sha256=crophash,neighbourhood_sha256=hashlib.sha256(bytes(occ)).hexdigest(),
        occupied_neighbour_count=sum(occ),target_occupancy=int(crops[k][tuple(TARGET-lo)])))
    cell_ledgers.append(dict(state=f'Y{k}',offset_order=OFF,occupancy=occ,cells_by_dimension=layers,boundary_columns_f2=[b1,b2]))

# Re-run the frozen probe; check it agrees with direct matrix ranks.
probe=PHASE/'code/attach_along_witness.py'
proc=subprocess.run([sys.executable,'-X','utf8',str(probe),str(PROJECT)],capture_output=True,text=True,encoding='utf-8')
(E/'FROZEN_PROBE_REPLAY.log').write_text(proc.stdout+proc.stderr,encoding='utf-8')
if proc.returncode:raise RuntimeError('Frozen probe failed')
assert "(106, 111, 162) ['(1,0,1)', '(1,0,0) SAFE', '(1,0,0) SAFE'" in proc.stdout
np.savez_compressed(E/'REPLAYED_FIRST_THREE_STATES.npz',crop_lo=lo,coords=coords,witness_states=stored[:3],crop_states=np.stack(crops[:3]),target_zyx=TARGET)
(E/'ATTACHMENT_CELL_LEDGER.json').write_text(json.dumps(cell_ledgers,indent=2),encoding='utf-8')
identifiers=dict(case='BraTS2021_00110',system='C00-F1',grid='G26',block='k2xy@106,111,162',coordinate_order='z,y,x',
    crop_lo=lo.tolist(),crop_shape=list(base.shape),window_voxels=180,
    bundle_file=str((fd/'C00-F1__G26.npz').relative_to(PROJECT)),bundle_sha256=sha(fd/'C00-F1__G26.npz'),
    ctx_file=str((fd/'ctx.npz').relative_to(PROJECT)),ctx_sha256=sha(fd/'ctx.npz'),
    coords_sha256=hashlib.sha256(np.ascontiguousarray(coords,dtype='<i8').tobytes()).hexdigest(),
    state_hash_encoding='C-order contiguous uint8 bytes, 180 entries in stored coords order; no header',
    matrix_code='code/close_evidence.py::cells_from_neighbours/columns/rank_f2',matrix_code_sha256=sha(Path(__file__)),
    frozen_probe=str(probe.relative_to(PROJECT))+'::attach_betti',frozen_probe_sha256=sha(probe),
    stored_states_match_replay=True,existing_probe_agrees=True,boundary_squared_zero=True,states=records)
(E/'ATTACHMENT_RECALCULATION.json').write_text(json.dumps(identifiers,indent=2),encoding='utf-8')
with (E/'ATTACHMENT_RECALCULATION.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f);w.writerow(['state','stored_key','completed_moves','beta0','beta1','beta2','SAFE','n0','n1','n2','rank_d1','rank_d2','state_sha256'])
    for r in records:w.writerow([r['state'],r['bundle_array'],r['completed_moves'],*r['beta'],int(r['safe']),*r['cells'],*r['boundary_ranks_f2'],r['witness_state_sha256']])

# Appendix is added after all original tables and Figure S1. Only the table
# range in the title changes; original prose, tables and figures stay intact.
texrows=[]
for r in records:
    label=r'$Y_0=X_T$' if r['completed_moves']==0 else r'$Y_'+str(r['completed_moves'])+r'$'
    texrows.append(label+' & '+r['last_move']+' & '+r'$('+','.join(map(str,r['beta']))+r')$'+' & '+('yes' if r['safe'] else 'no')+' & '+','.join(map(str,r['cells']))+' & '+','.join(map(str,r['boundary_ranks_f2']))+r' \\')
appendix=r'''
\clearpage
\begingroup
\small
\section*{Supplementary Table S4. Attachment of the filling voxel along the witness}
\addcontentsline{toc}{section}{Supplementary Table S4. Attachment of the filling voxel along the witness}
\label{tab-S4-attachment}

For BraTS2021\_00110, block \texttt{k2xy@106,111,162}, C00-F1 at G26,
let \(v=(106,111,162)\), with coordinates in \((z,y,x)\) order.
The states \(Y_0=X_T\), \(Y_1\) and \(Y_2\) are the initial state and the states
after the first and second removals of the four-flip witness in Table S1(a).
The removals act on \((106,110,162)\) and \((106,110,161)\), respectively.
The target voxel \(v\) is background in all three states and is filled at step 3.
We recompute \(A_{Y_k}(v)=K(Y_k\setminus\{v\})\cap C_v\) from the stored states.

\begin{center}
\small
\begin{tabular}{llllll}
\toprule
State & Last operation & \(\beta(A_{Y_k}(v))\) & SAFE & \((n_0,n_1,n_2)\) & \((r_1,r_2)\)\\
\midrule
__ROWS__
\bottomrule
\end{tabular}
\end{center}

Here \(n_j\) is the number of attachment cells of dimension \(j\),
and \(r_j\) is the rank of the cellular boundary map \(\partial_j\) over
\(\mathbb F_2\). The Betti vector is obtained by
\[
\beta_0=n_0-r_1,\qquad
\beta_1=n_1-r_1-r_2,\qquad
\beta_2=n_2-r_2.
\]
All three chain complexes satisfy \(\partial_1\partial_2=0\).
The computation uses the closed contact cells from the 26-neighbourhood of
\(v\), including their faces. It agrees with the attachment computation used
in the existing witness probe; no new search is performed. The filling voxel
is already SAFE after REM 1 and remains SAFE after REM 2.

\paragraph{Relation to Table S1.}
The SAFE column of Table S1 tests the \emph{voxel flipped in that row}
in the state \emph{immediately before that move}: \(v_i\) in \(Y_{i-1}\)
for display steps \(i=1,\ldots,4\). Table S4 instead tests the \emph{same fixed
voxel} \(v=(106,111,162)\) in each of \(Y_0,Y_1,Y_2\).
Only its \(Y_2\) row corresponds to the SAFE check for ADD 3 in Table S1(a).
The SAFE checks for REM 1 and REM 2 in Table S1 concern different voxels and
do not determine when the later filling voxel first becomes SAFE.

\paragraph{Stored-state and code identifiers.}
The original witness bundle is
\texttt{FIGURE\_DATA/shard\_1/FIG/}\allowbreak
\texttt{BraTS2021\_00110\_\_k2xy\_106\_111\_162/}\allowbreak
\texttt{C00-F1\_\_G26.npz}, arrays \texttt{wit\_states[0,:]},
\texttt{wit\_states[1,:]} and \texttt{wit\_states[2,:]} in the order of
\texttt{coords} (180 voxels). The context is \texttt{ctx.npz}, with crop origin
\((96,101,152)\) and shape \((21,22,22)\).
The stored path is matched to \texttt{PRODUCTION/shard\_1/PATHS.csv},
lines 211--214 (header is line 1); its steps 0--3 are display steps 1--4.
Replay matches all five stored witness states exactly.

The accompanying verification records contain \texttt{ATTACHMENT\_RECALCULATION.json},
\texttt{ATTACHMENT\_RECALCULATION.csv},
\texttt{ATTACHMENT\_CELL\_LEDGER.json}.
The internal replay array \texttt{REPLAYED\_FIRST\_THREE\_STATES.npz}
is retained with the source witness for verification.
It records SHA-256 digests of the context, witness bundle, path inputs,
coordinate array and each state, with the byte encoding specified in the JSON.
The boundary-rank calculation is identified by
\texttt{code/close\_evidence.py}, functions
\texttt{cells\_from\_neighbours}, \texttt{columns} and \texttt{rank\_f2}.
Its SHA-256 is \texttt{__CODEHASH__}.
The existing probe is
\texttt{MANUSCRIPT/PHASE\_S/code/attach\_along\_witness.py::attach\_betti},
SHA-256 \texttt{__PROBEHASH__}.

\begin{center}
\small
\begin{tabular}{lp{0.78\linewidth}}
\toprule
State & SHA-256 of the stored 180-entry occupancy array\\
\midrule
__HASHROWS__
\bottomrule
\end{tabular}
\end{center}
\endgroup
'''
def breakhash(h):return r'\allowbreak{}'.join(h[i:i+16] for i in range(0,len(h),16))
appendix=appendix.replace('__ROWS__','\n'.join(texrows)).replace('__CODEHASH__',breakhash(sha(Path(__file__)))).replace('__PROBEHASH__',breakhash(sha(probe)))
appendix=appendix.replace('__HASHROWS__','\n'.join(r['state']+r' & \texttt{'+breakhash(r['witness_state_sha256'])+r'} \\' for r in records))
(OR/'Table_S4_attachment_evidence.tex').write_text(appendix,encoding='utf-8')
oldor=(CLOSURE/'Online_Resource_3_R4H_candidate.tex').read_text(encoding='utf-8')
assert oldor.count(r'\end{document}')==1
newor=oldor.replace(r'\end{document}',appendix+'\n'+r'\end{document}')
title_old='Supplementary Tables S1–S3 and Supplementary Figure S1'
title_new='Supplementary Tables S1–S4 and Supplementary Figure S1'
assert newor.count(title_old)==1
newor=newor.replace(title_old,title_new,1)
assert newor.replace(appendix+'\n','').replace(title_new,title_old,1)==oldor
(OR/'Online_Resource_3_NARRATIVE_CLOSURE.tex').write_text(newor,encoding='utf-8')
(D/'OR3_added_Table_S4.diff').write_text(''.join(difflib.unified_diff(oldor.splitlines(True),newor.splitlines(True),fromfile='OR3_R4H_candidate',tofile='OR3_NARRATIVE_CLOSURE')),encoding='utf-8')
shutil.copy2(CLOSURE/'FIG_B1_R1_supp_control_00312.pdf',OR/'FIG_B1_R1_supp_control_00312.pdf')
# OR1/OR2 compile only as unchanged copies.
shutil.copy2(FROZEN/'Online_Resource_1.tex',OR/'Online_Resource_1_UNCHANGED.tex')
shutil.copy2(CLOSURE/'Online_Resource_2_R4H_candidate.tex',OR/'Online_Resource_2_UNCHANGED.tex')

for name in ('main_R4I_article.tex','main_R4I_sn-jnl_iicol.tex'):
    old=(SRC/name).read_text(encoding='utf-8');new=old
    # Narrow changes to the two SAFE-timing citations, keeping all prose intact.
    a='and becomes SAFE after the first flip, remaining SAFE after the second (Online Resource 3, Table S1).'
    b='and becomes SAFE after the first flip, remaining SAFE after the second (Online Resource 3, Table S4).'
    c='that voxel is already SAFE after step 1 (Online Resource 3, Table S1).'
    d='that voxel is already SAFE after step 1 (Online Resource 3, Table S4).'
    assert new.count(a)==new.count(c)==1,(name,'citation count')
    new=new.replace(a,b,1).replace(c,d,1)
    assert new.replace(b,a,1).replace(d,c,1)==old
    (MAN/name).write_text(new,encoding='utf-8')
    (D/(name+'.diff')).write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='APPROVED_NARRATIVE/'+name,tofile='EVIDENCE_CLOSED/'+name)),encoding='utf-8')
for p in SRC.iterdir():
    if p.suffix in ('.pdf','.bib','.cls','.bst') and not p.name.startswith('R4I_NARRATIVE'):shutil.copy2(p,MAN/p.name)

# 14 systems x two horizons = 28 existing result rows. Read the OR3 table and
# compare it with the unchanged effective-result CSV, without altering either.
part=oldor.split(r'\section*{Supplementary Table S2.')[1].split(r'\section*{Supplementary Table S3.')[0]
matches=re.findall(r'^(C\d\d(?:P|H)?-F\d) & (\d+ / \d+ / \d+ / \d+) & (\d+ / \d+ / \d+ / \d+)',part,re.M)
assert len(matches)==14
rows=[]
for system,g14,g26 in matches:
    for grid,s in (('G14',g14),('G26',g26)):rows.append(dict(system=system,grid=grid,**dict(zip(('FOUND','EXH','DL','CAP'),map(int,s.split(' / '))))))
with (PHASE/'logs/recomputed_table_S2_effective.csv').open(encoding='utf-8-sig',newline='') as f:
    expected={(r['system'],r['grid']):{k:int(r[k]) for k in ('FOUND','EXH','DL','CAP')} for r in csv.DictReader(f)}
assert all({k:r[k] for k in ('FOUND','EXH','DL','CAP')}==expected[(r['system'],r['grid'])] for r in rows)
(E/'PRESERVED_28_RESULT_ROWS.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
after={p:sha(PROJECT/p) for p in before}
assert before==after,'Protected input was modified'
check=dict(stored_witness_states_matched=5,betti_vectors=[r['beta'] for r in records],boundary_matrix_crosscheck=True,
    result_rows_preserved=28,or3_original_body_unchanged_except_title_range=True,or1_or2_sources_byte_identical=True,
    manuscript_changes_only_two_citations=True,figures_byte_identical={p.name:sha(p)==sha(SRC/p.name) for p in MAN.glob('F?_*.pdf')},
    protected_files_unchanged=len(before))
(ROOT/'EVIDENCE_CLOSURE_VERIFICATION.json').write_text(json.dumps(check,indent=2),encoding='utf-8')
print(json.dumps(check,indent=2))
