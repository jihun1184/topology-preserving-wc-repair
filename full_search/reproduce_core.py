"""Run the preserved bounded-search engine for BraTS2021_00110.

This adapter validates a single-case release contract instead of the historical
3,512-file cohort bootstrap. It does not implement or alter BFS, SAFE, admissible
moves, endpoint tests, state cap, or discovery order. See README for this scope.

Example: python -B full_search/reproduce_core.py --raw-wt /private/00110.nii.gz
         --out /private/results --work-dir /private/cache
Input masks and generated private working files must not be redistributed.
"""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json,platform,sys,time,types
import nibabel as nib
import numpy as np

HERE=Path(__file__).resolve().parent
FROZEN=HERE/'frozen'
def require(ok,message):
    if not ok:raise ValueError(message)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while chunk:=f.read(1024*1024):h.update(chunk)
    return h.hexdigest()
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module
def verify_sources():
    pins=json.loads((HERE/'FROZEN_SOURCE_SHA256.json').read_text())
    for rel,digest in pins.items():
        p=FROZEN/rel
        require(p.is_file() and sha(p)==digest,f'Frozen file missing or modified: {rel}')
    return len(pins)
def make_table(p0,path,expected):
    """Generate the exact pinned table with the original attachment classifier."""
    if not path.exists():
        rep={0:0}
        for k,closure in enumerate(p0.CL):
            for attachment,mask in list(rep.items()):
                rep.setdefault(attachment|closure,mask|(1<<k))
        keys=np.array(sorted(rep),dtype=np.uint32)
        values=np.array([7 if p0.classify(p0.attachment(rep[int(a)]))=='contractible' else 0 for a in keys],dtype=np.uint8)
        def lut(offset):
            result=np.zeros(8192,dtype=np.uint32)
            for i in range(1,8192):
                bit=i&-i
                result[i]=result[i^bit]|p0.CL[offset+bit.bit_length()-1]
            return result
        low,high=lut(0),lut(13)
        temporary=path.with_suffix('.tmp')
        with temporary.open('wb') as stream:
            for begin in range(0,1<<26,1<<20):
                masks=np.arange(begin,begin+(1<<20),dtype=np.uint32)
                attachment=low[masks&8191]|high[masks>>13]
                stream.write(values[np.searchsorted(keys,attachment)].tobytes())
        require(sha(temporary)==expected,'Generated SAFE table does not match frozen digest')
        temporary.replace(path)
    require(sha(path)==expected,'Cached SAFE table digest mismatch')

def main(a):
    start=time.perf_counter()
    require(sys.flags.optimize==0,'Run normal Python, not python -O')
    require(not (a.out/'results.json').exists(),'Output already exists; select a new --out directory')
    # Work/output directories may contain private data; refuse locations in source.
    repo=HERE.parent.resolve()
    for dest in (a.out,a.work_dir):
        require(not dest.resolve().is_relative_to(repo),'Choose output/work directories outside the release source tree')
    source_count=verify_sources()
    contract=json.loads((HERE/'core_case.json').read_text())
    expected=json.loads((HERE/'expected_results.json').read_text())
    require(len(contract['systems'])==len(expected)==4,'Expected exactly four core-search specifications')
    require(contract['systems']==[[r['row']['system'],r['row']['grid']] for r in expected],
        'Contract/expected search order mismatch')
    source=a.raw_wt if a.raw_wt is not None else a.segmentation
    require(source.is_file(),'Provide an authorized local BraTS segmentation or binary WT NIfTI')
    image=nib.load(str(source))
    labels=np.asanyarray(image.dataobj)
    require(labels.ndim==3 and np.isfinite(labels).all(),'Input must be a finite 3D volume')
    allowed=(0,1) if a.raw_wt is not None else (0,1,2,4)
    require(np.isin(labels,allowed).all(),f'Unexpected labels; expected {allowed}')
    raw_xyz=(labels!=0).astype(np.uint8) if a.raw_wt is not None else np.isin(labels,(1,2,4)).astype(np.uint8)
    raw=np.ascontiguousarray(raw_xyz.transpose(2,1,0))
    require(list(raw.shape)==contract['raw_shape_zyx'],'Wrong input shape/orientation')
    p0=load('p0_release',FROZEN/'P0_topology_audit/code/p0_audit.py')
    require(p0.mask_sha(raw)==contract['raw_mask_sha256'],'RAW mask digest mismatch: wrong case, release, labels or orientation')
    require(int(raw.sum())==contract['raw_foreground_voxels'],'RAW foreground count mismatch')
    # Materialize only the user-supplied case privately so original RunA.reconstruct
    # can load/replay it unchanged. The NIfTI container digest need not match the
    # historical gzip bytes; the fully specified mask-array digest must match.
    rawdir=a.work_dir/'materialization/raw'
    rawdir.mkdir(parents=True,exist_ok=True)
    rawpath=rawdir/f"{contract['case_id']}_WT_RAW.nii.gz"
    nib.save(nib.Nifti1Image(raw_xyz,image.affine),str(rawpath))
    table=a.work_dir/'gate_table_2p26.bin'
    make_table(p0,table,contract['gate_table_sha256'])
    vendor=FROZEN/'vendor'
    sys.path.insert(0,str(vendor))
    repair=load('r_plus_repair_production',vendor/'r_plus_repair_production.py')
    p1=load('p1_release',FROZEN/'P1_simple_constrained/P1_simple_constrained_rplus.py')
    p1.setup(vendor,FROZEN/'P0_topology_audit/code/p0_audit.py',table,None)
    p2b=load('p2b_release',FROZEN/'P2B_local_edit_search/P2B_search_R1.py')
    search=load('search_release',FROZEN/'P2AB_search_spec/p2ab_search.py')
    ablation=load('ablation_release',FROZEN/'P2AB_constraint_ablation/CODE/p2ab_ablate.py')
    p2b.env.update(P0=p0,P1=p1,R=repair,TBL=p0.load_table(table),
        BETTI=load('betti_release',FROZEN/'cubical_betti.py').cubical_betti,
        cfg=json.loads((FROZEN/'P2B_local_edit_search/P2B_CONFIG_R1.json').read_text()))
    p2b.NMFLAG=np.zeros(1<<26,dtype=np.uint8)
    # Same wiring as the historical ablation bootstrap, with release-scoped
    # source/mask pins replacing the unrelated whole-cohort materialization gate.
    search.P2B=p2b
    ablation.M=search
    search.allowed=ablation.cell_allowed
    args=types.SimpleNamespace(runa=FROZEN/'P2A_removal_addition/P2A_PRODUCTION/RUN_A',
        p1=FROZEN/'P1_simple_constrained',p4=a.work_dir/'materialization')
    runa=p2b.RunA(args)
    cid,block=contract['case_id'],contract['block']
    require(('C1',cid,block) in runa.targets(('C1',)),'Target does not match frozen C1 records')
    xp,reconstructed_raw,hash_ok=runa.reconstruct(cid)
    require(hash_ok and p0.mask_sha(xp[1:-1,1:-1,1:-1])==contract['xt_mask_sha256'],'XT reconstruction mismatch')
    require(np.array_equal(raw,reconstructed_raw),'Private input materialization changed array values')
    b=p2b.parse_bk(block)
    cache=dict(Xp=xp,raw=raw,C0=p2b.global_crit(xp[1:-1,1:-1,1:-1]))
    results=[]
    a.out.mkdir(parents=True,exist_ok=True)
    for (system,grid),wanted in zip(contract['systems'],expected):
        row,path=ablation.run_one(runa,system,grid,'C1',cid,block,b,cache,contract['nmax'],None)
        require(row['status']!='ERROR',str(row))
        differences={k:[row[k],v] for k,v in wanted['row'].items() if row[k]!=v}
        require(not differences,f'Result mismatch: {differences}')
        require(path==wanted['path'],'Witness path/SAFE/path-rule flags do not match frozen result')
        results.append(dict(row=row,path=path))
        print(system,grid,row['status'],row['n_states'],row['level_counts'],flush=True)
    record=dict(validation='PASS',bootstrap_profile='selected-case-v1 (not whole-cohort P4 bootstrap)',
        python=sys.version,numpy=np.__version__,nibabel=nib.__version__,platform=platform.platform(),
        frozen_files_verified=source_count,input_kind='binary WT' if a.raw_wt is not None else 'BraTS segmentation labels1,2,4',
        supplied_input_file_sha256=sha(source),raw_mask_sha256=p0.mask_sha(raw),
        xt_mask_sha256=p0.mask_sha(xp[1:-1,1:-1,1:-1]),gate_table_sha256=sha(table),
        search_engine_unmodified=True,original_cohort_bootstrap_executed=False,
        total_seconds=time.perf_counter()-start,results=results)
    (a.out/'results.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('PASS: all four original-engine results and witnesses match the frozen expected records')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--raw-wt',type=Path)
    group.add_argument('--segmentation',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    try:main(parser.parse_args())
    except Exception as exc:
        print(f'REPRODUCTION FAILED: {type(exc).__name__}: {exc}',file=sys.stderr)
        raise SystemExit(1) from exc
