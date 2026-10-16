# F1/F2 completeness — SPEC R0 (preregistration, Priority A)

작성 2026-10-08. 이 spec은 **어떤 completion 결과도 보기 전에** 고정한다(개발 중 differential 테스트는 scratch에서만 수행했고 결과에 포함하지 않는다; §8).
기존 산출물(Stage 1a/1b, restricted, fidelity F0/F1/F2, verify)은 read-only 입력이며 덮어쓰지 않는다. 새 결과는 `P2AB_fidelity_completion/` 아래에만 쓴다.
태그: [DEF] 정의 · [FIN] 재계산·탐색으로 얻은 사실 · [OBS] 관찰 · [HEUR] 가설 · [CAUSAL]·[EXTRAP] 금지.

## 0. 이 spec이 다루는 것과 다루지 않는 것
- 다룬다: 동결된 F1(endpoint fidelity)·F2(pathwise fidelity) 탐색에서 `CAP_EXCEEDED`였던 (target, tier) 행의 **판정**(PURE_FOUND 또는 부재 certificate).
- 다루지 않는다: locality/flip-once ablation(Priority B), 모든 최단 경로 열거(Priority C), 새로운 transition system, move set 확장, window·L·tier 값 변경. 모든 (R, L)은 동결 config(`P2B_CONFIG_R1.json`)와 같다: A (1, 4), B (2, 6), C (3, 8).
- 핸드오프 문서와의 차이(정정): (i) F2의 미결은 CAP 448만이 아니다. 최종 tier C에서 `NO_PURE_DEPTH_LIMITED` 128(C1 110, C2 18)도 L=8 이내만 배제한 상태이며, F2 최종 분류는 FOUND 15 · EXHAUSTED 89 · DEPTH_LIMITED 128 · CAP 448(합 680)이다. (ii) 동결 routing에서 tier B가 CAP이면 정지하므로(`B→C`는 소진·B_RESOLVED일 때만), F1의 tier B CAP 608 target은 tier C를 **한 번도 시도하지 않았다**. 따라서 tier B를 판정해도 그 뒤의 tier C 행이 새로 생길 수 있다(§4 stage S2).

## 1. 입력 (read-only; SHA-256)
| 파일 | SHA-256 (앞 16자…) |
|---|---|
| `P2AB_fidelity_search/F1/BLOCKS.csv` | `82f0273004f9acee…` |
| `P2AB_fidelity_search/F1/PATHS.csv` | `f7afef401dccf6fb…` |
| `P2AB_fidelity_search/F2/BLOCKS.csv` | `f8a771a29b289945…` |
| `P2AB_fidelity_search/F2/PATHS.csv` | `91eff40120187b79…` |
| `P2AB_search_spec/p2ab_search.py` (동결 러너) | `5f32741b78e4024e…` |
| `P2AB_search_spec/P2AB_SEARCH_SPEC_R0.md` | `20026cb11e8872e5…` |
| `P2B_local_edit_search/P2B_search_R1.py` | `ed683826f6bbf3c0…` |
| `P2B_local_edit_search/P2B_CONFIG_R1.json` | `3b2fdcf2e5de1de9…` |
전체 SHA는 `P2AB_COMPLETION_PRE_RUN_SHA256.txt`(freeze manifest)에 기록한다. 동결 러너·P2B 러너·spec은 **수정 없이 import**한다.

## 2. CAP 원인 분석 [FIN] (`profile/CAP_PROFILE_F{1,2}.csv`, `p2ab_cap_profile.py`)
- CAP는 구현 결함이 아니라 상태공간의 지수적 성장이다. reversible 단일 voxel flip 공간은 window가 |W|=180(R=2)·448(R=3)이고 BFS level 크기가 level당 수십 배 증가한다(예: F1 tier B 한 행 `1|97|4705|152308|842889`).
- 직전 두 level 비율의 감쇠를 외삽한 총 상태 수 추정 [HEUR](실측 검증: 추정 9.9×10⁵ → 실제 1.007×10⁶, 6.9×10⁶ → 7.1×10⁶): F1 tier B CAP 608행 중 moderate(10⁶–10⁸) 437, heavy(10⁸–10¹⁰) 171. F2 tier C CAP 329행은 heavy 191·pathological(>10¹⁰) 22를 포함한다. 단순 N_max 증가는 heavy 행에서 불가능하다(메모리 O(상태 수)).
- `d_cert`(CAP 행에서 마지막으로 완전 생성된 level = `len(level_counts)−2`)는 **이미 확정된 부재 범위**다: F1 tier B CAP 608행은 d_cert 3 (152), 4 (285), 5 (171); 즉 "R=2 window에서 ≤d_cert flip 이내에 PURE∧NRFR=0 state 없음"은 현재도 [FIN]이지만 F1 found의 lstar 분포(C1 4–5, C2 2–5)와 비교하면 L=6까지의 공백이 남는다. 이 공백이 이번 spec의 대상이다.
- duplicate-state ratio·branching factor 전체는 동결 CSV에 기록되어 있지 않다(`n_states`·`level_counts`만 있음). 새 엔진은 `n_children_eval`, `n_children_surv`를 기록한다.

## 3. 엔진 FS — 정의와 정당성
### 3.1 문제 [DEF]
시스템 S ∈ {F1, F2}(동결 spec §3), window W=W_R(B), 깊이 제한 L, 시작 X_T. 블록 b의 *status*는 b의 voxel 값만의 함수다(critical 여부). **Bad(X)** := {B if B critical in X} ∪ {b : b critical in X ∧ b가 X_T에서 critical이 아님}. 동결 PURE(X) ⇔ Bad(X)=∅. **success(X)** := PURE(X) (F2) / PURE(X) ∧ NRFR(X)=0 (F1). 문제: 깊이 ≤L의 허용 flip 열로 도달 가능한 success state가 있는가, 있다면 BFS(FIFO, 부모 순서→voxel (z,y,x) 순서) 정의의 canonical 첫 state와 경로는?

### 3.2 알고리즘
1. **저장 level 0..L−2**: 동결 BFS와 동일 코드 경로(정확한 byte-key dedup, 같은 `allowed()`, 같은 PURE·NRFR 검사, 첫 success에서 중단, cap=`nstore`).
2. **level L−2의 각 parent**에 대해 level L−1 child를 생성하되 **저장·dedup하지 않는다**(chunk 내 중복만 첫 출현 순서로 제거). 각 child c = X⊕v에 대해:
   (a) success(c)이면 lstar=L−1 (parent 순서, voxel 순서의 첫 것).
   (b) 아니면 **마지막 flip 후보** C(c) := {w ∈ allowed(c) : w가 Bad(c)의 모든 block에 속함} (Bad(c)=∅이면 제약 없음); F1이면 NRFR(c)≥2는 제외, NRFR(c)=1이면 w는 RAW 전경 voxel의 ADD(복원)여야 한다. 각 후보 w에 대해 success(c⊕w)를 직접 검사하여 lstar=L의 첫 것을 기록한다.
   (c) level L−1에서 success가 하나라도 있으면 L−1이 lstar이다(깊이 L 후보는 버린다).
3. **Prefilter** (child 생성 직후, 순수 필요조건): parent의 B가 critical이고 v가 B의 voxel이 아니면 c의 B는 여전히 critical이다(§3.3 L1). 이때 c에서 *flip 가능한 B voxel이 하나도 없으면* c는 success의 조상이 아니므로 버린다. flip 가능 = 동결 gate `SAFE`(parent 이웃 mask에 v의 비트를 xor하여 얻은 c의 이웃 mask를 동결 gate table로 조회, 사용한 mask를 동결 audit의 NMFLAG에 기록) ∧ (F2: ¬(bit=1∧RAW=1)). F1은 추가로 NRFR(c)≥2이면 버린다. B-flip(v∈B)과 B가 이미 해소된 parent의 child에는 B 조건을 적용하지 않는다.
4. 판정: success 발견 → `PURE_FOUND`(lstar, 경로; 동결 `replay_witness`로 전체 volume 재생). level≤L−2 저장 중 `nstore` 초과 → `CAP_EXCEEDED`. 저장 level 중 frontier가 비면 `NO_PURE_EXHAUSTED`. 그 밖에 lookahead가 `max_children` 안에 끝나면 `NO_PURE_COMPLETE`, 아니면 `CAP_BUDGET`.

### 3.3 정당성 [DEF]·증명
**L1 (블록 국소성).** block status는 그 block의 voxel 값만의 함수이므로 voxel v의 flip은 v를 포함하는 block의 status만 바꾼다. 특히 v∉B이면 B의 status는 불변이다.
**L2 (마지막 flip의 필요조건).** success(X⊕w)이면 (i) w∈allowed(X); (ii) Bad(X)의 모든 block이 w를 포함한다 — Bad block은 Bad(X⊕w)=∅에서 더 이상 bad이어서는 안 되는데, 그 block의 status는 voxel 값이 바뀔 때만 변하므로 L1에 의해 w가 그 block에 속한다; (iii) F1이면 NRFR(X)−NRFR(X⊕w)∈{0,±1}이고 NRFR(X⊕w)=0이므로 NRFR(X)≤1이고, NRFR(X)=1이면 w는 RAW 전경 voxel의 ADD이다. 따라서 §3.2(b)의 후보 집합은 success를 주는 w를 모두 포함하고, 각 후보를 직접 검사하므로 **정확**하다(상위집합 + 정확한 검사).
**L3 (prefilter 무손실).** 저장 level L−2의 parent X에서 B가 critical이고 v∉B이면 child c의 B는 critical(L1)이므로 success(c)는 거짓이고, c에서 한 번의 flip w로 success에 도달하려면 B가 변해야 하므로 w∈B이고 w∈allowed(c)여야 한다. allowed(c)의 SAFE 부분은 w의 이웃 mask만의 함수(동결 gate 정의)이고 그 mask는 parent의 mask ⊕ {v의 비트}이다(v가 w의 26-이웃이 아니면 변화 없음). 따라서 flip 가능한 B voxel이 없으면 c의 어떤 자손도 깊이 ≤L에서 success가 아니다. NRFR 조건은 L2(iii)이다.
**L4 (순서 동일성).** 동결 BFS는 level ≤L−2는 동일 코드로 처리한다. level L−1의 state 목록은 (parent 순서, voxel 순서)의 첫 출현 순서이고, 동결 BFS에서 깊이 ≤L−2 state와 같은 child(역방향 간선 등)는 이미 방문되어 제외되지만 그런 child의 후손이 success라면 그 success state는 깊이 ≤L−2에서 먼저 발견되므로(최단성) 이 지점에 도달하지 않는다. 중복 child를 제거하지 않아도 첫 출현이 동결 목록과 같은 위치이므로 첫 success의 (깊이, 순서)는 같다. 따라서 canonical 경로는 동결 BFS가 완주했을 때의 경로와 **동일**하다. (경험적 확인: §6 V1.)
**한계 [DEF].** 엔진은 level L−1·L을 저장하지 않으므로 `NO_PURE_EXHAUSTED`(frontier가 L 이전에 소진)와 `NO_PURE_DEPTH_LIMITED`를 level L−1·L에서는 구분하지 못한다. 결과 이름 `NO_PURE_COMPLETE`는 "window W_R(B) 안에서 ≤L flip 이내에 success state가 없다"만을 뜻한다. 동결 행의 `NO_PURE_EXHAUSTED`는 저장 level 중 소진된 경우에만 재현된다.

## 4. 실행 대상과 단계
- **S1** : 동결 F1 BLOCKS에서 `status=CAP_EXCEEDED ∧ tier=B`(608행)와 동결 F2 BLOCKS에서 같은 조건(119행). 엔진 `nstore=10,000,000`, `max_children=3×10⁹`(2 프로세스 병렬, 프로세스당 메모리 ≈2.5 GB). **tail 재실행**: S1에서 `CAP_EXCEEDED`(저장 한도) 또는 `CAP_BUDGET`인 행은 단일 프로세스로 `nstore=20,000,000`, `max_children=10¹⁰`에서 한 번 더 실행하며, 두 실행 결과를 모두 기록하고 마지막 실행 결과를 판정으로 쓴다.
- **S2** : (a) 동결 F1 `CAP_EXCEEDED ∧ tier=C`(45행), 동결 F2 `CAP_EXCEEDED ∧ tier=C`(329행); (b) S1에서 `NO_PURE_COMPLETE`이고 동결 routing(`B→C` on 소진)에 따라 tier C로 가야 하는 새 행. S2는 S1 결과와 비용 실측을 본 뒤 **별도 spec 부록**으로 범위를 고정한다(여기서는 포함 여부만 예약).
- 동결 행 중 `PURE_FOUND`/`NO_PURE_*`인 행은 다시 계산하지 않는다(회귀 대조용 표본만 §6 V1).
- target 순서·샤드: S1 행 목록 `P2AB_COMPLETION_ROWS_S1.csv`(동결 CSV의 정렬 후 추정 비용 오름차순, 2 shard round-robin). 결과는 shard별 write-once, 중단 행은 (system,case,block,tier) journal로 재개한다. 결과는 shard 배치와 무관해야 한다(§6 V4).

## 5. 보고 규칙과 허용·금지 문장
- 시스템(F1/F2), tier, cohort(C1/C2)를 합치지 않는다. 동결 결과와 completion 결과는 서로 다른 표로 보고하고, 합친 "최종" 표에서는 출처 열을 둔다.
- **허용 [FIN]**: "F1, (R=2, L=6): target ⟨n⟩개에서 window W_2(|W|=180) 안 ≤6 flip 이내에 PURE∧NRFR=0 state가 없다(`NO_PURE_COMPLETE`, 엔진 FS, 저장 level ≤4 + lookahead)." "F1/F2 PURE_FOUND 경로 ⟨표⟩; 전체 volume 독립 replay 통과." "`CAP_EXCEEDED`/`CAP_BUDGET` 행은 미결(부재 주장 없음)."
- **금지**: 위 문장을 L>6·R>2·다중 voxel move로 일반화; 엔진 결과를 동결 결과와 같은 "BFS"라고만 쓰는 것(알고리즘 이름을 밝힌다); F1과 F2 합산; "RAW 전경 제거가 위상 보존에 필수/불필요"; CAP 행의 부재 주장; "case SOLVED".

## 6. 검증 계획 (결과 확인 전에 고정)
- **V1 회귀 대조**: 동결 비-CAP 행 중 F1 found 27 + 무작위 음성 60, F2 found 15 + 무작위 음성 60(seed 20261008)을 엔진으로 재계산하여 status 클래스(FOUND / NO_PURE)·lstar·경로(type,z,y,x,nmask)·accounting 열이 동결 행과 **0 diff**. 통과하지 못하면 S1을 실행하지 않는다.
- **V2 동결 BFS 대조**: S1 결과 중 저장 상태 수가 동결 BFS로 완주 가능한 행(총 상태 ≤1.5×10⁷, 동결 `bfs` 수정 없이 N_max 상향)에서 층화 표본(음성 ≥60, found 전부)을 동결 BFS로 재계산하여 status·lstar·경로 일치. 불일치 1건이라도 있으면 해당 spec 결과를 보고하지 않고 원인을 기록한다.
- **V3 독립 구현 대조**: 독립 validator(`p2ab_search_validate.py qind`, per-voxel gate, set-based BFS)로 작은 행(≤20,000 state)을 대조하고, 모든 새 witness를 validator `witness`로 전체 volume replay.
- **V4 샤드 배치 불변성**: S1에서 40행을 1 shard와 2 shard로 재실행하여 동일.
- **V5 포함 관계**: F2 found ⇒ 같은 tier에서 F1 found, F1에서 lstar_F1 ≤ lstar_F2.
- **V6 gate audit**: 엔진이 직접 조회한 이웃 mask는 NMFLAG에 기록되어 동결 `audit_gate()`에서 frozen verdict와 대조된다(불일치=ERROR).

## 7. 환경
Python 3.13.16 + numpy 2.5.3 + scipy 1.18.1 + nibabel 5.4.2(전용 venv). 동결 실행의 `META.json`은 3.13.15이며 3.13.15 인터프리터를 이 환경에서 받을 수 없었으므로 **이 환경 변경을 명시**한다(P2B spec R1은 3.13.16을 기록). 환경 승인 조건: gate table 재생성 SHA `59407fef…4830` 일치(통과), 동결 F0/F1/F2 44행(F0 8, F1 16, F2 20)을 새 환경에서 재계산하여 동결 행과 0 diff(통과, scratch). gate table·P4 manifest·vendor·betti의 동결 pin은 동결 러너의 provenance gate가 매 실행 시작 시 검사한다. `PYTHONHASHSEED=0`, `python -P -s -B`.

## 8. 사전 공개: spec 이전의 개발 이력
엔진은 이 spec 이전에 scratch에서 개발·테스트했다: (i) 첫 prefilter 구현은 "B가 이미 해소된 parent"를 구분하지 않아 F2 한 행(`BraTS2021_01282 k2xy@85,134,113`, tier B)에서 동결 경로와 불일치했고, differential 테스트가 이를 잡아 §3.2(3)의 마지막 문장 조건으로 수정했다; 수정 전 엔진 산출물은 어떤 결과에도 사용하지 않았다. (ii) 수정 후 differential: 동결 비-CAP 162행 0 diff, 경계 테스트 155건 0 mismatch(found 행의 L=lstar−1…lstar+2, 작은 CAP 행의 L=6), 무작위 target 60개 × 2 시스템 × (R,L) 7격자 differential 806건(found 26, 음성 780; 동결 BFS cap 초과 34건 제외) 0 mismatch. 이 숫자는 개발 확인이며 V1–V5 사전 등록 검증을 대체하지 않는다. 위 테스트의 코드는 결과 폴더에 `dev/`로 보관한다.
