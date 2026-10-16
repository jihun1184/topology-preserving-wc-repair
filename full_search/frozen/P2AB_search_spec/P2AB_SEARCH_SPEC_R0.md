# P2A-restricted exhaustive search · Fidelity search — SPEC R0 (preregistration)

작성 2026-10-08, **어떤 탐색 결과도 보기 전에 고정**한다. 이 spec은 Task 2(P2A-restricted search)와 Task 3(fidelity search)를 함께 정의한다. 결과는 `P2AB_restricted_search/`와 `P2AB_fidelity_search/`에 별도로 쓴다. Stage 1a/1b 산출물은 변경하지 않는다. P2B Run A 산출물은 읽기 전용 입력이다.
태그: [DEF] 정의, [FIN] 재계산·탐색에서 얻은 사실, [OBS] 관찰된 패턴, [HEUR] 가설, [CAUSAL]·[EXTRAP] 금지.

## 0. 사전에 확인된 정의의 귀결 [DEF] — 이 spec이 transition system을 둘로 나누는 이유
X_T는 P2A-FO의 LOCAL_STUCK 최종 상태이므로 X_T에는 critical block의 voxel 안에 있고 SAFE이며 flip-once가 허용하는 단일 flip이 **없다**. 따라서 flip-once를 **P2A history(RAW 기준)로 적용**하면 L1/FO graph는 X_T에서 transition이 0개이고, 탐색은 root 한 state에서 자명하게 끝난다(Stage 1b S7이 222개 X_T 전부에서 n_L1_moves=0을 재계산). 이 탐색은 Stage 1a §4가 남긴 질문(P2A가 flip-once로만 제외한 in-block SAFE flip으로 시작하는 PURE 경로가 있는가)을 닫을 수 없다. 그래서 flip-once의 적용 범위를 둘로 사전 등록하고 **절대 섞지 않는다**.
- `HIST`: flip-once = RAW 규칙. ADD(v)는 RAW(v)=0일 때만, REM(v)는 RAW(v)=1일 때만 허용(P2A-FO의 `allowed_add/allowed_rem`과 동일; P2A history와 이후 경로를 모두 포함). **사전 예측: 680/680 target에서 transition 0, visited state 1.** 이 예측이 어긋나면 구현 오류로 보고 결과를 보고하지 않는다.
- `PATH`: flip-once = X_T 이후 탐색 경로 안에서만 적용. voxel v는 X(v)=X_T(v)(아직 flip되지 않음)일 때만 flip 가능. P2A가 X_T 이전에 한 flip의 복원·제거도 허용된다(P2A-RZ에 가까운 완화). 정보가 있는 시스템은 이쪽이다.

## 1. 공통 설정 (P2B와 동일한 동결 구조)
- **Target 집합**: P2B target 680 block(C1 614 / 213 case, C2 66 / 39 case). 합산 열 없음. P2B Run A에서 PURE witness가 없던 C1 5 block(`CAP_EXCEEDED`)도 target에 포함되며 별도 행으로 보고한다.
- **시작 상태·provenance**: P2B runner(SHA `ed683826…c62b`)를 수정 없이 import한다. `load_env`(P2A frozen 13 pin + P4 manifest 3,512 파일), `RunA.reconstruct`(P4 RAW → Run A edit replay → X_T, `mask_sha` 대조), P2B `Local`(window + 1-voxel halo 기하), `crit_mask`, 벡터화 gate 조회 + `audit_gate`. 환경은 Python 3.13.15, numpy 2.5.3, scipy 1.18.1, nibabel 5.4.2이며 이 기기에서 P2B Q0 680/680·X_T SHA 222/222 재현을 사전에 확인했다. Run B는 쓰지 않는다.
- **Window (허용 voxel 집합)**: W_R(B) := Ω ∩ {v : dist_∞(v, bbox(B)) ≤ R}, tier R ∈ {1,2,3}(P2B Tier A/B/C의 R). flip은 W 안에서만 일어난다. W 밖은 X_T로 고정된 halo이며 gate·critical 판정은 W ∪ halo에서 정확하다(P2B spec §4.2). 조건 "critical block 안"은 현재 state의 W ∪ halo 안에 완전히 포함된 critical block만 센다; W의 voxel을 포함하는 block은 모두 이 crop 안에 있으므로 이 제한은 판정을 바꾸지 않는다.
- **State / canonicalization**: state = window occupancy bitmask(bit i = window voxel i, 축 (z,y,x) 사전순; 최하위 비트가 첫 voxel). 대칭성 축소·다른 canonicalization 없음. 중복 판정은 정확한 byte key의 집합이다.
- **Gate (SAFE)**: frozen P0 gate table(`bits&7==7`); 사용된 모든 distinct neighbourhood bitmask를 frozen `P0.attachment`+`P0.verdict`로 재평가하는 감사를 block·tier마다 수행(불일치 = ERROR).
- **PURE (success, P2B와 동일)**: state에서 target block B가 critical이 아니고, X_T에는 없던 in-domain critical block이 없다.
- **cap**: N_max = 10⁶(root 포함 unique state). N_max+1번째 state를 삽입해야 하는 순간 `CAP_EXCEEDED`(P2B R1.1 규칙 동일); CAP는 어떤 부재 주장도 허용하지 않는다.
- **Depth**: restricted(§2) 시스템은 **depth 제한 없음**. 모든 transition이 아직 flip되지 않은 voxel을 flip하므로(HIST/PATH 모두 voxel당 최대 1회) 그래프는 비순환이고 depth ≤ |W|, 탐색은 frontier가 비어서 종료한다(또는 CAP). fidelity(§3) 시스템은 P2B와 같은 **depth 제한 L**(tier A 4, B 6, C 8)을 둔다.

## 2. Task 2 — P2A-restricted exhaustive search (4개 시스템, 서로 결과를 섞지 않는다)
현재 state X에서 voxel v ∈ W의 flip이 **허용**되는 조건(모든 시스템 공통 부분):
 (a) **locality**: v ∈ ∪_{b∈Crit(X)} Vox(b) (현재 state의 critical block 중 하나의 voxel); (b) **SAFE**(v, X); (c) **flip-once**(HIST 또는 PATH, 위).
- **A_HIST** = (a)+(b)+(c:HIST). **A_PATH** = (a)+(b)+(c:PATH). (P2A-L1 candidate graph; phase·tie 규칙 없음.)
- **B_HIST**, **B_PATH** = A + **phase 규칙(FO policy, tie relaxed)**: ADD(X(v)=0)는 위 조건만으로 허용; **REM(X(v)=1)은 fixed_point(X)일 때만** 허용. fixed_point(X) := Crit(X)의 voxel 중 X=0, SAFE, 해당 시스템의 flip-once가 허용하는 voxel이 없다. 즉 허용된 SAFE ADD가 하나라도 있으면 REM은 열리지 않는다. P2A-FO의 FIFO queue 순서, 사전식 최소 voxel, "첫 deferred block" 같은 tie 규칙은 history(queue 상태)가 필요하므로 **relaxed**(어떤 block·voxel이든 허용). 따라서 B는 실제 P2A-FO 정책의 **상태 기반 완화(superset)**이며 deterministic P2A-FO 재현이 아니다. 완화로 인해 B에서 PURE가 없다면 deterministic FO에도 없다(B ⊇ FO). B_PATH에서 PURE가 있어도 deterministic FO가 그 경로를 따른다는 뜻이 아니다.
- 포함 관계 [DEF]: A_HIST ⊆ A_PATH, B_x ⊆ A_x, B_HIST ⊆ B_PATH (transition 집합 기준). 탐색 결과에서 이 순서와 모순되는 PURE_FOUND 패턴이 있으면 구현 오류.
- **Tier 처리**: 각 target에 대해 R=1→2→3 순으로, **PURE_FOUND가 아닌 한 다음 tier로** 진행한다(window 포함 관계로 PURE_FOUND at R ⇒ PURE_FOUND at R'>R). 모든 tier 결과를 행으로 기록한다.
- **종료 상태**: `PURE_FOUND`(최단 depth의 첫 PURE; BFS FIFO 순서) / `NO_PURE_EXHAUSTED`(frontier 소진, PURE 없음; 깊이 제한 없음) / `CAP_EXCEEDED` / `ERROR`. 부가 열: `n_states`(visited), `n_moves_root`(X_T에서의 허용 move 수), `level_counts`, `n_bres`(B가 해소된 state 수; PURE 아님 포함), `max_depth`, `frontier_exhausted`(0/1), `depth_limit`("none").
- **허용되는 부재 문장**: "NO_PURE_EXHAUSTED at (system, R)" = "이 시스템의 window W_R(B) 안 state space에서 frontier가 소진될 때까지 탐색했고 PURE state가 없다(깊이 제한 없음)". 이것은 restricted state space 내 no path이며 same-grid impossibility나 P2A 일반에 대한 정리가 아니다.
- **사전 예측(구현 점검)**: C-H: A_HIST, B_HIST는 680/680에서 `n_moves_root=0`, `n_states=1`, `NO_PURE_EXHAUSTED`. C-T: A_PATH의 root move 수는 case의 tabu-excluded SAFE in-block flip(window 안)과 voxel 단위로 일치.

## 3. Task 3 — Fidelity search (P2B의 move set, 같은 tier·cap·routing)
P2B의 move set(**locality·flip-once 없음**, 모든 window voxel의 reversible SAFE ADD/REM), 같은 BFS(FIFO, exact dedup), 같은 tier(A: R=1,L=4; B: R=2,L=6; C: R=3,L=8), N_max 10⁶, 같은 routing(`P2B_CONFIG_R1.json`: A→B on {소진 NO_PATH, B_RESOLVED, CAP}, B→C on {소진 NO_PATH, B_RESOLVED}, 그 외 stop)을 쓴다. **fidelity 조건만 추가**한다.
- **RAW 전경 순제거 NRFR(X)** := |{v ∈ Ω : RAW(v)=1 ∧ X_T(v)=1 ∧ X(v)=0}| (P2B `net_raw_fg_removed`와 같은 정의; 기준선은 X_T이며 **P2A가 X_T까지 이미 제거한 원본 voxel은 세지 않는다**).
- **F1 (endpoint fidelity)**: 성공 = PURE(X) ∧ NRFR(X)=0. transition 제한 없음(remove-then-restore 허용). PURE이지만 NRFR>0인 state는 expand를 계속한다.
- **F2 (pathwise fidelity)**: 성공 = PURE(X). transition 제한: RAW(v)=1인 voxel의 **REM 금지**(어느 step에서도). RAW(v)=0인 voxel(P2A가 추가했거나 경로가 추가한 voxel)의 REM과 모든 ADD는 허용. [DEF] F2 경로는 NRFR≡0이므로 F2 ⊆ F1(F2의 PURE_FOUND ⇒ 같은 tier에서 F1의 PURE_FOUND, lstar_F1 ≤ lstar_F2).
- **F0 (control)**: fidelity 조건 없음 = P2B 재현. 680 target의 status·lstar·`n_states`·level_counts가 P2B Run A BLOCKS와 **행 단위로 일치**해야 한다(`n_bres`·`min_depth_bres` 포함).
- 부가 열(F1): `n_pure_seen`, `n_pure_unfaithful_seen`(NRFR>0인 PURE state), `min_nrfr_pure_seen`, 해당 NRFR 히스토그램. 탐색이 소진·depth 제한으로 끝나면 이 값은 탐색한 전 공간에 대한 것이며 PURE_FOUND로 일찍 끝나면 "발견한 범위"이다(과소계상 가능; 그렇게 표기).
- **허용되는 부재 문장**: "F1/F2 NO_PURE at (R, L): window W_R 안에서 ≤ L회 flip 이내에 …를 만족하는 PURE state가 없다(탐색 소진, 깊이 제한 L 존재)". 소진 종료는 두 종류를 구분한다: `NO_PURE_EXHAUSTED`(frontier가 L 이전에 소진)와 `NO_PURE_DEPTH_LIMITED`(depth L에서 frontier가 남은 채 종료). P2B의 `NO_PATH_EXHAUSTIVE`는 둘을 합친 이름이므로 비교 시 매핑한다.

## 4. Witness 기록과 독립 검증
- PURE_FOUND 경로는 PATHS에 step별로 기록(type, voxel, nmask). 모든 witness를 **전체 volume에서 독립 replay**한다(P1 gate verdict SAFE, step마다 `global_crit` 재계산, 각 step이 해당 시스템의 규칙을 만족하는지 재판정: locality·flip-once·phase / fidelity, endpoint PURE, NRFR, accounting).
- **독립 BFS 교차검증(Q-ind)**: 고정 seed(20261008)로 target 표본을 뽑아(시스템당 ≥ 60, ≤ 20,000 state인 것만) per-state 전역 재계산 pure-Python BFS와 vectorized 탐색의 (status, n_states, lstar)를 대조한다.
- **tabu 일치(C-T 확장)**: 222개 X_T 전부에서 "in-crit SAFE voxel 수 == STUCK의 n_safe_add_excluded + n_safe_rem_excluded(block별)" 재계산.
- **포함 관계 검사**: A_HIST ⊆ A_PATH ⊇ B_PATH ⊇ B_HIST (PURE_FOUND 방향), F2 ⇒ F1, F0 ⇒ P2B.
- 모든 shard는 write-once, 중단된 shard(`shard_meta.json` 없음)는 `_incomplete_*`로 옮기고 같은 입력으로 다시 실행한다(결과는 shard 배치와 무관해야 하며 검증한다).

## 5. 보고 규칙
cohort(C1/C2)는 합산하지 않는다. 시스템(A_HIST/A_PATH/B_HIST/B_PATH/F0/F1/F2)의 결과는 서로 다른 표로 보고하고 합치지 않는다.
쓸 수 있는 문장 [FIN]: "시스템 S, tier R에서 target ⟨n⟩개가 PURE_FOUND(최단 depth ⟨…⟩), ⟨m⟩개가 NO_PURE_EXHAUSTED(frontier 소진, 깊이 제한 없음/L), ⟨c⟩개가 CAP_EXCEEDED(미결)." "PURE_FOUND 경로의 RAW 순제거 수는 ⟨…⟩."
쓰지 않는 문장: "P2A-admissible PURE 경로가 없다/있다"(A_HIST 이외 시스템의 결과를 P2A-FO 일반으로 확대), "locality가 P2A 실패의 원인"(인과), "RAW 전경 제거가 위상 보존에 필수/불필요"(window·L·move set 한정을 벗어나는 일반화), "case SOLVED"(block 단위 결과), HIST와 PATH 또는 F1과 F2의 결과를 합친 문장, CAP_EXCEEDED에서의 부재 주장.

## 6. 한계
1. 단일 voxel flip move set(simple pair 등은 도달 못 함), 단일 데이터셋(BraTS WT), block별 독립 탐색(case 합성 미검증), 임상 평가 아님. 2. B 시스템은 tie 규칙을 완화한 상태 기반 FO이다. 3. PATH 시스템은 P2A-FO의 정의가 아니라 flip-once 적용 범위를 바꾼 변형이다(P2A-RZ에 가깝다). 4. window 안으로 제한된 state space이며 R=3 이상의 결과·window 밖 이동의 필요성은 외삽하지 않는다.
