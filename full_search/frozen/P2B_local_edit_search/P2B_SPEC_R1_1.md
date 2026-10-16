# TP-WC P2B — Specification R1.1 (bounded local SAFE-edit search on the original grid)

작성일 2026-10-07 · 상태: **SPEC R1.1 (freeze 대상; 사용자의 freeze 승인 전에는 production search 실행 없음).** `P2B_SPEC_R1.md`에 작은 patch를 더한 문서이며 R1·R0는 덮어쓰지 않고 보존한다. scientific gate, BFS, move set, PURE 판정, tier 값은 R1과 동일하다. 이 문서는 결과를 보기 전에 cohort, 시작 상태, window, halo, move set, state dedup, 성공 판정, 종료 상태, tier routing, 보고 규칙, validation contract를 고정한다.
태그: **[DEF/THM]** 정의·정리, **[FIN]** finite computation, **[EMP]** 데이터 관찰, **[HEUR]** 설계 판단·미검증 가설.
용어: **finite certificate** = 지정한 (start state, window, depth, move set, 성공 판정)에 대해 exhaustive 탐색이 끝났다는 사실. **obstruction** = 증명된 불가능성. 이 문서의 어떤 결과도 same-grid impossibility theorem으로 확대하지 않는다. LOCAL_STUCK, TABU_BLOCKED, GATE_ALL_UNSAFE, "no path"를 obstruction이라고 부르지 않는다.

## 0-A. R1 → R1.1 patch (사용자 승인 2026-10-07; 이 세 가지만 바뀌었다)
1. **Tier B의 `B_RESOLVED_WITH_NEW_CRIT`도 Tier C로 승격**한다. 승격 조건은 "B에서 PURE_FOUND가 아니고 탐색이 cap 없이 소진되었음"이며 `NO_PATH_EXHAUSTIVE`와 `B_RESOLVED_WITH_NEW_CRIT`가 모두 해당한다. `CAP_EXCEEDED`는 소진이 아니므로 C로 자동 승격하지 않고 computationally unresolved로 남는다. 변경은 config의 `routing.B.when`뿐이다.
2. **exhaustive-without-PURE의 의미를 명시**한다. 소진된 탐색에서 PURE가 없으면 해당 (X_T,R,L,move set,PURE)에 대해 "PURE 경로 없음"은 finite fact이다. enum 이름은 바꾸지 않는다: `B_RESOLVED_WITH_NEW_CRIT`는 diagnostic outcome, `NO_PATH_EXHAUSTIVE`는 B를 해소하는 endpoint조차 없었던 outcome이다. **분석 flag `PURE_NO_PATH_EXHAUSTIVE = (status ∈ {NO_PATH_EXHAUSTIVE, B_RESOLVED_WITH_NEW_CRIT})`**이며, "tier T에서 PURE no-path certified" 수는 두 outcome을 모두 포함해 센다(하나를 빼면 과소계상).
3. **Cap counting rule을 명시**한다(§4.5): N_max=10⁶은 root state를 포함해 depth 0..L에서 발견된 unique state 전체를 센다. depth-L state도 count에 포함하나 expand하지 않는다. 정확히 N_max개까지 발견한 뒤 queue가 정상 소진되면 exhaustive이며, **N_max+1번째 unique state를 삽입해야 하는 순간 CAP_EXCEEDED**로 종료한다(구현은 R1과 동일; 경계 시험 `CAP_BOUNDARY` 추가).
추가: validator assertion — `B_RESOLVED_WITH_NEW_CRIT`(및 `NO_PATH_EXHAUSTIVE`) 행은 PURE 미발견(경로 행 없음, lstar=−1), queue 소진, cap 미도달(n_states ≤ N_max, complete=1)을 만족해야 하고, `B_RESOLVED_WITH_NEW_CRIT`는 B-resolved endpoint가 ≥1개 관찰되어야 한다(n_bres ≥ 1, 최초 depth ∈ [1,L]). P4 pinned ZIP을 scratch에 재추출하고 ZIP SHA → manifest 3,512 → RAW → Run A replay → X_T SHA chain을 최종 SHA로 다시 통과시킨다. 재검증 범위: runner·validator·config·spec SHA가 바뀌었으므로 최종 SHA로 전체 qualification 단계를 다시 실행했다(freeze report R1.1).

## 0. R0 → R1 변경 요약

| 항목 | R0 | R1 |
|---|---|---|
| 종료 상태 이름 | PATH_PURE / NO_PURE_PATH_CERT(+UNRESOLVABLE·EXCHANGE_*) / CAP_EXCEEDED / ERROR | **PURE_FOUND / B_RESOLVED_WITH_NEW_CRIT / NO_PATH_EXHAUSTIVE / CAP_EXCEEDED / ERROR** (§4.6) |
| Finite certificate | NO_PURE_PATH_CERT (교환 라벨 포함) | **NO_PATH_EXHAUSTIVE만** certificate로 부른다 (§6) |
| 1-voxel halo | 미명시 | window 밖 1-voxel halo를 frozen X_T에서 가져와 gate와 critical 판정에 사용 (§4.2) |
| Q0 | Q2(X_T 복원)만 | **680/680 start-state qualification**, 하나라도 실패하면 production 금지 (§3.2) |
| Tier routing | B의 NO_PURE_PATH_CERT만 C로 | A에서 PURE_FOUND가 아니면 B로; B에서 PURE_FOUND가 아니고 소진된 경우(NO_PATH_EXHAUSTIVE, B_RESOLVED_WITH_NEW_CRIT) C로 (§5; R1에서는 NO_PATH_EXHAUSTIVE만이었으나 R1.1에서 승인되어 변경) |
| Cap 의미 | "distinct states > N_max" | depth-L leaf state 포함 admitted state 수, FIFO 발견 순서 rank 기준 (§4.5) |
| Validation | Q1–Q9 | Q0–Q9 + F1(결함 주입) + G(승인 gate) + T4 (§8) |
| 환경 | Python 3.13.15 | **Python 3.13.16**, numpy 2.5.3, scipy 1.18.1, nibabel 5.4.2 (P1 기록과 동일; R0의 3.13.15는 오기) |

## 1. 목적과 cohort

**질문.** P2A-FO 최종 상태 X_T에서도 critical인 RAW k=2 shell|shell block(C1)과 trajectory-induced k=2 shell|shell block(C2)에 대해, **block 밖 voxel과 reversible SAFE ADD/REM을 허용한 bounded exhaustive 탐색**에서 위상을 보존하며 해당 block을 해소하고 새 critical block을 남기지 않는 경로가 있는가.

**사전 근거 [FIN; `P2A_PARENT_ATTRIBUTION_REPORT.md`].** C1·C2 block은 어떤 REM의 parent 대상이 된 적이 없고 최종 상태에서 모두 `GATE_ALL_UNSAFE`이다. 이는 B의 4 voxel 각각의 단일 flip이 UNSAFE라는 depth-1 사실이며 B 밖 voxel은 평가된 적이 없다. REM 1,562개 중 1,154개(73.9%)가 immediate exchange였으므로 "B가 사라졌다"만으로는 성공이 아니다(§4.4).

| 역할 | cohort | blocks | cases |
|---|---|---|---|
| **Primary** | C1: P1 terminal RAW k=2 shell\|shell 중 P2A-FO 최종에도 critical | **614** | 213 |
| **Comparison** | C2: P1 terminal FINAL k=2 shell\|shell 중 P2A-FO 최종에도 critical | **66** | 39 |

총 **680 block rows**. C1과 C2를 합산한 단일 obstruction cohort로 보고하지 않는다(합계 열 없음). 정의: `P1_TERMINAL_STUCK.csv`의 행이 frozen P2A validator의 `cohort_of`로 C1/C2이고, 같은 case·block이 Run A `P2A_RUN_STUCK.csv`에 존재.

## 2. 입력, provenance, 환경

고정 입력(SHA-256; 불일치 시 exit 3, 과학 출력 없음):
```
bc146d6bfa7d7eda6e0129aefaaf4349985bd48376c1e6e9635d9acacf5450d6  P2A_RUN_CASES.csv   (Run A)
ba801ac55a0f436d5cd06e541bba7f9f048fa30585c888d06cadc78927cfac11  P2A_RUN_EDITS.csv   (Run A)
58fa33dedf19cdb196fb83991d55e0b9699b126d7007a6b6b4cef4b40fb3fa5e  P2A_RUN_STUCK.csv   (Run A)
8097a8e2539f8276985c3259c8b3d093214f47692bbb8eea8bcbd7ec23767b1d  P2A_RUN_META.json   (Run A)
944836057b2f6ec06083698192773c313bdea1401b5125b64ed440662efe4b10  P2A_frozen_R2.py    (P2A manifest 모듈, import)
```
`P2A_frozen_R2.FROZEN`이 추가로 검증하는 13개 pin(P1 코드·CSV, vendor R+ 3파일, `p0_audit.py`, gate table `59407fef…4830`, `cubical_betti.py`, P4 manifest `86d1d964…a99e`, audit CSV)과 P4 manifest 3,512개 항목 전체, 사용 case의 RAW 파일 manifest 포함 여부를 같은 방식으로 시작 시 검증한다. Run B와 parent-attribution 산출물은 입력이 아니다. P4 materialization은 계산 입력 provenance로만 쓰며 `canonical_lineage_ready=false`, `p3_administrative_freeze_attached=false`를 승격하지 않는다. gate table은 regenerated-and-pinned P2 artifact이다.

환경: 전용 venv(Python 3.13.16, numpy 2.5.3, scipy 1.18.1, nibabel 5.4.2), base 환경 미변경. 코드는 `python -P -s -B`로 실행한다. frozen P0/P1/vendor 함수는 import하며 재구현하지 않는다(벡터화한 gate 조회는 §4.3의 감사로 frozen verdict와 대조한다).

## 3. 시작 상태와 Q0

### 3.1 시작 상태 X_T [DEF]
case마다 **pinned P4 RAW → Run A ordered edit replay → reconstructed X_T**로 만든다(저장된 mask 파일을 가정하지 않는다). reconstructed X_T의 mask SHA(`P0.mask_sha`)가 Run A case 행의 `mask_sha`와 **정확히 일치**해야 한다. Run B는 복원·분석에 쓰지 않는다. 680개 block 각각의 탐색은 해당 case의 X_T에서 **독립적으로** 시작한다(§5.3).

### 3.2 Q0: start-state qualification (production 전 필수)
680/680 block 모두에서:
1. reconstructed X_T SHA = Run A `mask_sha`;
2. target BlockKey가 X_T의 vendor `_critical_keys_vectorized`에서 실제로 critical;
3. identity: cohort(C1/C2)가 `cohort_of`와 일치, k=2, P1 terminal `cand_classes=shell|shell`, Run A STUCK 행 `cause=GATE_ALL_UNSAFE`, 그리고 case의 재계산 critical 집합 = Run A STUCK 블록 집합;
4. block의 4 voxel에 대한 단일 ADD/REM 후보가 frozen `P1.gate`에서 **전부 UNSAFE**.
하나라도 실패하면 production search를 시작하지 않는다. 집계 일치(614/66 blocks, 213/39 cases, 222 distinct cases)도 확인한다. validator는 같은 항목을 독립 구현(명시적 loop gate, 자체 replay)으로 재계산해 runner의 Q0 표와 행 단위로 대조한다.

## 4. 탐색 문제의 정의 [DEF/THM]

### 4.1 State
Ω = 원본 volume(1-voxel zero padding은 gate 평가용이며 **flip 대상이 아니다**). window W ⊂ Ω와 depth L이 주어지면 state는 W 안 occupancy의 **정확한 정수 bitmask**이다: bit i = window voxel i(축 (z,y,x) 사전순, 최하위 비트가 사전순 첫 voxel)의 occupancy. W 밖은 X_T 값으로 고정된다.

### 4.2 Window와 halo
bbox(B) = B voxel들의 축정렬 bounding box(k=2에서 2×2×1 방향 plane별). **W_R(B) := Ω ∩ {v : dist_∞(v, bbox(B)) ≤ R}** (Ω 경계에서 clip, 항상 직육면체). k=2에서 R=1: 최대 48 voxel, R=2: 180, R=3: 448. **Halo** = W 바깥 1-voxel 층(padding 포함); 값은 **frozen X_T**에서 가져오고 탐색 중 변하지 않는다. window voxel의 26-이웃은 W ∪ halo 안에 있고, window voxel을 포함하는 block은 W ∪ halo ∩ Ω 안에 있으므로(block의 L∞ 직경 = 1) gate와 critical 판정은 이 crop에서 **정확하다**(근사 아님). 서로 다른 block의 window는 겹쳐도 합치지 않는다.

### 4.3 Move set (동결)
**단일 voxel flip**(v ∈ W, ADD: 0→1 또는 REM: 1→0). 현재 state에서 frozen simple gate가 SAFE일 때만 허용: A(v)=K(X∖v)∩C_v contractible ⇔ (h0,h1,h2)=(1,0,0), 구현은 26-이웃 bitmask(`P0.OFFS` 비트 순서) → 2²⁶ gate table(`bits&7==7`).
- reversible(flip-once 없음); **RAW 전경 voxel 제거 허용**; target block 밖 voxel이라도 window 안이면 허용; **다중 voxel 동시 flip 불허**.
- 모든 accepted transition은 같은 frozen gate를 통과한다. **Gate 감사 [FIN]**: search가 사용한 모든 distinct neighbourhood bitmask를 frozen `P0.attachment`+`P0.verdict`로 재평가해 table 조회와 SAFE 판정이 일치함을 block·tier마다 확인한다(불일치 = ERROR). table의 모든 2²⁶ 항목은 `bits&7 ∈ {0,7}`임을 qualification에서 확인한다.
- Theorem Z (P2A spec §3.2)에 의해 경로 위 모든 state S에 대해 K(X_T) ≃ K(S) ≃ K(X₀)(RAW)이다. flip-once나 window와 무관하다. 모든 PURE 경로의 endpoint는 validator가 from-scratch Betti (β0,β1,β2)로 RAW = X_T = endpoint임을 확인한다.

### 4.4 성공 판정: PURE (end-state, 동결)
C0 = X_T의 전역 critical block 집합, CY = endpoint Y의 전역 critical 집합.
**PURE(Y) ⇔ B ∉ CY ∧ CY ∖ C0 = ∅.**
- 다른 pre-existing critical block이 함께 해소되는 것은 허용되며 기록한다. 경로 **중간**에 임시 critical block이 생기는 것은 허용된다(판정은 endpoint에만 의존).
- **B ∉ CY ∧ CY∖C0 ≠ ∅**("B만 사라지고 새 critical이 endpoint에 있음")는 **diagnostic**이며 PURE 성공이 아니다.
- **국소성 보조정리 [THM].** window 안 voxel을 하나도 포함하지 않는 block은 X_T와 모든 state S에서 같은 값을 가지므로 critical 여부가 같다. 따라서 CY∖C0와 B∈CY는 window+halo crop의 in-domain block만 비교하면 정확히 결정된다. (구현은 crop에서 벡터화하며 vendor 전역 recount와의 차이 0을 qualification에서 확인한다.)
- ELI5: 방 하나(window)만 가구를 옮길 수 있고 방 밖 복도(halo)는 그대로라면, "새로 어지러워진 곳"은 그 방과 맞닿은 칸에서만 생길 수 있어 그곳만 확인하면 된다.

### 4.5 탐색 방법, dedup, cap (동결)
**Exhaustive FIFO BFS, depth ≤ L**(depth = accepted SAFE voxel flip 수). state 정체는 window bitmask 전체(정확 비교; hash 충돌에 의존하지 않음). 이미 발견된 mask는 다시 확장하지 않으므로 어떤 state도 경로에 두 번 나타나지 않는다. 확장 순서는 결정적: 부모는 FIFO 순서, 한 부모의 SAFE voxel은 window 사전순. 따라서 찾은 경로는 (R,L,move set) 안의 **최단 경로**이며 동률은 이 순서로 깨진다. 조기 종료는 **PURE state를 처음 발견했을 때만**이며(그 이후 state는 발견되지 않은 것으로 센다) 가지치기·대칭 제거·우선순위 정렬은 없다. depth-L state는 발견·판정하되 확장하지 않는다.
**Cap.** N_max = 1,000,000 = root state를 포함해 depth 0..L에서 발견된 unique state 전체(depth-L state도 count에 포함하나 expand하지 않음; block·tier마다). 정확히 N_max개까지 발견한 뒤 queue가 정상 소진되면 exhaustive 판정이 가능하고, **N_max+1번째 unique state를 삽입해야 하는 순간 CAP_EXCEEDED**로 중단한다. rank ≤ N_max에서 PURE가 먼저 발견되었으면 PURE_FOUND이다.

### 4.6 종료 상태 (block·tier 단위)
| status | 조건 | 허용되는 주장 |
|---|---|---|
| **PURE_FOUND** | PURE state 발견(최단 경로 길이 L*, flip 열 기록) | (R,L,move set,X_T)에서 PURE 해소 경로가 존재 [FIN]; Theorem Z로 위상 보존 [THM] |
| **B_RESOLVED_WITH_NEW_CRIT** | BFS가 cap 없이 소진되었고 PURE state는 없으나 **B ∉ Crit인 state가 ≥1개** | diagnostic outcome. 소진된 탐색이므로 "이 (X_T,R,L,move set,PURE)에는 PURE 경로가 없다"는 finite fact가 성립한다(분석 flag `PURE_NO_PATH_EXHAUSTIVE`, §6). outcome 이름에는 certificate라는 말을 붙이지 않는다 |
| **NO_PATH_EXHAUSTIVE** | BFS가 cap 없이 depth ≤ L의 모든 state를 소진했고 **B ∉ Crit인 state가 하나도 없음** | **finite certificate**(§6) |
| **CAP_EXCEEDED** | cap을 먼저 넘음 | 미결(unknown/computationally unresolved). certificate 아님 |
| **ERROR** | provenance·replay·accounting·gate 감사 위반, 예외 | 별도 집계. 해당 block은 상위 tier로 가지 않는다 |

## 5. Tier ladder, routing, 독립성

### 5.1 Tier와 routing (config `P2B_CONFIG_R1.json`, SHA는 freeze report에 기록)
| tier | R | L | 대상 |
|---|---|---|---|
| A | 1 | 4 | 모든 block |
| **B (headline)** | **2** | **6** | A에서 PURE_FOUND가 아닌 block (NO_PATH_EXHAUSTIVE, B_RESOLVED_WITH_NEW_CRIT, CAP_EXCEEDED) |
| C (sensitivity) | 3 | 8 | **B에서 PURE_FOUND가 아니고 소진된 block: NO_PATH_EXHAUSTIVE 또는 B_RESOLVED_WITH_NEW_CRIT** |

N_max = 1,000,000 (tier·block마다). A에서 PURE_FOUND이면 상위 tier를 실행하지 않는다. B의 CAP_EXCEEDED는 소진이 아니므로 C로 자동 승격하지 않고 computationally unresolved로 남는다. **Tier-B headline 집계는 A에서 이미 PURE_FOUND인 block도 "PURE reachable by B"로 포함**하며 어느 tier에서 찾았는지 별도 표기한다. ERROR는 승격하지 않는다. (routing은 config 값이며 변경하려면 config SHA가 바뀐다.)

### 5.2 정리 M (monotonicity) [DEF/THM]
시작 state와 성공 판정이 고정되어 있고 (R,L) ≤ (R',L')이면 (R,L)에서 도달 가능한 state 집합은 (R',L')의 부분집합이다. 따라서 **cap이 걸리지 않는 한** (R,L)에서 PURE_FOUND이면 (R',L')에서도 PURE_FOUND이고 L*는 증가하지 않으며, (R',L')에서 PURE가 없으면 더 작은 (R,L)에도 없다; 소진된 두 탐색에서 n_states, B-resolved 존재, 최초 B-resolved depth는 단조이다. *증명:* (R,L)의 경로는 voxel이 W_R ⊂ W_R'이고 길이 ≤ L ≤ L'이며 각 step의 gate는 현재 state에만 의존한다. 성공 판정이 window 상대가 아닌 전역 end-state 함수이므로 성립한다. ∎ **cap 아래에서는 상태(status) 단조성을 주장하지 않는다**(큰 tier가 L*의 depth에 닿기 전에 cap될 수 있음). 이것이 "A의 PURE_FOUND는 B를 실행하지 않고 PURE reachable로 센다"는 규칙의 이유이다. ELI5: 더 큰 방과 더 긴 시간을 주면 이전에 가능했던 길은 여전히 가능하다 — 단, 시간이 다 떨어지기 전에 그 길을 찾아낼 수 있다면.

### 5.3 Block 독립성
680 search는 각각 해당 case의 원래 X_T에서 독립 시작한다. 한 block의 witness를 다른 block의 시작 mask에 적용하지 않고 witness를 합성하지 않는다. 따라서 **block-level PURE_FOUND는 case-level simultaneous repair나 SOLVED를 의미하지 않는다**(같은 case의 block 경로는 서로 간섭할 수 있다).

## 6. Certificate 보고 규칙
**Headline 해석 (고정).** Tier별 `PURE_FOUND`는 bounded topology-preserving repair witness의 존재를 뜻한다. exhaustive search에서 PURE를 찾지 못한 경우에만 해당 (X_T,R,L,move set)에 대한 finite no-PURE certificate를 주장하며, `CAP_EXCEEDED`는 어떠한 부재 주장도 허용하지 않는다. "PURE no-path certified at tier T"의 수는 `NO_PATH_EXHAUSTIVE`와 `B_RESOLVED_WITH_NEW_CRIT`(둘 다 소진·PURE 없음)를 모두 포함해 센다.

NO_PATH_EXHAUSTIVE는 항상 *"start X_T, move set = 단일 voxel SAFE ADD/REM(reversible), window R, depth L, N_max에서 소진, 성공 판정 PURE(및 B 해소 state 자체가 없음)에 대해 exhaustive 탐색에서 해소 경로가 없다"*의 형태로만 쓴다. "불가능하다", "obstruction이다", "same-grid에서 풀 수 없다"는 쓰지 않는다. CAP_EXCEEDED는 unresolved로 별도 집계하며 NO_PATH_EXHAUSTIVE와 합치지 않는다. 한계: 단일 voxel SAFE flip이 만드는 경로에 한정되며 두 voxel이 함께만 simple한 flip(simple pair)은 도달하지 못할 수 있다.

## 7. 보고 지표 (사전 고정; cohort별, C1·C2 합계 열 없음)
1. tier별·headline status 분포(PURE_FOUND / B_RESOLVED_WITH_NEW_CRIT / NO_PATH_EXHAUSTIVE / CAP_EXCEEDED / ERROR); headline = tier-B 시점 PURE reachable(A 포함) + B 최종 status; 어느 tier에서 찾았는지.
2. PURE_FOUND의 L* 분포, path flip 수(ADD/REM), **net 변경 voxel 분류**: RAW 전경 제거 수(`net_raw_fg_removed`), P1/P2A가 추가한 voxel 되돌림 수, 새 배경 voxel 추가 수, P2A가 제거한 voxel 복원 수, B 밖 voxel 수, 해소된 pre-existing critical block 수(B 포함)와 B 이외 수.
3. 탐색량: n_states, level counts, 시간(진단), B-resolved state 수와 최초 depth.
4. tier C 민감도: tier B의 NO_PATH_EXHAUSTIVE에서 C의 결과 분포.
5. C1과 C2를 같은 표에 나란히 두되 합계 열 없음; 차이는 기술 통계로만 비교하고 검정하지 않는다.
6. 필수 한계 표기: 단일 데이터셋(BraTS 876 case WT mask), 이 (R,L,move set), 임상 평가 아님. 경로가 RAW 전경 voxel을 제거할 수 있으므로 제거 수를 반드시 보고하고 위상 보존을 분할 충실도·임상 타당성과 구분한다.

산출물: `P2B_RUN_BLOCKS.csv`(block·tier당 1행), `P2B_RUN_PATHS.csv`(PURE_FOUND 경로의 flip 열과 각 step의 nmask), `P2B_RUN_META.json`, validator의 `P2B_VALIDATION_REPORT.json`. 모든 출력은 write-once.

## 8. Validation contract (production 전 동결)

| ID | 검사 | 합격 기준 |
|---|---|---|
| Q0 | §3.2, runner + 독립 validator 행 단위 대조 | 680/680; 614/66, 213/39, 222 cases |
| Q1 | provenance gate (§2) + 음성: pin 불일치 시 exit 3 | 일치 / exit 3 |
| Q3 | **기존 single-flip null control**: window = B의 voxel(R=0), L=1을 680 block에 적용 | 전부 NO_PATH_EXHAUSTIVE, n_states = 1 (Run A가 B 안 단일 flip이 전부 UNSAFE임을 이미 검증) |
| Q4 | (i) 7³ solid + 대각 cavity 2개: 수기 추론 L*=2, PURE (REM 후 ADD), (ii) L=1 해소 없음, (iii) 무작위 소형 volume ≥200개(경계 clip, R∈{0,1}, L∈{3,4}, ≥25개의 non-PURE 결과 포함)에서 runner = 독립 순수 Python BFS | status, n_states, level counts, B-resolved 통계, 경로가 모두 일치 |
| Q5 | tier 단조성(정리 M): 정렬된 C5 block 50개를 A,B,C로 독립 탐색 | 위반 0 (cap-limited 비교는 별도 집계) |
| Q6 | cap 주입(N_max=25): cap 도달 block은 certificate가 아니라 CAP_EXCEEDED; 독립 BFS와 동일 | 일치 |
| Q7 | **witness 독립 replay**: 모든 PURE_FOUND 경로를 X_T에서 명시적 loop gate로 재생 | 매 step SAFE(table 경로와 attachment 경로 일치), voxel ∈ W, 길이 = L* ≤ L, state 반복 없음, 기록된 nmask 일치, endpoint에서 전역 critical 재계산으로 PURE 재판정, accounting 필드 일치, **from-scratch Betti RAW = X_T = endpoint**, **tier monotonicity**(witness가 더 큰 tier의 window·depth 범위에서도 admissible) |
| Q7b | NO_PATH/B_RESOLVED 행: validator가 결정적 균등 간격 표본을 독립 BFS로 재탐색 | status·n_states·level counts 정확 일치 |
| Q6b | cap 경계: N_max = n(정확히 n개 state를 가진 소진 탐색)이면 NO_PATH_EXHAUSTIVE, N_max = n−1이면 CAP_EXCEEDED(독립 BFS와 동일) | 일치 |
| R | routing: `next_tier` + frozen config가 §5.1 표와 일치(A: NO_PATH/B_RESOLVED/CAP→B, B: NO_PATH/B_RESOLVED→C, CAP·PURE·ERROR→stop) | 일치 |
| Q8 | validator 음성 대조: 경로 변조 8종(unsafe step, window 밖, 길이 초과, 반복 state, 끝 절단, type 변조, lstar 변조, accounting 변조)과 파일 변조 12종(n_states, status, lstar, accounting, voxel 이동, step 삭제, nmask, type, tier 행 삭제, 소진 아닌 B_RESOLVED 행, tier B의 B_RESOLVED가 C로 승격되지 않은 행, meta SHA) | 각각 기대 check id로 FAIL; 무변조 PASS |
| F1 | runner 결함 주입: (a) 사용 중인 gate table 항목 변조 → gate 감사가 검출, (b) dedup 해제 → state 수가 독립 BFS와 불일치 | 검출 |
| G | production 승인 gate: `search --cohort c1c2`를 승인 파일 없이/SHA 불일치로 실행 | exit 3, 출력 없음 |
| Q9 | timing dry-run: **C1/C2가 아닌 C5 block** 10개(정렬 순 처음 10개) | 실행 가능; 값(R,L,N_max)은 이 결과로 바꾸지 않으며 status는 분석에 쓰지 않음 |
| T4 | PYTHONHASHSEED 0 vs 1 재실행(C5 dry-run, 2 shard)과 shard 배치 불변성 | 시간 열을 제외한 BLOCKS/PATHS 일치 |

**실행 형태.** production = **one prespecified analysis run (Run A, PYTHONHASHSEED=0) + one determinism replicate (Run B, seed 1; 과학 수치에 쓰지 않음)**. block은 정렬된 (case_id, block) 순서의 `k::n` write-once shard로 나눠 실행하고(결과는 shard 배치와 무관함을 T4로 확인) merge한다. 기존 shard 출력의 재실행은 거부한다. A와 B가 다르면 T4 FAIL로 기록하고 재실행하지 않는다.
**승인 gate.** production은 `P2B_FREEZE_APPROVAL.json`(runner, config, spec, validator의 SHA-256)이 현재 파일과 일치할 때만 시작한다. spec·runner·validator·config SHA와 qualification 결과는 freeze report에 기록한다. 승인 후 §3–§5의 값과 이 절의 기준을 결과를 보고 바꾸지 않는다.

## 9. 결과 해석 규칙 (사전 고정)

| 결과 | 쓸 수 있는 문장 | 쓰지 않는 문장 |
|---|---|---|
| PURE_FOUND | "tier T에서 길이 L*의 위상 보존 경로가 존재하며 B를 해소하고 endpoint에 새 critical block이 없다"(B 밖·RAW 전경 제거 수 병기) | "이 case가 SOLVED된다", "임상적으로 안전" |
| NO_PATH_EXHAUSTIVE | §6의 finite certificate 문장 | 불가능성, obstruction, "same-grid에서 불가" |
| B_RESOLVED_WITH_NEW_CRIT | "B가 사라지는 state는 있으나 endpoint에 새 critical block이 남는다" | "교환 구조가 원인이다"(P2A 관찰과 같은 종류인지는 별도 분석), "해소됨" |
| CAP_EXCEEDED | "미결" | 어떤 방향의 결론 |
| C1 대 C2 | 기술 통계 | 인과("trajectory가 만들어서"), 검정 |

## 10. 범위 경계
하지 않는 것: **P2A-RZ**(B-local reversible zigzag; P2B 이후 flip-once 효과를 분리하는 ablation 후보), RAW 시작 ablation, 다중 voxel simple set move, 경로 합성과 case 단위 SOLVED 재평가, P2C(lattice refinement), C3/C4/C5 대상 production 탐색(C5는 qualification fixture로만 사용).

## 11. 결정 기록 (사용자 승인 2026-10-07)
1. B_RESOLVED_WITH_NEW_CRIT를 tier C로 승격: **승인**(R1.1에 반영).
2. B_RESOLVED_WITH_NEW_CRIT 이름에 certificate를 붙이지 않되, PURE no-path 집계에서는 두 소진 outcome을 모두 포함: **승인**(§0-A, §6).
3. Cap은 depth-L state를 포함하며 N_max+1번째 unique state에서 CAP_EXCEEDED; tier A의 CAP_EXCEEDED는 B로: **승인**(§4.5, §5.1).
4. P4 pinned ZIP scratch 재추출과 chain 재검증: **승인**.
남은 승인 대기 사항: freeze 승인 자체(R1.1 SHA 확인 후 `P2B_FREEZE_APPROVAL.json` 작성 → production).

## 12. 5인 리뷰어 자체 검토
1. **논리 검증.** PURE를 전역 end-state 함수로 정의해야 정리 M이 성립한다. cap 아래에서는 status 단조성이 깨질 수 있으며 §5.2에 명시했다. 국소성 보조정리는 block 직경 1과 halo 1에 의존한다. 미해결: 단일 voxel move set이 simple pair를 놓치는 한계(§6). B_RESOLVED_WITH_NEW_CRIT의 routing은 승인된 대로 C까지 간다(소진된 탐색만; cap은 승격하지 않음).
2. **데이터 정확성.** cohort 수치(614/66, 213/39)는 R0·attribution 보고서와 Q0에서 재계산한다. 이 spec 작성 시 어떤 block의 search 결과도 보지 않았다. 환경 버전 오기(3.13.15)를 3.13.16으로 정정했다.
3. **인용 적절성.** 새 인용 없음. Theorem Z의 문헌 번호는 P2A spec과 같이 미확인 상태로 두며 확인 전에 새 문헌 번호를 쓰지 않는다.
4. **언어.** 내부 spec이므로 한계 서술을 포함한다. 출판본에서는 방어적 어조를 제거한다.
5. **윤리 경계.** 경로가 RAW 전경(분할 mask)을 제거할 수 있으므로 제거 수를 반드시 보고하고 위상 보존과 임상 타당성을 구분한다. 임상 평가가 아니다.
