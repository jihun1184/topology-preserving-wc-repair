// P0 exhaustive consistency check over all 2^26 occupancy patterns of the 26-neighbourhood of a voxel p.
// Three independent decisions on "is toggling p topology-preserving in the (26,6) closed-cube model?":
//   A) attachment homology: A = K ∩ C(p) is contractible  (h0=1,h1=0,h2=0)       [definition used by the gate]
//   B) digital simple point (Bertrand-Malandain): T26(p,X)=1 and T6bar(p,Xbar)=1   [literature criterion, independent code path]
//   C) greedy elementary collapse of A to a single vertex (free pairs)              [Mazo et al. collapse/free-pair certificate]
// Output: summary on stdout, per-pattern table (uint8: bit0=A safe, bit1=B safe, bit2=C safe, bits4-5 = class) to argv[1].
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <map>
#include <array>
#include <algorithm>
#include <random>
#include <functional>
using u32 = uint32_t; using u8 = uint8_t;

static int offs[26][3];            // neighbour offsets (dz,dy,dx), C order, centre skipped  -> bit k
static u32 closureMask[26];        // 27-bit mask of faces g of C with g ⊆ f(k)
static int dimOf[27];              // dimension of face id (base-3 code, digit 2 = full interval)
static int code[27][3];

static bool subsetFace(int g, int f) {      // g ⊆ f ?  (per axis: f full, or equal endpoint)
  for (int a = 0; a < 3; ++a) if (code[f][a] != 2 && code[g][a] != code[f][a]) return false;
  return true;
}
struct Info { int V, E, S, h0, h1, h2; bool contractible, collapsible, orderIndep; };

static Info analyse(u32 A /*27-bit, bit26 (the cube) must be 0*/, std::mt19937& rng) {
  Info r{}; std::vector<int> verts, edges;
  for (int g = 0; g < 26; ++g) if (A >> g & 1) { if (dimOf[g] == 0) r.V++; else if (dimOf[g] == 1) r.E++; else r.S++; }
  // h0: components of the 1-skeleton (vertices + edges); an edge connects its two endpoint vertices
  int par[27]; for (int i = 0; i < 27; ++i) par[i] = i;
  std::function<int(int)> fnd = [&](int x) { return par[x] == x ? x : par[x] = fnd(par[x]); };
  int comps = r.V;
  for (int g = 0; g < 26; ++g) if ((A >> g & 1) && dimOf[g] == 1) {
    int ax = -1; for (int a = 0; a < 3; ++a) if (code[g][a] == 2) ax = a;
    int v0 = g - code[g][ax] * 0, c[3] = {code[g][0], code[g][1], code[g][2]};
    c[ax] = 0; int u = 9 * c[0] + 3 * c[1] + c[2]; c[ax] = 1; int w = 9 * c[0] + 3 * c[1] + c[2]; (void)v0;
    if (!(A >> u & 1) || !(A >> w & 1)) { fprintf(stderr, "closure violated\n"); exit(2); }
    int ru = fnd(u), rw = fnd(w); if (ru != rw) { par[ru] = rw; --comps; }
  }
  r.h0 = comps; r.h2 = (r.S == 6); r.h1 = r.h0 + r.h2 - (r.V - r.E + r.S);
  if (r.h1 < 0) { fprintf(stderr, "negative h1\n"); exit(2); }
  r.contractible = (r.h0 == 1 && r.h1 == 0 && r.h2 == 0);
  // greedy collapse, several random orders; collapsible iff ends in exactly one vertex; also check order independence
  bool first = true, res0 = false, indep = true;
  for (int trial = 0; trial < 6; ++trial) {
    u32 S = A;
    for (;;) {
      std::vector<std::pair<int,int>> pairs;
      for (int x = 0; x < 26; ++x) if (S >> x & 1) {
        int cover = -1, cnt = 0;
        for (int y = 0; y < 26; ++y) if (y != x && (S >> y & 1) && subsetFace(x, y)) { ++cnt; cover = y; }
        if (cnt == 1 && dimOf[cover] == dimOf[x] + 1) pairs.push_back({x, cover});
      }
      if (pairs.empty()) break;
      auto pr = trial == 0 ? pairs.front() : pairs[rng() % pairs.size()];
      S &= ~(1u << pr.first); S &= ~(1u << pr.second);
    }
    bool ok = (__builtin_popcount(S) == 1) && dimOf[__builtin_ctz(S)] == 0;
    if (first) { res0 = ok; first = false; } else if (ok != res0) indep = false;
  }
  r.collapsible = res0; r.orderIndep = indep; return r;
}

int main(int argc, char** argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s table.bin\n", argv[0]); return 1; }
  for (int id = 0; id < 27; ++id) { int t = id; for (int a = 2; a >= 0; --a) { code[id][a] = t % 3; t /= 3; }
    dimOf[id] = (code[id][0] == 2) + (code[id][1] == 2) + (code[id][2] == 2); }
  int k = 0; for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
    if (!dz && !dy && !dx) continue; offs[k][0] = dz; offs[k][1] = dy; offs[k][2] = dx; ++k; }
  for (int n = 0; n < 26; ++n) {              // face of C shared with neighbour n: offset -1 -> low(0), +1 -> high(1), 0 -> full(2)
    int f = 0; for (int a = 0; a < 3; ++a) f = f * 3 + (offs[n][a] == 0 ? 2 : (offs[n][a] < 0 ? 0 : 1));
    closureMask[n] = 0; for (int g = 0; g < 26; ++g) if (subsetFace(g, f)) closureMask[n] |= 1u << g; }
  // 13-bit chunk tables for A(mask)
  std::vector<u32> lo(1 << 13), hi(1 << 13);
  for (u32 m = 0; m < (1u << 13); ++m) { u32 a = 0, b = 0; for (int i = 0; i < 13; ++i) if (m >> i & 1) { a |= closureMask[i]; b |= closureMask[13 + i]; } lo[m] = a; hi[m] = b; }
  // BM tables: adjacency among the 26 neighbours
  u32 adj26[26] = {}, adj6[26] = {}, N18 = 0, F6 = 0;
  for (int i = 0; i < 26; ++i) {
    int nz = (offs[i][0] != 0) + (offs[i][1] != 0) + (offs[i][2] != 0);
    if (nz <= 2) N18 |= 1u << i; if (nz == 1) F6 |= 1u << i;
    for (int j = 0; j < 26; ++j) if (i != j) { int d = 0, mx = 0; for (int a = 0; a < 3; ++a) { int e = abs(offs[i][a] - offs[j][a]); d += e; mx = std::max(mx, e); }
      if (mx <= 1) adj26[i] |= 1u << j; if (d == 1) adj6[i] |= 1u << j; } }
  auto flood = [&](u32 set, const u32* adj) { // returns list of connected component masks of 'set'
    std::array<u32, 26> out; int n = 0; u32 rest = set;
    while (rest) { u32 comp = rest & (~rest + 1), prev = 0;
      while (comp != prev) { prev = comp; for (u32 t = comp; t; t &= t - 1) comp |= adj[__builtin_ctz(t)] & set; }
      out[n++] = comp; rest &= ~comp; }
    return std::vector<u32>(out.begin(), out.begin() + n); };

  std::mt19937 rng(12345); std::map<u32, Info> memo; std::vector<u8> table(1u << 26);
  uint64_t cnt[3] = {0, 0, 0}, dis = 0, nonIndep = 0, bothAC = 0; std::map<std::array<int,3>, uint64_t> klass; std::map<std::array<int,3>, uint64_t> klassSafeBM;
  for (u32 m = 0; m < (1u << 26); ++m) {
    u32 A = lo[m & 0x1FFF] | hi[m >> 13];
    auto it = memo.find(A); if (it == memo.end()) it = memo.emplace(A, analyse(A, rng)).first; const Info& I = it->second;
    // B) Bertrand–Malandain: T26 of occupied neighbours; T6bar = # of 6-components of background∩N18 that touch a face-neighbour of p
    int T26 = (int)flood(m, adj26).size();
    u32 bg = (~m) & ((1u << 26) - 1) & N18; int T6b = 0; for (u32 c : flood(bg, adj6)) if (c & F6) ++T6b;
    bool B = (T26 == 1 && T6b == 1);
    u8 v = (I.contractible ? 1 : 0) | (B ? 2 : 0) | (I.collapsible ? 4 : 0); table[m] = v;
    cnt[0] += I.contractible; cnt[1] += B; cnt[2] += I.collapsible;
    if ((I.contractible != B) || (I.contractible != I.collapsible)) ++dis;
    klass[{I.h0, I.h1, I.h2}]++; if (B) klassSafeBM[{I.h0, I.h1, I.h2}]++;
  }
  for (auto& kv : memo) if (!kv.second.orderIndep) ++nonIndep;
  printf("patterns=%u distinct_attachments=%zu\n", 1u << 26, memo.size());
  printf("safe_by_A_homology=%llu safe_by_BM=%llu safe_by_collapse=%llu disagreements=%llu collapse_order_dependent_attachments=%llu\n",
         (unsigned long long)cnt[0], (unsigned long long)cnt[1], (unsigned long long)cnt[2], (unsigned long long)dis, (unsigned long long)nonIndep);
  printf("class (h0,h1,h2): patterns | of which BM-safe\n");
  for (auto& kv : klass) printf("  (%d,%d,%d): %llu | %llu\n", kv.first[0], kv.first[1], kv.first[2], (unsigned long long)kv.second, (unsigned long long)klassSafeBM[kv.first]);
  FILE* f = fopen(argv[1], "wb"); if (!f || fwrite(table.data(), 1, table.size(), f) != table.size()) { fprintf(stderr, "write failed\n"); return 3; } fclose(f);
  return dis ? 4 : 0;
}
