// ============================================================================
// sex_bias_sim_v7.cpp
//
// Standalone C++ port of the SPPARKS "rxn/diff/custom" sex-bias tumour model
// (G. Giuliani; SPPARKS_Sex_Bias, branch Sex_bias_post-Zihai_with_CD4s) for
// the AR / IFNGR manuscript. Model and rates are exactly those of
// rsc_out_2026-07-03/in.full (= the Jul-2026 "sex_bias_model_rates" table):
//
//   app_style rxn/diff/custom 600 600 4 2 1 1 1 1 10000 20.0 0.015 0.1 0.0
//     0.0016667 1 1 0.6 0.6 0 0.03 0.2 0 14 1.0 0 4 4 2.0 1000 1000 1000 0
//     0.25 1 100 0.5     (60 reactions, seed from command line)
//
// TIME UNIT: 1 native unit = 600 s = 10 min (IFNG step 0.015 units = 9 s).
//   1 day = 144 units; day 14 -> day 21 = 1008 units (in.full "run 1008").
//   Rate constants below are per 10 min.
//
// Changes from sex_bias_sim_v6.cpp:
//   * fixed the stray "return bias_factor; }" that stopped v6 compiling;
//   * removed the non-model tuning options (--rates paper, --hill-rescale,
//     --pmax-multiplier, --carrying-capacity, --ifng-decay, v7 "modes");
//   * added --no-sequestration (AR+ cells keep everything except the
//     internalisation of bound IFNG) for the sequestration test;
//   * optional --te-death-weekly (rates-table value; in.full uses 1/day);
//   * command line in days; output at exact requested times without
//     consuming random numbers (adding outputs never changes a trajectory);
//   * speed-ups that leave trajectories bit-identical to v6/SPPARKS order
//     (zero-rate reactions skipped, empty sites skipped, O(n) tree reload).
//
// Female runs: --sex female applies the female rate changes (AR+ recruitment
// and Tp->AR+ PE set to 0, AR- Tp/Tm recruitment x1/0.63, Tp->Tm x2). Use the
// d14_M_4_madeF init (male tumour, AR+ cells relabelled AR-) to isolate AR,
// or d14_F_4 for the real female slide.
//
// Build: g++ -O3 -std=c++17 sex_bias_sim_v7.cpp -o sex_bias_sim_v7
//   (do not add -march=native / -ffast-math: they change floating-point
//    rounding and therefore the stochastic trajectory for a given seed)
// ============================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdint>

// ============================================================================
// COMPILE-TIME CONSTANTS matching in.full config
// ============================================================================

static const int    LATCOL           = 600;             // arg[1]
static const int    LATROW           = 600;             // arg[2]
static const int    NSITES           = LATCOL * LATROW; // 360000
static const int    MAX_SPECIES      = 18;              // from app_rxn_diff.h line 27
static const int    N_TRACKED        = 11;              // slots in ifng_tracker
static const int    IFNG_COL         = 200;             // = ceil(600/3)
static const int    IFNG_ROW         = 200;
static const int    IFNG_NSITES      = IFNG_COL * IFNG_ROW;
static const int    IFNG_LATTICE_CELL_CONVERTER = 3;
static const int    LATT_CONST       = 10;              // um
static const int    IFNG_LATT_CONST  = LATT_CONST * IFNG_LATTICE_CELL_CONVERTER; // 30 um
static const int    CP_MAX           = 12;              // cancer bias potential grid
static const int    CP_MAY           = 12;
static const int    CANCER_POTENTIAL_CONVERTER = 50;    // 600/12
static const int    N_REACTIONS      = 60;              // rxns 0..59
static const int    NEIGH_LIM        = 3;               // cancer division search depth
static const int    WIDTH_LIM        = 3 + 2*(NEIGH_LIM-1); // = 7
static const double TIME_CONVERSION  = 600.0;           // seconds per ifng_time_step unit

// Species IDs (add_species order in in.full)
enum : int {
  canID   = 0,  // A cancer
  TpID    = 1,  // B progenitor CD8
  TmID    = 2,  // C mature CD8
  TeID    = 3,  // D exhausted CD8
  ARTpID  = 4,  // E "GONE" (AR+ progenitor) - always 0
  ARTmID  = 5,  // F AR+ mature CD8
  ARTeID  = 6,  // G AR+ exhausted CD8
  recID   = 7,  // I recruitment marker
  errID   = 8,  // J error signal
  deadID  = 9,  // K dead cells
  latEID  = 10, // L lattice-edge marker
  boundID = 11, // M IFNG-diffusable region marker
  ThID    = 12, // O AR- CD4+ helper
  ARThID  = 13, // P AR+ CD4+ helper
  TrID    = 14, // Q AR- Treg
  ARTrID  = 15  // R AR+ Treg
};

// Reaction styles
enum { LOCAL = 0, NBR = 1 };

// Sex
enum Sex { MALE = 0, FEMALE = 1 };

// ============================================================================
// SECTION 1: RANDOM NUMBER GENERATORS
// Port of random_mars.cpp + random_park.cpp
// ============================================================================

// --- Park/Miller (random_park.cpp) ---
// Same 16807/2147483647 stream used by SPPARKS.
class RandomPark {
 public:
  int seed;

  static const int IA = 16807;
  static const int IM = 2147483647;
  static const int IQ = 127773;
  static const int IR = 2836;
  static constexpr double AM = 1.0 / 2147483647.0;

  explicit RandomPark(int iseed) { seed = iseed; }
  explicit RandomPark(double rseed) {
    // Same as SPPARKS ctor: seed = (int)(rseed * IM); if (seed==0) seed=1.
    seed = static_cast<int>(rseed * IM);
    if (seed == 0) seed = 1;
  }

  // Same as SPPARKS reset(rseed, offset, warmup):
  //   seed = (int)fmod(rseed*IM + offset, IM)
  //   if seed < 0: seed = -seed
  //   if seed == 0: seed = 1
  //   warmup calls
  void reset(double rseed, int offset, int warmup) {
    seed = static_cast<int>(std::fmod(rseed * (double)IM + (double)offset,
                                      (double)IM));
    if (seed < 0) seed = -seed;
    if (seed == 0) seed = 1;
    for (int i = 0; i < warmup; i++) uniform();
  }

  double uniform() {
    int k = seed / IQ;
    seed = IA * (seed - k * IQ) - IR * k;
    if (seed < 0) seed += IM;
    return AM * seed;
  }

  int irandom(int n) {
    int i = (int)(uniform() * n) + 1;
    if (i > n) i = n;
    return i;
  }
};

// --- Marsaglia (random_mars.cpp) ---
// Same algorithm and init as SPPARKS. Used to seed the two Park/Miller RNGs.
class RanMars {
  std::vector<double> u;  // size 97+1 to match SPPARKS's 1-indexed array
  double c, cd, cm;
  int i97, j97;
  bool initflag;
 public:
  RanMars() : c(0), cd(0), cm(0), i97(0), j97(0), initflag(false) {}

  void init(int seed) {
    initflag = true;
    while (seed > 900000000) seed -= 900000000;
    u.assign(97 + 1, 0.0);
    int ij = (seed - 1) / 30082;
    int kl = (seed - 1) - 30082 * ij;
    int i = (ij / 177) % 177 + 2;
    int j = ij % 177 + 2;
    int k = (kl / 169) % 178 + 1;
    int l = kl % 169;
    for (int ii = 1; ii <= 97; ii++) {
      double s = 0.0, t = 0.5;
      for (int jj = 1; jj <= 24; jj++) {
        int m = ((i * j) % 179) * k % 179;
        i = j;
        j = k;
        k = m;
        l = (53 * l + 1) % 169;
        if ((l * m) % 64 >= 32) s = s + t;
        t = 0.5 * t;
      }
      u[ii] = s;
    }
    c  = 362436.0  / 16777216.0;
    cd = 7654321.0 / 16777216.0;
    cm = 16777213.0 / 16777216.0;
    i97 = 97;
    j97 = 33;
    uniform(); // SPPARKS does one warmup call at end of init
  }

  double uniform() {
    double uni = u[i97] - u[j97];
    if (uni < 0.0) uni += 1.0;
    u[i97] = uni;
    i97--; if (i97 == 0) i97 = 97;
    j97--; if (j97 == 0) j97 = 97;
    c -= cd;
    if (c < 0.0) c += cm;
    uni -= c;
    if (uni < 0.0) uni += 1.0;
    return uni;
  }
};

// ============================================================================
// SECTION 2: SOLVE TREE (binary heap KMC event picker)
// Port of solve_tree.cpp
// ============================================================================

class SolveTree {
 public:
  RandomPark* random;    // owned; separate stream from ranapp
  std::vector<double> tree;
  int ntotal, offset;    // ntotal = 2*nround-1; offset = nround-1
  double sum;
  int num_active;

  SolveTree() : random(nullptr), ntotal(0), offset(0), sum(0.0), num_active(0) {}
  ~SolveTree() { delete random; }

  void init(int n, const double* propensity) {
    // Same as solve_tree.cpp init():
    //   m such that 2^m >= n
    //   offset = 2^m - 1, ntotal = 2*(2^m) - 1
    int m = 0;
    long long nround = 1;
    while (nround < n) { nround *= 2; m++; }
    (void)m;
    offset = (int)(nround - 1);
    ntotal = (int)(2 * nround - 1);
    tree.assign(ntotal, 0.0);
    for (int i = offset; i < offset + n; i++) tree[i] = propensity[i - offset];
    sum_tree();
  }

  void sum_tree() {
    for (int parent = offset - 1; parent >= 0; parent--) {
      int c1 = 2*parent + 1;
      int c2 = 2*parent + 2;
      tree[parent] = tree[c1] + tree[c2];
    }
    sum = tree[0];
    num_active = 0;
    for (int i = offset; i < ntotal; i++) if (tree[i] > 0.0) num_active++;
  }

  void set(int i, double value) {
    if (tree[offset + i] > 0.0) num_active--;
    if (value > 0.0) num_active++;
    tree[offset + i] = value;
    int idx = i + offset;
    while (idx > 0) {
      int sibling = (idx % 2) ? idx + 1 : idx - 1;
      int parent = (idx - 1) / 2;
      tree[parent] = tree[idx] + tree[sibling];
      idx = parent;
    }
    sum = tree[0];
  }

  // Update sites given a list of newly-changed indices.
  // Matches solve_tree::update(int n, int* indices, double* propensity) called
  // by AppRxnDiff via solve->update(nsites, esites, propensity) at end of
  // site_event. Also called after full-refresh in iterate_kmc_global.
  // Rewrite every leaf then re-sum bottom-up. Gives exactly the same tree as
  // calling set() on leaves 0..n-1 in order (each parent's final value is
  // the sum of its two final children), in O(n) instead of O(n log n).
  void reload_all(int n, const double* propensity) {
    for (int i = 0; i < n; i++) tree[offset + i] = propensity[i];
    sum_tree();
  }

  void update(int n, const int* indices, const double* propensity) {
    for (int k = 0; k < n; k++) set(indices[k], propensity[indices[k]]);
  }

  // Sample: draw an event index and a time step. Matches solve_tree::event().
  //   r2 = random->uniform()
  //   m  = find(r2 * sum)
  //   *pdt = -1/sum * log(random->uniform())
  int event(double* pdt) {
    if (sum == 0.0) { *pdt = 0.0; return -1; }
    double r2 = random->uniform();
    int m = find(r2 * sum);
    *pdt = -1.0 / sum * std::log(random->uniform());
    return m;
  }

  int find(double value) {
    int i = 0;
    while (i < offset) {
      int lc = 2*i + 1;
      if (value <= tree[lc]) i = lc;
      else { value -= tree[lc]; i = lc + 1; }
    }
    return i - offset;
  }
};

// ============================================================================
// SECTION 3: GLOBAL STATE (mirrors AppRxnDiff / AppLattice member data)
// ============================================================================

// --- Parameters (all 37 args parsed identically to constructors) ---
struct Params {
  // From AppRxnDiff ctor (app_rxn_diff.cpp lines 74-102):
  int    lattice_width;             // 600
  int    lattice_height;            // 600
  int    max_cell_area_per_voxel;   // 4
  int    cancer_area_count;         // 2
  int    Tp_area_count;             // 1
  int    Tm_area_count;             // 1
  int    Te_area_count;             // 1
  int    max_T_cells;               // 10000 (from arg[9])
  bool   update_all_always;         // false (arg[22]=0)

  // From AppLattice ctor (app_lattice.cpp lines 128-166):
  double ifng_diffusion;                 // 20.0    arg[10]
  double ifng_time_step;                 // 0.015   arg[11]
  double ifng_production;                // 0.1     arg[12]
  double ifng_decay;                     // 0.0     arg[13]
  double bIFNG_decay;                    // 0.0016667 arg[14]
  bool   ifng_cancer_proliferation_limitation;  // true arg[15]
  bool   ifng_cytotoxicity_augmentation;        // true arg[16]
  double ifng_ART_uptake;                // 14      arg[23]
  int    ifng_feedback_hill_degree;      // 4       arg[27]
  double max_ifng_production;            // 2.0     arg[28]
  double half_max_ifng_feedback;         // 0.2     arg[21]
  int    receptors_per_cell[3];          // {1000, 1000, 1000} arg[29..31]
  bool   periodic_boundary_conditions;   // false   arg[19]
  bool   use_simple_decay_model;         // false   arg[32]=0
  bool   use_leaky_receptor_model;       // true    (derived from use_simple==false)
  bool   use_cell_dependent_decay_model; // false   (only true if use_simple_decay_model)
  double ifng_binding;                   // 1       arg[34]
  double ifng_unbinding;                 // 100     arg[35]
  double dissC;                          // 100     = unbinding/binding
  int    do_carrying_capacity;           // 1       arg[8] (kept as int for the formula)

  // From AppRxnDiffCustom ctor (app_rxn_diff_custom.cpp lines 34-53):
  double ifng_half_max_proliferation;                       // 0.6  arg[17]
  double ifng_half_max_cytotoxic_augmentation_threshold;    // 0.6  arg[18]
  double tcell_stick_fac;                                   // 0.03 arg[20]
  double modulate_ifng_factor;                              // 1.0  arg[24]
  bool   meanIFNG;                                          // false arg[25]
  int    hill_degree;                                       // 4    arg[26]
  double cancer_proliferation_minimum_factor;               // 0.25 arg[33]
  double ifng_half_max_diffusion;                           // 0.5  arg[36]

  // ---- v7 run options (none of these change the in.full model unless set) ----
  // Ablation used for the sequestration test: AR+ cells keep everything
  // else (no IFNG secretion, same binding) but do NOT internalize bound IFNG.
  bool no_sequestration = false;
  // Exhausted CD8 death rate in native units (per 10 min). in.full uses
  // 0.00694 (= 1/day). The Jul-2026 rates table lists 1/week (0.000992).
  double te_death_rate = 0.00694;

  void init_from_infull() {
    lattice_width = 600;
    lattice_height = 600;
    max_cell_area_per_voxel = 4;
    cancer_area_count = 2;
    Tp_area_count = 1;
    Tm_area_count = 1;
    Te_area_count = 1;
    do_carrying_capacity = 1;
    max_T_cells = 10000;

    ifng_diffusion = 20.0;
    ifng_time_step = 0.015;
    ifng_production = 0.1;
    ifng_decay = 0.0;
    bIFNG_decay = 0.0016667;

    ifng_cancer_proliferation_limitation = true;
    ifng_cytotoxicity_augmentation = true;
    ifng_half_max_proliferation = 0.6;
    ifng_half_max_cytotoxic_augmentation_threshold = 0.6;

    periodic_boundary_conditions = false;
    tcell_stick_fac = 0.03;
    half_max_ifng_feedback = 0.2;

    update_all_always = false;
    ifng_ART_uptake = 14.0;
    modulate_ifng_factor = 1.0;
    meanIFNG = false;
    hill_degree = 4;
    ifng_feedback_hill_degree = 4;
    max_ifng_production = 2.0;

    receptors_per_cell[0] = 1000;
    receptors_per_cell[1] = 1000;
    receptors_per_cell[2] = 1000;

    use_simple_decay_model = false;
    if (!use_simple_decay_model) use_leaky_receptor_model = true;
    use_cell_dependent_decay_model = false;
    if (use_simple_decay_model) use_cell_dependent_decay_model = true;

    cancer_proliferation_minimum_factor = 0.25;
    ifng_binding = 1.0;
    ifng_unbinding = 100.0;
    dissC = ifng_unbinding / ifng_binding;
    ifng_half_max_diffusion = 0.5;

  }
};

// --- Reaction table (populated from in.full's add_rxn lines) ---
struct Reaction {
  int localReactants[MAX_SPECIES];
  int nbrReactants[MAX_SPECIES];
  int localDeltaPop[MAX_SPECIES];
  int nbrDeltaPop[MAX_SPECIES];
  double rate;
  int rxnStyle; // LOCAL or NBR
};

// --- Full simulation state ---
struct Sim {
  // Reference to params (owned externally)
  Params p;
  Sex    sex;

  // Species populations: population[species][site]
  int population[MAX_SPECIES][NSITES];

  // Continuous fields (SPPARKS's ifng_field_discretized[0..1][site])
  // These are the dumped/reported values; the "real" state lives in
  // ifng_field (coarse) and ifng_tracker (per-cell).
  std::vector<double> ifng_field_discretized_free;  // d1: free IFNG per fine voxel
  std::vector<double> ifng_field_discretized_bound; // d2: bound-IFNG-fraction

  // Coarse IFNG field (200x200)
  std::vector<double> ifng_field;       // size IFNG_NSITES
  std::vector<double> ifng_field_copy;  // buffer

  // Per-cell tracker: ifng_tracker[site][kind] = vector of {bound, receptors}
  // kind indexes 0..10 in the order {canID, TpID, TmID, TeID, ARTpID, ARTmID,
  // ARTeID, ThID, ARThID, TrID, ARTrID}. This matches SPPARKS's cellIDs[11]
  // order (see app_lattice.cpp line 758 and app_rxn_diff_custom.cpp line 104).
  //
  // Layout: cell_bound[site * N_TRACKED + kind] = vector of pair (bound, cap)
  //
  // Stored as parallel int vectors (bound count halved/moved integer-count-wise
  // in some places) — but the source uses doubles, so we do too.
  std::vector<std::vector<std::array<double, 2>>> cell_bound; // size NSITES * N_TRACKED

  // Aggregated bound_ifng per fine voxel: [0]=frac, [1]=totalBound, [2]=totalRec
  std::vector<std::array<double, 3>> bound_ifng; // size NSITES

  // Cancer bias potential (12x12, per axis) and discretized cancer counts (12x12)
  // Sizes match app_rxn_diff.h declarations (200x200 alloc, only 12x12 used).
  double cancer_bias_potential[CP_MAX][CP_MAY][2];
  int    discretized_cancer[CP_MAX][CP_MAY];

  // Event list for KMC (mirrors SPPARKS's linked-list events per site)
  struct Event {
    int style, which, jpartner, next;
    double propensity;
  };
  std::vector<Event> events;
  std::vector<int>   firstevent; // size NSITES
  int freeevent;
  int nevents_active;

  // Per-site aggregate propensity (matches SPPARKS's propensity[i2site[i]])
  std::vector<double> propensity;
  // For simple single-proc case, i2site[i]=i.
  std::vector<int>    i2site;    // size NSITES (identity)
  std::vector<int>    esites;    // scratch for touched-site list
  std::vector<char>   echeck;    // dedup marker for touched sites

  // Neighbor list: neighbor[i][0..3]
  std::vector<std::array<int, 4>> neighbor;

  // Reactions (60 slots)
  std::vector<Reaction> rxn;
  std::vector<long long> rxn_count;

  // Timing / flags (mirror AppLattice)
  double time;
  double ifng_time_keeper;
  int    tcell_total_population;
  bool   proliferation_update_bool;
  bool   cancer_proliferation_update_bool;
  bool   mapped_population_to_ifng_tracker;

  // Diagnostic tracking of per-cell bound-receptor fractions observed
  // during IFNG updates. Written after each IFNG update. Used to check
  // whether the Hill function K values are in the right regime -- if
  // max_bound_frac stays << K, suppression can never activate.
  double diag_max_bound_frac = 0.0;
  double diag_p95_bound_frac = 0.0;
  double diag_mean_bound_frac = 0.0;
  long long diag_ncells = 0;

  // Solver + RNGs
  SolveTree*   solve;
  RandomPark*  ranapp;  // second Park/Miller (site_event randomness)

  // For cancer division (need to communicate across site_event helpers)
  int cancerPos1;
  std::vector<int> tcellPos;

  // For recruitment (need to communicate recSpot to touched-site update)
  int recSpot;

  Sim() : sex(MALE), time(0.0), ifng_time_keeper(0.0), tcell_total_population(0),
          proliferation_update_bool(false), cancer_proliferation_update_bool(false),
          mapped_population_to_ifng_tracker(false),
          solve(nullptr), ranapp(nullptr),
          cancerPos1(-1), recSpot(-1) {
    std::memset(population, 0, sizeof(population));
    std::memset(cancer_bias_potential, 0, sizeof(cancer_bias_potential));
    std::memset(discretized_cancer, 0, sizeof(discretized_cancer));
    freeevent = 0;
    nevents_active = 0;
  }

  ~Sim() {
    delete solve;
    delete ranapp;
  }

  // -------------------------------------------------------------------------
  // NEIGHBOR ENUMERATION for sq/4n periodic BC.
  // Order from create_sites.cpp offsets_2d (i=-1..1, j=-1..1), which for
  // cut=1 gives neighbors at (-1,0), (0,-1), (0,1), (1,0). For a site
  // (x,y) with periodic wrap, this yields:
  //   nbr[0] = (x-1,y)   -> j = i - 1
  //   nbr[1] = (x, y-1)  -> j = i - LATCOL
  //   nbr[2] = (x, y+1)  -> j = i + LATCOL
  //   nbr[3] = (x+1,y)   -> j = i + 1
  // Confirmed by grep of create_sites.cpp lines 1181-1210.
  // -------------------------------------------------------------------------
  void build_neighbors() {
    neighbor.assign(NSITES, {0,0,0,0});
    for (int y = 0; y < LATROW; y++) {
      for (int x = 0; x < LATCOL; x++) {
        int i = x + LATCOL * y;
        int xm = (x - 1 + LATCOL) % LATCOL;
        int xp = (x + 1) % LATCOL;
        int ym = (y - 1 + LATROW) % LATROW;
        int yp = (y + 1) % LATROW;
        neighbor[i][0] = xm + LATCOL * y;
        neighbor[i][1] = x  + LATCOL * ym;
        neighbor[i][2] = x  + LATCOL * yp;
        neighbor[i][3] = xp + LATCOL * y;
      }
    }
  }

  // -------------------------------------------------------------------------
  // clear_events(i): identical to AppRxnDiff::clear_events (line 2141)
  // -------------------------------------------------------------------------
  void clear_events(int i) {
    int idx = firstevent[i];
    while (idx >= 0) {
      int next = events[idx].next;
      events[idx].next = freeevent;
      freeevent = idx;
      nevents_active--;
      idx = next;
    }
    firstevent[i] = -1;
  }

  // -------------------------------------------------------------------------
  // add_event: identical to AppRxnDiff::add_event (line 2160)
  // -------------------------------------------------------------------------
  void add_event(int i, int rstyle, int which, double prop, int jpartner) {
    if (freeevent < 0 || freeevent >= (int)events.size()) {
      int old = (int)events.size();
      int add = 100000; // DELTAEVENT
      events.resize(old + add);
      for (int m = old; m < old + add; m++) {
        events[m].next = m + 1;
      }
      // Link into free list: SPPARKS sets freeevent = nevents (i.e. old size)
      // after realloc AND makes the last element's next point past-end. In our
      // model, if freeevent was -1 (empty pool) we bootstrap it now.
      if (freeevent < 0) freeevent = old;
      else {
        // Walk the free chain to its end and splice.
        int cur = freeevent;
        while (events[cur].next < old || events[cur].next >= (int)events.size()) {
          if (events[cur].next < 0) break;
          if (events[cur].next >= (int)events.size()) { events[cur].next = old; break; }
          cur = events[cur].next;
        }
        events[cur].next = old;
      }
      events.back().next = -1;
    }
    int slot = freeevent;
    int next = events[slot].next;
    events[slot].style = rstyle;
    events[slot].which = which;
    events[slot].jpartner = jpartner;
    events[slot].propensity = prop;
    events[slot].next = firstevent[i];
    firstevent[i] = slot;
    freeevent = next;
    nevents_active++;
  }

  // -------------------------------------------------------------------------
  // filling helper: total occupied area at a voxel.
  // Matches formula used throughout SPPARKS (e.g. line 472).
  // -------------------------------------------------------------------------
  int filling(int site) const {
    const int cancA = p.cancer_area_count;
    const int tpA = p.Tp_area_count;
    const int tmA = p.Tm_area_count;
    const int teA = p.Te_area_count;
    const int thA = tpA; // "int thA = Tp_area_count;" line 431
    const int trA = thA;
    return cancA * population[canID][site]
         + tpA   * population[TpID][site]
         + tmA   * population[TmID][site]
         + teA   * population[TeID][site]
         + tpA   * population[ARTpID][site]
         + tmA   * population[ARTmID][site]
         + teA   * population[ARTeID][site]
         + thA   * population[ThID][site]
         + thA   * population[ARThID][site]
         + trA   * population[TrID][site]
         + trA   * population[ARTrID][site];
  }

  // Non-push filling: only cancer cells (for the cancer division "popos"
  // decision at line 473).
  int non_push_filling(int site) const {
    return p.cancer_area_count * population[canID][site];
  }
};

// ============================================================================
// SECTION 4: REACTION TABLE
// Populated from in.full lines 56-128 (60 reactions total).
// Every reaction is defined even if rate=0, because `which` is used as a
// raw index into cd8map arrays throughout site_event.
// ============================================================================

// Helper: add a reaction to the table given local/nbr reactant and product
// species lists (multiplicities allowed via repeated species).
// Sets rxnStyle = LOCAL if both nbrReactants and nbrProducts empty, else NBR.
static void set_reaction(Reaction& rx, double rate,
                         std::initializer_list<int> localReact,
                         std::initializer_list<int> nbrReact,
                         std::initializer_list<int> localProd,
                         std::initializer_list<int> nbrProd) {
  std::memset(&rx, 0, sizeof(rx));
  rx.rate = rate;
  int nNbrR = 0, nNbrP = 0;
  for (int s : localReact) { rx.localReactants[s]++; rx.localDeltaPop[s]--; }
  for (int s : nbrReact)   { rx.nbrReactants[s]++;   rx.nbrDeltaPop[s]--; nNbrR++; }
  for (int s : localProd)  { rx.localDeltaPop[s]++; }
  for (int s : nbrProd)    { rx.nbrDeltaPop[s]++; nNbrP++; }
  rx.rxnStyle = (nNbrR == 0 && nNbrP == 0) ? LOCAL : NBR;
}

static void build_reaction_table(Sim& sim, Sex sex) {
  sim.rxn.assign(N_REACTIONS, Reaction{});
  sim.rxn_count.assign(N_REACTIONS, 0);

  // Male rates from in.full lines 56-128.
  // Note: many reactions are style NBR even if their nbr list is empty in the
  // rule, because add_rxn only picks LOCAL when BOTH nbrReactants AND
  // nbrProducts are empty (line 2290-2293 of app_rxn_diff.cpp). The rule
  // "local A nbr 0.0009652 local A nbr" has empty nbr reactants and empty nbr
  // products; SPPARKS's parser sees no `nbr` keyword neighbor species so both
  // nNbrReactant and nNbrProduct are 0, hence LOCAL. So which=0 is LOCAL.
  //
  // But "local B nbr 0.2375 local nbr B" (Tp motility) puts B in nbr products
  // -> nNbrProduct>0 -> NBR style. Correct.
  //
  // Recruitment "local I nbr 6.87204 local I nbr" has no nbr reactants or
  // products -> LOCAL style. Correct.

  set_reaction(sim.rxn[0],  0.0009652, {canID},   {},        {canID},           {}); // A -> A (cancer prolif)
  set_reaction(sim.rxn[1],  0.00463,   {TpID},    {},        {TpID, TpID},      {}); // B -> B B
  set_reaction(sim.rxn[2],  0.0,       {TmID},    {},        {TmID, TmID},      {}); // C -> C C
  set_reaction(sim.rxn[3],  0.0,       {TeID},    {},        {TeID, TeID},      {}); // D -> D D
  set_reaction(sim.rxn[4],  0.0,       {ARTpID},  {},        {ARTpID, ARTpID},  {}); // E -> E E
  set_reaction(sim.rxn[5],  0.0,       {ARTmID},  {},        {ARTmID, ARTmID},  {}); // F -> F F
  set_reaction(sim.rxn[6],  0.0,       {ARTeID},  {},        {ARTeID, ARTeID},  {}); // G -> G G

  set_reaction(sim.rxn[7],  6.87204,   {recID},   {},        {recID},           {}); // Tp recruit
  set_reaction(sim.rxn[8],  6.87204,   {recID},   {},        {recID},           {}); // Tm recruit
  set_reaction(sim.rxn[9],  0.0,       {recID},   {},        {recID},           {}); // Te recruit
  set_reaction(sim.rxn[10], 0.0,       {recID},   {},        {recID},           {}); // ARTp recruit
  set_reaction(sim.rxn[11], 1.52712,   {recID},   {},        {recID},           {}); // ARTm recruit
  set_reaction(sim.rxn[12], 0.0,       {recID},   {},        {recID},           {}); // ARTe recruit

  set_reaction(sim.rxn[13], 0.2375,    {TpID},    {},        {},                {TpID});   // Tp motility
  set_reaction(sim.rxn[14], 0.2375,    {TmID},    {},        {},                {TmID});   // Tm motility
  set_reaction(sim.rxn[15], 0.2375,    {TeID},    {},        {},                {TeID});   // Te motility
  set_reaction(sim.rxn[16], 0.0,       {ARTpID},  {},        {},                {ARTpID}); // ARTp motility
  set_reaction(sim.rxn[17], 0.2375,    {ARTmID},  {},        {},                {ARTmID}); // ARTm motility
  set_reaction(sim.rxn[18], 0.2375,    {ARTeID},  {},        {},                {ARTeID}); // ARTe motility

  set_reaction(sim.rxn[19], 0.00347,   {TpID},    {},        {TmID},            {}); // Tp -> Tm
  set_reaction(sim.rxn[20], 0.0,       {ARTpID},  {},        {ARTmID},          {}); // ARTp -> ARTm

  set_reaction(sim.rxn[21], 0.0106,    {canID, TmID},   {},  {canID, TeID},     {}); // Tm exh cancer local
  set_reaction(sim.rxn[22], 0.0106,    {canID},   {TmID},    {canID},           {TeID}); // Tm exh cancer nbr
  set_reaction(sim.rxn[23], 0.0106,    {canID, ARTmID}, {},  {canID, ARTeID},   {}); // ARTm exh cancer local
  set_reaction(sim.rxn[24], 0.0106,    {canID},   {ARTmID},  {canID},           {ARTeID}); // ARTm exh cancer nbr

  set_reaction(sim.rxn[25], sim.p.te_death_rate,   {TeID},    {},        {},                {}); // Te death
  set_reaction(sim.rxn[26], sim.p.te_death_rate,   {ARTeID},  {},        {},                {}); // ARTe death

  set_reaction(sim.rxn[27], 0.06,      {TmID, canID}, {},    {TmID},            {}); // Cancer lysis by Tm local
  set_reaction(sim.rxn[28], 0.06,      {TmID},    {canID},   {TmID},            {}); // Cancer lysis by Tm nbr
  set_reaction(sim.rxn[29], 0.0,       {TeID, canID}, {},    {TeID},            {}); // lysis by Te local
  set_reaction(sim.rxn[30], 0.0,       {TeID},    {canID},   {TeID},            {}); // lysis by Te nbr
  set_reaction(sim.rxn[31], 0.06,      {ARTmID, canID}, {},  {ARTmID},          {}); // lysis by ARTm local
  set_reaction(sim.rxn[32], 0.06,      {ARTmID},  {canID},   {ARTmID},          {}); // lysis by ARTm nbr
  set_reaction(sim.rxn[33], 0.0,       {ARTeID, canID}, {},  {ARTeID},          {});
  set_reaction(sim.rxn[34], 0.0,       {ARTeID},  {canID},   {ARTeID},          {});

  set_reaction(sim.rxn[35], 0.0,       {TpID},    {},        {},                {}); // Tp death
  set_reaction(sim.rxn[36], 0.0,       {TmID},    {},        {},                {}); // Tm death
  set_reaction(sim.rxn[37], 0.0,       {ARTpID},  {},        {},                {});
  set_reaction(sim.rxn[38], 0.0,       {ARTmID},  {},        {},                {});

  set_reaction(sim.rxn[39], 0.2375,    {ThID},    {},        {},                {ThID});   // Th motility
  set_reaction(sim.rxn[40], 0.2375,    {ARThID},  {},        {},                {ARThID}); // ARTh motility
  set_reaction(sim.rxn[41], 0.2375,    {TrID},    {},        {},                {TrID});   // Treg motility
  set_reaction(sim.rxn[42], 0.2375,    {ARTrID},  {},        {},                {ARTrID}); // ARTreg motility

  set_reaction(sim.rxn[43], 0.0,       {ThID, TpID}, {},     {ThID, TmID},      {});
  set_reaction(sim.rxn[44], 0.0,       {ThID, ARTpID}, {},   {ThID, ARTmID},    {});
  set_reaction(sim.rxn[45], 0.0,       {ThID}, {TpID},       {ThID}, {TmID});
  set_reaction(sim.rxn[46], 0.0,       {ThID}, {ARTpID},     {ThID}, {ARTmID});
  set_reaction(sim.rxn[47], 0.0,       {ARThID, TpID}, {},   {ARThID, TmID},    {});
  set_reaction(sim.rxn[48], 0.0,       {ARThID, ARTpID}, {}, {ARThID, ARTmID},  {});
  set_reaction(sim.rxn[49], 0.0,       {ARThID}, {TpID},     {ARThID}, {TmID});
  set_reaction(sim.rxn[50], 0.0,       {ARThID}, {ARTpID},   {ARThID}, {ARTmID});

  set_reaction(sim.rxn[51], 0.0106,    {TrID, TmID}, {},     {TrID, TeID},      {});
  set_reaction(sim.rxn[52], 0.0106,    {TrID}, {TmID},       {TrID}, {TeID});
  set_reaction(sim.rxn[53], 0.0106,    {TrID, ARTmID}, {},   {TrID, ARTeID},    {});
  set_reaction(sim.rxn[54], 0.0106,    {TrID}, {ARTmID},     {TrID}, {ARTeID});

  set_reaction(sim.rxn[55], 0.0106,    {ARTrID, TmID}, {},   {ARTrID, TeID},    {});
  set_reaction(sim.rxn[56], 0.0106,    {ARTrID}, {TmID},     {ARTrID}, {TeID});
  set_reaction(sim.rxn[57], 0.0106,    {ARTrID, ARTmID}, {}, {ARTrID, ARTeID},  {});
  set_reaction(sim.rxn[58], 0.0106,    {ARTrID}, {ARTmID},   {ARTrID}, {ARTeID});

  set_reaction(sim.rxn[59], 0.00347,   {TpID}, {},           {ARTmID},          {}); // Tp -> ARTm

  // Sex compensation (paper's female modifications) -- this is the
  // "ported" default behavior, matching in.full / real SPPARKS exactly.
  if (sex == FEMALE) {
    const double comp = 1.0 / 0.63;
    sim.rxn[7].rate  *= comp;       // Tp recruit
    sim.rxn[8].rate  *= comp;       // Tm recruit
    sim.rxn[11].rate  = 0.0;        // ARTm recruit removed
    sim.rxn[19].rate *= 2.0;        // Tp -> Tm accelerated
    sim.rxn[59].rate  = 0.0;        // Tp -> ARTm removed
  }

}

// ============================================================================
// SECTION 5: CUSTOM_MULTIPLIER
// Port of AppRxnDiffCustom::custom_multiplier (app_rxn_diff_custom.cpp).
//
// Only the leaky-receptor branch (!use_simple_decay_model) is implemented,
// since arg[32]=0 for this run. Other branches would only fire with different
// input args.
// ============================================================================

// Hill helper: pow(x, n) / (pow(K, n) + pow(x, n)), guarded against x=0
static inline double hill(double x, double K, int n) {
  double xn = std::pow(x, n);
  double Kn = std::pow(K, n);
  return xn / (Kn + xn + 1e-30);
}

// Wrap a coordinate periodically using floor-based math to match SPPARKS
// (which uses `target -= floor(double(target)/latCol)*latCol` everywhere).
static inline int wrap(int v, int L) {
  int r = v - (int)std::floor((double)v / (double)L) * L;
  // Guard against edge cases from floor rounding
  if (r < 0)  r += L;
  if (r >= L) r -= L;
  return r;
}

// Custom multiplier. Signature: i = local site, rstyle = LOCAL/NBR,
// which = reaction index, jpartner = neighbor site (or -1 for LOCAL).
static double custom_multiplier(Sim& sim, int i, int rstyle, int which, int jpartner) {
  (void)rstyle; // matches SPPARKS: rstyle isn't inspected inside these branches
  const Params& p = sim.p;

  const int latCol = LATCOL;
  const int maxA = p.max_cell_area_per_voxel;
  const int cancA = p.cancer_area_count;
  const int tpA = p.Tp_area_count;
  const int tmA = p.Tm_area_count;
  const int teA = p.Te_area_count;
  const int thA = tpA;
  const int trA = thA;

  const int first_motility_interaction     = 13;
  const int first_cd4_motility_interaction = 39;
  const int first_T_prolif_interaction     = 1;
  const int first_lysis_interaction        = 27;

  // -------------------------------------------------------------------------
  // which == 0 : cancer proliferation
  //   Ring search 1,3,5,7,... for a target with room for cancer.
  //   Returns proliferation_factor based on suppression via IFNG binding.
  //   In leaky branch: avg per-cancer-cell modifier at LOCAL site i (not target).
  //   Returns 0 if no target found within widthLim.
  //
  // Carrying-capacity multiplier (--carrying-capacity <K>):
  //   Multiplies proliferation_factor by (1 - N_cancer / K), clamped >= 0.
  //   Produces a concave-down (parabolic) growth curve as tumor size
  //   approaches K. When K=0 (default), this feature is disabled --
  //   validated default behavior is preserved.
  // -------------------------------------------------------------------------
  if (which == 0) {
    double proliferation_factor = 1.0;
    bool done = false;
    int neighLim = 3;
    int widthLim = 3 + 2*(neighLim - 1);
    int width = 1;
    int ix = i % latCol;
    int iy = i / latCol;
    while (!done) {
      int halfw = (width - 1) / 2;
      for (int ly = 0; ly < 2; ly++) {
        for (int lx = 0; lx < width; lx++) {
          int t1x = ix + lx - halfw;
          int t1y = iy + (2*ly - 1) * halfw;
          int t2x = ix + (2*ly - 1) * halfw;
          int t2y = iy + lx - halfw;
          t1x = wrap(t1x, latCol); t1y = wrap(t1y, latCol);
          t2x = wrap(t2x, latCol); t2y = wrap(t2y, latCol);
          int target1 = t1x + latCol * t1y;
          int target2 = t2x + latCol * t2y;

          // target1
          if (sim.population[canID][target1] * cancA < maxA) {
            if (p.ifng_cancer_proliferation_limitation) {
              if (!p.use_simple_decay_model) {
                int local_cancer_count = sim.population[canID][i];
                double average_modifier = 0.0;
                auto& list = sim.cell_bound[i * N_TRACKED + 0]; // canID slot=0
                for (int cell = 0; cell < local_cancer_count && cell < (int)list.size(); cell++) {
                  double bound = list[cell][0];
                  double rec   = list[cell][1];
                  double bf    = (rec > 0.0) ? (bound / rec) : 0.0;
                  double h     = hill(bf, p.ifng_half_max_proliferation, p.hill_degree);
                  double pf    = p.cancer_proliferation_minimum_factor
                               + (1.0 - p.cancer_proliferation_minimum_factor) * (1.0 - h);
                  average_modifier += pf / local_cancer_count;
                }
                proliferation_factor = average_modifier;
              }
            }
            return proliferation_factor;
          }
          // target2
          if (sim.population[canID][target2] * cancA < maxA) {
            if (p.ifng_cancer_proliferation_limitation) {
              if (!p.use_simple_decay_model) {
                int local_cancer_count = sim.population[canID][i];
                double average_modifier = 0.0;
                auto& list = sim.cell_bound[i * N_TRACKED + 0];
                for (int cell = 0; cell < local_cancer_count && cell < (int)list.size(); cell++) {
                  double bound = list[cell][0];
                  double rec   = list[cell][1];
                  double bf    = (rec > 0.0) ? (bound / rec) : 0.0;
                  double h     = hill(bf, p.ifng_half_max_proliferation, p.hill_degree);
                  double pf    = p.cancer_proliferation_minimum_factor
                               + (1.0 - p.cancer_proliferation_minimum_factor) * (1.0 - h);
                  average_modifier += pf / local_cancer_count;
                }
                proliferation_factor = average_modifier;
              }
            }
            return proliferation_factor;
          }
        }
      }
      if (width >= widthLim) { done = true; return 0.0; }
      width += 2;
    }
    return 0.0;
  }

  // -------------------------------------------------------------------------
  // which == 7..12 : T cell recruitment
  // -------------------------------------------------------------------------
  if (7 <= which && which <= 12) {
    double t_carry_cap = (double)p.max_T_cells;
    double t_tot_pop   = (double)sim.tcell_total_population;
    double deadTrec_fac = 1.0 - (t_tot_pop * p.do_carrying_capacity / t_carry_cap);
    if (deadTrec_fac < 0.0) deadTrec_fac = 0.0;
    return deadTrec_fac;
  }

  // -------------------------------------------------------------------------
  // which == 1..6 : T cell proliferation
  //   Room check across ALL cell species (matches line 241 of custom).
  // -------------------------------------------------------------------------
  if (first_T_prolif_interaction <= which && which <= 6) {
    int filling = sim.filling(i);
    int weights[] = {tpA, tmA, teA, tpA, tmA, teA};
    int weight = weights[which - first_T_prolif_interaction];
    if (filling + weight <= maxA) {
      double t_carry_cap = (double)p.max_T_cells;
      double t_tot_pop   = (double)sim.tcell_total_population;
      double deadTrec_fac = 1.0 - (t_tot_pop * p.do_carrying_capacity / t_carry_cap);
      if (deadTrec_fac < 0.0) deadTrec_fac = 0.0;
      return deadTrec_fac;
    }
    return 0.0;
  }

  // -------------------------------------------------------------------------
  // Motility (13..18 CD8, 39..42 CD4/Treg)
  // -------------------------------------------------------------------------
  if ((first_motility_interaction <= which && which <= 18)
      || (first_cd4_motility_interaction <= which && which <= 42)) {

    int nfilling = sim.filling(jpartner);

    // axis/direction from |i - jpartner|
    int axis = 0;
    int direction = 0;
    if (std::abs(i - jpartner) > 1) axis = 1;
    if (i < jpartner) direction = 1;

    // Weights indexed by CD8 index (13..18 -> 0..5) or CD4 index (39..42 -> 6..9)
    int weights[10] = {tpA, tmA, teA, tpA, tmA, teA, tpA, tpA, tpA, tpA};
    // ^ last 4: Th, ARTh, Tr, ARTr (all use Tp_area_count per source line 277)
    int weight;
    double adhesion_factor = 1.0;
    int cd8_index;
    if (which - first_motility_interaction < 6) {
      cd8_index = which - first_motility_interaction; // 0..5
      weight = weights[cd8_index];
    } else {
      cd8_index = which - first_cd4_motility_interaction + 6; // 6..9
      weight = weights[cd8_index];
    }

    // bias_factor from cancer_bias_potential (0..1) via 2*potential
    int bx = (i % latCol) / CANCER_POTENTIAL_CONVERTER;
    int by = (i / latCol) / CANCER_POTENTIAL_CONVERTER;
    if (bx >= CP_MAX) bx = CP_MAX - 1;
    if (by >= CP_MAY) by = CP_MAY - 1;
    double bias_factor = 2.0 * sim.cancer_bias_potential[bx][by][axis];
    if (direction == 0) bias_factor = 2.0 - bias_factor;

    // Boundary-crossing check for motility exists in source (lines 1157-1218
    // of site_event, and 291-334 of custom) BUT is gated by `&& false` in
    // custom (line 1157), so it's dead code. Do not implement.

    // Room check
    if (maxA < nfilling + weight) return 0.0;

    // ifng_diffusion_factor (per-cell average) in leaky branch
    double ifng_diffusion_factor = 1.0;
    if (which - first_motility_interaction < 6) {
      // CD8 -> adhesion_factor = tcell_stick_fac
      adhesion_factor = p.tcell_stick_fac;
    } else {
      // CD4/Treg -> adhesion_factor stays 1.0
    }

    if (!p.use_simple_decay_model) {
      // cd8map: species IDs at the "cd8_index"
      static const int cd8map[10] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID,
                                     ThID, ARThID, TrID, ARTrID};
      // tctracker_map: tracker slot at that index. Slot ordering matches
      // cellIDs[11] = {canID, TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID,
      // ARThID, TrID, ARTrID}, so:
      //   TpID -> slot 1, TmID -> slot 2, TeID -> slot 3,
      //   ARTpID -> slot 4, ARTmID -> slot 5, ARTeID -> slot 6,
      //   ThID -> slot 7, ARThID -> slot 8, TrID -> slot 9, ARTrID -> slot 10
      static const int tctracker_map[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
      int cellID = cd8map[cd8_index];
      int trackID = tctracker_map[cd8_index];
      int local_tc_count = sim.population[cellID][i];

      int denom = local_tc_count;
      double average_modifier = 0.0;
      auto& list = sim.cell_bound[i * N_TRACKED + trackID];
      for (int cell = 0; cell < local_tc_count && cell < (int)list.size(); cell++) {
        double bound = list[cell][0];
        double rec   = list[cell][1];
        double bf    = (rec > 0.0) ? (bound / rec) : 0.0;
        double h     = hill(bf, p.ifng_half_max_diffusion, p.hill_degree);
        double idf   = adhesion_factor + (1.0 - adhesion_factor) * h;
        if (denom > 0) average_modifier += idf / denom;
      }
      if (denom > 0) ifng_diffusion_factor = average_modifier;
    }

    // Adhesion applies if cancer at i or ±y/±x nbrs (source lines 377-397).
    int ix = i % latCol;
    int iy = i / latCol;
    bool near_cancer = (sim.population[canID][i] > 0);
    if (!near_cancer) {
      for (int k = 0; k < 2 && !near_cancer; k++) {
        int targy = wrap(iy - 1 + 2*k, latCol);
        int targ  = ix + latCol * targy;
        if (sim.population[canID][targ] > 0) near_cancer = true;
      }
      for (int h_ = 0; h_ < 2 && !near_cancer; h_++) {
        int targx = wrap(ix - 1 + 2*h_, latCol);
        int targ  = targx + latCol * iy;
        if (sim.population[canID][targ] > 0) near_cancer = true;
      }
    }

    if (near_cancer) {
      // Original ported behavior (matching real SPPARKS): multiply by
      // adhesion_factor once inside ifng_diffusion_factor's blend AND once
      // as an outer multiplier. This effectively squares the adhesion
      // suppression when bound-IFNG is low, keeping T-cells strongly
      // adhered to cancer they contact. The "paper equation" version
      // (removing the outer multiply) was tested earlier but produced
      // tumor collapse under paper-rate T-cell proliferation boost, so
      // reverted to unified ported behavior for all rate modes.
      return bias_factor * adhesion_factor * ifng_diffusion_factor;
    }
    return bias_factor;
  }

  // -------------------------------------------------------------------------
  // Lysis (27..34)
  // -------------------------------------------------------------------------
  if (first_lysis_interaction <= which && which <= 34) {
    bool local = (which % 2 == 1);
    // Boundary crossing check on non-local lysis
    if (!local && sim.population[latEID][i] >= 1 && sim.population[latEID][jpartner] >= 1) {
      bool cross = false;
      int cornCount = 0, cornCount2 = 0;
      int type = -1, type2 = -1;
      if      (i / latCol == 0)          { cornCount++; type = 0; }
      else if (i / latCol == latCol - 1) { cornCount++; type = 2; }
      if      (i % latCol == 0)          { cornCount++; type = 1; }
      else if (i % latCol == latCol - 1) { cornCount++; type = 3; }
      if      (jpartner / latCol == 0)          { cornCount2++; type2 = 0; }
      else if (jpartner / latCol == latCol - 1) { cornCount2++; type2 = 2; }
      if      (jpartner % latCol == 0)          { cornCount2++; type2 = 1; }
      else if (jpartner % latCol == latCol - 1) { cornCount2++; type2 = 3; }
      if ((cornCount == 2 && cornCount2 == 2)
          || (cornCount == 1 && cornCount2 == 1 && type != type2)) cross = true;
      if (cross) return 0.0;
    }
    // Per source, cd8_index = (which - first_lysis)/2. cd8map = {TmID, TeID,
    // ARTmID, ARTeID}. cellID = species doing the lysing (has bound tracker).
    int cd8_index = (which - first_lysis_interaction) / 2;
    static const int cd8map[4] = {TmID, TeID, ARTmID, ARTeID};
    int cellID = cd8map[cd8_index];

    int ix = i % latCol;
    int iy = i / latCol;
    (void)ix; (void)iy;

    if (p.ifng_cytotoxicity_augmentation) {
      double modifier = 1.0;
      if (!p.use_simple_decay_model) {
        int local_tc_count = sim.population[cellID][i];
        int denom = local_tc_count;
        double average_modifier = 0.0;
        // Which tracker slot? cellID -> slot via species-ID-in-cellIDs order:
        //   TmID(2)->2, TeID(3)->3, ARTmID(5)->5, ARTeID(6)->6
        int trackID = cellID; // for species IDs 0..6, tracker slot = species ID
        auto& list = sim.cell_bound[i * N_TRACKED + trackID];
        for (int cell = 0; cell < local_tc_count && cell < (int)list.size(); cell++) {
          double bound = list[cell][0];
          double rec   = list[cell][1];
          double bf    = (rec > 0.0) ? (bound / rec) : 0.0;
          double h     = hill(bf, p.ifng_half_max_cytotoxic_augmentation_threshold, p.hill_degree);
          if (denom > 0) average_modifier += h / denom;
        }
        modifier = average_modifier;
      }
      return modifier;
    } else {
      return 0.5;
    }
  }

  // -------------------------------------------------------------------------
  // Exhaustion / Th-differentiation (21..24, 43..58)
  //   Boundary crossing check, else default 1.0.
  // -------------------------------------------------------------------------
  if ((21 <= which && which <= 24) || (43 <= which && which <= 58)) {
    if (sim.population[latEID][i] >= 1 && sim.population[latEID][jpartner] >= 1) {
      bool cross = false;
      int cornCount = 0, cornCount2 = 0;
      int type = -1, type2 = -1;
      if      (i / latCol == 0)          { cornCount++; type = 0; }
      else if (i / latCol == latCol - 1) { cornCount++; type = 2; }
      if      (i % latCol == 0)          { cornCount++; type = 1; }
      else if (i % latCol == latCol - 1) { cornCount++; type = 3; }
      if      (jpartner / latCol == 0)          { cornCount2++; type2 = 0; }
      else if (jpartner / latCol == latCol - 1) { cornCount2++; type2 = 2; }
      if      (jpartner % latCol == 0)          { cornCount2++; type2 = 1; }
      else if (jpartner % latCol == latCol - 1) { cornCount2++; type2 = 3; }
      if ((cornCount == 2 && cornCount2 == 2)
          || (cornCount == 1 && cornCount2 == 1 && type != type2)) cross = true;
      if (cross) return 0.0;
    }
    return 1.0;
  }

  return 1.0;
}

// ============================================================================
// SECTION 6: SITE_PROPENSITY
// Port of AppRxnDiff::site_propensity (app_rxn_diff.cpp lines 216-331).
// Standard mass-action multiplicity + rate + custom_multiplier.
// ============================================================================

// Speed-ups versus v6 (results are bit-identical, verified on seed 8377):
//   * reactions with rate 0 are skipped (they can never produce an event;
//     custom_multiplier has no side effects);
//   * each reaction's reactant species are pre-listed, in the same ascending
//     species order as the original 18-species scan, so the multiplicity
//     product is formed in exactly the same order;
//   * for NBR reactions, a missing LOCAL reactant is detected once instead of
//     once per neighbour (it fails for all four neighbours either way).
struct ReactantTerm { int s; int nloc; int nnbr; };
static std::vector<std::vector<ReactantTerm>> g_terms; // per reaction
static std::vector<int> g_local_species; // species that appear as a LOCAL reactant

static void build_reactant_terms(const Sim& sim) {
  g_terms.assign(N_REACTIONS, {});
  for (int m = 0; m < N_REACTIONS; m++) {
    const Reaction& rx = sim.rxn[m];
    for (int s = 0; s < MAX_SPECIES; s++) {
      if (rx.localReactants[s] > 0 || rx.nbrReactants[s] > 0)
        g_terms[m].push_back({s, rx.localReactants[s], rx.nbrReactants[s]});
    }
  }
  g_local_species.clear();
  for (int s = 0; s < MAX_SPECIES; s++) {
    bool used = false;
    for (int m = 0; m < N_REACTIONS; m++)
      if (sim.rxn[m].rate != 0.0 && sim.rxn[m].localReactants[s] > 0) used = true;
    if (used) g_local_species.push_back(s);
  }
  for (int m = 0; m < N_REACTIONS; m++) {          // safety: assumption check
    bool has_local = false;
    for (int s = 0; s < MAX_SPECIES; s++) if (sim.rxn[m].localReactants[s] > 0) has_local = true;
    if (!has_local) { fprintf(stderr, "reaction %d has no local reactant\n", m); std::exit(1); }
  }
}

static double site_propensity(Sim& sim, int i) {
  sim.clear_events(i);
  double proball = 0.0;

  // Fast path: every reaction needs >=1 local reactant, so a site holding
  // none of those species has zero propensity and no events.
  bool any = false;
  for (int s : g_local_species) if (sim.population[s][i] > 0) { any = true; break; }
  if (!any) return 0.0;

  for (int m = 0; m < N_REACTIONS; m++) {
    const Reaction& rx = sim.rxn[m];
    if (rx.rate == 0.0) continue;
    const std::vector<ReactantTerm>& terms = g_terms[m];

    // local reactants must be present at i (both LOCAL and NBR styles)
    bool missing_local = false;
    for (const ReactantTerm& t : terms)
      if (t.nloc > 0 && sim.population[t.s][i] < t.nloc) { missing_local = true; break; }
    if (missing_local) continue;

    if (rx.rxnStyle == LOCAL) {
      double multiplicity = 1.0;
      for (const ReactantTerm& t : terms) {
        for (int k = 1; k <= t.nloc; k++)
          multiplicity *= (double)(sim.population[t.s][i] + 1 - k) / (double)k;
      }
      double my_rate = rx.rate * multiplicity * custom_multiplier(sim, i, LOCAL, m, -1);
      if (my_rate == 0.0) continue;
      sim.add_event(i, LOCAL, m, my_rate, -1);
      proball += my_rate;
    } else { // NBR
      for (int jj = 0; jj < 4; jj++) {
        int j = sim.neighbor[i][jj];
        double multiplicity = 1.0;
        bool flag = false;
        for (const ReactantTerm& t : terms) {
          for (int k = 1; k <= t.nloc; k++)
            multiplicity *= (double)(sim.population[t.s][i] + 1 - k) / (double)k;
          if (t.nnbr > 0) {
            if (sim.population[t.s][j] < t.nnbr) { flag = true; break; }
            for (int k = 1; k <= t.nnbr; k++)
              multiplicity *= (double)(sim.population[t.s][j] + 1 - k) / (double)k;
          }
        }
        if (flag) continue;
        double my_rate = rx.rate * multiplicity * custom_multiplier(sim, i, NBR, m, j);
        if (my_rate == 0.0) continue;
        sim.add_event(i, NBR, m, my_rate, j);
        proball += my_rate;
      }
    }
  }
  return proball;
}

// ============================================================================
// SECTION 7: UPDATE_IFNG_FIELD
// Port of AppLattice::update_ifng_field (app_lattice.cpp lines 732-1024).
// Only the leaky-receptor branch is implemented.
// ============================================================================

static void update_ifng_field(Sim& sim, double time_step) {
  Params& p = sim.p;
  const int latCol = LATCOL;
  const int cvt = IFNG_LATTICE_CELL_CONVERTER;

  // IFNGR receptor internalization rate for AR+ cells (line 751 of source)
  const double IFNGR_decay_rate = std::exp(-p.bIFNG_decay * TIME_CONVERSION * p.ifng_time_step);

  // Cell IDs in tracker slot order (source line 758)
  static const int cellIDs[N_TRACKED] = {
    canID, TpID, TmID, TeID, ARTpID, ARTmID, ARTeID,
    ThID, ARThID, TrID, ARTrID
  };
  // AR+ set for receptor internalization (source line 763, cd8IDs[1])
  auto is_ARplus = [](int species_id) {
    return species_id == ARTpID || species_id == ARTmID
        || species_id == ARTeID || species_id == ARThID;
  };

  // 1) Copy field into buffer
  for (int i = 0; i < IFNG_NSITES; i++) sim.ifng_field_copy[i] = sim.ifng_field[i];

  // 2) Loop through coarse voxels
  for (int cy = 0; cy < IFNG_ROW; cy++) {
    for (int cx = 0; cx < IFNG_COL; cx++) {
      // Representative fine voxel: (cx*3, cy*3)
      int rep_fine = cx * cvt + LATCOL * cvt * cy; // source uses this exact index

      if (sim.population[boundID][rep_fine] != 1) {
        // Non-diffusable region
        sim.ifng_field[cx + IFNG_COL * cy] = 0.0;
        continue;
      }

      // Laplacian (non-periodic BC unless flag set): 4 neighbors
      double nbrs[4] = {0.0, 0.0, 0.0, 0.0};
      int ncounter = 0;
      for (int dx = -1; dx < 3; dx += 2) {
        for (int dy = 0; dy < 2; dy++) {
          int nx = cx + dx * dy;
          int ny = cy + dx * (1 - dy);
          if (p.periodic_boundary_conditions) {
            nx = wrap(nx, IFNG_COL);
            ny = wrap(ny, IFNG_ROW);
            nbrs[ncounter] = sim.ifng_field_copy[nx + IFNG_COL * ny];
          } else {
            if (nx < 0 || nx >= IFNG_COL || ny < 0 || ny >= IFNG_ROW) {
              nbrs[ncounter] = 0.0;
            } else {
              int rep_nbr = nx * cvt + LATCOL * cvt * ny;
              if (sim.population[boundID][rep_nbr] == 0) nbrs[ncounter] = 0.0;
              else nbrs[ncounter] = sim.ifng_field_copy[nx + IFNG_COL * ny];
            }
          }
          ncounter++;
        }
      }
      double nbr_sum = nbrs[0] + nbrs[1] + nbrs[2] + nbrs[3];
      double self_val = sim.ifng_field_copy[cx + IFNG_COL * cy];
      double laplacian = (nbr_sum - 4.0 * self_val) / (double)(IFNG_LATT_CONST * IFNG_LATT_CONST);

      // 3) Sweep 3x3 fine voxels within this coarse voxel:
      //    - Apply AR+ receptor internalization (multiplies cell[0] by decay rate)
      //    - Sum total_receptors, total_bound
      //    - Sum total_production from AR- CD8 producers (Tp, Tm) or AR- Th
      double total_receptors = 0.0;
      double total_bound = 0.0;
      double total_production = 0.0;
      for (int yy = 0; yy < cvt; yy++) {
        for (int xx = 0; xx < cvt; xx++) {
          int latx = cx * cvt + xx;
          int laty = cy * cvt + yy;
          if (latx < 0 || latx >= latCol || laty < 0 || laty >= latCol) continue;
          int site = latx + latCol * laty;
          // Clear per-fine-voxel bound aggregate
          sim.bound_ifng[site][0] = 0.0;
          sim.bound_ifng[site][1] = 0.0;
          sim.bound_ifng[site][2] = 0.0;

          for (int arn = 0; arn < N_TRACKED; arn++) {
            int cid = cellIDs[arn];
            auto& list = sim.cell_bound[site * N_TRACKED + arn];
            for (auto& cell : list) {
              // AR+ receptor internalization
              if (is_ARplus(cid) && !p.no_sequestration) {
                cell[0] = cell[0] * IFNGR_decay_rate;
              }
              double cell_bound = cell[0];
              double cell_receptors = cell[1];
              total_receptors += cell_receptors;
              total_bound += cell_bound;
              // Producers: AR- Tp (arn=1), AR- Tm (arn=2), AR- Th (arn=7).
              // Excludes Te (exhausted, paper section 4.5.2 CD8+ T cells).
              // Excludes all AR+ (paper section 4.5.2: "AR+ T cells do not
              // secrete IFNG").
              if ((1 <= arn && arn <= 2) || arn == 7) {
                double bound_reading = (cell_receptors > 0.0) ? (cell_bound / cell_receptors) : 0.0;
                double fb_hill = hill(bound_reading, p.half_max_ifng_feedback, p.ifng_feedback_hill_degree);
                // k_i,prod = k_max*(0.25 + 0.75*Hill(g_i)), with k_max = 0.1*2.0 = 0.2 /s
                double feedback = (p.max_ifng_production - 0.5) * fb_hill + 0.5;
                total_production += p.ifng_production * feedback;
              }
            }
          }
        }
      }

      // 4) Diffusion+source update (source line 917-929)
      // total_consumption is 0 in the !use_simple_decay_model branch (source
      // only accumulates it in the simple-decay branch), so omit it.
      double source_term = total_production - p.ifng_decay * self_val;
      double new_val = sim.ifng_field[cx + IFNG_COL * cy]
                     + time_step * TIME_CONVERSION * (p.ifng_diffusion * laplacian + source_term);
      if (new_val < 1e-10) new_val = 0.0;
      sim.ifng_field[cx + IFNG_COL * cy] = new_val;

      // 5) Equilibrium binding solve (source line 952-954)
      double total_local_ifng = total_bound + new_val;
      double R = total_receptors;
      double K = p.dissC;
      double a = R + K - total_local_ifng;
      double disc = a * a + 4.0 * K * total_local_ifng;
      double final_free  = -0.5 * (a - std::sqrt(disc));
      double final_bound = total_local_ifng - final_free;
      sim.ifng_field[cx + IFNG_COL * cy] = final_free;

      // 6) Partition bound proportionally across all cells in coarse voxel
      //    (leaky branch: !use_simple_decay_model)
      if (!p.use_simple_decay_model) {
        double partition_frac = (total_receptors > 0.0) ? (final_bound / total_receptors) : 0.0;
        for (int yy = 0; yy < cvt; yy++) {
          for (int xx = 0; xx < cvt; xx++) {
            int latx = cx * cvt + xx;
            int laty = cy * cvt + yy;
            if (latx < 0 || latx >= latCol || laty < 0 || laty >= latCol) continue;
            int site = latx + latCol * laty;
            for (int arn = 0; arn < N_TRACKED; arn++) {
              auto& list = sim.cell_bound[site * N_TRACKED + arn];
              for (auto& cell : list) {
                cell[0] = partition_frac * cell[1];
              }
            }
          }
        }
        // 7) Recompute per-fine-voxel bound_ifng from tracker
        for (int yy = 0; yy < cvt; yy++) {
          for (int xx = 0; xx < cvt; xx++) {
            int latx = cx * cvt + xx;
            int laty = cy * cvt + yy;
            if (latx < 0 || latx >= latCol || laty < 0 || laty >= latCol) continue;
            int site = latx + latCol * laty;
            for (int arn = 0; arn < N_TRACKED; arn++) {
              auto& list = sim.cell_bound[site * N_TRACKED + arn];
              for (auto& cell : list) {
                sim.bound_ifng[site][2] += cell[1];
                sim.bound_ifng[site][1] += cell[0];
              }
            }
            if (sim.bound_ifng[site][2] > 1.0) {
              sim.bound_ifng[site][0] = sim.bound_ifng[site][1] / sim.bound_ifng[site][2];
            } else {
              sim.bound_ifng[site][0] = -1.0;
            }
          }
        }
      }
    }
  }

  // 8) Update dumped fields (ifng_field_discretized[0..1])
  //    per fine voxel, matching source line 1013-1021.
  for (int yy = 0; yy < LATROW; yy++) {
    for (int xx = 0; xx < LATCOL; xx++) {
      int coarse = (xx / cvt) + IFNG_COL * (yy / cvt);
      double per_fine = sim.ifng_field[coarse] / (double)(cvt * cvt);
      sim.ifng_field_discretized_free [xx + LATCOL * yy] = per_fine;
      sim.ifng_field_discretized_bound[xx + LATCOL * yy] = sim.bound_ifng[xx + LATCOL * yy][0];
    }
  }

}

// ============================================================================
// SECTION 8: SITE_EVENT
// Port of AppRxnDiff::site_event (app_rxn_diff.cpp lines 338-2134).
// Only the leaky-receptor branch (!use_leaky_receptor_model=true) is
// implemented. Tracker slot indexing follows cellIDs[11] order.
//
// Structure (matches source):
//   1. Draw event from event list at site i via ranapp uniform (line 345)
//   2. Apply localDeltaPop / nbrDeltaPop (lines 373-395)
//   3. Reaction-specific logic:
//      which==0        : cancer proliferation (lines 439-780)
//      which==7..12    : recruitment          (lines 1305-1384)
//      which==13..18,
//              39..42  : motility tracker move (lines 1156-1302)
//      which==1..6     : T-cell proliferation tracker duplicate (line 1598)
//      which==25,26,
//              35..38  : T-cell death tracker erase (line 1611)
//      which==27..34   : cancer lysis tracker erase (line 1625)
//      which==19,20,59 : Tp differentiation tracker move (line 1638)
//      which==43..50   : Th-induced differentiation (line 1659) [all zero rate]
//      which==21..24   : Tm exhaustion by cancer (line 1682)
//      which==51..58   : Tm exhaustion by Treg (line 1701)
//   4. Cancer bias potential recompute (lines 1813-1861) if which==0 or 27-34
//   5. proliferation_update_bool / cancer_proliferation_update_bool logic
//   6. Touched-site propensity refresh (lines 1878-2128) + solve->update()
// ============================================================================

// Cell IDs in tracker slot order (matches cellIDs[11] in source)
static const int TRK_CELL_ID[N_TRACKED] = {
  canID, TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID, ARThID, TrID, ARTrID
};
// Reverse map: species ID -> tracker slot (or -1 if not tracked)
static int species_to_track_slot(int species_id) {
  for (int k = 0; k < N_TRACKED; k++) {
    if (TRK_CELL_ID[k] == species_id) return k;
  }
  return -1;
}

// Receptors-per-cell map indexed by tracker slot (source line 106):
//   {2, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1}
// where 0 = AR- pool, 1 = AR+ pool, 2 = cancer pool.
static const int RECEPTORS_PER_CELL_ID[N_TRACKED] = {2, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1};

// Cancer division: choose target voxel from opos/popos ring search.
// Also handles the T-cell pushing when the chosen target has T cells to bump.
// Returns the chosen daughter voxel (was written to population[canID][choice]).
// tcellPos out-parameter: sites where T cells got pushed to (empty if no push).
static int cancer_division_place(Sim& sim, int i, double r1_ring, RandomPark* ran,
                                 std::vector<int>& tcellPos_out) {
  const int latCol = LATCOL;
  const Params& p = sim.p;
  const int maxA = p.max_cell_area_per_voxel;
  const int cancA = p.cancer_area_count;

  int ix = i % latCol;
  int iy = i / latCol;
  int width = 1;
  bool done = false;
  std::vector<int> opos;   // free targets
  std::vector<int> popos;  // pushable targets (only T cells blocking)

  while (!done) {
    int halfw = (width - 1) / 2;
    for (int ly = 0; ly < 2; ly++) {
      for (int lx = 0; lx < width; lx++) {
        int t1x = ix + lx - halfw;
        int t1y = iy + (2*ly - 1) * halfw;
        int t2x = ix + (2*ly - 1) * halfw;
        int t2y = iy + lx - halfw;
        t1x = wrap(t1x, latCol); t1y = wrap(t1y, latCol);
        t2x = wrap(t2x, latCol); t2y = wrap(t2y, latCol);
        int target1 = t1x + latCol * t1y;
        int target2 = t2x + latCol * t2y;

        int f1 = sim.filling(target1);
        int npf1 = sim.non_push_filling(target1);
        int f2 = sim.filling(target2);
        int npf2 = sim.non_push_filling(target2);

        if (sim.population[canID][target1] * cancA < maxA) {
          if (maxA - f1 >= cancA) opos.push_back(target1);
          else if (maxA - npf1 >= cancA) popos.push_back(target1);
        }
        if (sim.population[canID][target2] * cancA < maxA) {
          if (maxA - f2 >= cancA) opos.push_back(target2);
          else if (maxA - npf2 >= cancA) popos.push_back(target2);
        }
      }
    }
    if (!opos.empty() || !popos.empty() || width >= WIDTH_LIM) {
      std::sort(opos.begin(), opos.end());
      opos.erase(std::unique(opos.begin(), opos.end()), opos.end());
      std::sort(popos.begin(), popos.end());
      popos.erase(std::unique(popos.begin(), popos.end()), popos.end());
      done = true;
    } else {
      width += 2;
    }
  }

  if (!opos.empty()) {
    int idx = (int)(r1_ring * opos.size());
    if (idx >= (int)opos.size()) idx = (int)opos.size() - 1;
    int choice = opos[idx];
    sim.population[canID][choice]++;
    return choice;
  }

  if (!popos.empty()) {
    // Push T cells: choose a popos target, move T cells outward to make room.
    // Matches source lines 638-1155 (leaky branch subset).
    int idx = (int)(r1_ring * popos.size());
    if (idx >= (int)popos.size()) idx = (int)popos.size() - 1;
    int chosen = popos[idx];

    static const int cd8map[10] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID,
                                   ThID, ARThID, TrID, ARTrID};
    // Weights indexed by cd8map order (all take area 1 in this config)
    const int weights[10] = {p.Tp_area_count, p.Tm_area_count, p.Te_area_count,
                             p.Tp_area_count, p.Tm_area_count, p.Te_area_count,
                             p.Tp_area_count, p.Tp_area_count, p.Tp_area_count,
                             p.Tp_area_count};

    // Snapshot counts of T cells at chosen
    int icd8s[10];
    for (int k = 0; k < 10; k++) icd8s[k] = sim.population[cd8map[k]][chosen];
    int cd8s_out[10] = {0,0,0,0,0,0,0,0,0,0};

    int cfilling = sim.filling(chosen);
    int space_left = maxA - cfilling;
    int tcount = cancA - space_left;

    // Decide which T cells get bumped (source line 784)
    int cd8_index_range = 2 * (p.Tp_area_count + p.Tm_area_count + p.Te_area_count
                              + p.Tp_area_count); // 2*(tpA+tmA+teA+thA), thA=tpA
    // Source uses (2*tpA + 2*tmA + 2*teA + 2*thA) as the sampling range for
    // "which of the 10 cd8map slots to try", but 2*(tpA+tmA+teA+thA) is 8 for
    // area=1 across all. Actually the correct range is 10 (# slots). Let's
    // reproduce source semantics: it picks int(range * uniform()) as the slot
    // index. If range=8 with uniform in [0,1), you never pick slots 8 or 9.
    // That's how source line 786 reads:
    //   int to_try = int((2*tpA + 2*tmA + 2*teA + 2*thA)*ran);
    // For this config: 2*1+2*1+2*1+2*1 = 8. So only 8 slot indices possible;
    // Tregs (indices 8, 9) NEVER get pushed. This is source's actual behavior.
    int slot_range = cd8_index_range;
    std::vector<std::array<double,2>> cd8_tracker_out[10];

    int guard = 0;
    while (tcount > 0 && guard++ < 1000) {
      double ran_u = ran->uniform();
      int to_try = (int)(slot_range * ran_u);
      if (to_try < 0) to_try = 0;
      if (to_try >= 10) to_try = 9;
      if (icd8s[to_try] > 0) {
        icd8s[to_try]--;
        cd8s_out[to_try]++;
        tcount--;
        sim.population[cd8map[to_try]][chosen]--;

        // Move a random cell from tracker (leaky branch)
        int track_idx = (to_try < 6) ? cd8map[to_try] : (to_try + 1);
        // ^ For slots 0..5 (CD8), track_idx = species ID (1..6).
        //   For slots 6..9 (Th, ARTh, Tr, ARTr), tracker slots are 7..10.
        auto& list_src = sim.cell_bound[chosen * N_TRACKED + track_idx];
        if (!list_src.empty()) {
          int rc = (int)(ran->uniform() * list_src.size());
          if (rc >= (int)list_src.size()) rc = (int)list_src.size() - 1;
          if (rc < 0) rc = 0;
          cd8_tracker_out[to_try].push_back(list_src[rc]);
          list_src.erase(list_src.begin() + rc);
        }
      }
    }

    sim.population[canID][chosen]++; // cancer takes the space

    // Now push the collected T cells to a nearby free target voxel.
    // The source's logic here (lines 830-1155) is complex; a faithful
    // implementation searches nearest ring for a voxel with room and dumps
    // all outward T cells there. If none found in NEIGH_LIM=3 rings, the
    // T cells are lost (which is what source does when it hits the goto-stop
    // exit path with target2 having room).
    //
    // We implement a simplified faithful version: search rings around chosen
    // for the first target with room for all cd8s_out. If not found within
    // ring depth 3, dump what fits into the first partial-room voxel found.
    int cx = chosen % latCol;
    int cy = chosen / latCol;
    int need = 0;
    for (int k = 0; k < 10; k++) need += cd8s_out[k] * weights[k];

    int put_target = -1;
    int pw = 3;
    for (int lo = 0; lo < NEIGH_LIM && put_target < 0; lo++) {
      int hw = (pw - 1) / 2;
      for (int ly = 0; ly < 2 && put_target < 0; ly++) {
        for (int lx = 0; lx < pw && put_target < 0; lx++) {
          int t1x = wrap(cx + lx - hw, latCol);
          int t1y = wrap(cy + (2*ly - 1) * hw, latCol);
          int t2x = wrap(cx + (2*ly - 1) * hw, latCol);
          int t2y = wrap(cy + lx - hw, latCol);
          int tg1 = t1x + latCol * t1y;
          int tg2 = t2x + latCol * t2y;
          int f1 = sim.filling(tg1);
          int f2 = sim.filling(tg2);
          if (f1 + need <= maxA) { put_target = tg1; break; }
          if (f2 + need <= maxA) { put_target = tg2; break; }
        }
      }
      pw += 2;
    }

    if (put_target >= 0) {
      int ptx = put_target % latCol;
      int pty = put_target / latCol;
      for (int k = 0; k < 10; k++) {
        if (cd8s_out[k] > 0) {
          sim.population[cd8map[k]][put_target] += cd8s_out[k];
          int track_idx = (k < 6) ? cd8map[k] : (k + 1);
          for (auto& cell : cd8_tracker_out[k]) {
            sim.cell_bound[put_target * N_TRACKED + track_idx].push_back(cell);
          }
          cd8_tracker_out[k].clear();
        }
      }
      tcellPos_out.push_back(put_target);
    } else {
      // Cells are dropped (rare; source has similar corner behavior).
      for (int k = 0; k < 10; k++) cd8_tracker_out[k].clear();
    }

    // Update tracker for daughter cancer cell (leaky branch, line 588-611)
    // Pick a parent cancer cell weighted by cell[0] (absorbed IFNG)
    int ix2 = i % latCol;
    int iy2 = i / latCol;
    auto& parent_list = sim.cell_bound[i * N_TRACKED + 0]; // canID slot=0
    if (!parent_list.empty()) {
      double gross = 0.0;
      for (auto& c : parent_list) gross += c[0];
      double rv = ran->uniform();
      int random_cell = 0;
      if (gross > 0.0) {
        double target_frac = rv * gross;
        double acc = 0.0;
        for (size_t cc = 0; cc < parent_list.size(); cc++) {
          random_cell = (int)cc;
          acc += parent_list[cc][0];
          if (target_frac <= acc) break;
        }
      } else {
        random_cell = (int)(rv * parent_list.size());
        if (random_cell >= (int)parent_list.size()) random_cell = (int)parent_list.size() - 1;
        if (random_cell < 0) random_cell = 0;
      }
      double saved = parent_list[random_cell][0];
      parent_list[random_cell][0] = saved / 2.0;
      std::array<double, 2> new_cell = { saved / 2.0, parent_list[random_cell][1] };
      sim.cell_bound[chosen * N_TRACKED + 0].push_back(new_cell);
    }
    (void)ix2; (void)iy2;
    return chosen;
  }

  // No target found at all
  return -1;
}

// Recruitment: pick a random slot across all boundID>0 voxels weighted by
// available room. Matches source lines 1305-1344.
static int recruitment_place(Sim& sim, int which, RandomPark* ran) {
  const Params& p = sim.p;
  const int latCol = LATCOL;
  static const int cd8map[10] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID,
                                 ThID, ARThID, TrID, ARTrID};
  const int weights[10] = {p.Tp_area_count, p.Tm_area_count, p.Te_area_count,
                           p.Tp_area_count, p.Tm_area_count, p.Te_area_count,
                           p.Tp_area_count, p.Tp_area_count, p.Tp_area_count,
                           p.Tp_area_count};
  const int maxA = p.max_cell_area_per_voxel;

  int tIndex = which - 7;
  int recType = cd8map[tIndex];
  int weight = weights[tIndex];

  double ranSpot_uni = ran->uniform();
  long long eCount = 0;
  for (int px = 0; px < NSITES; px++) {
    if (sim.population[boundID][px] > 0) {
      int f = sim.filling(px);
      int left = (maxA - f) / weight;
      if (left > 0) eCount += left;
    }
  }
  if (eCount == 0) return -1;
  long long ranSpot = (long long)(eCount * ranSpot_uni);
  if (ranSpot >= eCount) ranSpot = eCount - 1;

  eCount = 0;
  int recSpot = -1;
  for (int px = 0; px < NSITES; px++) {
    if (sim.population[boundID][px] > 0) {
      int f = sim.filling(px);
      int left = (maxA - f) / weight;
      if (left > 0 && eCount <= ranSpot && ranSpot <= eCount + left) {
        sim.population[recType][px] += 1;
        recSpot = px;
        break;
      }
      eCount += left;
    }
  }
  if (recSpot < 0) return -1;

  // Leaky branch: add new cell to tracker at recSpot (bound=0)
  int track_idx = (tIndex < 6) ? recType : (7 + (tIndex - 6));
  // ^ tIndex 0..5 -> track_idx = recType (species ID 1..6)
  //   tIndex 6..9 -> track_idx = 7,8,9,10 (Th, ARTh, Tr, ARTr slots)
  int rpc_id = RECEPTORS_PER_CELL_ID[track_idx];
  std::array<double, 2> cell = {0.0, (double)p.receptors_per_cell[rpc_id]};
  sim.cell_bound[recSpot * N_TRACKED + track_idx].push_back(cell);
  return recSpot;
}

// Motility tracker update: move a cell from tracker[i][slot] to tracker[j][slot].
// Leaky branch: cell picked weighted by its bound[0].
static void motility_tracker_move(Sim& sim, int i, int j, int which, RandomPark* ran) {
  static const int cd8map[10] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID,
                                 ThID, ARThID, TrID, ARTrID};
  int track_idx;
  if (which <= 18) {
    track_idx = cd8map[which - 13]; // 13..18 -> species 1..6
  } else {
    track_idx = which - 32;         // 39..42 -> 7,8,9,10
  }
  auto& list = sim.cell_bound[i * N_TRACKED + track_idx];
  if (list.empty()) return;

  double gross = 0.0;
  for (auto& c : list) gross += c[0];
  double rv = ran->uniform();
  int random_cell = 0;
  if (gross > 0.0) {
    double target_frac = rv * gross;
    double acc = 0.0;
    for (size_t cc = 0; cc < list.size(); cc++) {
      random_cell = (int)cc;
      acc += list[cc][0];
      if (target_frac <= acc) break;
    }
  } else {
    random_cell = (int)(rv * list.size());
    if (random_cell >= (int)list.size()) random_cell = (int)list.size() - 1;
    if (random_cell < 0) random_cell = 0;
  }
  std::array<double, 2> cell = list[random_cell];
  list.erase(list.begin() + random_cell);
  sim.cell_bound[j * N_TRACKED + track_idx].push_back(cell);
}

// T-cell proliferation tracker update: duplicate a random cell at i, halve bound
static void proliferation_tracker_dup(Sim& sim, int i, int which, RandomPark* ran) {
  static const int cd8map6[6] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID};
  int track_idx = cd8map6[which - 1]; // 1..6 -> species 1..6

  auto& list = sim.cell_bound[i * N_TRACKED + track_idx];
  if (list.empty()) return;
  int rc = (int)(ran->uniform() * list.size());
  if (rc >= (int)list.size()) rc = (int)list.size() - 1;
  if (rc < 0) rc = 0;
  double saved = list[rc][0];
  std::array<double, 2> new_cell = { saved / 2.0, list[rc][1] };
  list[rc][0] = saved / 2.0;
  list.push_back(new_cell);
}

// T-cell death: erase random cell at i from tracker slot
static void death_tracker_erase(Sim& sim, int i, int which, RandomPark* ran) {
  static const int cd8map6[6] = {TeID, ARTeID, TpID, TmID, ARTpID, ARTmID};
  int cellid;
  if (which <= 26) cellid = cd8map6[which - 25]; // 25->TeID, 26->ARTeID
  else             cellid = cd8map6[which - 33]; // 35..38 -> Tp,Tm,ARTp,ARTm

  auto& list = sim.cell_bound[i * N_TRACKED + cellid];
  if (list.empty()) return;
  int rc = (int)(ran->uniform() * list.size());
  if (rc >= (int)list.size()) rc = (int)list.size() - 1;
  if (rc < 0) rc = 0;
  list.erase(list.begin() + rc);
}

// Cancer lysis: erase random cancer cell at i (local) or j (nbr)
static void lysis_tracker_erase(Sim& sim, int i, int j, int which, RandomPark* ran) {
  (void)which;
  bool local = (j == -1);
  int site = local ? i : j;
  auto& list = sim.cell_bound[site * N_TRACKED + 0]; // canID slot=0
  if (list.empty()) return;
  int rc = (int)(ran->uniform() * list.size());
  if (rc >= (int)list.size()) rc = (int)list.size() - 1;
  if (rc < 0) rc = 0;
  list.erase(list.begin() + rc);
}

// Tp differentiation (19, 20, 59): move a Tp/ARTp cell at i to Tm/ARTm slot
static void differentiation_tracker_move(Sim& sim, int i, int which, RandomPark* ran) {
  static const int src_map[2]  = {TpID, ARTpID};
  static const int dst_map[2]  = {TmID, ARTmID};
  int cellid, target;
  if (which <= 20) {
    cellid = src_map[which - 19];
    target = dst_map[which - 19];
  } else { // which == 59: Tp -> ARTm
    cellid = src_map[0];
    target = dst_map[1];
  }
  auto& src_list = sim.cell_bound[i * N_TRACKED + cellid];
  if (src_list.empty()) return;
  int rc = (int)(ran->uniform() * src_list.size());
  if (rc >= (int)src_list.size()) rc = (int)src_list.size() - 1;
  if (rc < 0) rc = 0;
  std::array<double, 2> cell = src_list[rc];
  src_list.erase(src_list.begin() + rc);
  sim.cell_bound[i * N_TRACKED + target].push_back(cell);
}

// Tm exhaustion (21-24, 51-58): move Tm/ARTm cell at (i or j) to Te/ARTe
static void exhaustion_tracker_move(Sim& sim, int i, int j, int which, RandomPark* ran) {
  bool local = (j == -1);
  int site = local ? i : j;
  static const int src_map[2] = {TmID, ARTmID};
  static const int dst_map[2] = {TeID, ARTeID};
  int cellid, target;
  if (which <= 24) {
    int idx = (which - 21) / 2;   // 21,22->0; 23,24->1
    cellid = src_map[idx];
    target = dst_map[idx];
  } else {
    int arP = ((which - 51) / 2) % 2; // 51,52->0; 53,54->1; 55,56->0; 57,58->1
    cellid = src_map[arP];
    target = dst_map[arP];
  }
  auto& src_list = sim.cell_bound[site * N_TRACKED + cellid];
  if (src_list.empty()) return;
  int rc = (int)(ran->uniform() * src_list.size());
  if (rc >= (int)src_list.size()) rc = (int)src_list.size() - 1;
  if (rc < 0) rc = 0;
  std::array<double, 2> cell = src_list[rc];
  src_list.erase(src_list.begin() + rc);
  sim.cell_bound[site * N_TRACKED + target].push_back(cell);
}

// Cancer bias potential recompute (source lines 1813-1861)
// PRESERVES SOURCE BUGS:
//   - 0 < nx (strict) but 0 <= ny (non-strict)
//   - average_derivative is NOT reset between axis iterations (it mutates)
static void recompute_cancer_bias_potential(Sim& sim) {
  for (int cpx = 0; cpx < CP_MAX; cpx++) {
    for (int cpy = 0; cpy < CP_MAY; cpy++) {
      sim.cancer_bias_potential[cpx][cpy][0] = 0.0;
      sim.cancer_bias_potential[cpx][cpy][1] = 0.0;
      sim.discretized_cancer[cpx][cpy] = 0;
    }
  }
  int max_cancer = 0;
  for (int popy = 0; popy < LATROW; popy++) {
    for (int popx = 0; popx < LATCOL; popx++) {
      int pcx = popx / CANCER_POTENTIAL_CONVERTER;
      int pcy = popy / CANCER_POTENTIAL_CONVERTER;
      sim.discretized_cancer[pcx][pcy] += sim.population[canID][popx + LATCOL * popy];
      if (sim.discretized_cancer[pcx][pcy] > max_cancer)
        max_cancer = sim.discretized_cancer[pcx][pcy];
    }
  }
  if (max_cancer == 0) return;
  for (int popy = 0; popy < LATROW; popy++) {
    for (int popx = 0; popx < LATCOL; popx++) {
      int pcx = popx / CANCER_POTENTIAL_CONVERTER;
      int pcy = popy / CANCER_POTENTIAL_CONVERTER;
      int number_derivatives = 0;
      double average_derivative = 0.0;
      for (int axis = 0; axis < 2; axis++) {
        for (int ns = 0; ns < 2; ns++) {
          int nx = pcx + (1 - axis) * (2*ns - 1);
          int ny = pcy + axis * (2*ns - 1);
          if (0 < nx && nx < CP_MAX && 0 <= ny && ny < CP_MAY) {
            number_derivatives++;
            int lr = 2*ns - 1;
            average_derivative += (double)lr
              * (double)(sim.discretized_cancer[nx][ny] - sim.discretized_cancer[pcx][pcy])
              / (double)max_cancer;
          }
        }
        if (number_derivatives > 0) {
          average_derivative = average_derivative / (double)number_derivatives;
        }
        sim.cancer_bias_potential[pcx][pcy][axis] = (average_derivative + 1.0) / 2.0;
      }
    }
  }
}

static void site_event(Sim& sim, int i, RandomPark* ran) {
  // 1) Choose event from linked list (source lines 345-359)
  int isite_prop_index = i; // single-proc: i2site[i] = i
  double threshold = ran->uniform() * sim.propensity[isite_prop_index];
  double proball = 0.0;
  int ievent = sim.firstevent[i];
  while (ievent >= 0) {
    proball += sim.events[ievent].propensity;
    if (proball >= threshold) break;
    ievent = sim.events[ievent].next;
  }
  if (ievent < 0) return; // no events; shouldn't happen if propensity>0

  int rstyle = sim.events[ievent].style;
  int which  = sim.events[ievent].which;
  int j      = sim.events[ievent].jpartner;

  // 2) Apply population deltas (source lines 373-395)
  for (int s = 0; s < MAX_SPECIES; s++) {
    sim.population[s][i] += sim.rxn[which].localDeltaPop[s];
  }
  if (rstyle == NBR) {
    for (int s = 0; s < MAX_SPECIES; s++) {
      sim.population[s][j] += sim.rxn[which].nbrDeltaPop[s];
    }
  }

  // Track cancer proliferation destination and tcellPos (for touched-site updates)
  sim.tcellPos.clear();
  int cancerPos1 = -1;

  // 3) Reaction-specific logic
  if (which == 0) {
    // Cancer proliferation: ring search + place daughter
    double r1 = ran->uniform();
    cancerPos1 = cancer_division_place(sim, i, r1, ran, sim.tcellPos);
    sim.cancerPos1 = cancerPos1;
  }
  else if (7 <= which && which <= 12) {
    // Recruitment
    sim.recSpot = recruitment_place(sim, which, ran);
  }
  else if ((13 <= which && which <= 18) || (39 <= which && which <= 42)) {
    // Motility: move cell in tracker from i to j
    if (rstyle == NBR && j >= 0) {
      motility_tracker_move(sim, i, j, which, ran);
    }
  }

  // Tracker updates for all other events (leaky branch, source lines 1597-1723)
  if (1 <= which && which <= 6) {
    proliferation_tracker_dup(sim, i, which, ran);
  } else if (which == 25 || which == 26 || (35 <= which && which <= 38)) {
    death_tracker_erase(sim, i, which, ran);
  } else if (27 <= which && which <= 34) {
    lysis_tracker_erase(sim, i, j, which, ran);
  } else if (which == 19 || which == 20 || which == 59) {
    differentiation_tracker_move(sim, i, which, ran);
  } else if ((21 <= which && which <= 24) || (51 <= which && which <= 58)) {
    exhaustion_tracker_move(sim, i, j, which, ran);
  }
  // (43..50 Th-induced differentiation, all rate 0, tracker move similar to
  // 19/20 but never fires in this config)

  // 4) proliferation_update_bool triggers
  if ((1 <= which && which <= 12) || which == 25 || which == 26
      || (35 <= which && which <= 38)) {
    sim.proliferation_update_bool = true;
  }
  if (which == 0 || (27 <= which && which <= 34)) {
    sim.proliferation_update_bool = true;
    // Cancer bias potential recompute
    recompute_cancer_bias_potential(sim);
  }
  if (sim.proliferation_update_bool) {
    // Recompute tcell_total_population (source line 1866-1872)
    long long tot = 0;
    for (int ic = 0; ic < NSITES; ic++) {
      tot += sim.population[TpID][ic] + sim.population[TmID][ic] + sim.population[TeID][ic]
           + sim.population[ARTpID][ic] + sim.population[ARTmID][ic] + sim.population[ARTeID][ic];
    }
    if (tot > sim.p.max_T_cells) tot = sim.p.max_T_cells;
    sim.tcell_total_population = (int)tot;
  }

  sim.rxn_count[which]++;

  // 5) Touched-site propensity refresh (source lines 1878-2128)
  std::vector<int> esites_v;
  auto add_touched = [&](int m) {
    if (m < 0 || m >= NSITES) return;
    if (!sim.echeck[m]) {
      sim.propensity[m] = site_propensity(sim, m);
      esites_v.push_back(m);
      sim.echeck[m] = 1;
    }
  };

  if (sim.p.update_all_always) {
    for (int nx_ = 0; nx_ < NSITES; nx_++) add_touched(nx_);
  } else {
    add_touched(i);
    for (int nn = 0; nn < 4; nn++) add_touched(sim.neighbor[i][nn]);
    if (rstyle == NBR) {
      for (int nn = 0; nn < 4; nn++) add_touched(sim.neighbor[j][nn]);
    }

    // Cancer/recruit ring updates (source lines 1922-2029)
    if (which == 0 || (7 <= which && which <= 12) || (27 <= which && which <= 34)) {
      std::vector<int> cSites;

      if (which == 0 && cancerPos1 >= 0) {
        cSites.push_back(cancerPos1);
        for (int t : sim.tcellPos) cSites.push_back(t);
        // Expand rings around cancerPos1 (depth 3) and around each tcellPos (depth 1)
        for (int ci = 0; ci < (int)sim.tcellPos.size() + 1; ci++) {
          int cx, cy;
          int width = 3;
          if (ci == 0) {
            cx = cancerPos1 % LATCOL;
            cy = cancerPos1 / LATCOL;
          } else {
            cx = sim.tcellPos[ci - 1] % LATCOL;
            cy = sim.tcellPos[ci - 1] / LATCOL;
          }
          for (int lo = 0; lo < NEIGH_LIM; lo++) {
            int hw = (width - 1) / 2;
            for (int ly = 0; ly < 2; ly++) {
              for (int lx = 0; lx < width; lx++) {
                int t1x = wrap(cx + lx - hw, LATCOL);
                int t1y = wrap(cy + (2*ly - 1) * hw, LATCOL);
                int t2x = wrap(cx + (2*ly - 1) * hw, LATCOL);
                int t2y = wrap(cy + lx - hw, LATCOL);
                cSites.push_back(t1x + LATCOL * t1y);
                cSites.push_back(t2x + LATCOL * t2y);
              }
            }
            width += 2;
            if (ci > 0) break;
          }
        }
      } else if (27 <= which && which <= 34) {
        int cancPos = (which % 2 == 0) ? j : i;
        int cx = cancPos % LATCOL;
        int cy = cancPos / LATCOL;
        int width = 3;
        for (int lo = 0; lo < NEIGH_LIM; lo++) {
          int hw = (width - 1) / 2;
          for (int ly = 0; ly < 2; ly++) {
            for (int lx = 0; lx < width; lx++) {
              int t1x = wrap(cx + lx - hw, LATCOL);
              int t1y = wrap(cy + (2*ly - 1) * hw, LATCOL);
              int t2x = wrap(cx + (2*ly - 1) * hw, LATCOL);
              int t2y = wrap(cy + lx - hw, LATCOL);
              cSites.push_back(t1x + LATCOL * t1y);
              cSites.push_back(t2x + LATCOL * t2y);
            }
          }
          width += 2;
        }
      } else if (7 <= which && which <= 12 && sim.recSpot >= 0) {
        int rx = sim.recSpot % LATCOL;
        int ry = sim.recSpot / LATCOL;
        for (int ty = 0; ty < 3; ty++) {
          for (int tx = 0; tx < 3; tx++) {
            int tgx = wrap(rx + tx - 1, LATCOL);
            int tgy = wrap(ry + ty - 1, LATCOL);
            cSites.push_back(tgx + LATCOL * tgy);
          }
        }
      }

      std::sort(cSites.begin(), cSites.end());
      cSites.erase(std::unique(cSites.begin(), cSites.end()), cSites.end());
      if (!sim.proliferation_update_bool) {
        for (int m : cSites) add_touched(m);
      }
    }

    // Full T-cell voxel refresh
    if (sim.proliferation_update_bool) {
      for (int la = 0; la < NSITES; la++) {
        if (sim.population[TpID][la] > 0 || sim.population[TmID][la] > 0
            || sim.population[TeID][la] > 0 || sim.population[ARTpID][la] > 0
            || sim.population[ARTmID][la] > 0 || sim.population[ARTeID][la] > 0
            || sim.population[recID][la] > 0) {
          add_touched(la);
        }
      }
      sim.proliferation_update_bool = false;
    }
    // Full cancer voxel refresh
    if (sim.cancer_proliferation_update_bool) {
      for (int la = 0; la < NSITES; la++) {
        if (sim.population[canID][la] > 0) add_touched(la);
      }
      sim.cancer_proliferation_update_bool = false;
    }
  }

  // Push touched sites to solver
  if (!esites_v.empty()) {
    sim.solve->update((int)esites_v.size(), esites_v.data(), sim.propensity.data());
  }
  // Clear echeck
  for (int m : esites_v) sim.echeck[m] = 0;
}

// ============================================================================
// SECTION 9: INITIALIZATION
// ============================================================================

// Load the init file (SPPARKS-style Values section).
// File format:
//   line 1: "Site file written by dump sites 1 command"
//   line 2: (blank)
//   line 3: "id i1 i2 i3 i4 i5 i6 i7 d1 i8 i9 i10 i11 i12 d2 i13 i14 values"
//   line 4: (blank)
//   line 5: "Values"
//   line 6: (blank)
//   lines 7+: "<id> <int fields...> <d1> <more ints> <d2> <more ints>"
static void load_init_file(Sim& sim, const std::string& path) {
  std::ifstream fin(path);
  if (!fin) {
    fprintf(stderr, "Cannot open init file: %s\n", path.c_str());
    std::exit(1);
  }
  std::string line;
  // Read header column line to find field order.
  // File has: id i1 i2 i3 i4 i5 i6 i7 d1 i8 i9 i10 i11 i12 d2 i13 i14 values
  //
  // We rely on the fixed order from d14_M_4_SPPARKS_init.
  bool in_values = false;
  int lines_read = 0;
  while (std::getline(fin, line)) {
    if (line.find("Values") != std::string::npos && !in_values && lines_read > 2) {
      in_values = true;
      continue;
    }
    lines_read++;
    if (!in_values) continue;
    if (line.empty()) continue;

    std::istringstream iss(line);
    long long id;
    if (!(iss >> id)) continue;
    int site = (int)(id - 1);
    if (site < 0 || site >= NSITES) continue;

    int i1, i2, i3, i4, i5, i6, i7;
    double d1;
    int i8, i9, i10, i11, i12;
    double d2;
    int i13, i14;
    iss >> i1 >> i2 >> i3 >> i4 >> i5 >> i6 >> i7
        >> d1 >> i8 >> i9 >> i10 >> i11 >> i12
        >> d2 >> i13 >> i14;
    if (!iss) continue;

    sim.population[canID][site]   = i1;  // A
    sim.population[TpID][site]    = i2;  // B
    sim.population[TmID][site]    = i3;  // C
    sim.population[TeID][site]    = i4;  // D
    sim.population[ARTpID][site]  = i5;  // E (GONE)
    sim.population[ARTmID][site]  = i6;  // F
    sim.population[ARTeID][site]  = i7;  // G
    sim.ifng_field_discretized_free [site] = d1;
    sim.population[recID][site]   = i8;  // I
    sim.population[errID][site]   = i9;  // J
    sim.population[deadID][site]  = i10; // K
    sim.population[latEID][site]  = i11; // L
    sim.population[boundID][site] = i12; // M
    sim.ifng_field_discretized_bound[site] = d2;
    sim.population[ThID][site]    = i13; // O = Th
    sim.population[ARThID][site]  = i14; // P = ARTh
    // TrID (14), ARTrID (15), macS (16), macT (17) stay 0
  }
}

// Initialize per-cell tracker from population + init d2 field.
// Runs once, matches app_rxn_diff_custom.cpp line 103.
static void initialize_tracker_from_populations(Sim& sim) {
  sim.cell_bound.assign((size_t)NSITES * N_TRACKED, {});
  for (int y = 0; y < LATROW; y++) {
    for (int x = 0; x < LATCOL; x++) {
      int loc = x + LATCOL * y;
      double bound_frac = sim.ifng_field_discretized_bound[loc];
      if (bound_frac < 0.0) bound_frac = 0.0; // guard: init has 0 everywhere
      for (int arn = 0; arn < N_TRACKED; arn++) {
        int cid = TRK_CELL_ID[arn];
        int pop = sim.population[cid][loc];
        int rpc = sim.p.receptors_per_cell[RECEPTORS_PER_CELL_ID[arn]];
        std::array<double, 2> cell = { bound_frac * (double)rpc, (double)rpc };
        auto& list = sim.cell_bound[(size_t)loc * N_TRACKED + arn];
        for (int c = 0; c < pop; c++) list.push_back(cell);
      }
    }
  }
}

// Initialize coarse ifng_field from the fine ifng_field_discretized values.
// Actual real SPPARKS starts with ifng_field=0 everywhere (line 217 of ctor),
// so this just zeros it — but we keep the code path in case someone loads a
// warm-start init file.
static void initialize_ifng_field(Sim& sim) {
  sim.ifng_field.assign(IFNG_NSITES, 0.0);
  sim.ifng_field_copy.assign(IFNG_NSITES, 0.0);
  // Aggregate the fine free-IFNG field into the coarse field (typically all 0)
  for (int y = 0; y < LATROW; y++) {
    for (int x = 0; x < LATCOL; x++) {
      int cx = x / IFNG_LATTICE_CELL_CONVERTER;
      int cy = y / IFNG_LATTICE_CELL_CONVERTER;
      // Source stores the DIVIDED value in ifng_field_discretized[0][loc],
      // so aggregating back means summing per-fine-voxel values (which are
      // 1/9 of the coarse). Here init is all 0 anyway, so nothing happens.
      sim.ifng_field[cx + IFNG_COL * cy] += sim.ifng_field_discretized_free[x + LATCOL * y];
    }
  }
}

// Full setup routine: compute initial propensities, populate SolveTree,
// initialize bias potential.
static void initial_setup(Sim& sim) {
  sim.i2site.assign(NSITES, 0);
  for (int k = 0; k < NSITES; k++) sim.i2site[k] = k;
  sim.esites.reserve(NSITES);
  sim.echeck.assign(NSITES, 0);
  sim.firstevent.assign(NSITES, -1);
  sim.events.reserve(1000000);
  sim.events.resize(100000);
  for (size_t m = 0; m < sim.events.size(); m++) sim.events[m].next = (int)m + 1;
  sim.events.back().next = -1;
  sim.freeevent = 0;
  sim.nevents_active = 0;

  sim.propensity.assign(NSITES, 0.0);
  sim.bound_ifng.assign(NSITES, {0.0, 0.0, 0.0});
  // Set d2 = -1 marker for all voxels initially without cells (matches source
  // convention). For voxels with cells, bound_ifng[0] = 0 will be overwritten
  // by update_ifng_field on the first step.

  // Initialize cancer bias potential to 0.5 (neutral). Real SPPARKS leaves
  // this uninitialized (undefined behavior); we choose 0.5 which gives
  // bias_factor=1 for all directions. Once the first cancer event fires,
  // recompute_cancer_bias_potential() overwrites it.
  for (int a = 0; a < CP_MAX; a++)
    for (int b = 0; b < CP_MAY; b++)
      for (int c = 0; c < 2; c++)
        sim.cancer_bias_potential[a][b][c] = 0.5;

  // Compute initial propensities. This also triggers the first-call
  // tracker initialization inside custom_multiplier if we were using SPPARKS's
  // lazy init — but we've already done it explicitly above.
  for (int k = 0; k < NSITES; k++) {
    sim.propensity[k] = site_propensity(sim, k);
  }
  sim.solve->init(NSITES, sim.propensity.data());
}

// ============================================================================
// SECTION 10: MAIN KMC LOOP
// Port of AppLattice::iterate_kmc_global (app_lattice.cpp lines 653-727)
// ============================================================================

// Output hook: called with the state exactly at each requested time point.
// Writing output never draws random numbers or alters the event sequence,
// so adding/removing output times does not change the trajectory.
struct OutputPlan {
  std::vector<double> report_times;    // native units, ascending
  std::vector<double> snapshot_times;  // native units, ascending
  size_t ri = 0, si = 0;
};
static void write_report_row(std::FILE* fp, Sim& sim, double t);
static void write_snapshot(Sim& sim, double t, const std::string& prefix);

static void emit_outputs_upto(Sim& sim, double t_next, OutputPlan& plan,
                              std::FILE* fp, const std::string& prefix) {
  // Emit every output whose time is < t_next (state is unchanged until then).
  while (plan.ri < plan.report_times.size() && plan.report_times[plan.ri] < t_next) {
    write_report_row(fp, sim, plan.report_times[plan.ri]);
    plan.ri++;
  }
  while (plan.si < plan.snapshot_times.size() && plan.snapshot_times[plan.si] < t_next) {
    write_snapshot(sim, plan.snapshot_times[plan.si], prefix);
    plan.si++;
  }
}

// Same loop as AppLattice::iterate_kmc_global (app_lattice.cpp 653-727).
static void iterate_kmc(Sim& sim, double stoptime, OutputPlan& plan,
                        std::FILE* fp, const std::string& prefix) {
  int done = 0;
  int ifng_update_count = 0;
  long long naccept = 0;
  double dt_step;
  double next_progress = 0.0;

  while (!done) {
    int isite = sim.solve->event(&dt_step);

    if (dt_step + sim.ifng_time_keeper > sim.p.ifng_time_step) {
      if (sim.time <= stoptime) {
        double t_new = sim.time + sim.p.ifng_time_step - sim.ifng_time_keeper;
        emit_outputs_upto(sim, t_new, plan, fp, prefix);
        update_ifng_field(sim, sim.p.ifng_time_step);
        sim.proliferation_update_bool = sim.p.ifng_cytotoxicity_augmentation;
        sim.cancer_proliferation_update_bool = sim.p.ifng_cancer_proliferation_limitation;
        sim.time = t_new;
        sim.ifng_time_keeper = 0.0;
        ifng_update_count++;

        for (int n = 0; n < NSITES; n++) sim.propensity[n] = site_propensity(sim, n);
        sim.solve->reload_all(NSITES, sim.propensity.data());
      } else {
        done = 1;
      }
    } else if (isite >= 0) {
      if (sim.time <= stoptime) {
        emit_outputs_upto(sim, sim.time + dt_step, plan, fp, prefix);
        sim.ifng_time_keeper += dt_step;
        sim.time += dt_step;
        site_event(sim, isite, sim.ranapp);
        naccept++;
      } else {
        done = 1;
      }
    } else {
      done = 1;
    }
    if (sim.time >= next_progress) {
      fprintf(stderr, "  t = %7.1f units (%.2f days)  cancer events so far: %lld\n",
              sim.time, sim.time / 144.0, (long long)sim.rxn_count[0]);
      next_progress += 24.0; // every 4 simulated hours
    }
  }
  // Anything requested at or before stoptime that is still pending
  emit_outputs_upto(sim, stoptime + 1e-9, plan, fp, prefix);
  fprintf(stderr, "run done: t=%.4f units  naccept=%lld  ifng_updates=%d\n",
          sim.time, naccept, ifng_update_count);
}

// ============================================================================
// SECTION 11: OUTPUT
// ============================================================================
static void write_report_header(std::FILE* fp) {
  fprintf(fp, "time_units,time_days,cancer,Tp,Tm,Te,ARpos_PE,ARpos_Te,CD4,ARpos_CD4,"
              "CD8_total,free_ifng,bound_ifng,mean_boundfrac_cancer,mean_boundfrac_Tcell\n");
}

static void write_report_row(std::FILE* fp, Sim& sim, double t) {
  long long n[MAX_SPECIES] = {0};
  for (int k = 0; k < NSITES; k++)
    for (int s = 0; s < MAX_SPECIES; s++) n[s] += sim.population[s][k];
  double free_ifng = 0.0;
  for (int c = 0; c < IFNG_NSITES; c++) free_ifng += sim.ifng_field[c];
  double bound = 0.0, bf_c = 0.0, bf_t = 0.0;
  long long nc = 0, nt = 0;
  for (int site = 0; site < NSITES; site++) {
    for (int arn = 0; arn < N_TRACKED; arn++) {
      for (auto& cell : sim.cell_bound[(size_t)site * N_TRACKED + arn]) {
        bound += cell[0];
        double f = cell[1] > 0 ? cell[0] / cell[1] : 0.0;
        if (arn == 0) { bf_c += f; nc++; } else { bf_t += f; nt++; }
      }
    }
  }
  long long cd8 = n[TpID] + n[TmID] + n[TeID] + n[ARTpID] + n[ARTmID] + n[ARTeID];
  fprintf(fp, "%.4f,%.6f,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%.6g,%.6g,%.6g,%.6g\n",
          t, t / 144.0, n[canID], n[TpID], n[TmID], n[TeID], n[ARTmID], n[ARTeID],
          n[ThID], n[ARThID], cd8, free_ifng, bound,
          nc ? bf_c / nc : 0.0, nt ? bf_t / nt : 0.0);
  fflush(fp);
}

// Spatial snapshot: one row per occupied or IFNG-carrying fine site.
//   x y cancer Tp Tm Te ARpos_PE ARpos_Te CD4 ARpos_CD4 free_ifng_per_site bound_frac
static void write_snapshot(Sim& sim, double t, const std::string& prefix) {
  char name[512];
  std::snprintf(name, sizeof(name), "%s_snap_t%04d.txt", prefix.c_str(), (int)std::lround(t));
  std::FILE* f = std::fopen(name, "w");
  if (!f) { fprintf(stderr, "cannot write %s\n", name); return; }
  fprintf(f, "# time_units %.4f time_days %.4f\n", t, t / 144.0);
  fprintf(f, "# x y cancer Tp Tm Te ARpos_PE ARpos_Te CD4 ARpos_CD4 free_ifng bound_frac\n");
  for (int y = 0; y < LATROW; y++) {
    for (int x = 0; x < LATCOL; x++) {
      int s = x + LATCOL * y;
      int c[8] = {sim.population[canID][s], sim.population[TpID][s], sim.population[TmID][s],
                  sim.population[TeID][s], sim.population[ARTmID][s], sim.population[ARTeID][s],
                  sim.population[ThID][s], sim.population[ARThID][s]};
      double fi = sim.ifng_field[(x / 3) + IFNG_COL * (y / 3)] / 9.0;
      bool any = fi > 0.0;
      for (int k = 0; k < 8; k++) any = any || c[k] > 0;
      if (!any) continue;
      fprintf(f, "%d %d %d %d %d %d %d %d %d %d %.5g %.4g\n", x, y,
              c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], fi,
              sim.bound_ifng[s][0]);
    }
  }
  std::fclose(f);
  fprintf(stderr, "  wrote %s\n", name);
}

// ============================================================================
// SECTION 12: MAIN
// ============================================================================
static std::vector<double> parse_list(const std::string& csv) {
  std::vector<double> v;
  std::stringstream ss(csv);
  std::string tok;
  while (std::getline(ss, tok, ',')) if (!tok.empty()) v.push_back(std::atof(tok.c_str()));
  std::sort(v.begin(), v.end());
  return v;
}

static void usage(const char* p) {
  fprintf(stderr,
    "Usage: %s --init <SPPARKS_init file> --sex male|female [options]\n"
    "  --days D               simulated days (default 7 = day 14 -> day 21; 1 day = 144 units)\n"
    "  --seed N               RNG seed (default 8377)\n"
    "  --out PREFIX           output prefix (default run): writes PREFIX_timecourse.csv\n"
    "  --report-hours H       time-course sampling interval in hours (default 6)\n"
    "  --snapshot-days LIST   spatial snapshots, e.g. 0,3.5,7 (default none)\n"
    "  --no-sequestration     AR+ cells do not internalize bound IFNG (ablation)\n"
    "  --te-death-weekly      exhausted CD8 death 1/week (rates table) instead of\n"
    "                         1/day (in.full default)\n"
    "Sex 'female' applies the female rate changes (no AR+ recruitment or Tp->AR+PE,\n"
    "AR- recruitment x1/0.63, Tp->Tm x2). Pair it with d14_M_4_madeF_SPPARKS_init\n"
    "(or d14_F_4_SPPARKS_init for the real female slide).\n", p);
}

int main(int argc, char** argv) {
  std::string init_path, prefix = "run", snapCSV;
  double days = 7.0, report_hours = 6.0;
  int seed = 8377;
  Sex sex = MALE;
  bool sex_given = false, no_seq = false, te_weekly = false;

  for (int a = 1; a < argc; a++) {
    std::string f = argv[a];
    auto need = [&](void) -> std::string {
      if (a + 1 >= argc) { usage(argv[0]); std::exit(1); }
      return std::string(argv[++a]);
    };
    if (f == "--init") init_path = need();
    else if (f == "--sex") { std::string s = need(); sex = (s == "female") ? FEMALE : MALE; sex_given = (s == "female" || s == "male"); }
    else if (f == "--days") days = std::atof(need().c_str());
    else if (f == "--seed") seed = std::atoi(need().c_str());
    else if (f == "--out") prefix = need();
    else if (f == "--report-hours") report_hours = std::atof(need().c_str());
    else if (f == "--snapshot-days") snapCSV = need();
    else if (f == "--no-sequestration") no_seq = true;
    else if (f == "--te-death-weekly") te_weekly = true;
    else { fprintf(stderr, "unknown option %s\n", f.c_str()); usage(argv[0]); return 1; }
  }
  if (init_path.empty() || !sex_given) { usage(argv[0]); return 1; }

  const double UNITS_PER_DAY = 144.0;           // 1 native unit = 600 s
  double stoptime = days * UNITS_PER_DAY;

  OutputPlan plan;
  double dt_rep = report_hours * 6.0;           // hours -> units
  for (double t = 0.0; t <= stoptime + 1e-9; t += dt_rep) plan.report_times.push_back(t);
  if (plan.report_times.back() < stoptime - 1e-9) plan.report_times.push_back(stoptime);
  for (double d : parse_list(snapCSV)) plan.snapshot_times.push_back(d * UNITS_PER_DAY);

  fprintf(stderr, "sex_bias_sim_v7: init=%s sex=%s days=%g (%.0f units) seed=%d%s%s\n",
          init_path.c_str(), sex == FEMALE ? "female" : "male", days, stoptime, seed,
          no_seq ? " NO-SEQUESTRATION" : "", te_weekly ? " Te-death=1/week" : "");

  Sim* simp = new Sim();
  Sim& sim = *simp;
  sim.p.init_from_infull();
  sim.p.no_sequestration = no_seq;
  if (te_weekly) sim.p.te_death_rate = 1.0 / (7.0 * UNITS_PER_DAY);
  sim.sex = sex;

  sim.ifng_field_discretized_free .assign(NSITES, 0.0);
  sim.ifng_field_discretized_bound.assign(NSITES, 0.0);
  sim.build_neighbors();
  build_reaction_table(sim, sex);
  build_reactant_terms(sim);

  // RNG initialisation identical to SPPARKS (RanMars seeds two Park-Miller streams)
  RanMars mars;
  mars.init(seed);
  sim.solve = new SolveTree();
  sim.solve->random = new RandomPark(mars.uniform());
  sim.solve->random->reset(mars.uniform(), 0, 100);
  sim.ranapp = new RandomPark(mars.uniform());
  sim.ranapp->reset(mars.uniform(), 0, 100);

  load_init_file(sim, init_path);
  initialize_ifng_field(sim);
  initialize_tracker_from_populations(sim);
  sim.mapped_population_to_ifng_tracker = true;
  initial_setup(sim);

  std::string csv = prefix + "_timecourse.csv";
  std::FILE* fout = std::fopen(csv.c_str(), "w");
  if (!fout) { fprintf(stderr, "cannot write %s\n", csv.c_str()); return 1; }
  write_report_header(fout);

  iterate_kmc(sim, stoptime, plan, fout, prefix);
  std::fclose(fout);
  fprintf(stderr, "Done. Time course -> %s\n", csv.c_str());
  return 0;
}
