// ============================================================================
// ifng_effective_diffusion.cpp
//
// Simplified model from thesis Ch. 4.3 (Fig 4.8) used to isolate how
// IFNG-IFNGR binding slows the spread of IFNG into the tumour nest.
//
// Same IFNG physics as sex_bias_sim_v7 / SPPARKS app_lattice.cpp:
//   * 600x600 cell lattice (10 um), IFNG solved on 200x200 pixels (30 um),
//     explicit finite differences, D = 20 um^2/s, time step 9 s,
//     non-periodic (zero) boundary, no free-IFNG decay;
//   * after each diffusion step, free IFNG and IFNGR in each IFNG pixel are
//     put at binding equilibrium (K = 100 molecules per IFNG pixel) and the
//     bound IFNG is shared among the cells in proportion to receptors;
//   * secretion by an AR- T cell = 0.1/s * (1.5*Hill(g;0.2,4) + 0.5), i.e.
//     k_max*(0.25 + 0.75 Hill) with k_max = 0.2/s, g = its bound fraction.
//
// Simplified geometry (as in the thesis):
//   * one cancer cell in every cell site (9 per IFNG pixel, 1000 IFNGR each,
//     R0 = 9000 per IFNG pixel), no proliferation, no killing, no motion;
//   * one fixed IFNG source (an AR- mature T cell, 1000 IFNGR) at the centre
//     secreting at 20x the single-cell rate, representing a small immobile
//     T-cell cluster.
//
// A second run with no receptors ("bare" diffusion, D = 20 um^2/s) is used to
// calibrate how the front position relates to D for this source geometry.
//
// Output:
//   <prefix>_profiles.csv : time_h, x_um, free_ifng_per_site   (bound run)
//   <prefix>_front.csv    : run, time_h, front_um, front_um2
//
// Build: g++ -O3 -std=c++17 ifng_effective_diffusion.cpp -o ifng_effective_diffusion
// Run:   ./ifng_effective_diffusion [prefix] [days] [front_fraction]
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

static const int    N        = 200;       // IFNG pixels per side
static const double A        = 30.0;      // IFNG pixel size (um)
static const double D        = 20.0;      // um^2/s
static const double DT       = 9.0;       // s  (0.015 native units x 600 s)
static const double KD       = 100.0;     // molecules per IFNG pixel (unbinding/binding = 100/1)
static const double REC      = 1000.0;    // IFNGR per cell
static const double K_PROD   = 0.1;       // ifng_production (per s)
static const double MAX_FB   = 2.0;       // max_ifng_production
static const double G_FB     = 0.2;       // half_max_ifng_feedback
static const int    HILL_N   = 4;
static const double SOURCE_X = 20.0;      // source strength multiplier

static double hill(double x, double K, int n) {
  double xn = std::pow(x, n), Kn = std::pow(K, n);
  return xn / (Kn + xn + 1e-30);
}

struct Result { std::vector<double> t_h, front; };

// with_receptors=false -> bare diffusion calibration run
static Result run(bool with_receptors, double days, double frac,
                  std::FILE* profiles, double stop_front_um) {
  const int c = N / 2;                         // source pixel (centre)
  std::vector<double> f(N * N, 0.0), fc(N * N, 0.0), B(N * N, 0.0), R(N * N, 0.0);
  if (with_receptors) {
    for (int k = 0; k < N * N; k++) R[k] = 9.0 * REC;   // 9 cancer cells / pixel
    R[c + N * c] += REC;                                // the source T cell
  }
  Result res;
  const long nsteps = (long)std::llround(days * 86400.0 / DT);
  const long front_every = (long)std::llround(1800.0 / DT);   // every 0.5 h
  const long prof_every  = (long)std::llround(86400.0 / DT);  // every 24 h
  const double lapfac = D / (A * A);

  for (long step = 1; step <= nsteps; step++) {
    fc = f;
    for (int y = 0; y < N; y++) {
      for (int x = 0; x < N; x++) {
        int k = x + N * y;
        double nb = 0.0;
        if (x > 0)     nb += fc[k - 1];
        if (x < N - 1) nb += fc[k + 1];
        if (y > 0)     nb += fc[k - N];
        if (y < N - 1) nb += fc[k + N];
        double src = 0.0;
        if (x == c && y == c) {
          double g = (R[k] > 0.0) ? B[k] / R[k] : 0.0;   // source's own bound fraction
          src = SOURCE_X * K_PROD * ((MAX_FB - 0.5) * hill(g, G_FB, HILL_N) + 0.5);
        }
        double fn = f[k] + DT * (lapfac * (nb - 4.0 * fc[k]) + src);
        if (fn < 1e-10) fn = 0.0;
        if (R[k] > 0.0) {                               // equilibrium binding
          double tot = B[k] + fn;
          double a = R[k] + KD - tot;
          double ffree = -0.5 * (a - std::sqrt(a * a + 4.0 * KD * tot));
          B[k] = tot - ffree;
          fn = ffree;
        }
        f[k] = fn;
      }
    }

    double t_h = step * DT / 3600.0;
    if (step % front_every == 0) {
      // front along the row through the source: furthest x where
      // free IFNG >= frac * max (linear interpolation between pixels)
      const double* row = &f[N * c];
      double fmax = *std::max_element(row, row + N);
      double thr = frac * fmax, front = 0.0;
      for (int x = c; x < N - 1; x++) {
        if (row[x] >= thr && row[x + 1] < thr) {
          double s = (row[x] - thr) / (row[x] - row[x + 1]);
          front = (x - c + s) * A;
          break;
        }
      }
      res.t_h.push_back(t_h);
      res.front.push_back(front);
      if (front >= stop_front_um) break;
    }
    if (profiles && step % prof_every == 0) {
      for (int x = 0; x < N; x++)
        fprintf(profiles, "%.1f,%.1f,%.6g\n", t_h, (x - c) * A, f[x + N * c] / 9.0);
    }
  }
  return res;
}

int main(int argc, char** argv) {
  std::string prefix = argc > 1 ? argv[1] : "ifng_diffusion";
  double days = argc > 2 ? std::atof(argv[2]) : 7.0;
  double frac = argc > 3 ? std::atof(argv[3]) : 0.01;

  std::FILE* prof = std::fopen((prefix + "_profiles.csv").c_str(), "w");
  std::fprintf(prof, "time_h,x_um,free_ifng_per_site\n");
  std::fprintf(stderr, "bound run (R0 = 9000 per pixel, K = 100) for %g days...\n", days);
  Result rb = run(true, days, frac, prof, 2700.0);
  std::fclose(prof);

  std::fprintf(stderr, "bare-diffusion calibration run...\n");
  Result r0 = run(false, days, frac, nullptr, 2700.0);   // stops near the box edge

  std::FILE* fr = std::fopen((prefix + "_front.csv").c_str(), "w");
  std::fprintf(fr, "run,time_h,front_um,front_um2\n");
  for (size_t i = 0; i < rb.t_h.size(); i++)
    std::fprintf(fr, "bound,%.2f,%.3f,%.3f\n", rb.t_h[i], rb.front[i], rb.front[i] * rb.front[i]);
  for (size_t i = 0; i < r0.t_h.size(); i++)
    std::fprintf(fr, "bare,%.2f,%.3f,%.3f\n", r0.t_h[i], r0.front[i], r0.front[i] * r0.front[i]);
  std::fclose(fr);
  std::fprintf(stderr, "wrote %s_profiles.csv and %s_front.csv\n", prefix.c_str(), prefix.c_str());
  return 0;
}
