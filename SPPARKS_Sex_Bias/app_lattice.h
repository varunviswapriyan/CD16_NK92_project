/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef SPK_APP_LATTICE_H
#define SPK_APP_LATTICE_H

#include "stdio.h"
#include "app.h"

//Giuseppe import:
#include <iostream>
#include <vector>
#include <array>
//////////////////
//
#include <stack>

namespace SPPARKS_NS {

class AppLattice : public App {
  friend class CommLattice;

 public:
  int delpropensity;           // # of sites away needed to compute propensity
  int delevent;                // # of sites away affected by an event

  int nmax;                    // max # of sites per-site arrays can store
  int maxneigh;                // max neighbors of any site in entire system
  int *numneigh;               // # of neighbors of each site
  int **neighbor;              // local indices of neighbors of each site
  
  // GIUSEPPE'S CODE //////////////////////////////////////////////////////
  float  dt_diff; // how often to update dead_mem
  int diffCount;  // int count of updates to dead_mem              //Giuseppe's diffusion constants
  int twindow; // how long of a time to measure dead cells
  int **population; // variable dimension array two pointers allow you to construct the 2D array
  double **ifng_field_discretized; // the "population" array but for ifng
  bool mapped_population_to_ifng_tracker;
  std::vector<std::vector <double>> ***ifng_tracker; // [x][y][cellid][cell][bound, receptors] for advanced decay model: [recently absorbed, cell absorption rate]
  std::array<double, 3> **bound_ifng; // bound fraction, total bound, total receptors

  bool use_simple_decay_model;
  bool use_cell_dependent_decay_model;
  bool use_leaky_receptor_model;

//  double **bound_ifng;  // bound ifng per lattice site
//  int **receptors;  // number of receptors in each ifng lattice site
  int *receptors_per_cell; // number of receptors per cell in the system (maybe make an int[] later)
  std::vector<std::vector<int>> nDead_mem; // keep dead cancer counts stored
  int nDead; // number of dead cells to consider now
  int nDeadi; // approximate number of dead cells used at initial time? Maybe
  int nDeadCounti;
  bool proliferation_update_bool; // boolean that declares when next to update all CD8+ T cell proliferation rates across all cells because dead cell number update
  bool cancer_proliferation_update_bool; // update all cancer cells' proliferation rates
  bool ifng_cancer_proliferation_limitation; // if ifng limits cancer prolif
  bool ifng_cytotoxicity_augmentation; // if ifng augments t cell cytotoxicity
  double ifng_diffusion;
  double ifng_decay;      // linear decay rate of ifng
  double ifng_production;
  double ifng_binding;
  double ifng_unbinding;
  double dissC;
  double ifng_uptake;
  double ifng_ART_uptake;
  double ifng_feedback_saturation_threshold;
  double ifng_saturation_threshold;
  double bIFNG_decay;
  bool update_bIFNG_decay;

  double arn_cd8_consumption_rate; // rate at which ar- CD8+ T cells consume free IFNG
  double arp_cd8_consumption_rate; // rate at which ar+ CD8+ T cells consume free IFNG
  double cancer_consumption_rate;  // rate at which cancer cells consume free IFNG

//  [x][y][id]<vector><vector>
//  [x][y][id]<>

  int ifng_feedback_hill_degree;
  double max_ifng_production;
  double half_max_ifng_feedback;

  double ifng_time_step;
  double ifng_time_keeper;

  double ifng_lattice_cell_converter;
  int latt_const;
  int ifng_latt_const;

  double total_ifng;

  bool periodic_boundary_conditions;
  double **ifng_field;
  double **ifng_field_copy;
//  std::array<std::array<double, 200>, 200> ifng_field; // CHANGE THIS WHEN CHANGING SIZE OF FIELD!!!
//  std::array<std::array<double, 200>, 200> ifng_field_copy;

//  int ***src_field;
//  std::array<std::array<std::array<int, 3>, 200>, 200> src_field;

  int ifng_col;
  int ifng_row;

    //int tcell_total_population; // total t cell population


  class CommLattice *comm;

  AppLattice(class SPPARKS *, int, char **);
  virtual ~AppLattice();
  void input(char *, int, char **);
  void init();
  void setup();
  void iterate();
  virtual bigint memory_usage();

  void grow(int);
  void add_site(tagint, double, double, double);
  void add_ghost(tagint, double, double, double, int, int);
  void add_neighbors(int, int, char **);
  void add_values(int, char **);
  void add_value(int, int, int, char *);
  void print_connectivity();

  // pure virtual functions, must be defined in child class

  virtual void grow_app() = 0;
  virtual double site_energy(int) = 0;
  virtual void site_event_rejection(int, class RandomPark *) = 0;
  virtual double site_propensity(int) = 0;
  virtual void site_event(int, class RandomPark *) = 0;

  // virtual functions, may be overridden by child class

  virtual void input_app(char *, int, char **);
  virtual void init_app() {}
  virtual void setup_app() {}
  virtual void setup_end_app() {}
  virtual void *extract_app(char *) {return NULL;}

  virtual void push_new_site(int, int *, int, std::stack<int>*);
  virtual void push_connected_neighbors(int, int *, int, std::stack<int>*);
  virtual void connected_ghosts(int, int *, class Cluster *, int);

  virtual void app_update(double) {}

 protected:
  int me,nprocs;

  bigint naccept,nattempt;    // number of accepted and attempted events
  int nsweeps;                // number of sweeps performed
  double temperature,t_inverse;  // temperature settings
  double dt_sweep;            // rKMC time for nglobal attemped events
  double dt_rkmc;             // rKMC time for one pass thru all sectors
  double dt_kmc;              // KMC time for one pass thru all sectors
  double dt_step;             // KMC time for single global KMC step

  int allow_kmc;               // 1 if app supports KMC
  int allow_rejection;         // 1 if app supports rejection KMC
  int allow_masking;           // 1 if app supports rKMC masking
  int allow_app_update;        // 1 if app provides app_update()
  int numrandom;               // # of RN used by rejection routine

  int sweepflag;               // set if rejection KMC solver
  int sectorflag;              // 1 if partition my domain into sectors
  int nsector;                 // 1,2,4,8 = # of sectors
  int nsector_user;            // 0 if default, else 2,4,8
  int ncolors;                 // # of colors, depends on lattice
  int bothflag;                // 1 if both sectors and colors
  int app_update_only;         // 1 if skip KMC and rKMC updates

  class RandomPark *ranapp;    // RN generator for KMC and rejection KMC
  class RandomPark *ranstrict; // RN generator for per-site strict rKMC
  int *siteseeds;              // per-site seeds for ransite
  int *sitelist;               // randomized list of site indices

  bool Lmask;                  // masking on/off
  char *mask;                  // size of nlocal + nghost sites

  double Ladapt;               // adaptive sector time increments for KMC
  double tstop;                // requested time increment in sector
  double nstop;                // requested events per site in sector
  
  // Giuseppe added:
  int *esites;
  int latCol;
  //////////////////
                               // arrays for owned + ghost sites
  int *owner;                  // proc who owns the site
  int *index;                  // index of site on owning proc

  double *propensity;          // probabilities for each owned site
  int *i2site;                 // mapping of owned lattice to site index

  struct Set {                 // subset of lattice sites I own
    int nlocal;                // # of owned sites in set
    int nselect;               // # of selections from set for rKMC
    int nloop;                 // # of loops over set for rKMC
    int nborder;               // # of sites with non-set site as neighbor
    int *border;               // lattice index for each border site
    int *bsites;               // list of border sites to pass to solver
    class Solve *solve;        // KMC solver
    double *propensity;        // propensities for set sites
    int *site2i;               // map from set sites to lattice index
    int *i2site;               // map from lattice index to set sites
  };
  Set *set;                    // list of subsets
  int nset;                    // # of subsets of lattice sites

  struct Site {
    int id,proc,index;
    double x,y,z;
  };

  void iterate_kmc_global(double);
  void update_ifng_field(double);
  void iterate_kmc_sector(double);
  virtual void iterate_rejection(double);
  void iterate_app_update_only(double,double);

  typedef void (AppLattice::*FnPtrSweep)(int, int *);
  FnPtrSweep sweep;                         // ptr to< sweep functions
  void sweep_nomask_nostrict(int, int *);
  void sweep_mask_nostrict(int, int *);
  void sweep_nomask_strict(int, int *);
  void sweep_mask_strict(int, int *);

  void ghosts_from_connectivity();
  void connectivity_within_cutoff();

  void create_set(int, int, int, class Solve *);
  class Solve *free_set(int);
  int id2color(int);
  int find_border_sites(int);
  void boundary_clear_mask(int);

  void stats(char *);
  void stats_header(char *);

  void set_sector(int, char **);
  void set_sweep(int, char **);
  void set_temperature(int, char **);
  void set_app_update_only(int, char **);
  void bounds(char *, int, int, int &, int &);
};

}

#endif

/* ERROR/WARNING messages:

E: App needs a KMC or rejection KMC solver

You must define either a solver or sweep option.

E: App cannot use both a KMC and rejection KMC solver

You cannot define both a solver and sweep option.

E: KMC events are not implemented in app

Not every application supports KMC solvers.

E: Rejection events are not implemented in app

Self-explanatory.

E: Mask logic not implemented in app

Not every application supports masking.

E: Cannot use KMC solver in parallel with no sectors

Self-explanatory.

E: Cannot use random rejection KMC in parallel with no sectors

Self-explanatory.

E: Cannot use raster rejection KMC in parallel with no sectors

Self-explanatory.

E: Cannot use color/strict rejection KMC with sectors

Self-explanatory.

E: App did not set dt_sweep

Internal SPPARKS error.

E: Invalid number of sectors

Self-explanatory.

E: Cannot color without a lattice definition of sites

UNDOCUMENTED

E: Cannot color without contiguous site IDs

UNDOCUMENTED

E: Cannot color this combination of lattice and app

Coloring is not supported on this lattice for the neighbor
dependencies of this application.

E: Choice of sector stop led to no rKMC events

Self-explanatory.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: App does not permit app_update_only yes

This app does not have it's own update method

E: Per-processor system is too big

UNDOCUMENTED

*/
