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

#ifdef APP_CLASS
AppStyle(rxn/diff,AppRxnDiff)

#else

#ifndef SPK_APP_RXN_DIFF_H
#define SPK_APP_RXN_DIFF_H

#include "app_lattice.h"

namespace SPPARKS_NS {

#define MAX_SPECIES 18 // was 10, need to make 12 for all our species

class AppRxnDiff : public AppLattice {
//  friend class DiagRxnDiff;

 public:
  AppRxnDiff(class SPPARKS *, int, char **);
  ~AppRxnDiff();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);
  
  /////////////////////////////////// Giuseppe edit ///////////////////////

  int tcell_total_population;
  double cancer_potential_converter;
  std::array<std::array<int, 12>, 12> discretized_cancer; // CHANGE THIS WHEN CHANGING SIZE OF FIELD!!!
  std::array<std::array<std::array<double, 2>, 200>, 200> cancer_bias_potential;

  ////////////////////////////////////////////////////////////////////////

 protected:
  int engstyle;
 // int **population;      // each lattice site has an array of populations
  int firsttime;

  int *esites;
  int *echeck;


  int nspecies;                  // # of unique species
  char **sname;                  // ID of each species

  int nreactions;                // # of user defined reactions
  char **rname;                  // ID of each reaction
  int **localReactants;          // local reactants;
                                 // localReactants[I][J] = local number of reactants for reaction I of species J
  int **nbrReactants;            // neighbor reactants
  int **localDeltaPop;           // local population change if even is carried out;
                                 // localDeltaPop[I][J] = change in local population of species J if reaction I is carried out
  int **nbrDeltaPop;             // neighbor population change if even is carried out
  double *rate;                  // rate[I] = input rate for reaction I
  int *rxnStyle;                 // rxnStyle[I] = LOCAL or NBR for reaction I

  int *rxn_count;               // statistics

  ////////////////////////////////////////////// GIUSEPPE GIULIANI EDITS /////////////////////////////////////////////////////
  int max_cell_area_per_voxel;
  int lattice_width;
  int lattice_height;
  int cancer_area_count;
  int Tp_area_count;
  int Tm_area_count;
  int Te_area_count;
  int max_T_cells; // max number of CD8+ T cells
  int cpMax; // lattice x size for cancer potential for motility bias
  int cpMay; // lattice y size for cancer potential for motility bias
  bool update_all_always; // whether or not to update all lattice sites each Gillespie move
  //int rxn_mem; // for remembering previous interactions

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct Event {           // one event for an owned site
    int style;             // reaction style = LOCAL, NBR
    int which;             // which reaction
    int jpartner;          // which J neighbor of I are part of event
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);

  void add_event(int, int, int, double, int);
  void add_rxn(int, char **);
  void add_species(int, char **);

  int find_reaction(char *);
  int find_species(char *);

  virtual double custom_multiplier(int, int, int, int);

  void grow_reactions(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:


E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.


*/
