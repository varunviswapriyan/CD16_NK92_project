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
AppStyle(rxn/diff/custom,AppRxnDiffCustom)

#else

#ifndef SPK_APP_RXN_DIFF_CUSTOM_H
#define SPK_APP_RXN_DIFF_CUSTOM_H

#include "app_rxn_diff.h"

namespace SPPARKS_NS {

class AppRxnDiffCustom : public AppRxnDiff {
 public:
  AppRxnDiffCustom(class SPPARKS *, int, char **);
  ~AppRxnDiffCustom() {}

 protected:
  // Giuseppe defined? Are there any that aren't?
  int max_cell_area_per_voxel;
  int cancer_area_count;
  int Tp_area_count;
  int Tm_area_count;
  int Te_area_count;
  int do_carrying_capacity;
  int tcell_carrying_capacity;
  int lattice_width;
  int lattice_height;
  bool do_ifng_cancer_proliferation_limit;
  bool do_ifng_cytotoxicity_augmentation;
  double ifng_half_max_proliferation;
  double ifng_half_max_cytotoxic_augmentation_threshold;
    double ifng_half_max_diffusion;
  double tcell_stick_fac;
  bool meanIFNG;
  double modulate_ifng_factor;
  int hill_degree;
    double cancer_proliferation_minimum_factor;


  double custom_multiplier(int, int, int, int);

};

}

#endif
#endif
