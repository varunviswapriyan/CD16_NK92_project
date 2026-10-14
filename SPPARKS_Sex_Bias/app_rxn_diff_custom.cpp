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

#include <iostream>

#include "math.h"
#include "string.h"
#include "app_rxn_diff_custom.h"
#include "random_park.h"
#include "error.h"
#include <random>
#include <iostream>

using namespace SPPARKS_NS;

enum{LOCAL,NBR};

/* ---------------------------------------------------------------------- */
AppRxnDiffCustom::AppRxnDiffCustom(SPPARKS *spk, int narg, char **arg) : 
  AppRxnDiff(spk,narg,arg)
{
  // parse arguments for RxnDiffCustom class only, not children
  if (strcmp(style,"rxn/diff/custom") != 0) return;
  if (narg != 37) error->all(FLERR,"Illegal app_style command"); // gdag, (narg != 37) diffusion half max
  lattice_width = atoi(arg[1]);
  lattice_height = atoi(arg[2]);
  max_cell_area_per_voxel = atoi(arg[3]);
  cancer_area_count = atoi(arg[4]);
  Tp_area_count = atoi(arg[5]);
  Tm_area_count = atoi(arg[6]);
  Te_area_count = atoi(arg[7]);
  do_carrying_capacity = atoi(arg[8]);
  tcell_carrying_capacity = atoi(arg[9]);
  do_ifng_cancer_proliferation_limit = atoi(arg[15]) == 1;
  do_ifng_cytotoxicity_augmentation = atoi(arg[16]) == 1;
  ifng_half_max_proliferation = atof(arg[17]);
  ifng_half_max_cytotoxic_augmentation_threshold = atof(arg[18]);
  tcell_stick_fac = atof(arg[20]);
  modulate_ifng_factor = atof(arg[24]);
  meanIFNG = atoi(arg[25]) == 1;
  hill_degree = atoi(arg[26]);
    cancer_proliferation_minimum_factor = atof(arg[33]);
    ifng_half_max_diffusion = atof(arg[36]);
  /*std::cout<< max_Tcell_per_voxel << std::endl;
  std::cout<< max_cancer_per_voxel << std::endl;
  std::cout<< tcell_stick_fac << std::endl;
  std::cout<< mac_activation_inhib_fac << std::endl;
  std::cout<< pd1_fac << std::endl;*/
}

/* ----------------------------------------------------------------------
 * Custom rate multiplier function
------------------------------------------------------------------------- */
double AppRxnDiffCustom::custom_multiplier(int i, int rstyle, int which, int jpartner)
{
	int latCol = lattice_width; //
	int latRow = lattice_height;
//    hill_degree = 4;

    //data IDs to access the populations within the population array
    int canID = 0;
    int TpID = 1;
    int TmID = 2;
    int TeID = 3;
    int ARTpID = 4;
    int ARTmID = 5;
    int ARTeID = 6;
    int ThID = 12; // helper CD4+ T cells
    int ARThID = 13; // AR+ helper CD4+ T cells
    int TrID = 14; // regulatory T cells
    int ARTrID = 15; // AR+ regulatory T cells
    int first_motility_interaction = 13;
    int first_cd4_motility_interaction = 39;
    int first_T_prolif_interaction = 1;
    int first_lysis_interaction = 27;
    int ifngID = 0;
    int boundifngID = 1;
    int recID = 7;
    int latEID = 10;
    int deadID = 9;
    int errID = 8;

    //area units of different cell types
    int maxA = max_cell_area_per_voxel;
    int cancA = cancer_area_count;
    int tpA = Tp_area_count;
    int tmA = Tm_area_count;
    int teA = Te_area_count;
    int thA = Tp_area_count;
    int trA = thA;


    if(!mapped_population_to_ifng_tracker){ // initialize ifng tracker if the first pass of ifng_tracker
        int cellIDs[11] = {canID, TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID, ARThID, TrID, ARTrID};

        int receptors_per_cell_id[11] = {2, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1};
        for (int x = 0; x < latCol; x++) {
            bound_ifng[x] = new std::array<double, 3>[latCol];
            ifng_tracker[x] = new std::vector<std::vector <double>>*[latCol];
            for (int y = 0; y < latCol; y++) {
                bound_ifng[x][y][0] = 0;
                bound_ifng[x][y][1] = 0;
                bound_ifng[x][y][2] = 0;
                int location = x+latCol*y;
                ifng_field[int(x/ifng_lattice_cell_converter)][int(y/ifng_lattice_cell_converter)] += ifng_field_discretized[ifngID][location];
                ifng_tracker[x][y] = new std::vector<std::vector <double>>[11]; // for 11 species to track
                for(int id=0; id<11; id++){
                    // ifng_tracker[x][y][id] = new std::vector<std::vector <double>>;
                    for(int p=0; p < population[cellIDs[id]][location]; p++){
                        double number_bound = ifng_field_discretized[boundifngID][location]*receptors_per_cell[receptors_per_cell_id[id]];
                        ifng_tracker[x][y][id].push_back(std::vector <double>{number_bound, (double) receptors_per_cell[receptors_per_cell_id[id]]});
                        bound_ifng[x][y][2] += receptors_per_cell[receptors_per_cell_id[id]];
                    }
                }
            }
        }
        mapped_population_to_ifng_tracker = true;
    }


    if(which == 0){ // need to return the average modifier per cells in pixel
        double proliferation_factor = 1.0;
        bool done = false;
        int neighLim = 3; // number of layers of neighbors that cancer can jump to
        int widthLim = 3 + 2*(neighLim-1);
        int width = 1;
        int ix = i%latCol;
        int iy = i/latCol;
        while(!done){
            int halfw = int((width-1) / 2);
            for(int ly=0; ly<2; ly++){
                for(int lx=0; lx<width; lx++){
                    int target1x = ix + lx - halfw;
                    int target1y = iy + (2*ly - 1)*halfw;
                    int target2x = ix + (2*ly-1)*halfw;
                    int target2y = iy + lx - halfw;
                    target1x -= floor(double(target1x) / latCol) * latCol;
                    target1y -= floor(double(target1y) / latCol) * latCol;
                    target2x -= floor(double(target2x)/latCol)*latCol;
                    target2y -= floor(double(target2y)/latCol)*latCol;
                    int target1 = target1x + latCol * target1y;
                    int target2 = target2x + latCol*target2y;
                    if(population[canID][target1]*cancA < maxA){
                        if(do_ifng_cancer_proliferation_limit) {
//                            double ifng_value = total_ifng/pow(ifng_col, 2);
//                            double ifng_reading = bound_ifng[x][y]/(double)receptors[x][y];
//                            if(!meanIFNG) {
//                                ifng_value = ifng_field[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength
//                            }
//                            double bound_fraction = bound_ifng[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]/(double)receptors[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength

                            if(!use_simple_decay_model) {
                                int local_cancer_count = population[canID][i];
                                double average_modifier = 0;
                                for(int cell=0; cell < population[canID][i]; cell++){
                                    double bound_fraction = ifng_tracker[ix][iy][canID][cell][0]/ifng_tracker[ix][iy][canID][cell][1]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                                    double hill_function = pow(bound_fraction, hill_degree) / (pow(ifng_half_max_proliferation, hill_degree)+pow(bound_fraction, hill_degree));
                                    double proliferation_limit_factor = cancer_proliferation_minimum_factor + (1-cancer_proliferation_minimum_factor)*(1 - hill_function); //(1 - (ifng_value / ifng_max_proliferation_threshold));
                                    average_modifier += proliferation_limit_factor/local_cancer_count;
                                }

                                //                            proliferation_factor = 0.5 * (1 + proliferation_limit_factor); // makes it so that the most the proliferation rate can dive is to half
                                proliferation_factor = average_modifier; // makes it so that the most the proliferation rate can dive is to half
                            }else {
                                if(!use_cell_dependent_decay_model){
                                    double local_ifng = ifng_field[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                                    double hill_function = pow(local_ifng, hill_degree) / (pow(ifng_half_max_proliferation, hill_degree)+pow(local_ifng, hill_degree));
                                    double proliferation_limit_factor = cancer_proliferation_minimum_factor + (1-cancer_proliferation_minimum_factor)*(1 - hill_function); //(1 - (ifng_value / ifng_max_proliferation_threshold));
                                    proliferation_factor = proliferation_limit_factor;
                                }else { // using the same mode of local ifng sensing right now as the cell should always sense true concentration
                                    double local_ifng = ifng_field[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                                    double hill_function = pow(local_ifng, hill_degree) / (pow(ifng_half_max_proliferation, hill_degree)+pow(local_ifng, hill_degree));
                                    double proliferation_limit_factor = cancer_proliferation_minimum_factor + (1-cancer_proliferation_minimum_factor)*(1 - hill_function); //(1 - (ifng_value / ifng_max_proliferation_threshold));
                                    proliferation_factor = proliferation_limit_factor;
                                }
                            }
                        }
                        return proliferation_factor; // should be 1.0
                        done = true;
                    }
                    if(population[canID][target2]*cancA < maxA){
                        if(do_ifng_cancer_proliferation_limit) {
                            //   double bound_fraction = bound_ifng[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]/(double)receptors[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength

                            //    ifng_value *= modulate_ifng_factor;
                            if(!use_simple_decay_model) {
                                int local_cancer_count = population[canID][i];
                                double average_modifier = 0;
                                for(int cell=0; cell < population[canID][i]; cell++){
                                    double bound_fraction = ifng_tracker[ix][iy][canID][cell][0]/ifng_tracker[ix][iy][canID][cell][1]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                                    double hill_function = pow(bound_fraction, hill_degree) / (pow(ifng_half_max_proliferation, hill_degree)+pow(bound_fraction, hill_degree));
                                    double proliferation_limit_factor = cancer_proliferation_minimum_factor + (1-cancer_proliferation_minimum_factor)*(1 - hill_function); //(1 - (ifng_value / ifng_max_proliferation_threshold));
                                    average_modifier += proliferation_limit_factor/local_cancer_count;
                                }

                                //                            proliferation_factor = 0.5 * (1 + proliferation_limit_factor); // makes it so that the most the proliferation rate can dive is to half
                                proliferation_factor = average_modifier; // makes it so that the most the proliferation rate can dive is to half
                            }else {
                                if(!use_cell_dependent_decay_model){
                                    double local_ifng = ifng_field[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                                    double hill_function = pow(local_ifng, hill_degree) / (pow(ifng_half_max_proliferation, hill_degree)+pow(local_ifng, hill_degree));
                                    double proliferation_limit_factor = cancer_proliferation_minimum_factor + (1-cancer_proliferation_minimum_factor)*(1 - hill_function); //(1 - (ifng_value / ifng_max_proliferation_threshold));
                                    proliferation_factor = proliferation_limit_factor;
                                }else { // using the same mode of local ifng sensing right now as the cell should always sense true concentration
                                    double local_ifng = ifng_field[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                                    double hill_function = pow(local_ifng, hill_degree) / (pow(ifng_half_max_proliferation, hill_degree)+pow(local_ifng, hill_degree));
                                    double proliferation_limit_factor = cancer_proliferation_minimum_factor + (1-cancer_proliferation_minimum_factor)*(1 - hill_function); //(1 - (ifng_value / ifng_max_proliferation_threshold));
                                    proliferation_factor = proliferation_limit_factor;
                                }
                            }
                        }
                        return proliferation_factor;
                        done = true;
                    }
                }
            }
            if(width >= widthLim){
                return 0.0;
                done =  true;
            }else{width += 2;}
        }
    }else if(7 <= which && which <= 12){ // T cell recruitment
        int t_carry_cap = tcell_carrying_capacity;
        int t_tot_pop = tcell_total_population;
        double deadTrec_fac = (1-(double(t_tot_pop)*do_carrying_capacity/double(t_carry_cap)));
        if(deadTrec_fac < 0){
            deadTrec_fac = 0;
        }
        return deadTrec_fac;
    }else if(first_T_prolif_interaction <= which && which <=6){ // T cell proliferation
        int filling = tpA*population[TpID][i] + tmA*population[TmID][i] + teA*population[TeID][i] + tpA*population[ARTpID][i] + tmA*population[ARTmID][i] + teA*population[ARTeID][i] + thA*population[ThID][i]+ thA*population[ARThID][i] + trA*population[TrID][i] + trA*population[ARTrID][i] + cancA*population[canID][i];
        int weights[] = {tpA, tmA, teA, tpA, tmA, teA}; // cells which may proliferate
        int weight = weights[which - first_T_prolif_interaction];
        if(filling + weight <= maxA){
          int t_carry_cap = tcell_carrying_capacity;
          int t_tot_pop = tcell_total_population;
          double deadTrec_fac = 1-(double(t_tot_pop)*do_carrying_capacity/double(t_carry_cap));
            if(deadTrec_fac < 0){
                deadTrec_fac = 0;
            }
          return deadTrec_fac;
        }else{
          return 0;
        }

    }else if ((first_motility_interaction <= which && which <= 18) || (first_cd4_motility_interaction <= which && which <= 42)) {// T cell motility
//        int filling = tpA*population[TpID][i] + tmA*population[TmID][i] + teA*population[TeID][i] + tpA*population[ARTpID][i] + tmA*population[ARTmID][i] + teA*population[ARTeID][i] + cancA*population[canID][i];
        int nfilling =
                tpA * population[TpID][jpartner] + tmA * population[TmID][jpartner] + teA * population[TeID][jpartner] +
                tpA * population[ARTpID][jpartner] + tmA * population[ARTmID][jpartner] +
                teA * population[ARTeID][jpartner] + thA*population[ThID][jpartner]+ thA*population[ARThID][jpartner] +
                    trA*population[TrID][jpartner] + trA*population[ARTrID][jpartner] + cancA * population[canID][jpartner];

        // find direction of motility
        int axis = 0; // left/right
        int direction = 0;
        if (abs(i - jpartner) > 1) {
            axis = 1; // up/down
        }
        if(i < jpartner){ // is it forward?
            direction = 1;
        }

        int t_carry_cap = tcell_carrying_capacity;
        int t_tot_pop = tcell_total_population;

        int weights[] = {tpA, tmA, teA, tpA, tmA, teA, tpA, tpA, tpA, tpA}; // last 4 are Th, AR+Th, Treg, AR+ Treg
        int weight = 0;
        double adhesion_factor = 1;
        if(which - first_motility_interaction < 6) { // CD8s
            weight = weights[which - first_motility_interaction];
        }else { // CD4s
            weight = weights[which - first_cd4_motility_interaction + 6];
        }

        double bias_factor = 2*cancer_bias_potential[int((i%latCol)/cancer_potential_converter)][int(int(i/latCol)/cancer_potential_converter)][axis]; // scaled 0-1, 0 being complete bias to left, 1 bias to right, 0.5 is no bias, so 2* makes it diffuse normally when potential is 1/2 (no bias)
        if(direction == 0){
            bias_factor = 2 - bias_factor; // bias factor is 0-2, 1 gives normal diffusion. >1 gives bias to direction
        }

        if(population[latEID][i] >= 1 && population[latEID][jpartner] >= 1){ // check if hopping over boundary and block
            bool cross = false;
            int cornCount = 0;
            int cornCount2 = 0;
            int type = -1;
            int type2 = -1;
            if( i/latCol == 0){ // top
                cornCount++;
                type = 0;
            }else if(i/latCol == latCol - 1){ // bottom
                cornCount++;
                type = 2;
            }
            if(i%latCol == 0){ // left
                cornCount++;
                type = 1;
            }else if(i%latCol == latCol - 1){ // right
                cornCount++;
                type = 3;
            }

            if(jpartner/latCol == 0){ // top
                cornCount2++;
                type2 = 0;
            }else if(jpartner/latCol == latCol - 1){ // bottom
                cornCount2++;
                type2 = 2;
            }
            if(jpartner%latCol == 0){ // left
                cornCount2++;
                type2 = 1;
            }else if(jpartner%latCol == latCol - 1){ // right
                cornCount2++;
                type2 = 3;
            }

            if((cornCount == 2 && cornCount2 == 2) || (cornCount == 1 && cornCount2 == 1 && type != type2)){
                cross = true;
            }

            if(cross){
                return 0;
            }
        }

        double local_ifng = ifng_field[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength
        double hill_function = pow(local_ifng, hill_degree) / (pow(ifng_half_max_proliferation, hill_degree)+pow(local_ifng, hill_degree));
        double proliferation_limit_factor = cancer_proliferation_minimum_factor + (1-cancer_proliferation_minimum_factor)*(1 - hill_function); //(1 - (ifng_value / ifng_max_proliferation_threshold));

        if(maxA < nfilling + weight){
            return 0.0;
        }else{
            ///////////////////////////////////////////////////////////////////////////////////////////// diffusion IFNG dependence here
            double ifng_diffusion_factor = 1;
            int cd8_index = 0;
            if(which - first_motility_interaction < 6) { // CD8s
                cd8_index = int(which - first_motility_interaction);
                adhesion_factor = tcell_stick_fac; // make it so CD4s don't stick
            }else { // CD4s
                cd8_index = int(which - first_cd4_motility_interaction + 6);
            }

            if(!use_simple_decay_model) {
                int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID, ARThID, TrID, ARTrID};
                int tctracker_map[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
                int cellID = cd8map[cd8_index]; // cell ID of cell type
                int trackID = tctracker_map[cd8_index];

                int ix = i%latCol;
                int iy = i/latCol;

                int local_tc_count = population[cellID][i];

                int denom = local_tc_count; //*target_can_count;
                double average_modifier = 0;
                for(int cell=0; cell < population[cellID][i]; cell++){
                    double bound_fraction = ifng_tracker[ix][iy][trackID][cell][0]/ifng_tracker[ix][iy][trackID][cell][1]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                    double hill_function = pow(bound_fraction, hill_degree) / (pow(ifng_half_max_diffusion, hill_degree)+pow(bound_fraction, hill_degree));
                    double ifng_diffusion_factor = adhesion_factor+(1-adhesion_factor)*hill_function;
                    average_modifier += ifng_diffusion_factor/denom;
                }
                ifng_diffusion_factor = average_modifier;
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            int ix = i%latCol;
            int iy = i/latCol;
            if(population[canID][i] > 0){
                return bias_factor*adhesion_factor*ifng_diffusion_factor;
            }
            for(int k=0;k < 2; k++){  // check if any cancer cells surrounding
                int targx = ix;
                int targy = iy - 1 + 2*k;
                targy -= floor(double(targy)/latCol)*latCol;
                int targ = targx + targy*latCol;
                if(population[canID][targ] > 0){
                    return bias_factor*adhesion_factor*ifng_diffusion_factor;
                }
            }
            for(int h=0; h < 2; h++){
                 int targx = ix - 1 + 2*h;
                 int targy = iy;
                 targx -= floor(double(targx)/latCol)*latCol;
                 int targ = targx + targy*latCol;
                 if(population[canID][targ] > 0){
                     return bias_factor*adhesion_factor*ifng_diffusion_factor;
                 }
            }
           return bias_factor;
        }
    /*}else if(false){ // (NO MORE Tn) Tn -> Tc
      int macFactor = 1; // 0.5 without anti-pd1, 1
      int max_filling = 4;
      int ix = i%latCol;
      int iy = i/latCol;
      bool canPresent = false;
      bool macPresent = false;
      for(int k=0; k<3; k++){
        for(int h=0; h<3; h++){
          int targx = ix + (h-1);
          int targy = iy + (k-1);
          targx -= floor(double(targx)/latCol)*latCol;
          targy -= floor(double(targy)/latCol)*latCol;
          int targ = targx + targy*latCol;
          if(population[1][targ] > 0){
            canPresent = true;
          }else if(population[7][targ] > 0){
            macPresent = true;
          }
        }
      }
      if(macPresent && canPresent){
        double mac_fac = 1; //mac_activation_inhib_fac/pd1_fac;
        if(mac_fac > 1){
          mac_fac = 1;
        }
        return mac_fac;
      }else if(canPresent){
        return 1.0;
      }else{
        return 0.0;
      }*/
    }else if(first_lysis_interaction <= which && which <= 34){ // update T cell killing
        bool local = which % 2 == 1;
        if(population[latEID][i] >= 1 && population[latEID][jpartner] >= 1 && !local){ // check if lysing over boundary and block
            bool cross = false;
            int cornCount = 0;
            int cornCount2 = 0;
            int type = -1;
            int type2 = -1;
            if( i/latCol == 0){ // top
                cornCount++;
                type = 0;
            }else if(i/latCol == latCol - 1){ // bottom
                cornCount++;
                type = 2;
            }
            if(i%latCol == 0){ // left
                cornCount++;
                type = 1;
            }else if(i%latCol == latCol - 1){ // right
                cornCount++;
                type = 3;
            }

            if(jpartner/latCol == 0){ // top
                cornCount2++;
                type2 = 0;
            }else if(jpartner/latCol == latCol - 1){ // bottom
                cornCount2++;
                type2 = 2;
            }
            if(jpartner%latCol == 0){ // left
                cornCount2++;
                type2 = 1;
            }else if(jpartner%latCol == latCol - 1){ // right
                cornCount2++;
                type2 = 3;
            }

            if((cornCount == 2 && cornCount2 == 2) || (cornCount == 1 && cornCount2 == 1 && type != type2)){
                cross = true;
            }

            if(cross){
                return 0;
            }
        }

        int cd8_index = int((which - first_lysis_interaction)/2);
        int cd8map[] = {TmID, TeID,ARTmID, ARTeID};
        int cellID = cd8map[cd8_index]; // cell ID of cell type

        int ix = i%latCol;
        int iy = i/latCol;

        if(do_ifng_cytotoxicity_augmentation) {
            int local_tc_count = population[cellID][i];
//            int target_can_count = 0;
//            if(local){
//                target_can_count = population[canID][i];
//            }else{
//                target_can_count = population[canID][jpartner];
//            }
            double modifier = 1;
            if(!use_simple_decay_model) {
                int denom = local_tc_count; //*target_can_count;
                double average_modifier = 0;
                for(int cell=0; cell < population[cellID][i]; cell++){
                    double bound_fraction = ifng_tracker[ix][iy][cellID][cell][0]/ifng_tracker[ix][iy][cellID][cell][1]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                    double hill_function = pow(bound_fraction, hill_degree) / (pow(ifng_half_max_cytotoxic_augmentation_threshold, hill_degree)+pow(bound_fraction, hill_degree));
                    double cyto_augmentation_fac = hill_function; //(1 - (ifng_value / ifng_max_proliferation_threshold));
                    average_modifier += cyto_augmentation_fac/denom;
                }
                modifier = average_modifier;
            }else { // if using simple model
                double local_ifng = ifng_field[int((i % latCol) / ifng_lattice_cell_converter)][int((int(i / latCol)) / ifng_lattice_cell_converter)]; // total_ifng/pow(ifng_col, 2);  // find local field strength
                double hill_function = pow(local_ifng, hill_degree) / (pow(ifng_half_max_cytotoxic_augmentation_threshold, hill_degree)+pow(local_ifng, hill_degree));
                modifier = hill_function;
            }

            return modifier;
        }else{
            return 0.5;
        }
    }else if((21 <= which && which <= 24) || (43 <= which && which <= 58)){ // update exhaustion of T cells as well as activation of progenitor cells, blocked if across a barrier
        if(population[latEID][i] >= 1 && population[latEID][jpartner] >= 1){ // check if lysing over boundary and block
            bool cross = false;
            int cornCount = 0;
            int cornCount2 = 0;
            int type = -1;
            int type2 = -1;
            if( i/latCol == 0){ // top
                cornCount++;
                type = 0;
            }else if(i/latCol == latCol - 1){ // bottom
                cornCount++;
                type = 2;
            }
            if(i%latCol == 0){ // left
                cornCount++;
                type = 1;
            }else if(i%latCol == latCol - 1){ // right
                cornCount++;
                type = 3;
            }

            if(jpartner/latCol == 0){ // top
                cornCount2++;
                type2 = 0;
            }else if(jpartner/latCol == latCol - 1){ // bottom
                cornCount2++;
                type2 = 2;
            }
            if(jpartner%latCol == 0){ // left
                cornCount2++;
                type2 = 1;
            }else if(jpartner%latCol == latCol - 1){ // right
                cornCount2++;
                type2 = 3;
            }

            if((cornCount == 2 && cornCount2 == 2) || (cornCount == 1 && cornCount2 == 1 && type != type2)){
                cross = true;
            }

            if(cross){
                return 0;
            }
        }
    }else if(false){ // Cells floating in from exterior
      return 0.0; /////////////////// CHANGE
      int filling = tpA*population[TpID][i] + tmA*population[TmID][i] + teA*population[TeID][i] + tpA*population[ARTpID][i] + tmA*population[ARTmID][i] + teA*population[ARTeID][i] + cancA*population[canID][i];

      if(filling >= maxA){
        return 0.0;
      }else{
        return 1.0;
      }
    }
//    else if(false){ // Macrophage motility
//      int filling = tcA*population[cd8aID][i] + teA*population[cd8eID][i] + cancA*(population[canID][i] + macA*population[macID][i]);
//      if(false && which == 13){ // proliferation
//        if(maxA - filling < macA){
//          return 0.0;
//        }else{return 1.0;}
//      }else if(which == 10){
//        int nfilling = tcA*population[cd8aID][jpartner] + teA*population[cd8eID][jpartner] + cancA*(population[canID][jpartner] + macA*population[macID][jpartner]);
//
//        if(maxA - nfilling < macA){
//          return 0.0;
//        }else{return 1.0;}
//      }
//    }
//    if(which == 7 || which == 11 || which == 13 || which == 15){ // apply pd1 ratio to these rates so we don't have to in the input file directly
//      //return pd1_fac;
//    }
//    if(which == 3 || which == 14){
//      //return 1/pd1_fac;
//    }
    return 1.0;
}

//void AppRxnDiffCustom::stats(char *strtmp)
//{
  //char big[8],format[64];
  //strcpy(big,BIGINT_FORMAT);

  /*int rxn_count_b2B, rxn_count_p2P;
  MPI_Allreduce(&rxn_count[17],&rxn_count_b2B,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&rxn_count[18],&rxn_count_p2P,1,MPI_INT,MPI_SUM,world);
  food_left = num_nutrient - rxn_count_b2B - rxn_count_p2P;

  bigint naccept_all;
  MPI_Allreduce(&naccept,&naccept_all,1,MPI_SPK_BIGINT,MPI_SUM,world);

  sprintf(format,"%%10g %%10%s %%10d %%10d",&big[1]);
  sprintf(strtmp,format,time,naccept_all,nsweeps,food_left);

}*/

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

/*void AppRxnDiffCustom::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s %10s %10s %10s","Time","Naccept","Nsweeps","Nfood");
}*/
/* ----------------------------------------------------------------------
  	}
  }else if (which==2 || which == 14){
    int il2_conc = population[6][i];

    if ( bact_fill_ahead >= max_Tcell_per_voxel || tcell_fill_ahead >= max_cancer_per_voxel) {
      return 0.0;*/
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

#include <iostream>

#include "math.h"
#include "string.h"
#include "app_rxn_diff_custom.h"
