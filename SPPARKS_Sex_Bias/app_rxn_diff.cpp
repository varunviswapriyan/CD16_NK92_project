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
#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_rxn_diff.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
// Giuseppe import:
#include <random>
#include <iostream>
#include <algorithm>

using namespace SPPARKS_NS;

enum{LOCAL,NBR};

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppRxnDiff::AppRxnDiff(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = MAX_SPECIES;
  ndouble = 2; // gdag edit - this was 0, but I need an ifng field and a bound ifng field
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  create_arrays();

  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  nevents = 0;
  maxevent = 0;
  firstevent = NULL;

  // species lists
  nspecies = 0;
  sname = NULL;

  // reaction lists
  nreactions = 0;
  rname = NULL;
  localReactants = NULL;
  nbrReactants = NULL;
  localDeltaPop = NULL;
  nbrDeltaPop = NULL;
  rate = NULL;
  rxnStyle = NULL;

  rxn_count = NULL; // count of how many times each reaction takes place

  /////////////////////////////////////// GIUSEPPE GIULIANI EDITS ///////////////////////////////////
  //max_Tcell_per_voxel = atoi(arg[1]);
  //max_cancer_per_voxel = atoi(arg[2]);
  lattice_width = atoi(arg[1]);
  lattice_height = atoi(arg[2]);
  max_cell_area_per_voxel = atoi(arg[3]);
  cancer_area_count = atoi(arg[4]);
  Tp_area_count = atoi(arg[5]);
  Tm_area_count = atoi(arg[6]);
  Te_area_count = atof(arg[7]);
  max_T_cells = atof(arg[9]);
  update_all_always = atoi(arg[22]) == 1;


  int latCol = lattice_width;

  int canID = 0;
  int TpID = 1;
  int TmID = 2;
  int TeID = 3;
  int ARTpID = 4;
  int ARTmID = 5;
  int ARTeID = 6;
  int ifngID = 0;
  int recID = 7;
  int latEID = 10;
  int deadID = 9;
  int errID = 8;
  
  tcell_total_population = 0;
  cpMax = 12;
  cpMay = 12;
  cancer_potential_converter = 50;
  /*
  for(int ic=0; ic < latCol*latCol; ic++){
    tcell_total_population += population[cd8aID][ic] + population[cd8eID][ic];
  }*/

  // parse arguments for RxnDiff class only, not children
  if (strcmp(style,"rxn/diff") != 0) return;
  if (narg != 1) error->all(FLERR,"Illegal app_style command");
}

/* ---------------------------------------------------------------------- */

AppRxnDiff::~AppRxnDiff()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->destroy(firstevent);

  memory->destroy(sname);
  memory->destroy(rname);
  memory->destroy(localReactants);
  memory->destroy(nbrReactants);
  memory->destroy(localDeltaPop);
  memory->destroy(nbrDeltaPop);
  memory->destroy(rate);
  memory->destroy(rxnStyle);
  memory->destroy(rxn_count);
}

/* ---------------------------------------------------------------------- */

void AppRxnDiff::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"add_rxn") == 0) add_rxn(narg,arg);
  else if (strcmp(command,"add_species") == 0) add_species(narg,arg);
  else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppRxnDiff::grow_app() // gdag: called by app_lattice.cpp grow()
{
  population = iarray;
  ifng_field_discretized = darray; // double field, gdag added for ifng field
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppRxnDiff::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for 2 sites and their 1st neighbors, must be large enough for entire lattice for t cell/cancer cell updates!
    esites = new int[360000];//new int[118];//2 + 2*maxneigh]; gdag edit, for 2 cancer cell sites, their neighbors, their displaced t cells and their neighbors: 2: 25+16+4*9, 3: 49+33+4*9
  }

  // zero reaction counts
  delete [] rxn_count;
  memory->create(rxn_count,nreactions,"rxn/diff:rxn_count");
  for (int m = 0; m < nreactions; m++) rxn_count[m] = 0;

  // site validity
  int flag = 0;
  for (int ispecies = 0; ispecies < MAX_SPECIES; ispecies++) {
    for (int i = 0; i < nlocal; i++) {
      if (population[ispecies][i] < 0) flag = 1;
    }
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ---------------------------------------------------------------------- */

void AppRxnDiff::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppRxnDiff::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppRxnDiff::site_propensity(int i)
{
  int j,k,m;

  // valid single, double, triple events are in tabulated lists
  // propensity for each event is input by user

  clear_events(i);

  double proball = 0.0;
  double my_rate, multiplicity;
  int ispecies, flag;

  // populate event list for this site
  for (m = 0; m < nreactions; m++) {

    // Local reactions
    if (rxnStyle[m] == LOCAL){
      multiplicity = 1.0;
      flag = 0;

      // check that reactants are present
      for (ispecies=0; ispecies < MAX_SPECIES; ispecies++) {
        if (localReactants[m][ispecies]) {
          if (population[ispecies][i] < localReactants[m][ispecies]){
            flag = 1;
            break;
          }

          // calculate multiplicity:
          // population[ispecies][i] choose localReactants[m][ispecies]
          for (k = 1; k<=localReactants[m][ispecies]; k++)
            multiplicity *= (double)(population[ispecies][i]+1-k)/(double)k;
        }
      }
      if (flag) continue;

      my_rate = rate[m] * multiplicity * custom_multiplier(i,LOCAL,m,-1);
      if (my_rate == 0.0) continue;
      add_event(i, LOCAL, m, my_rate, -1);
      proball += my_rate;

    // Neighbor reactions
    } else { // rxnStyle[m] == NBR
      for (int jj = 0; jj < numneigh[i]; jj++) { // through e.g. 4 nearest neighbors
        j = neighbor[i][jj]; // location of neighbor jj?
        multiplicity = 1.0; 
        flag = 0;

        // check that reactants are present
        for (ispecies=0; ispecies < MAX_SPECIES; ispecies++) {
          // local
          if (localReactants[m][ispecies]) {
            if (population[ispecies][i] < localReactants[m][ispecies]){
              flag = 1;
              break;
            }

            // calculate multiplicity:
            // population[ispecies][i] choose localReactants[m][ispecies]
            for (k = 1; k<=localReactants[m][ispecies]; k++)
              multiplicity *= (double)(population[ispecies][i]+1-k)/(double)k;
          }

          // neighbor
          if (nbrReactants[m][ispecies]) {
            if (population[ispecies][j] < nbrReactants[m][ispecies]){
              flag = 1;
              break;
            }

            // calculate multiplicity:
            // population[ispecies][j] choose nbrReactants[m][ispecies]
            for (k = 1; k<=nbrReactants[m][ispecies]; k++)
              multiplicity *= (double)(population[ispecies][j]+1-k)/(double)k;
          }
        }
        if (flag) continue;

        my_rate = rate[m] * multiplicity * custom_multiplier(i,NBR,m,j);
        if (my_rate == 0.0) continue;
        add_event(i, NBR, m, my_rate, j);
        proball += my_rate;
      }
    }
  }
//////////////////////// gdag edit /////////////////////////
  /*std::vector<int> poss;
  bool problo = false;
  int ievent = firstevent[i];
  
  while(ievent >= 0){
  	int wh = events[ievent].which;
	for(int sp=0; sp < MAX_SPECIES; sp++){
		if(localReactants[wh][sp]){
             		if (population[sp][i] < localReactants[wh][sp]){
				problo = true;
				poss.push_back(ievent);
				break;
             		}
		}
		if(nbrReactants[wh][sp]){
			if (population[sp][j] < nbrReactants[wh][sp]){
				problo = true;
				poss.push_back(ievent);
				break;
			}
		}
	}
	ievent = events[ievent].next;
  }
  
  poss.clear();
  */
  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppRxnDiff::site_event(int i, class RandomPark *random)
{
  int j,m,n;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i]; // gdag- firstevent[i] refers to the first event in the chain of possible events for site i
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform the event

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  
  //bool pro = (which == 5 && population[0][i] == 0)||(which == 6 && population[5][i] == 0)||(which == 7 && population[2][i] == 0)||(which == 9 && population[0][i] == 0) || (which == 3 && population[5][i] == 0);
  bool prob = false;
  bool overage = false;
  int countermi = 0;
  for (int ispecies=0; ispecies < MAX_SPECIES; ispecies++){
    population[ispecies][i] += localDeltaPop[which][ispecies];
    if(population[ispecies][i] < 0){
      prob = true;
    }
    if(2*population[0][i] + population[1][i] + population[2][i] + population[3][i] + population[4][i] + population[5][i] + population[6][i] > 4){
        overage = true;
        countermi += 1;
    }
  }
  if (rstyle == NBR){
    for (int ispecies=0; ispecies < MAX_SPECIES; ispecies++){
      population[ispecies][j] += nbrDeltaPop[which][ispecies];
      if(population[ispecies][j] < 0){
        prob = true;
      }
      if(2*population[0][j] + population[1][j] + population[2][j] + population[3][j] + population[4][j] + population[5][j] + population[6][j] > 4){
          overage = true;
          countermi += 1;
      }
    }
  }

    //////////////////////////////////////////////////// GIUSEPPE GIULIANI EDITS ///////////////////////////////////////////
    int latCol = lattice_width; //
    int latRow = lattice_height;

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
    int ifngID = 0;
    int recID = 7;
    int latEID = 10;
    int deadID = 9;
    int errID = 8;
    int boundID = 11;

    bool update_ifng_distribs = false;

    if(prob){
        population[errID][0] = i;
    }

    //area units of different cell types
    int maxA = max_cell_area_per_voxel;
    int cancA = cancer_area_count;
    int tpA = Tp_area_count;
    int tmA = Tm_area_count;
    int teA = Te_area_count;
    int thA = Tp_area_count;
    int trA = thA;

    int cancerPos[] = {0,0};
    std::vector<int> tcellPos;
    int tcCount = 0;
    int neighLim = 3; // number of layers of neighbors that cancer can jump to
    int widthLim = 3 + 2*(neighLim-1);
    if(which == 0){
    //    int latCol = 100; ///////// SHOULD BE 100
        /*std::random_device rd;
        std::mt19937 mt(rd());
          std::uniform_real_distribution<double> dist(0,1.0);
        double r1 = dist(mt); // neeeeeeeeeeeeedd to put this back!!!!!*/
        bool tcellThere = false;
        double r1 = random->uniform();
        cancerPos[0] = i;
        std::vector<int> opos;
        std::vector<int> popos;
        int ix = i%latCol;
        int iy = i/latCol;
        int width = 1;
        bool done = false;

        int way = -1;
        while(!done){
            int halfw = int((width-1)/2);
            for(int ly=0; ly < 2; ly++){
                for(int lx=0; lx < width; lx++){
                    int target1x = ix + lx - halfw;
                    int target1y = iy + (2*ly - 1)*halfw;
                    int target2x = ix + (2*ly - 1)*halfw;
                    int target2y = iy + lx - halfw;
                    target1x -= floor(double(target1x)/latCol)*latCol;
                    target1y -= floor(double(target1y)/latCol)*latCol;
                    target2x -= floor(double(target2x)/latCol)*latCol;
                    target2y -= floor(double(target2y)/latCol)*latCol;
                    int target1 = target1x + latCol*target1y;
                    int target2 = target2x + latCol*target2y;
                    int tcellCount1 = population[TpID][target1] + population[TmID][target1] + population[TeID][target1] + population[ARTpID][target1] + population[ARTmID][target1] + population[ARTeID][target1] + population[ThID][target1]+ population[ARThID][target1] + population[TrID][target1] + population[ARTrID][target1];
                    int tcellCount2 = population[TpID][target2] + population[TmID][target2] + population[TeID][target2] + population[ARTpID][target2] + population[ARTmID][target2] + population[ARTeID][target2] + population[ThID][target2]+ population[ARThID][target2] + population[TrID][target2] + population[ARTrID][target2];
                    int filling1 = tpA*population[TpID][target1] + tmA*population[TmID][target1] + teA*population[TeID][target1] + tpA*population[ARTpID][target1] + tmA*population[ARTmID][target1] + teA*population[ARTeID][target1] + cancA*population[canID][target1] + thA*population[ThID][target1]+ thA*population[ARThID][target1] + trA*population[TrID][target1] + trA*population[ARTrID][target1];
                    int non_push_filling1 = cancA*population[canID][target1];
                    int filling2 = tpA*population[TpID][target2] + tmA*population[TmID][target2] + teA*population[TeID][target2] + tpA*population[ARTpID][target2] + tmA*population[ARTmID][target2] + teA*population[ARTeID][target2] + cancA*population[canID][target2] + thA*population[ThID][target2]+ thA*population[ARThID][target2] + trA*population[TrID][target2] + trA*population[ARTrID][target2];
                    int non_push_filling2 = cancA*population[canID][target2];
                    if(cancA*population[canID][target1] < maxA){  // is target 1 filled by cancer cells already?
                        if(maxA - filling1 >= cancA){ // is the room left in the site enough for a cancer cell?
                            opos.push_back(target1); // if so, add this site to the free open site list
                        }else if(maxA - non_push_filling1 >= cancA){  // if not, add this site to the list of possible sites where T cells need bumped
                            popos.push_back(target1);
                        }
                    }
                    if(cancA*population[canID][target2] < maxA){
                        if(maxA - filling2 >= cancA){
                            opos.push_back(target2);
                        }else if(maxA - non_push_filling2 >= cancA){
                            popos.push_back(target2);
                        }
                    }
                }
            }
            if(opos.size() > 0 || popos.size() > 0 || width >= widthLim){
                done =  true;
                if(width >= widthLim && !(opos.size() > 0 || popos.size() > 0)){
                    population[errID][7] = 1000;
                }
                std::sort(opos.begin(), opos.end());
                auto first = std::unique(opos.begin(), opos.end());
                opos.erase(first, opos.end());

                std::sort(popos.begin(), popos.end());
                auto sec = std::unique(popos.begin(), popos.end());
                popos.erase(sec, popos.end());
            }else{width += 2;}
        }
    //	tpopos = opos;
    //	tpopos.insert(tpopos.end(), popos.begin(), popos.end());
        //if(opos.size() > 0){
        if(opos.size() > 0) {
            int choice = opos[int(r1*opos.size())];
            //int choice = opos[int(r1*oCount)];
            population[canID][choice]++;
            cancerPos[1] = choice;

            int ilatx = i%latCol;
            int ilaty = int(i/latCol);
            int clatx = choice%latCol;
            int claty = int(choice/latCol);

            int local_ifngx = int(ilatx/ifng_lattice_cell_converter);
            int local_ifngy = int(ilaty/ifng_lattice_cell_converter);
            int chosen_ifngx = int(clatx/ifng_lattice_cell_converter);
            int chosen_ifngy = int(claty/ifng_lattice_cell_converter);

            if(!use_simple_decay_model) {
                // update the tracker and receptors which we do not need to do with simple model
                if(!use_leaky_receptor_model){
                    double random_value =random->uniform();

                    double gross_fraction = 0;
                    for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                        gross_fraction += ifng_tracker[ilatx][ilaty][canID][cc][0]/ifng_tracker[ilatx][ilaty][canID][cc][1];
                    }
                    double chosen_fraction = random_value*gross_fraction;
                    double bound_receptor_fraction_counter = 0;
                    int chosen_cell_index_count = 0;
                    for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                        chosen_cell_index_count = cc;
                        bound_receptor_fraction_counter += ifng_tracker[ilatx][ilaty][canID][cc][0]/ifng_tracker[ilatx][ilaty][canID][cc][1];
                        if(chosen_fraction <= bound_receptor_fraction_counter){
                            break;
                        }
                    }
                    int random_cell = chosen_cell_index_count;
                    double saved_bound = ifng_tracker[ilatx][ilaty][canID][random_cell][0];
                    std::vector<double> new_cancer_cell = {saved_bound/2, ifng_tracker[ilatx][ilaty][canID][random_cell][1]};
                    ifng_tracker[ilatx][ilaty][canID][random_cell][0] = saved_bound/2;
                    bool same_pixel = (local_ifngx == chosen_ifngx && local_ifngy == chosen_ifngy);
                    bool same_pixel_update = same_pixel && ifng_field[local_ifngx][local_ifngy] > 0;
                    bool two_pixel = (local_ifngx != chosen_ifngx || local_ifngy != chosen_ifngy);
                    // bind local free ifng to the old and new cancer cells
                    if(same_pixel_update){
                        double partition = ifng_field[local_ifngx][local_ifngy]/2;
                        if(partition < ifng_tracker[ilatx][ilaty][canID][random_cell][1] - ifng_tracker[ilatx][ilaty][canID][random_cell][0]){
                            ifng_tracker[ilatx][ilaty][canID][random_cell][0] += partition;
                            new_cancer_cell[0] += partition;
                            ifng_field[local_ifngx][local_ifngy] = 0;
                        }else{
                            double ifng_leftover = 2*(partition - (ifng_tracker[ilatx][ilaty][canID][random_cell][1] - ifng_tracker[ilatx][ilaty][canID][random_cell][0]));
                            ifng_tracker[ilatx][ilaty][canID][random_cell][0] = ifng_tracker[ilatx][ilaty][canID][random_cell][1];
                            new_cancer_cell[0] = new_cancer_cell[1];
                            ifng_field[local_ifngx][local_ifngy] = ifng_leftover;
                        }
                    }else if(!same_pixel){
                        if(ifng_field[local_ifngx][local_ifngy] > 0){
                            double ifng_needed = ifng_tracker[ilatx][ilaty][canID][random_cell][1] - ifng_tracker[ilatx][ilaty][canID][random_cell][0];
                            if(ifng_field[local_ifngx][local_ifngy] >= ifng_needed){
                                ifng_field[local_ifngx][local_ifngy] -= ifng_needed;
                                ifng_tracker[ilatx][ilaty][canID][random_cell][0] = ifng_tracker[ilatx][ilaty][canID][random_cell][1];
                            }else{
                                ifng_tracker[ilatx][ilaty][canID][random_cell][0] += ifng_field[local_ifngx][local_ifngy];
                                ifng_field[local_ifngx][local_ifngy] = 0;
                            }
                        }
                        if(ifng_field[chosen_ifngx][chosen_ifngy] > 0){
                            double ifng_needed = new_cancer_cell[1] - new_cancer_cell[0];
                            if(ifng_field[chosen_ifngx][chosen_ifngy] >= ifng_needed){
                                ifng_field[chosen_ifngx][chosen_ifngy] -= ifng_needed;
                                new_cancer_cell[0] = new_cancer_cell[1];
                            }else{
                                new_cancer_cell[0] += ifng_field[chosen_ifngx][chosen_ifngy];
                                ifng_field[chosen_ifngx][chosen_ifngy] = 0;
                            }
                        }
                    }

                    ifng_tracker[clatx][claty][canID].push_back(new_cancer_cell); // update ifng tracker with the new cell and distribute the bound IFNG
                }else { // leaky receptor model
                    double random_value = random->uniform();

                    double gross_fraction = 0;
                    for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                        gross_fraction += ifng_tracker[ilatx][ilaty][canID][cc][0];// add up all the recently absorbed ifng (per one time step)
                    }
                    double chosen_fraction = random_value*gross_fraction;
                    double absorbed_ifng_fraction_counter = 0;
                    int chosen_cell_index_count = 0;
                    for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                        chosen_cell_index_count = cc;
                        absorbed_ifng_fraction_counter += ifng_tracker[ilatx][ilaty][canID][cc][0];
                        if(chosen_fraction <= absorbed_ifng_fraction_counter){
                            break;
                        }
                    }
                    int random_cell = chosen_cell_index_count;
                    double saved_absorption = ifng_tracker[ilatx][ilaty][canID][random_cell][0];
                    ifng_tracker[ilatx][ilaty][canID][random_cell][0] = saved_absorption/2;
                    std::vector<double> new_cancer_cell = {saved_absorption/2, ifng_tracker[ilatx][ilaty][canID][random_cell][1]}; // new cancer cell has same proliferative state

                    ifng_tracker[clatx][claty][canID].push_back(new_cancer_cell); // update ifng tracker with the new cell and distribute the bound IFNG
                }
            }else if(use_cell_dependent_decay_model){
                double random_value = random->uniform();

                double gross_fraction = 0;
                for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                    gross_fraction += ifng_tracker[ilatx][ilaty][canID][cc][0];// add up all the recently absorbed ifng (per one time step)
                }
                double chosen_fraction = random_value*gross_fraction;
                double absorbed_ifng_fraction_counter = 0;
                int chosen_cell_index_count = 0;
                for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                    chosen_cell_index_count = cc;
                    absorbed_ifng_fraction_counter += ifng_tracker[ilatx][ilaty][canID][cc][0];
                    if(chosen_fraction <= absorbed_ifng_fraction_counter){
                        break;
                    }
                }
                int random_cell = chosen_cell_index_count;
                double saved_absorption = ifng_tracker[ilatx][ilaty][canID][random_cell][0];
                std::vector<double> new_cancer_cell = {saved_absorption, ifng_tracker[ilatx][ilaty][canID][random_cell][1]}; // new cancer cell has same proliferative state

                ifng_tracker[clatx][claty][canID].push_back(new_cancer_cell); // update ifng tracker with the new cell and distribute the bound IFNG
            }
            opos.clear();
            popos.clear();
            way = 1;
        }else if(popos.size() > 0) {
            way = 2;
            int chosen = popos[int(r1*popos.size())];
            //int chosen = popos[int(r1*pCount)];
            cancerPos[1] = chosen;
            // int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID};
            int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID, ARThID, TrID, ARTrID};
            int tctracker_map[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
            int icd8s[] = {population[cd8map[0]][chosen], population[cd8map[1]][chosen], population[cd8map[2]][chosen], population[cd8map[3]][chosen], population[cd8map[4]][chosen], population[cd8map[5]][chosen], population[cd8map[6]][chosen], population[cd8map[7]][chosen], population[cd8map[8]][chosen], population[cd8map[9]][chosen]};// activated, exhausted
            int cd8s[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            int cfilling = tpA*population[cd8map[0]][chosen] + tmA*population[cd8map[1]][chosen] + teA*population[cd8map[2]][chosen] + tpA*population[cd8map[3]][chosen] + tmA*population[cd8map[4]][chosen] + teA*population[cd8map[5]][chosen] + cancA*population[canID][chosen] + thA*population[cd8map[6]][chosen] + thA*population[cd8map[7]][chosen] + trA*population[cd8map[8]][chosen] + trA*population[cd8map[9]][chosen];
            int space_left = maxA - cfilling;
            int tcount = cancA - space_left;
            //            std::vector <std::vector <double>>[6] cd8_tracker;
            std::vector <std::vector <double>> cd8_tracker[10];
            double total_free_receptors = 0;
            // double total_bound_receptors = 0;

            //        for(int h=0; h<6; h++){
            //            if(icd8s[h] < 0){
            //                //population[6][5] = 100;
            //            }
            //        }
            // update ifng tracker
            int ilatx = i%latCol;
            int ilaty = int(i/latCol);
            int clatx = chosen%latCol;
            int claty = int(chosen/latCol);

            int local_ifngx = int(ilatx/ifng_lattice_cell_converter);
            int local_ifngy = int(ilaty/ifng_lattice_cell_converter);
            int chosen_ifngx = int(clatx/ifng_lattice_cell_converter);
            int chosen_ifngy = int(claty/ifng_lattice_cell_converter);

            if(!use_simple_decay_model){
                if(!use_leaky_receptor_model) {
                    double random_value =random->uniform();

                    double gross_fraction;
                    for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                        gross_fraction += ifng_tracker[ilatx][ilaty][canID][cc][0]/ifng_tracker[ilatx][ilaty][canID][cc][1];
                    }
                    double chosen_fraction = random_value*gross_fraction;
                    double bound_receptor_fraction_counter = 0;
                    int chosen_cell_index_count = 0;
                    for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                        chosen_cell_index_count = cc;
                        bound_receptor_fraction_counter += ifng_tracker[ilatx][ilaty][canID][cc][0]/ifng_tracker[ilatx][ilaty][canID][cc][1];
                        if(chosen_fraction <= bound_receptor_fraction_counter){
                            break;
                        }
                    }
                    int random_cell = chosen_cell_index_count;
                    double saved_bound = ifng_tracker[ilatx][ilaty][canID][random_cell][0];
                    std::vector <double> new_cancer_cell = std::vector<double>{saved_bound/2, ifng_tracker[ilatx][ilaty][canID][random_cell][1]};
                    ifng_tracker[ilatx][ilaty][canID][random_cell][0] = saved_bound/2;
                    bool same_pixel = (local_ifngx == chosen_ifngx && local_ifngy == chosen_ifngy);
                    bool same_pixel_update = same_pixel && ifng_field[local_ifngx][local_ifngy] > 0;
                    bool two_pixel = (local_ifngx != chosen_ifngx || local_ifngy != chosen_ifngy);
                    // bind local free ifng to the old and new cancer cells
                    if(same_pixel_update){
                        double partition = ifng_field[local_ifngx][local_ifngy]/2;
                        if(partition <= ifng_tracker[ilatx][ilaty][canID][random_cell][1] - ifng_tracker[ilatx][ilaty][canID][random_cell][0]){
                            ifng_tracker[ilatx][ilaty][canID][random_cell][0] += partition;
                            new_cancer_cell[0] += partition;
                            ifng_field[local_ifngx][local_ifngy] = 0;
                        }else{
                            double ifng_leftover = 2*(partition - (ifng_tracker[ilatx][ilaty][canID][random_cell][1] - ifng_tracker[ilatx][ilaty][canID][random_cell][0]));
                            ifng_tracker[ilatx][ilaty][canID][random_cell][0] = ifng_tracker[ilatx][ilaty][canID][random_cell][1];
                            new_cancer_cell[0] = new_cancer_cell[1];
                            ifng_field[local_ifngx][local_ifngy] = ifng_leftover;
                        }
                    }else if(!same_pixel){
                        if(ifng_field[local_ifngx][local_ifngy] > 0){
                            double ifng_needed = ifng_tracker[ilatx][ilaty][canID][random_cell][1] - ifng_tracker[ilatx][ilaty][canID][random_cell][0];
                            if(ifng_field[local_ifngx][local_ifngy] >= ifng_needed){
                                ifng_field[local_ifngx][local_ifngy] -= ifng_needed;
                                ifng_tracker[ilatx][ilaty][canID][random_cell][0] = ifng_tracker[ilatx][ilaty][canID][random_cell][1];
                            }else{
                                ifng_tracker[ilatx][ilaty][canID][random_cell][0] += ifng_field[local_ifngx][local_ifngy];
                                ifng_field[local_ifngx][local_ifngy] = 0;
                            }
                        }
                        if(ifng_field[chosen_ifngx][chosen_ifngy] > 0){
                            double ifng_needed = new_cancer_cell[1] - new_cancer_cell[0];
                            if(ifng_field[chosen_ifngx][chosen_ifngy] >= ifng_needed){
                                ifng_field[chosen_ifngx][chosen_ifngy] -= ifng_needed;
                                new_cancer_cell[0] = new_cancer_cell[1];
                            }else{
                                new_cancer_cell[0] += ifng_field[chosen_ifngx][chosen_ifngy];
                                ifng_field[chosen_ifngx][chosen_ifngy] = 0;
                            }
                        }
                    }

                    ifng_tracker[clatx][claty][canID].push_back(new_cancer_cell); // update ifng tracker with the new cell and distribute the bound IFNG
                }else { // use leaky receptor model
                    double random_value = random->uniform();

                    double gross_fraction = 0;
                    for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                        gross_fraction += ifng_tracker[ilatx][ilaty][canID][cc][0];// add up all the recently bound ifng (per one time step)
                    }
                    double chosen_fraction = random_value*gross_fraction;
                    double absorbed_ifng_fraction_counter = 0;
                    int chosen_cell_index_count = 0;
                    for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                        chosen_cell_index_count = cc;
                        absorbed_ifng_fraction_counter += ifng_tracker[ilatx][ilaty][canID][cc][0];
                        if(chosen_fraction <= absorbed_ifng_fraction_counter){
                            break;
                        }
                    }
                    int random_cell = chosen_cell_index_count;
                    double saved_absorption = ifng_tracker[ilatx][ilaty][canID][random_cell][0];
                    ifng_tracker[ilatx][ilaty][canID][random_cell][0] = saved_absorption/2;
                    std::vector<double> new_cancer_cell = {saved_absorption/2, ifng_tracker[ilatx][ilaty][canID][random_cell][1]}; // new cancer cell has same proliferative state

                    ifng_tracker[clatx][claty][canID].push_back(new_cancer_cell); // update ifng tracker with the new cell and distribute the bound IFNG
                }
            }else if(use_cell_dependent_decay_model){
                double random_value = random->uniform();

                double gross_fraction = 0;
                for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                    gross_fraction += ifng_tracker[ilatx][ilaty][canID][cc][0];// add up all the recently absorbed ifng (per one time step)
                }
                double chosen_fraction = random_value*gross_fraction;
                double absorbed_ifng_fraction_counter = 0;
                int chosen_cell_index_count = 0;
                for(int cc=0; cc < ifng_tracker[ilatx][ilaty][canID].size(); cc++){
                    chosen_cell_index_count = cc;
                    absorbed_ifng_fraction_counter += ifng_tracker[ilatx][ilaty][canID][cc][0];
                    if(chosen_fraction <= absorbed_ifng_fraction_counter){
                        break;
                    }
                }
                int random_cell = chosen_cell_index_count;
                double saved_absorption = ifng_tracker[ilatx][ilaty][canID][random_cell][0];
                std::vector<double> new_cancer_cell = {saved_absorption, ifng_tracker[ilatx][ilaty][canID][random_cell][1]}; // new cancer cell has same proliferative state

                ifng_tracker[clatx][claty][canID].push_back(new_cancer_cell); // update ifng tracker with the new cell and distribute the bound IFNG
            }


            int tcCopy = tcount;
            while(tcCopy > 0){ // decide which t cells will be bumped
                double ran = random->uniform();
                int to_try = int((2*tpA + 2*tmA + 2*teA + 2*thA)*ran); // would be 2*ran but if teA=0, we only choose Tcs. This is a fix for now. Later should do something more general if necessary.
                if(icd8s[to_try] > 0){
                    icd8s[to_try]--;
                    cd8s[to_try]++;
                    tcCopy--;
                    population[cd8map[to_try]][chosen]--;
                    // remove cells from chosen ifng tracker site and save them in cd8 tracker
                    if(!use_simple_decay_model) {
                        if(!use_leaky_receptor_model) {
                            double random_cell = random->uniform();
                            int random_cell_int = int(random_cell*ifng_tracker[clatx][claty][tctracker_map[to_try]].size());
                            cd8_tracker[to_try].push_back(ifng_tracker[clatx][claty][tctracker_map[to_try]][random_cell]);

                            // total_bound_receptors += ifng_tracker[clatx][claty][cd8map[to_try]][random_cell][0];
                            total_free_receptors += ifng_tracker[clatx][claty][tctracker_map[to_try]][random_cell][1] - ifng_tracker[clatx][claty][tctracker_map[to_try]][random_cell][0];
                            ifng_tracker[clatx][claty][tctracker_map[to_try]].erase(ifng_tracker[clatx][claty][tctracker_map[to_try]].begin() + random_cell_int);
                        }else { // use leaky receptor model
                            double random_cell = random->uniform();
                            int random_cell_int = int(random_cell*ifng_tracker[clatx][claty][tctracker_map[to_try]].size());
                            cd8_tracker[to_try].push_back(ifng_tracker[clatx][claty][tctracker_map[to_try]][random_cell]);

                            ifng_tracker[clatx][claty][tctracker_map[to_try]].erase(ifng_tracker[clatx][claty][tctracker_map[to_try]].begin() + random_cell_int);
                        }
                    }else if(use_cell_dependent_decay_model) {
                        double random_cell = random->uniform();
                        int random_cell_int = int(random_cell*ifng_tracker[clatx][claty][tctracker_map[to_try]].size());
                        cd8_tracker[cd8map[to_try]-1].push_back(ifng_tracker[clatx][claty][tctracker_map[to_try]][random_cell]);

                        ifng_tracker[clatx][claty][tctracker_map[to_try]].erase(ifng_tracker[clatx][claty][tctracker_map[to_try]].begin() + random_cell_int);
                    }
                    if(population[cd8map[to_try]][chosen] < 0){
                        population[errID][chosen] = 47;
                    }
                }
            }
            population[canID][chosen]++;

            bool done = false;
            int width = 3;
            int cx = chosen%latCol;
            int cy = chosen/latCol;
            //if(chosen < 0){// || chosen > latCol*latCol){
            //population[6][chosen] = 3;
            //}
            while(!done){
                for(int ly=0; ly < 2; ly++){
                    int halfw = int((width-1)/2);
                    for(int lx=0; lx < width; lx++){
                        int target1x = cx + lx - halfw;
                        int target1y = cy + (2*ly - 1)*halfw;
                        int target2x = cx + (2*ly - 1)*halfw;
                        int target2y = cy + lx - halfw;
                        //periodic boundary conditions:
                        target1x -= floor(double(target1x)/latCol)*latCol;
                        target1y -= floor(double(target1y)/latCol)*latCol;
                        target2x -= floor(double(target2x)/latCol)*latCol;
                        target2y -= floor(double(target2y)/latCol)*latCol;
                        int target1 = target1x + latCol*target1y;
                        int target2 = target2x + latCol*target2y;
                        int voxHad1 = tpA*population[cd8map[0]][target1] + tmA*population[cd8map[1]][target1]
                                + teA*population[cd8map[2]][target1] + tpA*population[cd8map[3]][target1]
                                + tmA*population[cd8map[4]][target1] + teA*population[cd8map[5]][target1]
                        + cancA*population[canID][target1]+ thA*population[cd8map[6]][target1]+
                            thA*population[cd8map[7]][target1] + thA*population[cd8map[8]][target1] + thA*population[cd8map[9]][target1];
                        int voxHad2 = tpA*population[cd8map[0]][target2] + tmA*population[cd8map[1]][target2]
                                + teA*population[cd8map[2]][target2] + tpA*population[cd8map[3]][target2]
                                + tmA*population[cd8map[4]][target2] + teA*population[cd8map[5]][target2]
                                + cancA*population[canID][target2] + thA*population[cd8map[6]][target2]+
                                    thA*population[cd8map[7]][target2] + thA*population[cd8map[8]][target2] +
                                        thA*population[cd8map[9]][target2];
                        if(tcount <= maxA - voxHad1 && tcount > 0){
                            int ifngx = int(target1x/ifng_lattice_cell_converter);
                            int ifngy = int(target1y/ifng_lattice_cell_converter);

                            if(!use_simple_decay_model) {
                                if(!use_leaky_receptor_model) {
                                    bool update_binding = (ifngx != chosen_ifngx || ifngy != chosen_ifngy) && ifng_field[ifngx][ifngy] > 0;
                                    double free_leftover = ifng_field[ifngx][ifngy] - total_free_receptors; // gives how much ifng will be left over after filling all cells
                                    double leftover = ifng_field[ifngx][ifngy]; // gives how much more ifng needs to be distributed to receptors

                                    bool saturate = free_leftover >= 0;
                                    if(saturate) {
                                        ifng_field[ifngx][ifngy] = free_leftover;
                                    }else{
                                        ifng_field[ifngx][ifngy] = 0; //  2024-09-11, added to delete ifng when t cells eat it
                                    }
                                    bool began = false;
                                    while((saturate && !began) || (!saturate && leftover > 0) || (!update_binding && !began)){
                                        int filled_count = 0;
                                        double partition_fraction = leftover/total_free_receptors;
                                        began = true;
                                        for (int cd = 0; cd < 10; cd++) { // loop through CD8 types
                                            population[cd8map[cd]][target1] += cd8s[cd];
                                            if (update_binding) {
                                                for (int cc = 0; cc < cd8_tracker[cd8map[cd]-1].size(); cc++) { // loop over cells of type
                                                    //                                        std::vector<double> cell = cd8_tracker[cd8map[cd]][cc];
                                                    double free_receptors = cd8_tracker[cd8map[cd]-1][cc][1] - cd8_tracker[cd8map[cd]-1][cc][0];
                                                    double partition = partition_fraction*free_receptors;
                                                    if (saturate) { // if local ifng can fill all the receptors
                                                        cd8_tracker[cd8map[cd]-1][cc][0] = cd8_tracker[cd8map[cd]-1][cc][1];
                                                    } else { // if only a fraction of cell receptors can be filled
                                                        if(cd8_tracker[cd8map[cd]-1][cc][0] + partition < cd8_tracker[cd8map[cd]-1][cc][1]){ // if partition does not saturate cell
                                                            cd8_tracker[cd8map[cd]-1][cc][0] += partition;
                                                            total_free_receptors -= partition;
                                                            leftover -= partition;
                                                        }else if(cd8_tracker[cd8map[cd]-1][cc][1] != cd8_tracker[cd8map[cd]-1][cc][0]){  // if partition does saturate cell and cell is not already saturated
                                                            filled_count += 1;
                                                            total_free_receptors -= free_receptors;
                                                            cd8_tracker[cd8map[cd]-1][cc][0] = cd8_tracker[cd8map[cd]-1][cc][1];
                                                            double partition_left = partition - free_receptors;
                                                            leftover -= free_receptors;
                                                        }
                                                    }
                                                }
                                            }
                                            cd8s[cd] = 0;
                                        }
                                        if(filled_count == 0){
                                            leftover = 0;
                                        }
                                    }
                                    for (int cd = 0; cd < 8; cd++) { // loop through CD8 types
                                        for (int cc = 0; cc < cd8_tracker[cd8map[cd]-1].size(); cc++) { // loop over cells of type
                                            ifng_tracker[target1x][target1y][cd].push_back(cd8_tracker[cd][cc]);
                                        }
                                    }
                                }else { // leaky receptor model
                                    for (int cd = 0; cd < 8; cd++) { // loop through CD8 types
                                        population[cd8map[cd]][target1] += cd8s[cd];
                                        for (int cc = 0; cc < cd8_tracker[cd].size(); cc++) { // loop over cells of type
                                            ifng_tracker[target1x][target1y][tctracker_map[cd]].push_back(cd8_tracker[cd][cc]);
                                        }
                                    }
                                }
                            }else if(use_cell_dependent_decay_model) {
                                for (int cd = 0; cd < 8; cd++) { // loop through CD8 types
                                    for (int cc = 0; cc < cd8_tracker[cd].size(); cc++) { // loop over cells of type
                                        ifng_tracker[target1x][target1y][cd].push_back(cd8_tracker[cd][cc]);
                                    }
                                }
                            }
                            tcellPos.push_back(target1);
                            done = true;
                            tcount = 0;
                            goto stop;
                        }else if(tcount <= maxA - voxHad2 && tcount > 0){

                            int ifngx = int(target2x/ifng_lattice_cell_converter);
                            int ifngy = int(target2y/ifng_lattice_cell_converter);

                            if(!use_simple_decay_model) {
                                if(!use_leaky_receptor_model){
                                    bool update_binding = (ifngx != chosen_ifngx || ifngy != chosen_ifngy) && ifng_field[ifngx][ifngy] > 0;
                                    double free_leftover = ifng_field[ifngx][ifngy] - total_free_receptors; // gives how much ifng will be left over after filling all cells
                                    double leftover = ifng_field[ifngx][ifngy]; // gives how much more ifng needs to be distributed to receptors

                                    bool saturate = free_leftover >= 0;
                                    if(saturate) {
                                        ifng_field[ifngx][ifngy] = free_leftover;
                                    }else{
                                        ifng_field[ifngx][ifngy] = 0;
                                    }
                                    bool began = false;
                                    while((saturate && !began) || (!saturate && leftover > 0) || (!update_binding && !began)){
                                        int filled_count = 0;
                                        double partition_fraction = leftover/total_free_receptors;
                                        began = true;
                                        for (int cd = 0; cd < 6; cd++) { // loop through CD8 types
                                            population[cd8map[cd]][target2] += cd8s[cd];
                                            if (update_binding) {
                                                for (int cc = 0; cc < cd8_tracker[cd8map[cd]-1].size(); cc++) { // loop over cells of type
                                                    //                                        std::vector<double> cell = cd8_tracker[cd8map[cd]][cc];
                                                    double free_receptors = cd8_tracker[cd8map[cd]-1][cc][1] - cd8_tracker[cd8map[cd]-1][cc][0];
                                                    double partition = partition_fraction*free_receptors;
                                                    if (saturate) { // if local ifng can fill all the receptors
                                                        cd8_tracker[cd8map[cd]-1][cc][0] = cd8_tracker[cd8map[cd]-1][cc][1];
                                                    } else { // if only a fraction of cell receptors can be filled
                                                        if(cd8_tracker[cd8map[cd]-1][cc][0]+ partition < cd8_tracker[cd8map[cd]-1][cc][1]){ // if partition does not saturate cell
                                                            cd8_tracker[cd8map[cd]-1][cc][0] += partition;
                                                            total_free_receptors -= partition;
                                                            leftover -= partition;
                                                        }else if(cd8_tracker[cd8map[cd]-1][cc][1] != cd8_tracker[cd8map[cd]-1][cc][0]){  // if partition does saturate cell
                                                            filled_count += 1;
                                                            total_free_receptors -= free_receptors;
                                                            cd8_tracker[cd8map[cd]-1][cc][0] = cd8_tracker[cd8map[cd]-1][cc][1];
                                                            double partition_left = partition - free_receptors;
                                                            leftover -= free_receptors;
                                                        }
                                                    }
                                                }
                                            }
                                            cd8s[cd] = 0;
                                        }
                                        if(filled_count == 0){
                                            leftover = 0;
                                        }
                                    }
                                    for (int cd = 0; cd < 6; cd++) { // loop through CD8 types
                                        for (int cc = 0; cc < cd8_tracker[cd8map[cd]-1].size(); cc++) { // loop over cells of type
                                            ifng_tracker[target2x][target2y][cd8map[cd]].push_back(
                                                    cd8_tracker[cd8map[cd]-1][cc]);
                                        }
                                    }
                                }else {  // leaky receptor model
                                    for (int cd = 0; cd < 8; cd++) { // loop through CD8 types
                                        population[cd8map[cd]][target2] += cd8s[cd];
                                        for (int cc = 0; cc < cd8_tracker[cd].size(); cc++) { // loop over cells of type
                                            ifng_tracker[target2x][target2y][tctracker_map[cd]].push_back(cd8_tracker[cd][cc]);
                                        }
                                    }
                                }
                            }else if(use_cell_dependent_decay_model) {
                                for (int cd = 0; cd < 8; cd++) { // loop through CD8 types
                                    for (int cc = 0; cc < cd8_tracker[cd8map[cd]-1].size(); cc++) { // loop over cells of type
                                        ifng_tracker[target2x][target2y][tctracker_map[cd]].push_back(cd8_tracker[cd][cc]);
                                    }
                                }
                            }
                            tcellPos.push_back(target2);
                            //tcellPos[tcCount] = target2;
                            //tcCount++;
                            done = true;
                            tcount = 0;
                            goto stop;
                        }else if(tcount > 0){
                            if(maxA - voxHad1 > 0){
                                int left = maxA - voxHad1;
                                while(left > 0 && tcount > 0){
                                    double rando = random->uniform(); //dist(mt);
                                    int try_out = int((2*tpA + 2*tmA + 2*teA + 2*thA)*rando);
                                    if(cd8s[try_out] > 0){
                                        population[cd8map[try_out]][target1]++;
                                        if(!use_simple_decay_model) {
                                            if(!use_leaky_receptor_model) {
                                                // note that this fills the free receptors of one cell at a time so not all at once so the first to a single spot might get more than it should
                                                // choose a random cell with bound ifng to be transported here
                                                double random_bound_cell = random->uniform();
                                                int int_random_bound_cell = int(random_bound_cell*cd8_tracker[cd8map[try_out]-1].size());
                                                std::vector <double> cell = cd8_tracker[cd8map[try_out]-1][int_random_bound_cell];
                                                int ifngx = int(target1x/ifng_lattice_cell_converter);
                                                int ifngy = int(target1y/ifng_lattice_cell_converter);
                                                if(ifngx != chosen_ifngx || ifngy != chosen_ifngy) {
                                                    if (ifng_field[ifngx][ifngy] > 0 && cell[0] < cell[1]) {
                                                        double need_bound = cell[1] - cell[0];
                                                        if (ifng_field[ifngx][ifngy] > need_bound) {
                                                            cell[0] = cell[1];
                                                            ifng_field[ifngx][ifngy] -= need_bound;
                                                        } else {
                                                            cell[0] += ifng_field[ifngx][ifngy];
                                                            ifng_field[ifngx][ifngy] = 0;
                                                        }
                                                    }
                                                }
                                                ifng_tracker[target1x][target1y][cd8map[try_out]].push_back(cell);
                                                cd8_tracker[cd8map[try_out]-1].erase(cd8_tracker[cd8map[try_out]-1].begin() + int_random_bound_cell);
                                            }else {  // leaky receptor model
                                                // choose a random CD8+ T cell to be transported here
                                                double random_bound_cell = random->uniform();
                                                int int_random_bound_cell = int(random_bound_cell*cd8_tracker[try_out].size());
                                                std::vector <double> cell = cd8_tracker[try_out][int_random_bound_cell];
                                                ifng_tracker[target1x][target1y][tctracker_map[try_out]].push_back(cell);
                                                cd8_tracker[try_out].erase(cd8_tracker[try_out].begin() + int_random_bound_cell);
                                            }
                                        }else if(use_cell_dependent_decay_model) {
                                            // choose a random CD8+ T cell to be transported here
                                            double random_bound_cell = random->uniform();
                                            int int_random_bound_cell = int(random_bound_cell*cd8_tracker[cd8map[try_out]-1].size());
                                            std::vector <double> cell = cd8_tracker[cd8map[try_out]-1][int_random_bound_cell];
                                            ifng_tracker[target1x][target1y][cd8map[try_out]].push_back(cell);
                                            cd8_tracker[cd8map[try_out]-1].erase(cd8_tracker[cd8map[try_out]-1].begin() + int_random_bound_cell);
                                        }
                                        cd8s[try_out]--;
                                        left--;
                                        tcount--;
                                    }
                                }
                                tcellPos.push_back(target1);
                            }
                            if(maxA - voxHad2 > 0 && tcount > 0 && target1 != target2){
                                int left = maxA - voxHad2;
                                while(left > 0 && tcount > 0){
                                    double rando = random->uniform(); //dist(mt);
                                    int try_out = int((2*tpA + 2*tmA + 2*teA + 2*thA)*rando);
                                    if(cd8s[try_out] > 0){
                                        population[cd8map[try_out]][target2]++;

                                        if(!use_simple_decay_model) {
                                            if(!use_leaky_receptor_model) {
                                                // choose a random cell with bound ifng to be transported here
                                                double random_bound_cell = random->uniform();
                                                int int_random_bound_cell = int(random_bound_cell*cd8_tracker[cd8map[try_out]-1].size());
                                                std::vector <double> cell = cd8_tracker[cd8map[try_out]-1][int_random_bound_cell];
                                                int ifngx = int(target2x/ifng_lattice_cell_converter);
                                                int ifngy = int(target2y/ifng_lattice_cell_converter);
                                                if(ifngx != chosen_ifngx || ifngy != chosen_ifngy) { // if entering a new ifng pixel
                                                    if (ifng_field[ifngx][ifngy] > 0 && cell[0] < cell[1]) {
                                                        double need_bound = cell[1] - cell[0];
                                                        if (ifng_field[ifngx][ifngy] > need_bound) {
                                                            cell[0] = cell[1];
                                                            ifng_field[ifngx][ifngy] -= need_bound;
                                                        } else {
                                                            cell[0] += ifng_field[ifngx][ifngy];
                                                            ifng_field[ifngx][ifngy] = 0;
                                                        }
                                                    }
                                                }
                                                ifng_tracker[target2x][target2y][cd8map[try_out]].push_back(cell);
                                                cd8_tracker[cd8map[try_out]-1].erase(cd8_tracker[cd8map[try_out]-1].begin() + int_random_bound_cell);
                                            } else {  // leaky receptor model
                                                // choose a random CD8+ T cell to be transported here
                                                double random_bound_cell = random->uniform();
                                                int int_random_bound_cell = int(random_bound_cell*cd8_tracker[try_out].size());
                                                std::vector <double> cell = cd8_tracker[try_out][int_random_bound_cell];
                                                ifng_tracker[target2x][target2y][tctracker_map[try_out]].push_back(cell);
                                                cd8_tracker[try_out].erase(cd8_tracker[try_out].begin() + int_random_bound_cell);
                                            }
                                        }else if(use_cell_dependent_decay_model) {
                                            // choose a random CD8+ T cell to be transported here
                                            double random_bound_cell = random->uniform();
                                            int int_random_bound_cell = int(random_bound_cell*cd8_tracker[cd8map[try_out]-1].size());
                                            std::vector <double> cell = cd8_tracker[cd8map[try_out]-1][int_random_bound_cell];
                                            ifng_tracker[target2x][target2y][cd8map[try_out]].push_back(cell);
                                            cd8_tracker[cd8map[try_out]-1].erase(cd8_tracker[cd8map[try_out]-1].begin() + int_random_bound_cell);
                                        }
                                        cd8s[try_out]--;
                                        left--;
                                        tcount--;
                                    }
                                }
                                tcellPos.push_back(target2);
                            }
                            if(tcount <= 0){done=true; goto stop;}
                            for(int h=0; h<6; h++){
                                if(cd8s[h] < 0){
                                    population[errID][5] = 30000;
                                }
                            }
                        }
                        if(tcount <= 0){done=true; goto stop;}
                    }
                }
                width += 2;
                if(width > 200){done=true; goto stop;} // stops the search for space for a new T cell when cancer prolfierates if there is no space left. The T cell is just removed. This is only when the cancer is winning to a large extent (at the end of the simulation
            }
        }else{population[errID][4] = 20000; way = 3;}
        stop:
        popos.clear();
        opos.clear();

        /*for(int cd=0; cd < latCol*latCol; cd++){
            if(population[5][cd] < 0){// || population[0][i] + population[2][i] + population[5][i] < 0){
                population[6][cd] = 3;
            }else if(population[0][cd] < 0){
                population[6][cd] = 7;
            }else if(population[2][cd] < 0){
                population[6][cd] = 11;
            }
    //		if(population[0][i] + population[2][i] + population[5][i] < 0){
    //			population[6][i] = 7;
    //		}
        }*/
        //population[0][0] = population[0][0];
    //    for(int ch=0; ch<tcellPos.size(); ch++){ // some old debugging check
    //        if(tcellPos[ch] == 7659){
    //            tcellThere = true;
    //        }
    //    }
    }

    // motility
    if((13 <= which && which <= 18) || (39 <= which && which <= 42)){ // motility interactions
        if(population[latEID][i] >= 1 && population[latEID][j] >= 1 && false){
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

            if(j/latCol == 0){ // top
              cornCount2++;
              type2 = 0;
            }else if(j/latCol == latCol - 1){ // bottom
              cornCount2++;
              type2 = 2;
            }
            if(j%latCol == 0){ // left
              cornCount2++;
              type2 = 1;
            }else if(j%latCol == latCol - 1){ // right
              cornCount2++;
              type2 = 3;
            }

            if((cornCount == 2 && cornCount2 == 2) || (cornCount == 1 && cornCount2 == 1 && type != type2)){
              cross = true;
            }

            if(cross){
                double dissProb = 0.5;
                int tIndex = which - 13; // decide which cell type is moving
                int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID};
                int tType = cd8map[tIndex];

                //std::random_device rt;
                //std::mt19937 mta(rt());
                //std::uniform_real_distribution<double> distr(0, 1.0);
                double ranGen = random->uniform();//distr(mta);

                if(dissProb >= ranGen){
                    if(population[tType][j] > 0){
                      population[tType][j] -= 1;
                    }else{
                      population[errID][456] = 43;
                      std::cout << "problem with boundary" << std::endl;
                    }
                    proliferation_update_bool = true;
                }
            }
        }

        //update ifng_tracker
        if(!use_simple_decay_model) {
            int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID, ARThID, TrID, ARTrID};
            int cell_type = 0;
            int track_index = 0;
            if(which <= 18) {
                cell_type = cd8map[which - 13];
                track_index = cell_type;
            }else {
                cell_type = cd8map[which - 33];
                track_index = which - 32;
            }

            int ilatx = i%latCol;
            int ilaty = int(i/latCol);
            int jlatx = j%latCol;
            int jlaty = int(j/latCol);

            if(!use_leaky_receptor_model) {
                int i_ifngx = int(ilatx/ifng_lattice_cell_converter);
                int i_ifngy = int(ilaty/ifng_lattice_cell_converter);
                int j_ifngx = int(jlatx/ifng_lattice_cell_converter);
                int j_ifngy = int(jlaty/ifng_lattice_cell_converter);

                double randomv = random->uniform();
                int random_cell = int(randomv*ifng_tracker[ilatx][ilaty][track_index].size());
                std::vector <double> cell = ifng_tracker[ilatx][ilaty][track_index][random_cell];
                ifng_tracker[ilatx][ilaty][cell_type].erase(ifng_tracker[ilatx][ilaty][track_index].begin() + random_cell);

                bool same_pixel = (i_ifngx == j_ifngx && i_ifngy == j_ifngy);

                if(!same_pixel && ifng_field[j_ifngx][j_ifngy] > 0){
                    double ifng_needed = cell[1]-cell[0];
                    if(ifng_field[j_ifngx][j_ifngy] >= ifng_needed){
                        cell[0] = cell[1];
                        ifng_field[j_ifngx][j_ifngy] -= ifng_needed;
                    }else{
                        cell[0] += ifng_field[j_ifngx][j_ifngy];
                        ifng_field[j_ifngx][j_ifngy] = 0;
                    }
                }
                ifng_tracker[jlatx][jlaty][cell_type].push_back(cell);
            }else {  // leaky receptor model
                // added code which chooses cell which jumps based on receptors bound
                double random_value = random->uniform();

                double gross_fraction = 0;
                for(int cc=0; cc < ifng_tracker[ilatx][ilaty][track_index].size(); cc++){
                    gross_fraction += ifng_tracker[ilatx][ilaty][track_index][cc][0];// add up all the recently absorbed ifng (per one time step)
                }
                double chosen_fraction = random_value*gross_fraction;
                double absorbed_ifng_fraction_counter = 0;
                int chosen_cell_index_count = 0;
                for(int cc=0; cc < ifng_tracker[ilatx][ilaty][track_index].size(); cc++){
                    chosen_cell_index_count = cc;
                    absorbed_ifng_fraction_counter += ifng_tracker[ilatx][ilaty][track_index][cc][0];
                    if(chosen_fraction <= absorbed_ifng_fraction_counter){
                        break;
                    }
                }
                int random_cell = chosen_cell_index_count;
                // double randomv = random->uniform();
                // int random_cell = int(randomv*ifng_tracker[ilatx][ilaty][cell_type].size());
                std::vector <double> cell = ifng_tracker[ilatx][ilaty][track_index][random_cell];
                ifng_tracker[ilatx][ilaty][track_index].erase(ifng_tracker[ilatx][ilaty][track_index].begin() + random_cell);
                ifng_tracker[jlatx][jlaty][track_index].push_back(cell);
            }
        }else if(use_cell_dependent_decay_model) {
            int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID};
            int cell_type = cd8map[which - 13];

            int ilatx = i%latCol;
            int ilaty = int(i/latCol);
            int jlatx = j%latCol;
            int jlaty = int(j/latCol);

            double randomv = random->uniform();
            int random_cell = int(randomv*ifng_tracker[ilatx][ilaty][cell_type].size());
            std::vector <double> cell = ifng_tracker[ilatx][ilaty][cell_type][random_cell];
            ifng_tracker[ilatx][ilaty][cell_type].erase(ifng_tracker[ilatx][ilaty][cell_type].begin() + random_cell);
            ifng_tracker[jlatx][jlaty][cell_type].push_back(cell);
        }
    }

    int recSpot;
    if(7 <= which && which <= 12){ // Tp and Tm recruitment
      int tIndex = which - 7; // decide which cell type is moving
      // int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID};
        int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID, ARThID, TrID, ARTrID};

      int weights[] = {tpA, tmA, teA, tpA, tmA, teA, thA, thA, trA, trA};
      int receptor_map[] = {2, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1};

      int recType = cd8map[tIndex];
      int weight = weights[tIndex];

      double ranSpot = random->uniform(); //spDist(mtd);

      int eCount = 0;
      for(int px=0; px < latCol*latCol; px++){
          if(population[boundID][px] > 0) {
              int filling = tpA * population[TpID][px] + tmA * population[TmID][px] + teA * population[TeID][px] +
                            tpA * population[ARTpID][px] + tmA * population[ARTmID][px] + teA * population[ARTeID][px] +
                            cancA * population[canID][px] + thA*population[cd8map[6]][px] + thA*population[cd8map[7]][px] + trA*population[cd8map[8]][px] + trA*population[cd8map[9]][px];
              int left = int((maxA - filling) / weight);
              eCount += left;
          }
      }
      ranSpot = (int)(eCount * ranSpot);
      eCount = 0;

      for(int px=0; px < latCol*latCol; px++){
          if(population[boundID][px] > 0) {
              int filling = tpA * population[TpID][px] + tmA * population[TmID][px] + teA * population[TeID][px] +
                            tpA * population[ARTpID][px] + tmA * population[ARTmID][px] + teA * population[ARTeID][px] +
                            cancA * population[canID][px] + thA*population[cd8map[6]][px] + thA*population[cd8map[7]][px] + trA*population[cd8map[8]][px] + trA*population[cd8map[9]][px];
              int left = int((maxA - filling) / weight);
              if (left > 0 && eCount <= ranSpot && ranSpot <= eCount + left) {
                  population[recType][px] += 1;
                  recSpot = px;
                  break;
              }
              eCount += left;
          }
      }
        if(!use_simple_decay_model) {
            if(!use_leaky_receptor_model) {
                // bind ifng to new cell, cell type- recType, location = recSpot
                int ilatx = recSpot%latCol;
                int ilaty = int(recSpot/latCol);

                int i_ifngx = int(ilatx/ifng_lattice_cell_converter);
                int i_ifngy = int(ilaty/ifng_lattice_cell_converter);
                std::vector<double> cell = {0, (double) receptors_per_cell[receptor_map[recType]]};

                if(ifng_field[i_ifngx][i_ifngy] > 0){
                    double ifng_needed = cell[1];
                    if(ifng_field[i_ifngx][i_ifngy] >= ifng_needed){
                        cell[0] = cell[1];
                        ifng_field[i_ifngx][i_ifngy] -= ifng_needed;
                    }else{
                        cell[0] += ifng_field[i_ifngx][i_ifngy];
                        ifng_field[i_ifngx][i_ifngy] = 0;
                    }
                }
                ifng_tracker[ilatx][ilaty][recType].push_back(cell);
            }else {  // leaky receptor model
                // bind ifng to new cell, cell type- recType, location = recSpot
                int ilatx = recSpot%latCol;
                int ilaty = int(recSpot/latCol);

                std::vector<double> cell = {0, (double) receptors_per_cell[receptor_map[recType]]};

                ifng_tracker[ilatx][ilaty][recType].push_back(cell);
            }
        }else if(use_cell_dependent_decay_model) {
            // bind ifng to new cell, cell type- recType, location = recSpot
            int ilatx = recSpot%latCol;
            int ilaty = int(recSpot/latCol);

            std::vector<double> cell = {0, (double) receptors_per_cell[receptor_map[recType]]};

            ifng_tracker[ilatx][ilaty][recType].push_back(cell);
        }
    }

    // Boundary producing random T cells (according to current distribution of T cells) WAS OFF
    if(which == 35 && false){
        //std::random_device rdv;
        //std::mt19937 mtt(rdv());
        //std::uniform_real_distribution<double> distrib(0,1.0);
        double rand = random->uniform(); //distrib(mtt);

        int cd8map[] = {TpID, TmID, TeID}; // Activated, exhausted
        int cd8Count[] = {0,0, 0};
        int totcd8 = 0;
        for(int ii=0; ii < latCol*latCol; ii++){
            for(int ci=0; ci < 3; ci++){
                cd8Count[ci] += population[cd8map[ci]][ii];
                totcd8 += population[cd8map[ci]][ii];
            }
        }

        int decision = round(totcd8*rand);
        int label = 0;
        int w_type = -1;
        for(int ii=0; ii < 3; ii++){
            if(label <= decision && decision <= label + cd8Count[ii]){
                w_type = ii;
            }
            label += cd8Count[ii];
        }

        population[cd8map[w_type]][i] += 1;
    }
    //check to see if boundary labels are correct
    //int county = 0;
    //for(int m=0;m<10000;m++){
    //	if(population[7][m] >= 1){
    //		county += population[7][m];
    //	}
    //}
    //population[6][63] = county;

    std::vector<int> nLocs;
    //int latCol = 100; //  SHOULD BE 100
    bool problem = false;
    int wProblem[] = {0,0,0,0,0,0,0};
    for(int ii=0; ii < latCol*latCol; ii++){
        if(population[TpID][ii] < 0){// || population[0][i] + population[2][i] + population[5][i] < 0){
            //population[outSigID][ii] = 3;
            problem = true;
            wProblem[0] += 1;
            nLocs.push_back(ii);
            //std::cout << "negative " << << std::endl;
        }else if(population[TmID][ii] < 0){
            //population[outSigID][ii] = 7;
            problem = true;
            wProblem[1] += 1;
            nLocs.push_back(ii);
            //std::cout << "negative";
        }else if(population[TeID][ii] < 0){
            //population[outSigID][ii] = 7;
            problem = true;
            wProblem[2] += 1;
            nLocs.push_back(ii);
            //std::cout << "negative";
        }else if(population[ARTpID][ii] < 0){// || population[0][i] + population[2][i] + population[5][i] < 0){
            //population[outSigID][ii] = 3;
            problem = true;
            wProblem[3] += 1;
            nLocs.push_back(ii);
            //std::cout << "negative " << << std::endl;
        }else if(population[ARTmID][ii] < 0){
            //population[outSigID][ii] = 7;
            problem = true;
            wProblem[4] += 1;
            nLocs.push_back(ii);
            //std::cout << "negative";
        }else if(population[ARTeID][ii] < 0){
            //population[outSigID][ii] = 7;
            problem = true;
            wProblem[5] += 1;
            nLocs.push_back(ii);
            //std::cout << "negative";
        }else if(population[canID][ii] < 0){
            //population[outSigID][ii] = 11;
            problem = true;
            wProblem[6] += 1;
            nLocs.push_back(ii);
            //std::cout << "negative";
        }
    }
    nLocs.clear();

    if(!use_simple_decay_model) {
        if(!use_leaky_receptor_model) {
            if(1 <= which && which <= 6){ // proliferation, update of ifng receptor binding
                int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID, ARThID, TrID, ARTrID};

                int cellid = cd8map[which-1];
                int ilatx = i%latCol;
                int ilaty = int(i/latCol);

                int local_ifngx = int(ilatx/ifng_lattice_cell_converter);
                int local_ifngy = int(ilaty/ifng_lattice_cell_converter);

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
                double saved_bound = ifng_tracker[ilatx][ilaty][cellid][random_cell][0];
                std::vector <double> new_cd8 = std::vector<double>{saved_bound/2, ifng_tracker[ilatx][ilaty][cellid][random_cell][1]};
                ifng_tracker[ilatx][ilaty][cellid][random_cell][0] = saved_bound/2;
                bool bind_update = ifng_field[local_ifngx][local_ifngy] > 0;
                if(bind_update){
                    double ifng_needed = 2*(new_cd8[1] - new_cd8[0]);
                    double partition = ifng_field[local_ifngx][local_ifngy]/2;
                    if(ifng_field[local_ifngx][local_ifngy] >= ifng_needed){
                        new_cd8[0] = new_cd8[1];
                        ifng_tracker[ilatx][ilaty][cellid][random_cell][0] = ifng_tracker[ilatx][ilaty][cellid][random_cell][1];
                        ifng_field[local_ifngx][local_ifngy] -= ifng_needed;
                    }else{
                        ifng_field[local_ifngx][local_ifngy] = 0;
                        new_cd8[0] += partition;
                        ifng_tracker[ilatx][ilaty][cellid][random_cell][0] += partition;
                    }
                }
                ifng_tracker[ilatx][ilaty][cellid].push_back(new_cd8);
            }else if(which == 25 || which == 26 || (35 <= which && which <= 38)){ // T cell death, update of ifng receptor binding
                int cd8map[] = {TeID, ARTeID, TpID, TmID, ARTpID, ARTmID};
                int cellid;
                if(which <= 26){
                    cellid = cd8map[which-25];
                }else{
                    cellid = cd8map[which-33];
                }
                int ilatx = i%latCol;
                int ilaty = int(i/latCol);

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
                ifng_tracker[ilatx][ilaty][cellid].erase(ifng_tracker[ilatx][ilaty][cellid].begin() + random_cell);
            }else if(27 <= which && which <= 34){ // cancer cell death, update of ifng receptor binding
                bool local = j == -1;
                int ilatx = i%latCol;
                int ilaty = int(i/latCol);

                if(!local){
                    ilatx = j%latCol;
                    ilaty = int(j/latCol);
                }

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][canID].size());
                ifng_tracker[ilatx][ilaty][canID].erase(ifng_tracker[ilatx][ilaty][canID].begin() + random_cell);
            }else if(which == 19 || which == 20 || which == 59){ // Tp differentiation, update of ifng receptor binding
                int cd8map[] = {TpID, ARTpID};
                int target_cd8map[] = {TmID, ARTmID};
				int cellid;
				int target_cellid;
				if(which <=20){
                	target_cellid = target_cd8map[which - 19];
					cellid = cd8map[which - 19];
				}else if(which == 59){
					target_cellid = target_cd8map[1];
					cellid = cd8map[0];
				}

                int ilatx = i%latCol;
                int ilaty = int(i/latCol);

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
                std::vector<double> cell = ifng_tracker[ilatx][ilaty][cellid][random_cell];
                ifng_tracker[ilatx][ilaty][cellid].erase(ifng_tracker[ilatx][ilaty][cellid].begin() + random_cell);
                ifng_tracker[ilatx][ilaty][target_cellid].push_back(cell);
            }else if(43 <= which && which <= 50) { // Th induced differentiation of tp to tm
                int cd8map[] = {TpID, ARTpID};
                int target_cd8map[] = {TmID, ARTmID}; // Th 2 local, Th 2 nbr, AR+ Th 2 local, AR+ Th 2 nbr
                bool arP = int(which-43)%2; // is the interaction target AR+? 0 if AR-
                bool local = int((which-43)/2)%2 == 0; // is the interaction local
                int ilatx = int(i%latCol);
                int ilaty = int(i/latCol);

                if(!local) { // Th
                    ilatx = j%latCol;
                    ilaty = int(j/latCol);
                }

                int cellid = cd8map[arP];
                int target_cellid = target_cd8map[arP];

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
                std::vector<double> cell = ifng_tracker[ilatx][ilaty][cellid][random_cell];
                ifng_tracker[ilatx][ilaty][cellid].erase(ifng_tracker[ilatx][ilaty][cellid].begin() + random_cell);
                ifng_tracker[ilatx][ilaty][target_cellid].push_back(cell);

            }else if(21 <= which && which <= 24){  // Tm exhaustion, update of ifng receptor binding
                int cd8map[] = {TmID, ARTmID};
                int target_cd8map[] = {TeID, ARTeID};
                int cellid = cd8map[int((which - 21)/2)];
                int target_cellid = target_cd8map[int((which - 21)/2)];
                bool local = j == -1;

                int jlatx = i%latCol;
                int jlaty = int(i/latCol);
                if(!local){
                    jlatx = j%latCol;
                    jlaty = int(j/latCol);
                }

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[jlatx][jlaty][cellid].size());
                std::vector<double> cell = ifng_tracker[jlatx][jlaty][cellid][random_cell];
                ifng_tracker[jlatx][jlaty][cellid].erase(ifng_tracker[jlatx][jlaty][cellid].begin() + random_cell);
                ifng_tracker[jlatx][jlaty][target_cellid].push_back(cell);
            }
        }else {  // leaky receptor model
            if(1 <= which && which <= 6){ // proliferation, update of ifng receptor binding
                int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID};

                int cellid = cd8map[which-1];
                int ilatx = i%latCol;
                int ilaty = int(i/latCol);

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
                double saved_absorption = ifng_tracker[ilatx][ilaty][cellid][random_cell][0];
                std::vector <double> new_cd8 = std::vector<double>{saved_absorption/2, ifng_tracker[ilatx][ilaty][cellid][random_cell][1]};
                ifng_tracker[ilatx][ilaty][cellid][random_cell][0] = saved_absorption/2;
                ifng_tracker[ilatx][ilaty][cellid].push_back(new_cd8);
            }else if(which == 25 || which == 26 || (35 <= which && which <= 38)){ // T cell death, update of ifng receptor binding
                int cd8map[] = {TeID, ARTeID, TpID, TmID, ARTpID, ARTmID};
                int cellid;
                if(which <= 26){
                    cellid = cd8map[which-25];
                }else{
                    cellid = cd8map[which-33];
                }
                int ilatx = i%latCol;
                int ilaty = int(i/latCol);

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
                ifng_tracker[ilatx][ilaty][cellid].erase(ifng_tracker[ilatx][ilaty][cellid].begin() + random_cell);
            }else if(27 <= which && which <= 34){ // cancer cell death, update of ifng receptor binding
                bool local = j == -1;
                int ilatx = i%latCol;
                int ilaty = int(i/latCol);

                if(!local){
                    ilatx = j%latCol;
                    ilaty = int(j/latCol);
                }

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][canID].size());
                ifng_tracker[ilatx][ilaty][canID].erase(ifng_tracker[ilatx][ilaty][canID].begin() + random_cell);
            }else if((which == 19 || which == 20 || which == 59)){ // Tp differentiation, update of ifng receptor binding
                int cd8map[] = {TpID, ARTpID};
                int target_cd8map[] = {TmID, ARTmID};
				int cellid;
				int target_cellid;
				if(which <=20){
                	target_cellid = target_cd8map[which - 19];
					cellid = cd8map[which - 19];
				}else if(which == 59){
					target_cellid = target_cd8map[1];
					cellid = cd8map[0];
				}

                int ilatx = i%latCol;
                int ilaty = int(i/latCol);

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
                std::vector<double> cell = ifng_tracker[ilatx][ilaty][cellid][random_cell];
                ifng_tracker[ilatx][ilaty][cellid].erase(ifng_tracker[ilatx][ilaty][cellid].begin() + random_cell);
                ifng_tracker[ilatx][ilaty][target_cellid].push_back(cell);
            }else if(43 <= which && which <= 50) { // Th induced differentiation of tp to tm
                int cd8map[] = {TpID, ARTpID};
                int target_cd8map[] = {TmID, ARTmID}; // Th 2 local, Th 2 nbr, AR+ Th 2 local, AR+ Th 2 nbr
                int arP = int(which-43)%2; // is the interaction target AR+? 0 if AR-
                // bool local = int((which-43)/2)%2 == 0; // is the interaction local
                bool local = j == -1;
                int ilatx = int(i%latCol);
                int ilaty = int(i/latCol);

                if(!local) { // Th
                    ilatx = j%latCol;
                    ilaty = int(j/latCol);
                }

                int cellid = cd8map[arP];
                int target_cellid = target_cd8map[arP];

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
                std::vector<double> cell = ifng_tracker[ilatx][ilaty][cellid][random_cell];
                ifng_tracker[ilatx][ilaty][cellid].erase(ifng_tracker[ilatx][ilaty][cellid].begin() + random_cell);
                ifng_tracker[ilatx][ilaty][target_cellid].push_back(cell);

            }else if(21 <= which && which <= 24){  // Tm exhaustion, update of ifng receptor binding
                int cd8map[] = {TmID, ARTmID};
                int target_cd8map[] = {TeID, ARTeID};
                int cellid = cd8map[int((which - 21)/2)];
                int target_cellid = target_cd8map[int((which - 21)/2)];
                bool local = j == -1;

                int jlatx = i%latCol;
                int jlaty = int(i/latCol);
                if(!local){
                    jlatx = j%latCol;
                    jlaty = int(j/latCol);
                }

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[jlatx][jlaty][cellid].size());
                std::vector<double> cell = ifng_tracker[jlatx][jlaty][cellid][random_cell];
                ifng_tracker[jlatx][jlaty][cellid].erase(ifng_tracker[jlatx][jlaty][cellid].begin() + random_cell);
                ifng_tracker[jlatx][jlaty][target_cellid].push_back(cell);
            }else if(51 <= which && which <= 58){  // Tm exhaustion, update of ifng receptor binding
                int cd8map[] = {TmID, ARTmID};
                int target_cd8map[] = {TeID, ARTeID};

                int arP = int((which-51)/2)%2; // is the interaction target AR+? 0 if AR-
                // bool local = int(which-51)%2 == 0; // is the interaction local
                int cellid = cd8map[arP];
                int target_cellid = target_cd8map[arP];
                bool local = j == -1;

                int jlatx = i%latCol;
                int jlaty = int(i/latCol);
                if(!local){
                    jlatx = j%latCol;
                    jlaty = int(j/latCol);
                }

                double random_cell =random->uniform();
                random_cell = int(random_cell*ifng_tracker[jlatx][jlaty][cellid].size());
                std::vector<double> cell = ifng_tracker[jlatx][jlaty][cellid][random_cell];
                ifng_tracker[jlatx][jlaty][cellid].erase(ifng_tracker[jlatx][jlaty][cellid].begin() + random_cell);
                ifng_tracker[jlatx][jlaty][target_cellid].push_back(cell);
            }
        }
    }else if(use_cell_dependent_decay_model) {
        if(1 <= which && which <= 6){ // proliferation, update of ifng receptor binding
            int cd8map[] = {TpID, TmID, TeID, ARTpID, ARTmID, ARTeID};

            int cellid = cd8map[which-1];
            int ilatx = i%latCol;
            int ilaty = int(i/latCol);

            double random_cell =random->uniform();
            random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
            double saved_absorption = ifng_tracker[ilatx][ilaty][cellid][random_cell][0];
            std::vector <double> new_cd8 = std::vector<double>{saved_absorption, ifng_tracker[ilatx][ilaty][cellid][random_cell][1]};
            ifng_tracker[ilatx][ilaty][cellid].push_back(new_cd8);
        }else if(which == 25 || which == 26 || (35 <= which && which <= 38)){ // T cell death, update of ifng receptor binding
            int cd8map[] = {TeID, ARTeID, TpID, TmID, ARTpID, ARTmID};
            int cellid;
            if(which <= 26){
                cellid = cd8map[which-25];
            }else{
                cellid = cd8map[which-33];
            }
            int ilatx = i%latCol;
            int ilaty = int(i/latCol);

            double random_cell =random->uniform();
            random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
            ifng_tracker[ilatx][ilaty][cellid].erase(ifng_tracker[ilatx][ilaty][cellid].begin() + random_cell);
        }else if(27 <= which && which <= 34){ // cancer cell death, update of ifng receptor binding
            bool local = j == -1;
            int ilatx = i%latCol;
            int ilaty = int(i/latCol);

            if(!local){
                ilatx = j%latCol;
                ilaty = int(j/latCol);
            }

            double random_cell =random->uniform();
            random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][canID].size());
            ifng_tracker[ilatx][ilaty][canID].erase(ifng_tracker[ilatx][ilaty][canID].begin() + random_cell);
        }else if(which == 19 || which == 20 || which == 59){ // Tp differentiation, update of ifng receptor binding
            int cd8map[] = {TpID, ARTpID};
            int target_cd8map[] = {TmID, ARTmID};
			int cellid;
			int target_cellid;
			if(which <=20){
                target_cellid = target_cd8map[which - 19];
				cellid = cd8map[which - 19];
			}else if(which == 59){
				target_cellid = target_cd8map[1];
				cellid = cd8map[0];
			}

            int ilatx = i%latCol;
            int ilaty = int(i/latCol);

            double random_cell =random->uniform();
            random_cell = int(random_cell*ifng_tracker[ilatx][ilaty][cellid].size());
            std::vector<double> cell = ifng_tracker[ilatx][ilaty][cellid][random_cell];
            ifng_tracker[ilatx][ilaty][cellid].erase(ifng_tracker[ilatx][ilaty][cellid].begin() + random_cell);
            ifng_tracker[ilatx][ilaty][target_cellid].push_back(cell);
        }else if(21 <= which && which <= 24){  // Tm exhaustion, update of ifng receptor binding
            int cd8map[] = {TmID, ARTmID};
            int target_cd8map[] = {TeID, ARTeID};
            int cellid = cd8map[int((which - 21)/2)];
            int target_cellid = target_cd8map[int((which - 21)/2)];
            bool local = j == -1;

            int jlatx = i%latCol;
            int jlaty = int(i/latCol);
            if(!local){
                jlatx = j%latCol;
                jlaty = int(j/latCol);
            }

            double random_cell =random->uniform();
            random_cell = int(random_cell*ifng_tracker[jlatx][jlaty][cellid].size());
            std::vector<double> cell = ifng_tracker[jlatx][jlaty][cellid][random_cell];
            ifng_tracker[jlatx][jlaty][cellid].erase(ifng_tracker[jlatx][jlaty][cellid].begin() + random_cell);
            ifng_tracker[jlatx][jlaty][target_cellid].push_back(cell);
        }
    }

    if((1 <= which && which <= 12) || which == 25 || which == 26 || (35 <= which && which <= 38)){ // proliferation, recruitment, Te death, cells floating in (these update when the T cell # number changes?)
      proliferation_update_bool = true;
    }

    // bias the T cell random walks
    if((which == 0 || (27 <= which && which <= 34))){ // if a cancer cell proliferates or if a cancer cell is lysed, we need to update all T cells for biased motility
        proliferation_update_bool = true;
    //    std::array<std::array<double, 12>, 12> discretized_cancer; // discretize cancer pop to larger lattice constant

        for(int cpx=0; cpx < cpMax; cpx++){
            for(int cpy=0; cpy < cpMay; cpy++){
                for(int nb=0; nb<2; nb++) {
                    cancer_bias_potential[cpx][cpy][nb] = 0; // neighbors: 0-left/right, 1-above (subtracted in y)/below
                }
                discretized_cancer[cpx][cpy] = 0;
            }
        }

        int max_cancer = 0;
        for(int popx=0; popx < latCol; popx++){
            for(int popy=0; popy < latCol; popy++){
                discretized_cancer[int(popx/cancer_potential_converter)][int(popy/cancer_potential_converter)] += population[canID][popx + latCol*popy];
                int cancer_value = discretized_cancer[int(popx/cancer_potential_converter)][int(popy/cancer_potential_converter)];
                if(cancer_value > max_cancer){
                    max_cancer = cancer_value;
                }
            }
        }

        for(int popx=0; popx < latCol; popx++){
            for(int popy=0; popy < latCol; popy++){
                int pcx = int(popx/cancer_potential_converter);
                int pcy = int(popy/cancer_potential_converter);
                int number_derivatives = 0;
                double average_derivative = 0;

                for(int axis=0; axis<2; axis++){ // search neighbors
    //                int neighbor_cancer_pops[] = {-1, -1}; // left, right
                    for(int ns=0; ns<2; ns++){ // left then right
                        int nx = pcx + (1-axis)*(2*ns - 1);
                        int ny = pcy + axis*(2*ns - 1);
                        if(0 < nx && nx < cpMax && 0 <= ny && ny < cpMay){
                            number_derivatives += 1;
                            int left_or_right = 2*ns - 1;
                            average_derivative += left_or_right*(discretized_cancer[nx][ny] - discretized_cancer[pcx][pcy])/max_cancer;
    //                        neighbor_cancer_pops[ns] = ; // left then right
                        }
                    }
                    average_derivative = average_derivative/number_derivatives;
                    cancer_bias_potential[pcx][pcy][axis] = (average_derivative+1)/2; // normalized average derivative, scaled to 0->1
                }
            }
        }
    }



    if(proliferation_update_bool){ // reset the total count of cd8+ T cells
         tcell_total_population = 0;
         for(int ic=0; ic < latCol*latCol; ic++){
             tcell_total_population += population[TpID][ic] + population[TmID][ic]+population[TeID][ic] + population[ARTpID][ic] + population[ARTmID][ic]+population[ARTeID][ic];
         }
         if(tcell_total_population > max_T_cells){
             tcell_total_population = max_T_cells;
         }
      //proliferation_update_bool = true;
    }

    rxn_count[which]++;

    // compute propensity changes for participating sites and first neighs
    // ignore update of sites with isite < 0
    // use echeck[] to avoid resetting propensity of same site

    int nsites = 0;

    if(update_all_always){
        for(int nx=0; nx<latCol*latCol; nx++){
            int isite = i2site[nx];
            propensity[isite] = site_propensity(nx);
            esites[nsites++] = isite;
            echeck[isite] = 1;
        }
    }else{

        int isite = i2site[i];
        propensity[isite] = site_propensity(i);
        esites[nsites++] = isite;
        echeck[isite] = 1;

        for (n = 0; n < numneigh[i]; n++) {
            m = neighbor[i][n];
            isite = i2site[m];
            if (isite >= 0 && echeck[isite] == 0) {
                propensity[isite] = site_propensity(m);
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
        }

        if (rstyle == NBR) {
            for (n = 0; n < numneigh[j]; n++) {
                m = neighbor[j][n];
                isite = i2site[m];
                if (isite >= 0 && echeck[isite] == 0) {
                    propensity[isite] = site_propensity(m);
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
                }
            }
        }
        // cancerPos tcellPos
        std::vector<int> cSites;
        //cSites.push_back(9433);
        if (which == 0 || (7 <= which && which <= 12) || (27 <= which && which <=
                                                                         34)) { // cancer prolif, tp tm and te recruitment, local cancer death, nbr cancer death
            if (which == 0) { // for farther updates for cancer proliferation

                cSites.push_back(cancerPos[1]);
                //cSites[0] = cancerPos[1];
                for (int ts = 0; ts < tcellPos.size(); ts++) {
                    cSites.push_back(tcellPos[ts]);
                }
                //removed the following and replaced with lattice-wide update for testing purposes

                for (int ci = 0; ci < tcellPos.size() + 1; ci++) {
                    //for(int ci=0;ci<tcCount+1; ci++){
                    int cx;
                    int cy;
                    int width = 3;
                    if (ci == 0) { // deal with cancer cell site
                        cx = cancerPos[1] % latCol; // target cancer proliferation site coords
                        cy = cancerPos[1] / latCol;
                    } else {
                        cx = tcellPos[ci - 1] % latCol;
                        cy = tcellPos[ci - 1] / latCol;
                    }
                    for (int lo = 0; lo < neighLim; lo++) { // updates around cancer cells to a certain depth
                        for (int ly = 0; ly < 2; ly++) { //lx and ly take us around the neighbors
                            int halfw = int((width - 1) / 2);
                            for (int lx = 0; lx < width; lx++) {
                                int target1x = cx + lx - halfw;
                                int target1y = cy + (2 * ly - 1) * halfw;
                                int target2x = cx + (2 * ly - 1) * halfw;
                                int target2y = cy + lx - halfw;
                                //periodic boundary conditions:
                                target1x -= floor(double(target1x) / latCol) * latCol;
                                target1y -= floor(double(target1y) / latCol) * latCol;
                                target2x -= floor(double(target2x) / latCol) * latCol;
                                target2y -= floor(double(target2y) / latCol) * latCol;
                                int target1 = target1x + latCol * target1y;
                                int target2 = target2x + latCol * target2y;
                                cSites.push_back(target1);
                                cSites.push_back(target2);
                            }
                        }
                        width += 2;
                        if (ci > 0) {
                            break;
                        }
                    }
                }
                /*int IDsToUse[] = {canID, cd8aID, cd8eID, macID};
                for(int hw=0; hw < latCol*latCol; hw++){
                  for(int sp=0; sp < 4; sp++){
                    if(population[IDsToUse[sp]][hw] > 0){
                      cSites.push_back(hw);
                      break;
                    }
                  }
                }*/
            } else if (27 <= which && which <= 34) { // cancer lysis
                int cancPos = -1;
                int width = 3;
                if ((which % 2) == 0) { // nbr lysis (even)
                    cancPos = j;
                } else { // local lysis (even)
                    cancPos = i;
                }

                // pasted code below...
                int cx = cancPos % latCol;
                int cy = cancPos / latCol;
                for (int lo = 0; lo < neighLim; lo++) { // updates around cancer cells to a certain depth
                    for (int ly = 0; ly < 2; ly++) { // lx and ly take us around the neighbors
                        int halfw = int((width - 1) / 2);
                        for (int lx = 0; lx < width; lx++) {
                            int target1x = cx + lx - halfw;
                            int target1y = cy + (2 * ly - 1) * halfw;
                            int target2x = cx + (2 * ly - 1) * halfw;
                            int target2y = cy + lx - halfw;
                            //periodic boundary conditions:
                            target1x -= floor(double(target1x) / latCol) * latCol;
                            target1y -= floor(double(target1y) / latCol) * latCol;
                            target2x -= floor(double(target2x) / latCol) * latCol;
                            target2y -= floor(double(target2y) / latCol) * latCol;
                            int target1 = target1x + latCol * target1y;
                            int target2 = target2x + latCol * target2y;
                            cSites.push_back(target1);
                            cSites.push_back(target2);
                        }
                    }
                    width += 2;
                }

                // above

            } else if (7 <= which && which <= 12) {  // recruitment
                //proliferation_update_bool = true;
                int recX = recSpot % latCol;
                int recY = recSpot / latCol;
                for (int tx = 0; tx <
                                 3; tx++) { // update neighboring corner cancer sites why? cancer proliferation propensity doesn't care about surrounding cd8s
                    for (int ty = 0; ty < 3; ty++) {
                        int targx = recX + tx - 1;
                        int targy = recY + ty - 1;
                        targx -= floor(double(targx) / latCol) * latCol;
                        targy -= floor(double(targy) / latCol) * latCol;
                        int targ = targx + latCol * targy;
                        cSites.push_back(targ);
                    }
                }
            }

            sort(cSites.begin(), cSites.end());
            cSites.erase(unique(cSites.begin(), cSites.end()), cSites.end());

            if (!proliferation_update_bool) {
                for (n = 0; n < cSites.size(); n++) {
                    m = cSites[n];
                    isite = i2site[m];
                    if (isite >= 0 && echeck[isite] == 0) {
                        propensity[isite] = site_propensity(m);
                        esites[nsites++] = isite;
                        echeck[isite] = 1;
                    }
                }
            }

            //std::cout << nsites  << " " << 819 << std::endl;
            //for (m = 0; m < ncSites; m++) cSites[m] = 0;
            //std::cout << i  << " " << 821 << std::endl;
            //cSites.resize(0);
            //cSites.clear();
            //std::vector<int>().swap(cSites);
            //cSites.swap(std::vector<int>(temp));
        }

        // SHOULD NOT ALWAYS BE TRUE
        //proliferation_update_bool = true;

        if (proliferation_update_bool) {
            //cSites.clear();
            for (int la = 0; la < latCol * latCol; la++) {
                if (population[TpID][la] > 0 || population[TmID][la] > 0 || population[TeID][la] > 0 ||
                    population[ARTpID][la] > 0 || population[ARTmID][la] > 0 || population[ARTeID][la] > 0 ||
                    population[recID][la] >
                    0) {  // make sure to update all sites that recruit or proliferative T cells with the new rates (population threshold), this is also used to update all cytotoxic CD8s on the ifng field
                    cSites.push_back(la);
                }
            }
            proliferation_update_bool = false;

            //std::sort(cSites.begin(), cSites.end());
            //auto last = std::unique(cSites.begin(), cSites.end());
            //cSites.erase(last, cSites.end());
            //cSites.clear();
        }

        //std::vector<int> canSites;

        if (cancer_proliferation_update_bool) {
            for (int la = 0; la < latCol * latCol; la++) {
                if (population[canID][la] > 0) {  // make sure to update all sites that recruit or proliferative T cells with the new rates
                    cSites.push_back(la);
                }
            }
            cancer_proliferation_update_bool = false;
            //    sort(cSites.begin(), cSites.end());
            //    cSites.erase(unique(cSites.begin(), cSites.end()), cSites.end());
            //
            //    for (n = 0; n < canSites.size(); n++) {
            //        m = canSites[n];
            //        isite = i2site[m];
            //        if (isite >= 0 && echeck[isite] == 0) {
            //            propensity[isite] = site_propensity(m);
            //            esites[nsites++] = isite;
            //            echeck[isite] = 1;
            //        }
            //    }

        }

        sort(cSites.begin(), cSites.end());
        cSites.erase(unique(cSites.begin(), cSites.end()), cSites.end());

        for (int ch = 0; ch < cSites.size(); ch++) {
            int x = cSites[ch]; // use with break point to check what sites update
        }

        for (n = 0; n < cSites.size(); n++) {
            m = cSites[n];
            isite = i2site[m];
            if (isite >= 0 && echeck[isite] == 0) {
                propensity[isite] = site_propensity(m);
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
        }

        //canSites.clear();
        cSites.clear();
        tcellPos.clear();

        //if((which != 1 && which != 9) && )
    }


    ///////////////////// NEED TO PUT UPDATE FOR LONG-RANGE CANCER GROWTH HERE AND WHATEVER EFFECTS PAST NEAREST NEIGHBOR

  solve->update(nsites,esites,propensity); // gdag- updates nsites # of sites in esites array, do I need to add sites that are farther away that I update?
  //std::cout << i  << " " << 831 << std::endl;
  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
  //std::cout << i  << " " << 835 << std::endl;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppRxnDiff::clear_events(int i)// gdag- each set of events in the list of events for site i is a permutation loop. This closes the most recent set of events that refer to each other into a new loop. Maybe this is a way to reuse old interaction propensities without recacalculating for certain model types?
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppRxnDiff::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;//gdag- Could we have memory issues here with our customizations???
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ---------------------------------------------------------------------- */

void AppRxnDiff::add_rxn(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal reaction command");

  // store ID

  if (find_reaction(arg[0]) >= 0) {
    char *str = new char[128];
    sprintf(str,"Reaction ID %s already exists",arg[0]);
    error->all(FLERR,str);
  }

  int n = nreactions + 1;
  rname = (char **) memory->srealloc(rname,n*sizeof(char *),
					  "rxn/diff:rname");
  int nlen = strlen(arg[0]) + 1;
  rname[nreactions] = new char[nlen];
  strcpy(rname[nreactions],arg[0]);

  // grow reaction arrays
  memory->grow(localReactants,n,MAX_SPECIES,"rxn/diff:localReactants");
  memory->grow(nbrReactants,n,MAX_SPECIES,"rxn/diff:nbrReactants");
  memory->grow(localDeltaPop,n,MAX_SPECIES,"rxn/diff:localDeltaPop");
  memory->grow(nbrDeltaPop,n,MAX_SPECIES,"rxn/diff:nbrDeltaPop");
  memory->grow(rate,n,"rxn/diff:rate");
  memory->grow(rxnStyle,n,"rxn/diff:rxnStyle");

  for (int ispecies = 0; ispecies < MAX_SPECIES; ispecies++){
    localReactants[nreactions][ispecies] = 0;
    nbrReactants[nreactions][ispecies] = 0;
    localDeltaPop[nreactions][ispecies] = 0;
    nbrDeltaPop[nreactions][ispecies] = 0;
  }

  // find which arg is numeric reaction rate
  char c;
  int rateArg = 1;
  while (rateArg < narg) {
    c = arg[rateArg][0];
    if ((c >= '0' && c <= '9') || c == '+' || c == '-' || c == '.') break;
    rateArg++;
  }

  // find which args are the local and nbr reactants and products
  int localReactArg, nbrReactArg, localProdArg, nbrProdArg;
  localReactArg = nbrReactArg = localProdArg = nbrProdArg = 0;
  for (int iarg = 1; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"local") == 0){
      if (iarg<rateArg) localReactArg = iarg;
      else localProdArg = iarg;
    } else if (strcmp(arg[iarg],"nbr") == 0){
      if (iarg<rateArg) nbrReactArg = iarg;
      else nbrProdArg = iarg;
    }
  }

  // error checks
  if (rateArg == narg) error->all(FLERR,"Reaction has no numeric rate");
  if (!localReactArg || !nbrReactArg || !localProdArg || !nbrProdArg)
    error->all(FLERR,"One or more 'local' or 'nbr' keywords missing");
  if (localReactArg > nbrReactArg || nbrReactArg > rateArg || rateArg > localProdArg ||
      localProdArg > nbrProdArg)
    error->all(FLERR,"One or more 'local' or 'nbr' keywords out of order");

  // extract reactant and product species names
  int nLocalReactant = 0;
  for (int i = localReactArg+1; i < nbrReactArg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    localReactants[nreactions][ispecies]++;
    localDeltaPop[nreactions][ispecies]--;
    nLocalReactant++;
  }

  int nNbrReactant = 0;
  for (int i = nbrReactArg+1; i < rateArg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    nbrReactants[nreactions][ispecies]++;
    nbrDeltaPop[nreactions][ispecies]--;
    nNbrReactant++;
  }

  rate[nreactions] = atof(arg[rateArg]);

  int nLocalProduct = 0;
  for (int i = localProdArg+1; i < nbrProdArg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    localDeltaPop[nreactions][ispecies]++;
    nLocalProduct++;
  }
  
  int nNbrProduct = 0;
  for (int i = nbrProdArg+1; i < narg; i++) {
    int ispecies = find_species(arg[i]);
    if (ispecies == -1) error->all(FLERR,"Unknown species in reaction command");
    nbrDeltaPop[nreactions][ispecies]++;
    nNbrProduct++;
  }

  // Set the reaction style as local only or not
  rxnStyle[nreactions] = NBR;
  if(!nNbrReactant && !nNbrProduct){
    rxnStyle[nreactions] = LOCAL;
  }

  // additional error checking
  if(!nLocalReactant && !nLocalProduct && !nNbrReactant && !nNbrProduct)
    error->all(FLERR,"Reaction must have at least one reactant or product");
  if(!nLocalReactant && !nLocalProduct)
    error->all(FLERR,"Reaction cannot only act on neighbor");

  nreactions++;
}

/* ---------------------------------------------------------------------- */

void AppRxnDiff::add_species(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal species command");

  // grow species arrays

  int n = nspecies + narg;
  sname = (char **) memory->srealloc(sname,n*sizeof(char *),
					  "rxn/diff:sname");

  for (int iarg = 0; iarg < narg; iarg++) {
    if (find_species(arg[iarg]) >= 0) {
      char *str = new char[128];
      sprintf(str,"Species ID %s already exists",arg[iarg]);
      error->all(FLERR,str);
    }
    int nlen = strlen(arg[iarg]) + 1;
    sname[nspecies+iarg] = new char[nlen];
    strcpy(sname[nspecies+iarg],arg[iarg]);
  }
  nspecies += narg;
}

/* ----------------------------------------------------------------------
   return reaction index (0 to N-1) for a reaction ID
   return -1 if doesn't exist
------------------------------------------------------------------------- */

int AppRxnDiff::find_reaction(char *str)
{
  for (int i = 0; i < nreactions; i++)
    if (strcmp(str,rname[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   return species index (0 to N-1) for a species ID
   return -1 if doesn't exist
------------------------------------------------------------------------- */

int AppRxnDiff::find_species(char *str)
{
  for (int i = 0; i < nspecies; i++)
    if (strcmp(str,sname[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
 * Stub for custom rate multiplier function to be replaced in rxn/diff/custom
------------------------------------------------------------------------- */
double AppRxnDiff::custom_multiplier(int i, int rstyle, int which, int jpartner)
{
  return 1.0;
}

