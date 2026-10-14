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
  ndouble = 0;
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
  max_cell_area_per_voxel = atoi(arg[1]);
  cancer_area_count = atoi(arg[2]);
  cytotoxic_T_area_count = atoi(arg[3]);
  exhausted_T_area_count = atof(arg[4]);
  mac_area_count = atof(arg[5]);
  int latCol = 100;
  int macID = 6; 
  int cd8aID = 0; 
  int canID = 1; 
  int cd8eID = 2; 
  int recID = 4; 
  int latEID = 7; 
  int deadID = 3; 
  int outSigID = 5;
  
  tcell_total_population = 0;/*
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

void AppRxnDiff::grow_app()
{
  population = iarray;
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

    // esites must be large enough for 2 sites and their 1st neighbors
    esites = new int[10000];//new int[118];//2 + 2*maxneigh]; gdag edit, for 2 cancer cell sites, their neighbors, their displaced t cells and their neighbors: 2: 25+16+4*9, 3: 49+33+4*9
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
          // population[ispecies]i[i] choose localReactants[m][ispecies]
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
      for (int jj = 0; jj < numneigh[i]; jj++) {
        j = neighbor[i][jj];
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
  for (int ispecies=0; ispecies < MAX_SPECIES; ispecies++){
    population[ispecies][i] += localDeltaPop[which][ispecies];
    if(population[ispecies][i] < 0){
      prob = true;
    }
  }
  if (rstyle == NBR){
    for (int ispecies=0; ispecies < MAX_SPECIES; ispecies++){
      population[ispecies][j] += nbrDeltaPop[which][ispecies];
      if(population[ispecies][j] < 0){
        prob = true;
      }
    }
  }
  if(prob && which != 1 && which != 9 && population[5][877] == 0){
    population[5][877] = i + 10;
  }
  else if(prob && population[5][877] == 0){
    population[5][877] = 30;
  }
//////////////////////////////////////////////////// GIUSEPPE GIULIANI EDITS ///////////////////////////////////////////
int latCol = 100;

int macID = 6;
int cd8aID = 0;
int canID = 1;
int cd8eID = 2;
int recID = 4;
int latEID = 7;
int deadID = 3;
int outSigID = 5;

//area units of different cell types
int maxA = max_cell_area_per_voxel;
int tcA = cytotoxic_T_area_count;
int teA = exhausted_T_area_count;
int macA = mac_area_count;
int cancA = cancer_area_count;
//teA = 1; // we want exhausted cells to always be bumped in space in the same manner

int cancerPos[] = {0,0};
std::vector<int> tcellPos;
//int tcellPos[] = {0,0,0,0};
int tcCount = 0;
int neighLim = 3; // number of layers of neighbors that cancer can jump to
int widthLim = 3 + 2*(neighLim-1);
if(which == 0){
	int latCol = 100; ///////// SHOULD BE 100
	/*std::random_device rd;
    std::mt19937 mt(rd());
  	std::uniform_real_distribution<double> dist(0,1.0);
	double r1 = dist(mt); // neeeeeeeeeeeeedd to put this back!!!!!*/
    bool tcellThere = false;
    double r1 = random->uniform();
    cancerPos[0] = i;
	std::vector<int> opos;
	//int opos[16];
    //int oCount = 0;
	std::vector<int> popos;
	//int popos[16];
    //int pCount = 0;
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
                            int target2x = ix + (2*ly-1)*halfw;
                            int target2y = iy + lx - halfw;
			            	target1x -= floor(double(target1x)/latCol)*latCol;
                            target1y -= floor(double(target1y)/latCol)*latCol;
                            target2x -= floor(double(target2x)/latCol)*latCol;
                            target2y -= floor(double(target2y)/latCol)*latCol;
                            int target1 = target1x + latCol*target1y;
                            int target2 = target2x + latCol*target2y;
				            int tcellCount1 = population[cd8aID][target1] + population[cd8eID][target1];
				            int tcellCount2 = population[cd8aID][target2] + population[cd8eID][target2];
                            int filling1 = tcA*population[cd8aID][target1] + teA*population[cd8eID][target1] + cancA*population[canID][target1] + macA*population[macID][target1];
                            int filling2 = tcA*population[cd8aID][target2] + teA*population[cd8eID][target2] + cancA*population[canID][target2] + macA*population[macID][target2];
                            if(cancA*population[canID][target1] + macA*population[macID][target1] < maxA){
					            if(maxA - filling1 >= cancA){
						            opos.push_back(target1);
                                    //opos[oCount] = target1;
                                    //oCount++;
					            }else{
						            popos.push_back(target1);
						            //popos[pCount] = target1;
                                    //pCount++;
					            }
				            }
				            if(cancA*population[canID][target2] + macA*population[macID][target2]  < maxA){
					            if(maxA - filling2 >= 2){
                       	            opos.push_back(target2);
                       	            //opos[oCount] = target1;
                                    //oCount++;
					            }else{
                       	            popos.push_back(target2);
                       	            //popos[pCount] = target1;
                                    //pCount++;
					            }
                            }
                        }
        }
		if(opos.size() > 0 || popos.size() > 0 || width >= widthLim){
            done =  true;
	        if(width >= widthLim && !(opos.size() > 0 || popos.size() > 0)){
				population[outSigID][7] = 1000;
			}
            std::sort(opos.begin(), opos.end());
            auto first = std::unique(opos.begin(), opos.end());
            opos.erase(first, opos.end());
            
            std::sort(popos.begin(), popos.end());
            auto sec = std::unique(popos.begin(), popos.end());
            popos.erase(sec, popos.end()); 
//else if(opos.size() > 0){
//				population[4][5] = 3000;
//			}else if(popos.size() > 0){
//                              population[4][5] = 4000;
//			}
        }else{width += 2;}
    }
//	tpopos = opos;
//	tpopos.insert(tpopos.end(), popos.begin(), popos.end());
	//if(opos.size() > 0){
    if(opos.size() > 0){
		int choice = opos[int(r1*opos.size())];
		//int choice = opos[int(r1*oCount)];
		population[canID][choice]++;
		cancerPos[1] = choice;
		opos.clear();
		popos.clear();
		way = 1;
	//}else if(popos.size() > 0){
	/*}else if(popos.size() > 0){
		way = 2;
		int chosen = popos[int(r1*popos.size())];
        //int chosen = popos[int(r1*pCount)];
		cancerPos[1] = chosen;
		int cd8map[] = {5, 0, 2};
		opos.clear();
		popos.clear();
		way = 1;
	//}else if(popos.size() > 0){
	*/}else if(popos.size() > 0){
		way = 2;
		int chosen = popos[int(r1*popos.size())];
        //int chosen = popos[int(r1*pCount)];
		cancerPos[1] = chosen;
		int cd8map[] = {cd8aID, cd8eID};
		int icd8s[] = {population[cd8map[0]][chosen], population[cd8map[1]][chosen]};// activated, exhausted
        int cd8s[] = {0,0};
        int cfilling = tcA*population[cd8map[0]][chosen] + teA*population[cd8map[1]][chosen] + cancA*population[canID][chosen] + macA*population[macID][chosen];
        int space_left = maxA - cfilling;
        int tcount = cancA - space_left;
		for(int h=0; h<3; h++){
                	if(icd8s[h] < 0){
                        	//population[6][5] = 100;
                        }
        }

        int tcCopy = tcount;
        while(tcCopy > 0){ // decide which t cells will be bumped
          double ran = random->uniform();
          int to_try = int((teA+tcA)*ran); // would be 2*ran but if teA=0, we only choose Tcs. This is a fix for now. Later should do something more general if necessary.
          if(icd8s[to_try] > 0){
            icd8s[to_try]--;
            cd8s[to_try]++;
            tcCopy--;
            population[cd8map[to_try]][chosen]--;
            if(population[cd8map[to_try]][chosen] < 0){
              population[outSigID][chosen] = 47;
            }
          }
        }
		//int tcount = cd8s[0] + cd8s[1] + cd8s[2];
		//population[cd8map[0]][chosen] = 0;
		//population[cd8map[1]][chosen] = 0;
		//population[cd8map[2]][chosen] = 0;
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
					int target2x = cx + (2*ly-1)*halfw;
					int target2y = cy + lx - halfw;
					//periodic boundary conditions:
					target1x -= floor(double(target1x)/latCol)*latCol;
					target1y -= floor(double(target1y)/latCol)*latCol;
					target2x -= floor(double(target2x)/latCol)*latCol;
                    target2y -= floor(double(target2y)/latCol)*latCol;
					int target1 = target1x + latCol*target1y;
					int target2 = target2x + latCol*target2y;
					//if(target1 < 0 || target1 >latCol*latCol || target2 < 0 || target2 > latCol*latCol){
						//population[6][chosen] = 3;
					//}
					int voxHad1 = tcA*population[cd8map[0]][target1]
						+ teA*population[cd8map[1]][target1]
						 + cancA*population[canID][target1] + macA*population[macID][target1];
					int voxHad2 = tcA*population[cd8map[0]][target2] 
						+ teA*population[cd8map[1]][target2]
                                                + cancA*population[canID][target2] + macA*population[macID][target2];
					if(tcount <= maxA - voxHad1 && tcount > 0){
						for(int cd=0; cd < 2; cd++){
							population[cd8map[cd]][target1] += cd8s[cd];
							cd8s[cd] = 0;
						}
						tcellPos.push_back(target1);
                        //tcellPos[tcCount] = target1;
                        //tcCount++;
						done = true;
                        tcount = 0;
						goto stop;
					}else if(tcount <= maxA - voxHad2 && tcount > 0){
                        for(int cd=0; cd < 2; cd++){
                          population[cd8map[cd]][target2] += cd8s[cd];
       					  cd8s[cd] = 0;
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
                              int try_out = int((tcA+teA)*rando);
                              if(cd8s[try_out] > 0){
                                population[cd8map[try_out]][target1]++;
                                cd8s[try_out]--;
                                left--;
                                tcount--;
                              }
                            }
							tcellPos.push_back(target1);
                            //tcellPos[tcCount] = target1;
                            //tcCount++;
						}
						if(maxA - voxHad2 > 0 && tcount > 0 && target1 != target2){
							int left = maxA - voxHad2;
							while(left > 0 && tcount > 0){
                              double rando = random->uniform(); //dist(mt);
                              int try_out = int((tcA+teA)*rando);
                              if(cd8s[try_out] > 0){
                                population[cd8map[try_out]][target2]++;
                                cd8s[try_out]--;
                                left--;
                                tcount--;
                              }
                            }
                            tcellPos.push_back(target2);
							//tcellPos[tcCount] = target2;
                            //tcCount++;
						}
						if(tcount <= 0){done=true; goto stop;}
						for(int h=0; h<2; h++){
							if(cd8s[h] < 0){
								population[outSigID][5] = 30000;
							}
						}
					}
					if(tcount <= 0){done=true; goto stop;}
				}
			}
			width += 2;
            if(width > 200){done=true; goto stop;} // stops the search for space for a new T cell when cancer prolfierates if there is no space left. The T cell is just removed. This is only when the cancer is winning to a large extent (at the end of the simulation
		}
	}else{population[outSigID][4] = 20000;way = 3;}
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
    for(int ch=0; ch<tcellPos.size(); ch++){
      if(tcellPos[ch] == 7659){
        tcellThere = true;
      }
    }
  }

// T cells disappearing when crossing boundary, WAS OFF

//teA = exhausted_T_area_count; // put teA back to whatever it was preset

if(which == 4  || which == 5){
	if(population[latEID][i] >= 1 && population[latEID][j] >= 1){
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
			int tIndex = which - 4; // decide which cell type is moving
			int cd8map[] = {cd8aID,cd8eID};
			int tType = cd8map[tIndex];
			
            //std::random_device rt;
       		//std::mt19937 mta(rt());
        	//std::uniform_real_distribution<double> distr(0, 1.0);
            double ranGen = random->uniform();//distr(mta);
			
			if(dissProb >= ranGen){ // THIS WAS MADE FALSE TO CHECK EFFECTS OF T CELL EGRESS, NEED TO CHANGE BACK FOR FULL ICS
                if(population[tType][j] > 0){
				  population[tType][j] -= 1;
                }else{
                  population[outSigID][456] = 43;
                  std::cout << "problem with boundary" << std::endl;
                }
                proliferation_update_bool = true;
			}/*else{ // this condition is to employ bouncy reflective boundary conditions by reversing the hopping action
              population[tType][j] -= 1;
              population[tType][i] += 1;
            }*/
		}
	}
}

if(which == 10 and false){ // tam motility -> to put into place the reflective boundaries
  if(population[latEID][i] >= 1 && population[latEID][j] >= 1){
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
          population[macID][j] -= 1;
          population[macID][i] += 1;
        }
    }
}

int recSpot;
if(which == 1 || which == 9){ // Ta and macrophage recruitment
  int weight;
  int recType;
  if(which == 1){ // T cell
    weight = tcA;
    recType = cd8aID;
    proliferation_update_bool = true;
  }else if(which == 9){ // macrophage
    weight = macA;
    recType = macID;
  }
  //std::random_device ch;
  //std::mt19937 mtd(ch());
  //std::uniform_real_distribution<double> spDist(0, 1.0);
  double ranSpot = random->uniform(); //spDist(mtd);
  
  int eCount = 0;
  for(int px=0; px < latCol*latCol; px++){
    int filling = tcA*population[cd8aID][px] + teA*population[cd8eID][px] + cancA*population[canID][px] + macA*population[macID][px];
    int left = int((maxA - filling)/weight);
    eCount += left;
  }
  ranSpot = (int)(eCount * ranSpot);
  eCount = 0;
  
  for(int px=0; px < latCol*latCol; px++){
    int filling = tcA*population[cd8aID][px] + teA*population[cd8eID][px] + cancA*population[canID][px] + macA*population[macID][px];
    int left = int((maxA - filling)/weight);
    if(left > 0 && eCount <= ranSpot && ranSpot <= eCount + left){
      population[recType][px] += 1;
      recSpot = px;
      break;
    }
    eCount += left;
  }
}

// Boundary producing random T cells (according to current distribution of T cells) WAS OFF
if(which == 8 && false){
	//std::random_device rdv;
    //std::mt19937 mtt(rdv());
    //std::uniform_real_distribution<double> distrib(0,1.0);
    double rand = random->uniform(); //distrib(mtt);
	
	int cd8map[] = {cd8aID, cd8eID}; // Activated, exhausted
	int cd8Count[] = {0,0};
	int totcd8 = 0;
	for(int ii=0; ii < latCol*latCol; ii++){
		for(int ci=0; ci < 2; ci++){
			cd8Count[ci] += population[cd8map[ci]][ii];
			totcd8 += population[cd8map[ci]][ii];
		}
	}

	int decision = round(totcd8*rand);
	int label = 0;
	int w_type = -1;
	for(int ii=0; ii < 2; ii++){
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
int wProblem[] = {0,0,0};
for(int ii=0; ii < latCol*latCol; ii++){
                if(population[cd8aID][ii] < 0){// || population[0][i] + population[2][i] + population[5][i] < 0){
                    //population[outSigID][ii] = 3;
			        problem = true;
			        wProblem[0] += 1;
			        nLocs.push_back(ii);
                    //std::cout << "negative " << << std::endl;
                }else if(population[cd8eID][ii] < 0){
                    //population[outSigID][ii] = 7;
			        problem = true;
			        wProblem[1] += 1;
			        nLocs.push_back(ii);
                    //std::cout << "negative";
                }else if(population[canID][ii] < 0){
                    //population[outSigID][ii] = 11;
			        problem = true;
			        wProblem[2] += 1;
			        nLocs.push_back(ii);
                    //std::cout << "negative";
		        }
}
nLocs.clear();

if(which == 2 || which == 6 || which == 8){
  proliferation_update_bool = true;
}

if(proliferation_update_bool){//which == 1 || which == 2 || which == 4 || which == 5 || which == 6 || which == 8){ // reset the total count of cd8+ T cells
  tcell_total_population = 0;
  for(int ic=0; ic < latCol*latCol; ic++){
    tcell_total_population += population[cd8aID][ic] + population[cd8eID][ic];
  }
  //proliferation_update_bool = true;
}

  rxn_count[which]++;

  // compute propensity changes for participating sites and first neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;
  bool hit1 = false;
  bool hit = false;

  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;
  
  
  if(i == 3406){
    hit1 = true;
  }

  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;

      if(isite == 3406){
        hit1 = true;
      }
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

        if(isite == 3406){
          hit1 = true;
        }
      }
    }
  }
// cancerPos tcellPos
  std::vector<int> cSites;
  //cSites.push_back(9433);
  if(which == 0 || which == 1 || which == 9 || which == 10 || which == 12 || which == 14 || which == 3){ // cancer prolif, Tc recruitment, mac recruitment, mac motility, mac death, local cancer death, nbr cancer death
    //std::vector<int> cSites;
    //cSites.push_back(9433);
    if (which == 0) { // for farther updates for cancer proliferation
      //std::cout << nsites  << " " << 762 << std::endl;
      //int cSites[5*widthLim];
      //int ncSites = 1;
      cSites.push_back(cancerPos[1]);
      //cSites[0] = cancerPos[1];
      for(int ts=0; ts < tcellPos.size(); ts++){
        cSites.push_back(tcellPos[ts]);  
        //cSites[ncSites++] = tcellPos[ts];
      }
      //removed the following and replaced with lattice-wide update for testing purposes
      
      for(int ci=0;ci<tcellPos.size()+1; ci++){
      //for(int ci=0;ci<tcCount+1; ci++){
        int cx;
        int cy;
        int width = 3;
        if(ci == 0){
          cx = cancerPos[1]%latCol;
          cy = cancerPos[1]/latCol;
        }else{
	      cx = tcellPos[ci - 1]%latCol;
          cy = tcellPos[ci - 1]/latCol;
        }
        for(int lo=0; lo < neighLim; lo++){ // updates around cancer cells to a certain depth
          for(int ly=0; ly < 2; ly++){ //lx and ly take us around the neighbors
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
	          cSites.push_back(target1);
	          cSites.push_back(target2);
	          //cSites[ncSites] = target1;
              //ncSites++;
              //cSites[ncSites] = target2;
              //ncSites++;
    	    }
          }
	      width += 2;
          if(ci > 0){
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
    }else if(which == 3 || which == 14){
      int cancPos = -1;
      int width = 3;
      if(which == 3){
        cancPos = j;
      }else if(which == 14){
        cancPos = i;
      }
      
      // pasted code below...
      int cx = cancPos % latCol;
      int cy = cancPos / latCol;
      for(int lo=0; lo < neighLim; lo++){ // updates around cancer cells to a certain depth
        for(int ly=0; ly < 2; ly++){ //lx and ly take us around the neighbors
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
            cSites.push_back(target1);
            cSites.push_back(target2);
          }
        }
        width += 2;
      }

      // above
      
    }else if(which == 1){
      //proliferation_update_bool = true;
      int recX = recSpot % latCol;
      int recY = recSpot / latCol;
      for(int tx=0; tx < 3; tx++){
        for(int ty=0; ty < 3; ty++){
          int targx = recX + tx - 1;
          int targy = recY + ty - 1;
          targx -= floor(double(targx)/latCol)*latCol;
          targy -= floor(double(targy)/latCol)*latCol;
          int targ = targx + latCol*targy;
          cSites.push_back(targ);
        }
      }
      //cSites.push_back(9433);
    }else if(which == 9 || which == 10 || which == 12){ // macrophages
      std::vector<int> mChax;
      std::vector<int> mChay;
      if(which == 12){ // death
        mChax.push_back(i%latCol);
        mChay.push_back(i/latCol);
      }else if(which == 10){ // diffusion
        mChax.push_back(i%latCol);
        mChay.push_back(i/latCol);
        mChax.push_back(j%latCol);
        mChay.push_back(j/latCol);
      }else if(which == 9){ // recruitment
        mChax.push_back(recSpot%latCol);
        mChay.push_back(recSpot/latCol);
      }
      for(int k=0;k < mChax.size(); k++){
        for(int lx=0; lx < widthLim; lx++){
          for(int ly=0; ly < widthLim; ly++){
            int targx = mChax[k] - neighLim + lx;
            int targy = mChay[k] - neighLim + ly;
            targx -= floor(double(targx)/latCol)*latCol;
            targy -= floor(double(targy)/latCol)*latCol;
            int targ = targx + latCol*targy;
            cSites.push_back(targ);
          }
        }
      }
      //cSites.push_back(9433);
    }
    
    //if(which == 9 || which == 1 || which == 2 || which == 3 || which == 4 || which == 5 || which == 6 || which == 7){
    if(which == 0){  
      //cSites.push_back(9433);
    }
    //std::cout << i  << " " << 804 << std::endl;
    //std::cout << ncSites << std::endl;
    //std::sort(cSites.begin(), cSites.end());
    //auto last = std::unique(cSites.begin(), cSites.end());
    //cSites.erase(last, cSites.end());
    
    sort(cSites.begin(), cSites.end());
    cSites.erase(unique(cSites.begin(), cSites.end()), cSites.end());
    
    for(int i=0; i < cSites.size(); i++){
      if(cSites[i] == 3406){
        hit1 = true;
      }
    }
    if(!proliferation_update_bool){
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
 
int max_filling = 4;
  if(proliferation_update_bool){
    //cSites.clear();
    for(int la=0; la < latCol*latCol; la++){
      if(population[cd8aID][la] > 0 || population[recID][la] > 0){
        cSites.push_back(la);
        if(la == 3406){
          hit1 = true;
        }
      }
    }
    proliferation_update_bool = false;
    
    // gung-ho updates that always trigger for all cells in the lattice  
    /*int IDsToUse[] = {canID, cd8aID, cd8eID, macID};
    for(int hw=0; hw < latCol*latCol; hw++){
      for(int sp=0; sp < 4; sp++){
        if(population[IDsToUse[sp]][hw] > 0){
          cSites.push_back(hw);
          break;
        }
      }
    }*/

  
    //std::sort(cSites.begin(), cSites.end());
    //auto last = std::unique(cSites.begin(), cSites.end());
    //cSites.erase(last, cSites.end());
    
    sort(cSites.begin(), cSites.end());
    cSites.erase(unique(cSites.begin(), cSites.end()), cSites.end());
    
    for (n = 0; n < cSites.size(); n++) {
      m = cSites[n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
    }
    //cSites.clear();
  }
  //bool hit =  false;
  for(int i=0; i < cSites.size(); i++){
    if(cSites[i] == 3406){
      hit = true;
    }
  }
  cSites.clear();
  tcellPos.clear();
  
  //if((which != 1 && which != 9) && )
  

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

