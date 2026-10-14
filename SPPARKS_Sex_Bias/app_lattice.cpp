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

#include "spktype.h"
#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "solve.h"
#include "domain.h"
#include "lattice.h"
#include "random_mars.h"
#include "random_park.h"
#include "cluster.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

//Giuseppe imports
#include <iostream>
#include <algorithm>
#include <numeric>
#include <array>


using namespace SPPARKS_NS;

#define DELTA 32768

enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};

/* ---------------------------------------------------------------------- */

AppLattice::AppLattice(SPPARKS *spk, int narg, char **arg) : App(spk,narg,arg) {
    appclass = LATTICE;

    MPI_Comm_rank(world, &me);
    MPI_Comm_size(world, &nprocs);

    // default settings

    sectorflag = 0;
    nset = 0;
    set = NULL;
    nstop = 1.0;
    tstop = 0.0;

    sweepflag = NOSWEEP;
    ranapp = NULL;
    ranstrict = NULL;
    siteseeds = NULL;
    sitelist = NULL;
    Lmask = false;
    mask = NULL;

    allow_app_update = 0;

    temperature = 0.0;

    propensity = NULL;
    i2site = NULL;

    comm = NULL;
    sweep = NULL;

    nlocal = nghost = nmax = 0;
    owner = NULL;
    index = NULL;

    maxneigh = 0;
    numneigh = NULL;
    neighbor = NULL;

    dt_sweep = 0.0;
    naccept = nattempt = 0;
    nsweeps = 0;


    app_update_only = 0;

//  dt_diff = 6;   //Once an hour, update dead_mem		//Intialize Giuseppe's variables
//  diffCount = 0;
//  twindow = atof(arg[11]); // was 24 for a day in time steps of dtdiff, a day of lag
//  nDeadCounti = 0;
    proliferation_update_bool = false; // initiates update of all cd8+ T cell proliferation rates
    cancer_proliferation_update_bool = false;

    mapped_population_to_ifng_tracker = false;

    esites = new int[1000000]; // enough for 360000 sites or 600x600
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

//    std::vector<std::vector<int>> receptor_map(3, std::vector<int>(3)) = {{TpID, TmID, TeID}, {ARTpID, ARTmID, ARTeID}, {canID}};
//    std::array<std::array<int, 3>, 3> receptor_map;
//    int receptor_map[3][3] = {{TpID, TmID, TeID}, {ARTpID, ARTmID, ARTeID}, {canID, -1, -1}};

    /*tcell_population_total = 0;
    for(int ic=0; ic < latCol*latCol; ic++){
      tcell_population_total += population[cd8aID][ic] + population[cd8eID][ic];
    }*/

    //initialize ifng_field from double pointer to array?

    latCol = atoi(arg[1]); // get lattice size (in sites)
    ifng_diffusion = atof(arg[10]); // in um2/s
    ifng_time_step = atof(arg[11]); // in units of 10min
    ifng_production = atof(arg[12]); // rate of ifng production in AR- cells
    ifng_decay = atof(arg[13]); // rate of ifng decay
    ifng_ART_uptake = atof(arg[23]); // rate of ifng consumption in AR+ T cells
//    ifng_saturation_threshold = atof(arg[14]); // threshold at which ifng uptake is maximal
//    ifng_feedback_saturation_threshold = atof(arg[21]); // threshold at which ifng feedback is maximal

    ifng_feedback_hill_degree = atoi(arg[27]);  // degree of ifng feedback hill function
    max_ifng_production = atof(arg[28]); // maximum ifng feedback factor (was 2)
    half_max_ifng_feedback = atof(arg[21]); // ifng at which feedback augmentation is half max

    ifng_cancer_proliferation_limitation = atoi(arg[15]) == 1; // ifng limits cancer proliferation
    ifng_cytotoxicity_augmentation = atoi(arg[16]) == 1; // ifng augments T cell lysis rates

    receptors_per_cell = new int[3]{atoi(arg[29]), atoi(arg[30]), atoi(arg[31])}; // AR- CD8, AR+ CD8, cancer

    periodic_boundary_conditions = atoi(arg[19]) == 1; // on/off ifng periodic boundary conditions

    bIFNG_decay = atof(arg[14]);  // decay rate of bound IFNG s-1

    use_simple_decay_model = atoi(arg[32]) == 1;
    if(!use_simple_decay_model) {
      use_leaky_receptor_model = true;
    }

    ifng_binding = atof(arg[34]);
    ifng_unbinding = atof(arg[35]);
    dissC = ifng_unbinding/ifng_binding;

    use_cell_dependent_decay_model = false;
    if(use_simple_decay_model) {
      use_cell_dependent_decay_model = true;
      arn_cd8_consumption_rate = receptors_per_cell[0]; // taking the receptor count as input for degradation rates
      arp_cd8_consumption_rate = receptors_per_cell[1];
      cancer_consumption_rate = receptors_per_cell[2];
    }
    update_bIFNG_decay = bIFNG_decay > 0;

    ifng_time_keeper = 0.0;

    ifng_col = ceil((double)latCol/3); //200;
    ifng_row = ceil((double)latCol/3); //

    latt_const = 10; // 10 microns
    ifng_lattice_cell_converter = 3; //double(latCol/ifng_col);
    ifng_latt_const = int(ifng_lattice_cell_converter*latt_const);

    //initialize the ifng_field
    ifng_field = new double*[ifng_col];
    ifng_field_copy = new double*[ifng_col];
//    src_field = new int**[ifng_col];
//    bound_ifng = new double*[ifng_col];
//    receptors = new int*[ifng_col];
    ifng_tracker = new std::vector<std::vector <double>>**[latCol];
    bound_ifng = new std::array<double, 3>*[latCol];

//////////////////////////////////
//    for (int x = 0; x < latCol; x++) {
//        ifng_tracker[x] = new std::vector<std::vector <double>>*[latCol];
//        for (int y = 0; y < latCol; y++) {
//            int location = x+latCol*y;
//            ifng_tracker[x][y] = new std::vector<std::vector <double>>[7]; // for 7 species to track
//            for(int id=0; id<7; id++){
////                ifng_tracker[x][y][id] = new std::vector<std::vector <double>>;
//                for(int p=0; p < population[cellIDs[id]][location]; p++){
//                    ifng_tracker[x][y][id].push_back(std::vector <double>{0, (double) receptors_per_cell[receptors_per_cell_id[cellIDs[id]]]});
//                }
//            }
//        }
//    }
////////////////////////////////////
//    for (int i = 0; i < ifng_col; i++) ifng_field[i] = new double[ifng_row];


//    for (int i = 0; i < ifng_col; i++) ifng_field_copy[i] = new double[ifng_row];

//    for (int i = 0; i < ifng_col; i++) src_field[i] = new int*[ifng_row];

    for (int x = 0; x < ifng_col; x++){  // initialize fields
        ifng_field[x] = new double[ifng_row];
        ifng_field_copy[x] = new double[ifng_row];

//        src_field[x] = new int*[ifng_row];
//        bound_ifng[x] = new double[ifng_row];
//        receptors[x] = new int[ifng_row];
        for (int y = 0; y < ifng_row; y++) {
//            double ifng_value = population[ifngID][int(ifng_lattice_cell_converter * (y * latCol + x))];
            ifng_field[x][y] = 0;//0; // set ifng field to that given by init populations
            ifng_field_copy[x][y] = 0; //0;
        }
//            src_field[x][y] = new int[3]{0, 0, 0};
//            bound_ifng[x][y] = 0;
//            receptors[x][y] = 0;
//
//            for (int ar = 0; ar < 2; ar++) { // ar+ or not (ar-, ar+) # DON'T NEED TO COUNT CANCER CELLS ANYMORE
//                if (ar == 2) {
//                    src_field[int(x / ifng_lattice_cell_converter)][int(
//                            y / ifng_lattice_cell_converter)][ar] += population[canID][x + (latCol * y)]; // add cancer cells as ar+ (ifng taking but not producing)
//                } else {
//                    for (int t = 0; t < 3; t++) { // cell types (p, m, e)
//                        src_field[int(x / ifng_lattice_cell_converter)][int(y / ifng_lattice_cell_converter)][ar] += population[cd8IDs[ar][t]][x + (latCol * y)]; // src field becomes [x][y][ar-,ar+,caner]
//                    }
//                }
//            }
//
//            for(int lx=0; lx < ifng_lattice_cell_converter; lx++){ // find all receptors
//                for(int ly=0; ly < ifng_lattice_cell_converter; ly++) {
//                    for(int p=0; p < 3; p++){
//                        for(int p2=0; p<3; p++){ // through all species in each category
//                            int id = receptor_map[p][p2];
//                            if(id >=0){
//                                int location = ((ifng_lattice_cell_converter*x+lx)+latCol*(ifng_lattice_cell_converter*y+ly));
//                                receptors[x][y] += receptors_per_cell[p]*population[id][location];
//                            }
//                        }
//                    }
//                }
//            }
//        }

    }


//    std::array<std::array<double, ifng_row>, ifng_col> ifng_field_copy;


//    int cd8IDs[2][3] = {{TpID, TmID, TeID}, {ARTpID, ARTmID, ARTeID}};
//    //initialize source field by finding cd8s and cancer cells
//    for(int x=0; x < latCol; x++){ // population columns
//        for(int y=0; y < latCol; y++){ // population rows
//            for(int ar=0; ar < 3; ar++){ // ar+ or not (ar-, ar+)
//                if(ar == 2){
//                    src_field[int(x / 10)][int(y / 10)][1] += population[canID][x + (latCol * y)]; // add cancer cells as ar+ (ifng taking but not producing)
//                }else {
//                    for (int t = 0; t < 3; t++) { // cell types (p, m, e)
//                        src_field[int(x / 10)][int(y / 10)][ar] += population[cd8IDs[ar][t]][x + (latCol * y)]; // latCol is maybe not defined?
//                    }
//                }
//            }
//        }
//    }

}

/* ---------------------------------------------------------------------- */

AppLattice::~AppLattice()
{
  Solve *s;
  for (int i = 0; i < nset; i++) {
    s = free_set(i);
    delete s;
  }
  delete [] set;

  delete ranapp;
  delete ranstrict;
  memory->destroy(siteseeds);
  memory->destroy(sitelist);
  memory->destroy(mask);

  delete comm;

  memory->destroy(owner);
  memory->destroy(index);

  memory->destroy(numneigh);
  memory->destroy(neighbor);
}

/* ---------------------------------------------------------------------- */

void AppLattice::input(char *command, int narg, char **arg)
{
  if (strcmp(command,"sector") == 0) set_sector(narg,arg);
  else if (strcmp(command,"sweep") == 0) set_sweep(narg,arg);
  else if (strcmp(command,"temperature") == 0) set_temperature(narg,arg);
  else if (strcmp(command,"app_update_only") == 0) set_app_update_only(narg,arg);
  else input_app(command,narg,arg);
}

/* ---------------------------------------------------------------------- */

void AppLattice::init()
{
  // error checks

  if (solve == NULL && sweepflag == NOSWEEP)
    error->all(FLERR,"App needs a KMC or rejection KMC solver");
  if (solve && sweepflag != NOSWEEP)
    error->all(FLERR,"App cannot use both a KMC and rejection KMC solver");

  if (solve && allow_kmc == 0)
    error->all(FLERR,"KMC events are not implemented in app");
  if (sweepflag != NOSWEEP && allow_rejection == 0)
    error->all(FLERR,"Rejection events are not implemented in app");
  if (sweepflag != NOSWEEP && Lmask && allow_masking == 0)
    error->all(FLERR,"Mask logic not implemented in app");

  if (sweepflag != NOSWEEP && Lmask && allow_masking == 0)
    error->all(FLERR,"Mask logic not implemented in app");

  if (nprocs > 1 && sectorflag == 0 && solve)
    error->all(FLERR,"Cannot use KMC solver in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RANDOM)
    error->all(FLERR,"Cannot use random rejection KMC in parallel with no sectors");
  if (nprocs > 1 && sectorflag == 0 && sweepflag == RASTER)
    error->all(FLERR,"Cannot use raster rejection KMC in parallel with no sectors");
  if (sectorflag && sweepflag == COLOR_STRICT)
    error->all(FLERR,"Cannot use color/strict rejection KMC with sectors");

  // if sectors, set number of sectors

  int dimension = domain->dimension;

  nsector = 1;
  if (sectorflag) {
    if (nsector_user) nsector = nsector_user;
    else if (dimension == 1) nsector = 2;
    else if (dimension == 2) nsector = 4;
    else nsector = 8;

    if (dimension == 3) {
      if (nsector == 2 && (domain->procgrid[1] != 1 || 
			   domain->procgrid[2] != 1))
	error->all(FLERR,"Invalid number of sectors");
      if (nsector == 4 && domain->procgrid[2] != 1)
	error->all(FLERR,"Invalid number of sectors");
    }
    if (dimension == 2) {
      if (nsector == 2 && domain->procgrid[1] != 1)
	error->all(FLERR,"Invalid number of sectors");
      if (nsector == 8)
	error->all(FLERR,"Invalid number of sectors");
    }
    if (dimension == 1 && nsector != 2)
      error->all(FLERR,"Invalid number of sectors");
  }

  // if coloring, determine number of colors
  // setup test for create_set
  // check periodicity against lattice extent

  ncolors = 1;
  if (sweepflag == COLOR || sweepflag == COLOR_STRICT) {
    int delcolor = delevent + delpropensity;
    if (domain->lattice == NULL)
      error->all(FLERR,"Cannot color without a lattice definition of sites");
    if (contiguous_sites() == 0)
      error->all(FLERR,"Cannot color without contiguous site IDs");
    ncolors = domain->lattice->ncolors(delcolor);
    if (ncolors == 0)
      error->all(FLERR,"Cannot color this combination of lattice and app");
  }

  if (nsector > 1 && ncolors > 1) bothflag = 1;
  else bothflag = 0;

  // create sets based on sectors and coloring
  // set are either all sectors or all colors or both
  // for both, first sets are entire sections, remaining are colors in sectors
  // if new nset is same as old nset, pass each set's solver to create_set,
  //   so it can reuse solver and its RNG,
  //   so consecutive runs are identical to one continuous run

  int nsetold = nset;
  Solve **sold = new Solve*[nsetold];

  for (int i = 0; i < nset; i++) sold[i] = free_set(i);
  delete [] set;

  if (nsector == 1 && ncolors == 1) {
    nset = 1;
    set = new Set[nset];
    if (nset == nsetold) create_set(0,0,0,sold[0]);
    else create_set(0,0,0,NULL);
  } else if (nsector > 1 && ncolors == 1) {
    nset = nsector;
    set = new Set[nset];
    for (int i = 0; i < nset; i++) {
      if (nset == nsetold) create_set(i,i+1,0,sold[i]);
      else create_set(i,i+1,0,NULL);
    }
  } else if (ncolors > 1 && nsector == 1) {
    nset = ncolors;
    set = new Set[nset];
    for (int i = 0; i < nset; i++) {
      if (nset == nsetold) create_set(i,0,i+1,sold[i]);
      else create_set(i,0,i+1,NULL);
    }
  } else if (bothflag) {
    nset = nsector + ncolors*nsector;
    set = new Set[nset];
    int m = 0;
    for (int i = 0; i < nsector; i++) {
      if (nset == nsetold) create_set(m,i+1,0,sold[m]);
      else create_set(m,i+1,0,NULL);
      m++;
    }
    for (int i = 0; i < nsector; i++) 
      for (int j = 0; j < ncolors; j++) {
	if (nset == nsetold) create_set(m,i+1,j+1,sold[m]);
        else create_set(m,i+1,j+1,NULL);
	m++;
      }
  }

  if (nset != nsetold)
    for (int i = 0; i < nsetold; i++) delete sold[i];
  delete [] sold;

  // initialize mask array

  if (!Lmask && mask) {
    memory->destroy(mask);
    mask = NULL;
  }
  if (Lmask && mask == NULL) {
    memory->create(mask,nlocal+nghost,"app:mask");
    for (int i = 0; i < nlocal+nghost; i++) mask[i] = 0;
  }

  // setup RN generators, only on first init
  // ranapp is used for all options except sweep color/strict
  // setup ranapp so different on every proc
  // if color/strict, initialize per-lattice site seeds

  if (ranapp == NULL) {
    ranapp = new RandomPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    ranapp->reset(seed,me,100);
  }

  if (sweepflag != COLOR_STRICT) {
    delete ranstrict;
    memory->destroy(siteseeds);
    ranstrict = NULL;
    siteseeds = NULL;
  }

  if (sweepflag == COLOR_STRICT && ranstrict == NULL) {
    ranstrict = new RandomPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    memory->create(siteseeds,nlocal,"app:siteseeds");
    for (int i = 0; i < nlocal; i++) {
      ranstrict->tagreset(seed,id[i],100);
      siteseeds[i] = ranstrict->seed;
    }
  }

  // initialize comm, both for this proc's full domain and sectors
  // recall comm->init in case sectoring has changed

  if (comm == NULL) comm = new CommLattice(spk);
  comm->init(nsector,delpropensity,delevent,NULL);

  // set sweep function ptr

  if (sweepflag != NOSWEEP) {
    if (sweepflag != COLOR_STRICT && !Lmask)
      sweep = &AppLattice::sweep_nomask_nostrict;
    else if (sweepflag != COLOR_STRICT && Lmask)
      sweep = &AppLattice::sweep_mask_nostrict;
    else if (sweepflag == COLOR_STRICT && !Lmask)
      sweep = &AppLattice::sweep_nomask_strict;
    else if (sweepflag == COLOR_STRICT && Lmask)
      sweep = &AppLattice::sweep_mask_strict;
  } else sweep = NULL;

  // app-specific initialization, after general initialization

  init_app();

  // error checks that cannot be done until after init_app()

  if (sweepflag && dt_sweep == 0.0) 
    error->all(FLERR,"App did not set dt_sweep");

  // initialize output

  output->init(time);
}

/* ---------------------------------------------------------------------- */

void AppLattice::setup()
{
  // app-specific setup, before propensities are computed

  setup_app();

  // initialize propensities for KMC solver within each set
  // comm insures ghost sites are up to date

  if (solve) {
    comm->all();
    for (int i = 0; i < nset; i++) {
      for (int m = 0; m < set[i].nlocal; m++)
	set[i].propensity[m] = site_propensity(set[i].site2i[m]);
      set[i].solve->init(set[i].nlocal,set[i].propensity);
    }
  }

  // convert per-sector time increment info to KMC params

  if (sectorflag && solve) {
    if (tstop > 0.0) {
      Ladapt = false;
      dt_kmc = tstop;
    }

    if (nstop > 0.0) {
      Ladapt = true;
      double pmax = 0.0;
      for (int i = 0; i < nset; i++) {
	int ntmp = set[i].solve->get_num_active();
	if (ntmp > 0) {
	  double ptmp = set[i].solve->get_total_propensity();
	  ptmp /= ntmp;
	  pmax = MAX(ptmp,pmax);
	}
      }
      double pmaxall;
      MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
      if (pmaxall > 0.0) dt_kmc = nstop/pmaxall;
      else dt_kmc = stoptime-time;
    }

    dt_kmc = MIN(dt_kmc,stoptime-time);
  }

  // convert rejection info to rKMC params
  // nloop and nselect are set whether sectoring is used or not
  // if bothflag, list of active sets starts with nsector

  if (sweepflag != NOSWEEP) {
    int first = 0;
    if (bothflag) first = nsector;

    if (sweepflag == RANDOM) {
      if (nstop > 0.0) {
	for (int i = first; i < nset; i++) {
	  set[i].nloop = 0;
	  set[i].nselect = static_cast<int> (nstop*set[i].nlocal);
	}
      }
      if (tstop > 0.0) {
	double n = tstop / (dt_sweep/nglobal);
	for (int i = first; i < nset; i++) {
	  set[i].nloop = 0;
	  set[i].nselect = static_cast<int> (n/nglobal * set[i].nlocal);
	}
      }

    } else if (sweepflag == RASTER || 
	       sweepflag == COLOR || sweepflag == COLOR_STRICT) {
      int n;
      if (nstop > 0.0) n = static_cast<int> (nstop);
      if (tstop > 0.0) n = static_cast<int> (tstop/dt_sweep);
      for (int i = first; i < nset; i++) {
	set[i].nloop = n;
	set[i].nselect = n * set[i].nlocal;
      }
    }

    double nme = 0.0;
    for (int i = first; i < nset; i++) nme += set[i].nselect;
    double ntotal;
    MPI_Allreduce(&nme,&ntotal,1,MPI_DOUBLE,MPI_SUM,world);

    dt_rkmc = ntotal/nglobal * dt_sweep;
    if (dt_rkmc == 0.0)
      error->all(FLERR,"Choice of sector stop led to no rKMC events");
    dt_rkmc = MIN(dt_rkmc,stoptime-time);
  }

  // setup sitelist if sweepflag = RANDOM
  // do this every run since sector timestep could have changed

  if (sweepflag == RANDOM) {
    memory->destroy(sitelist);
    int n = 0;
    for (int i = 0; i < nset; i++) n = MAX(n,set[i].nselect);
    memory->create(sitelist,n,"app:sitelist");
  }

  // second stage of app-specific setup

  setup_end_app();

  // setup future output

  nextoutput = output->setup(time,first_run);
}

/* ---------------------------------------------------------------------- */

void AppLattice::iterate()
{
  timer->barrier_start(TIME_LOOP);

  if (solve) {
    if (sectorflag == 0) 
      iterate_kmc_global(stoptime);
    else if (allow_app_update && app_update_only)
      iterate_app_update_only(stoptime,dt_kmc);
    else
      iterate_kmc_sector(stoptime);
  } else {
    if (allow_app_update && app_update_only)
      iterate_app_update_only(stoptime,dt_rkmc);
    else
      iterate_rejection(stoptime);
  }

  timer->barrier_stop(TIME_LOOP);
}

/* ----------------------------------------------------------------------
   KMC solver on entire domain
   can only be invoked in serial
 ------------------------------------------------------------------------- */

void AppLattice::iterate_kmc_global(double stoptime)
{
  int isite, i;
  
  // global KMC runs with one set
  // save ptr to system solver

  Solve *hold_solve = solve; // solve is the string taken from the first argument after "app_style" in the input file (in.const...), for us this is rxn/diff/custom
  solve = set[0].solve; // this decides where the
  propensity = set[0].propensity;
  i2site = set[0].i2site;

  int done = 0;
  int ifng_update_count = 0;
  while (!done) {

      timer->stamp();
      isite = solve->event(&dt_step); // go to solve_tree.cpp function event()
      timer->stamp(TIME_SOLVE);

      ///////////////////////// Giuseppe IFNG Code ////////////////////////

	  // if the calculated time to the next gillespie event is longer than the time to the next IFNG step
      if(dt_step + ifng_time_keeper > ifng_time_step){
          if(time <= stoptime) { // if the simulation is not over (maximum time of simulation)
              AppLattice::update_ifng_field(ifng_time_step); // update the ifng field
              proliferation_update_bool = ifng_cytotoxicity_augmentation; // turn on propensity checks for T cells etc.
              cancer_proliferation_update_bool = ifng_cancer_proliferation_limitation; // turn on propensity checks for cancer cells etc.
              time += ifng_time_step-ifng_time_keeper; // update time
              ifng_time_keeper = 0.0; // zero ifng time keeper
              ifng_update_count += 1;

              // now update the site propensities (do all sites for now)
//          std::vector<int> update_sites;
              int number_of_sites = 0;
              for (int n = 0; n < latCol * latCol; n++) {
                  isite = i2site[n];
                  if (isite >= 0) {// && echeck[isite] == 0) {
                      propensity[isite] = site_propensity(n);
                      esites[number_of_sites++] = isite;
//                  echeck[isite] = 1;
//                  esites[nsites++] = isite;
//                  echeck[isite] = 1;
                  }
              }
              solve->update(number_of_sites, esites, propensity);
              timer->stamp(TIME_APP);
          }else {
              done = 1;
              time = stoptime;
          }
      }else if(isite >= 0){
          if (time <= stoptime) {
              ifng_time_keeper += dt_step;
              time += dt_step;
              site_event(isite,ranapp);
              naccept++;
              timer->stamp(TIME_APP);
          }else {
              done = 1;
              time = stoptime;
          }
      }else{
          done = 1;
          time = stoptime;
      }

    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }

  // restore system solver

  solve = hold_solve;
}

/* ----------------------------------------------------------------------
   a function which takes values and returns the
 ------------------------------------------------------------------------- */
void AppLattice::update_ifng_field(double time_step){
    int canID = 0;
    int TpID = 1;
    int TmID = 2;
    int TeID = 3;
    int ARTpID = 4;
    int ARTmID = 5;
    int ARTeID = 6;
    int ifngID = 0;
    int recID = 7;
    int deadID = 9;
    int ThID = 12; // All AR- CD4+ T cells // formerly helper CD4+ T cells
    int ARThID = 13; // All AR+ CD4+ T cells // formerly AR+ helper CD4+ T cells
    int TrID = 14; // regulatory T cells
    int ARTrID = 15; // AR+ regulatory T cells
    int boundID = 11;
    int boundifngID = 1;
    double time_conversion_factor = 600; // x600s

    double IFNGR_decay_rate = std::exp(-bIFNG_decay*time_conversion_factor*ifng_time_step); // gives the exponential dropoff of the bound IFNGR

//    double too_long_step = 0.15; // used to track errors
//
//    if(ifng_time_step > too_long_step){
//        population[deadID][31] += 1; // deadID output the counts where time went past what I need it to be
//    }
    int cellIDs[11] = {canID, TpID, TmID, TeID, ARTpID, ARTmID, ARTeID, ThID, ARThID, TrID, ARTrID};
    int trackIDs[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int srcIDs[7] = {TpID, TmID, TeID};

    int receptor_map[3][3] = {{TpID, TmID, TeID}, {ARTpID, ARTmID, ARTeID}, {canID, -1, -1}};
    int cd8IDs[2][4] = {{TpID, TmID, TeID, ThID}, {ARTpID, ARTmID, ARTeID, ARThID}};

// SRC_FIELD UPDATE IS TO BE REMOVED UPON COMPLETION OF BINDING/UNBINDING ADDITION /////////////////////////////////
// update source field by finding cd8s
//    for(int x=0; x < latCol; x++){  // population columns
//        for(int y=0; y < latCol; y++){  // population rows
//            if(population[boundID][int(x+latCol*y)] > 0) {
//                for (int ar = 0; ar < 2; ar++) {  // ar+ or not (ar-, ar+) # DON'T NEED TO COUNT CANCER CELLS ANYMORE
//                    if (ar == 2) {
//                        src_field[int(x / ifng_lattice_cell_converter)][int(
//                                y / ifng_lattice_cell_converter)][ar] += population[canID][x + (latCol * y)]; // add cancer cells as ar+ (ifng taking but not producing)
//                    } else {
//                        for (int t = 0; t < 3; t++) { // cell types (p, m, e)
//                            src_field[int(x / ifng_lattice_cell_converter)][int(y / ifng_lattice_cell_converter)][ar] += population[cd8IDs[ar][t]][x + (latCol * y)]; // src field becomes [x][y][ar-,ar+,caner]
//                        }
//                    }
//                }
//                for(int p=0; p < 3; p++){ // through each category (AR-, AR+, cancer)
//                    for(int p2=0; p<3; p++){ // through all species in each category
//                        int id = receptor_map[p][p2];
//                        if(id >=0){
//                            int location = x + (latCol * y);
//                            receptors[int(x / ifng_lattice_cell_converter)][int(
//                                    y / ifng_lattice_cell_converter)] += receptors_per_cell[p]*population[id][location];
//                        }
//                    }
//                }
//            }
//        }
//    }
    // SRC_FIELD UPDATE IS TO BE REMOVED UPON COMPLETION OF BINDING/UNBINDING ADDITION /////////////////////////////////

    for(int x=0; x < ifng_col; x++) {
        for (int y = 0; y < ifng_row; y++) {
            ifng_field_copy[x][y] = ifng_field[x][y];
        }
    }

//    ifng_field_copy = ifng_field;
    for(int x=0; x < ifng_col; x++){ // all ifng pixels
        for(int y=0; y < ifng_row; y++){
            // find neighbors
            if(population[boundID][int(x * ifng_lattice_cell_converter + (latCol * ifng_lattice_cell_converter * y))] == 1) { // boundaryID- check if IFNG is able to diffuse in this region (not filling)
              double neighbors_list[] = {0.0, 0.0, 0.0, 0.0};
              int ncounter = 0;
              for (int dx = -1; dx < 3; dx += 2) { // through all 4 neighboring ifng pixels
                for (int dy = 0; dy < 2; dy++) {
                  int nx = x + dx * dy;
                  int ny = y + dx * (1 - dy);

                  if (periodic_boundary_conditions) {
                    nx -= floor(double(nx) / ifng_col) * ifng_col;
                    ny -= floor(double(ny) / ifng_row) * ifng_row; // periodic boundary conditions
                    neighbors_list[ncounter] = ifng_field_copy[nx][ny];
                  } else {
                    if (nx < 0 || nx >= ifng_col || ny < 0 || ny >= ifng_row || population[boundID][int(nx * ifng_lattice_cell_converter + (latCol * ifng_lattice_cell_converter * ny))] == 0) {
                      neighbors_list[ncounter] = 0;
                    } else {
                      neighbors_list[ncounter] = ifng_field_copy[nx][ny];
                    }
                  }
                  ncounter++;
                }
              }
              double neighbor_contribution = accumulate(std::begin(neighbors_list), std::end(neighbors_list), 0.0,std::plus<double>()); // may need to switch to a different method which would be fine
              double laplacian = (neighbor_contribution - 4 * ifng_field_copy[x][y]) / (ifng_latt_const * ifng_latt_const); // using pixel size for ifng field

              // build source term
              //                double ifng_saturation_term = ifng_field_copy[x][y] / ifng_saturation_threshold;
              //                if (ifng_saturation_term > 1) {
              //                    ifng_saturation_term = 1;
              //                }

              //double max_ifng_augmentation = 2;
              //                double ifng_feedback_term = ((max_ifng_production-0.5)*pow(ifng_field_copy[x][y], ifng_feedback_hill_degree)/(pow(half_max_ifng_feedback, ifng_feedback_hill_degree)+pow(ifng_field_copy[x][y], ifng_feedback_hill_degree))) + 0.5; //(ifng_field_copy[x][y] / ifng_feedback_saturation_threshold) + 0.5;

              //                double ifng_reading = bound_ifng[x][y]/(double)receptors[x][y];
              //                double ifng_feedback_term = ((max_ifng_production-0.5)*pow(ifng_reading, ifng_feedback_hill_degree)/(pow(half_max_ifng_feedback, ifng_feedback_hill_degree)+pow(ifng_reading, ifng_feedback_hill_degree))) + 0.5; //(ifng_field_copy[x][y] / ifng_feedback_saturation_threshold) + 0.5;
              //                if (ifng_feedback_term > max_ifng_augmentation) {
              //                    ifng_feedback_term = max_ifng_augmentation;
              //                }

              double total_production = 0; // number of sources * ifng feedback factor
              double total_consumption = 0; // number of sinks * (their consumption rates?)
              double total_receptors = 0;
              double total_bound = 0;
              int total_free_receptors = 0;
              // go to each site and calculate a) IFNG production and b) laplacian for diffusion then update site concentration
              for(int xx=0; xx < ifng_lattice_cell_converter; xx++){ // find sources in local ifng lattice pixel
                for(int yy=0; yy < ifng_lattice_cell_converter; yy++){ // go through IFNG pixel in cell lattice coordinates
                  int latx = x*ifng_lattice_cell_converter+xx;
                  int laty = y*ifng_lattice_cell_converter+yy;

                  if(0 <= latx && latx < latCol && 0 <= laty && laty < latCol) {
                    bound_ifng[latx][laty][0] = 0; // set local bound fraction
                    bound_ifng[latx][laty][1] = 0; // set local number of bound receptors to 0
                    bound_ifng[latx][laty][2] = 0; // set local number of total receptors to 0
                    for(int arn=0; arn < 11; arn++) {  // go through each cell type to find total IFNG production
                      int current_cellID = cellIDs[arn]; // used by population but not ifng_tracker, remember!
                      if(use_simple_decay_model) {
                        if(!use_cell_dependent_decay_model) {
                          if(std::find(std::begin(cd8IDs[0]), std::end(cd8IDs[0]), current_cellID) != std::end(cd8IDs[0])) { // if AR- CD8
                            double ifng_feedback_term = ((max_ifng_production - 0.5) * pow(ifng_field_copy[x][y], ifng_feedback_hill_degree) / (pow(half_max_ifng_feedback, ifng_feedback_hill_degree) + pow(ifng_field_copy[x][y], ifng_feedback_hill_degree))) + 0.5;
                            total_production += ifng_production * population[current_cellID][latx+latCol*laty] * ifng_feedback_term;
                          }
                        } else { // find total production as well as total decay/consumption
                          for(auto & cell : ifng_tracker[latx][laty][current_cellID]){ // go through each arn cell type
                            double total_consumed = cell[0];
                            double cell_consumption_rate = cell[1];
                            total_consumption += cell_consumption_rate;
                            cell[0] = cell_consumption_rate*ifng_field_copy[x][y]; //record consumption at this point ()

                            double ifng_feedback_term = ((max_ifng_production - 0.5) * pow(ifng_field_copy[x][y], ifng_feedback_hill_degree) / (pow(half_max_ifng_feedback, ifng_feedback_hill_degree) + pow(ifng_field_copy[x][y], ifng_feedback_hill_degree))) + 0.5;

                            if(1 <= arn && arn <= 2) {  // each AR- CD8+ T cell type
                              total_production += ifng_production * ifng_feedback_term;
                            }
                          }
                        }
                      }else{ // if not using simple decay model
                        for(auto & cell : ifng_tracker[latx][laty][arn]){ // go through each arn cell type, find total local production of IFNG

                          if(std::find(std::begin(cd8IDs[1]), std::end(cd8IDs[1]), current_cellID) != std::end(cd8IDs[1])) { // if AR+
                            double number_bound = cell[0]; // endo/exocytosis decays bound IFNGR
                            cell[0] = number_bound*IFNGR_decay_rate;
                          }

                          double cell_bound = cell[0];
                          double cell_receptors = cell[1];
                          total_receptors += cell_receptors;
                          total_bound += cell_bound;
                          total_free_receptors += cell_receptors - cell_bound;
                          double bound_reading = cell_bound/cell_receptors;
						  // I could adjust IFNG production rate based on how long each cell has been in the IFNG pixel
					      // but a) cells which move in 10x10um2 cell pixel will not often move in and out of 30x30um2
					      // IFNG pixel between IFNG updates and when cell are recruited, they will likely not make much
                          // more IFNG than is expected. This is a small error which fixing would be a pain in the ass. Fuck it.
                          if(1 <= arn && arn <= 2 || arn == 7) {  // each AR- CD8+ T cell type (no more Te IFNG production) or AR- CD4+ T cell
                            double ifng_feedback_term = ((max_ifng_production - 0.5) *
                                                         pow(bound_reading, ifng_feedback_hill_degree) /
                                                         (pow(half_max_ifng_feedback,
                                                              ifng_feedback_hill_degree) +
                                                          pow(bound_reading, ifng_feedback_hill_degree))) +
                                                        0.5; //(ifng_field_copy[x][y] / ifng_feedback_saturation_threshold) + 0.5;
                            total_production += ifng_production * ifng_feedback_term;
                          }
                        }
                      }
                    }
                  }
                }
              }  // end looping through cell pixels within ifng pixel

              //                double source_term = ifng_production * ifng_feedback_term * src_field[x][y][0] - ifng_decay * ifng_field_copy[x][y] - ifng_ART_uptake * ifng_saturation_term * src_field[x][y][1];
              double source_term = total_production - (ifng_decay + total_consumption) * ifng_field_copy[x][y]; // - ifng_ART_uptake * ifng_saturation_term * src_field[x][y][1];

              //                receptors[x][y] = 0; // clearing src and receptors to redefine them in the future
              //                for (int cv=0; cv<3; cv++) { // need to make these update every time a cell moves/added/removed from system for more efficiency
              //                    src_field[x][y][cv] = 0; //{0, 0, 0};
              //                }

              //                if(time > 1){
              //                    source_term = -ifng_uptake * ifng_saturation_term * (src_field[x][y][0] + src_field[x][y][1]);
              //                }
              double test_value = ifng_field[x][y];

              ifng_field[x][y] += time_step * time_conversion_factor * (ifng_diffusion * laplacian + source_term);

              if (ifng_field[x][y] < 0.0000000001) {
                ifng_field[x][y] = 0.0;
              }

              // decay bound ifng if there is decay
              //                if(update_bIFNG_decay){ // MAYBE PUT THIS IN THE LOOP ABOVE????? ///////////////////////////////////
              //                    bIFNG_decay
              //                    for (int xx = 0; xx < ifng_lattice_cell_converter; xx++) { // search local ifng lattice pixel
              //                        for (int yy = 0; yy < ifng_lattice_cell_converter; yy++) { // go through IFNG pixel in cell lattice coordinates
              //                            int latx = x * ifng_lattice_cell_converter + xx;
              //                            int laty = y * ifng_lattice_cell_converter + yy;
              //                            if(0 <= latx && latx < latCol && 0 <= laty && laty < latCol){
              //                                for(int id=0; id < 7; id++) {  // go through each all cell types
              //                                    ifng_tracker[latx][laty]
              //                                }
              //                            }
              //                        }
              //                    }
              //                }

              // update ifng_tracker as receptors bind local free ifng
              double total_local_ifng = total_bound + ifng_field[x][y];
              double final_free = -0.5*(total_receptors+dissC-total_local_ifng-std::sqrt(pow(total_receptors+dissC-total_local_ifng, 2)+4*dissC*total_local_ifng));
              double final_bound = total_local_ifng - final_free;
              double to_bind = final_bound - total_bound;
              bool bind = to_bind > 0;
              ifng_field[x][y] = final_free; // set local IFNG value
              if(!use_simple_decay_model){ // if we need to bind receptors
                double partition_fraction = final_bound / total_receptors;
                for (int xx = 0; xx < ifng_lattice_cell_converter; xx++) { // find sources in local ifng lattice pixel
                  for (int yy = 0; yy < ifng_lattice_cell_converter; yy++) { // go through IFNG pixel in cell lattice coordinates
                    int latx = x * ifng_lattice_cell_converter + xx;
                    int laty = y * ifng_lattice_cell_converter + yy;
                    if(0 <= latx && latx < latCol && 0 <= laty && laty < latCol){
                      for(int id=0; id < 11; id++) {  // go through each all cell types
                        int current_cellID = cellIDs[id]; // ifng tracker does not use these cellIDs!! These are for population, ifng_tracker uses 0-10
                        for(int cell=0; cell<ifng_tracker[latx][laty][id].size(); cell++) { // go through each cell of type id
                          double cell_bound = ifng_tracker[latx][laty][id][cell][0];
                          double cell_receptors = ifng_tracker[latx][laty][id][cell][1];
                          double partition = partition_fraction*cell_receptors;
                          ifng_tracker[latx][laty][id][cell][0] = partition;
                        }
                      }
                    }
                  }
                }

                // update/calculate bound ifng fractions from ifng_tracker
                for (int xx = 0; xx < ifng_lattice_cell_converter; xx++) { // find sources in local ifng lattice pixel
                  for (int yy = 0; yy < ifng_lattice_cell_converter; yy++) { // go through IFNG pixel in cell lattice coordinates
                    int latx = x * ifng_lattice_cell_converter + xx;
                    int laty = y * ifng_lattice_cell_converter + yy;
                    if(0 <= latx && latx < latCol && 0 <= laty && laty < latCol) {
                      for(int id=0; id < 11; id++) {  // go through each all cell types
                        for(int cell=0; cell<ifng_tracker[latx][laty][id].size(); cell++) { // go through each cell of type id
                          bound_ifng[latx][laty][2] += ifng_tracker[latx][laty][id][cell][1]; // total number of receptors in the ifng pixel
                          bound_ifng[latx][laty][1] += ifng_tracker[latx][laty][id][cell][0]; // number of bound receptors in the ifng pixel
                        }
                      }
                      if(bound_ifng[latx][laty][2] > 1){ // number of IFNG receptors must be at least 1, 19-10-2024
                        bound_ifng[latx][laty][0] = bound_ifng[latx][laty][1]/bound_ifng[latx][laty][2]; // fraction of bound receptors in the ifng pixel
                      }else{
                        bound_ifng[latx][laty][0] = -1; // was 0, make pixels with no cells have a value of -1 to differentiate from those with cells but no bound IFNG 19-10-2024
                      }
                    }
                  }
                }
              }

            }else{
                ifng_field[x][y] = 0.0;
            }
        }
    }

    // update ifng population and bound ifng in the dump files.
    total_ifng = 0.0;
    for (int xx = 0; xx < latCol; xx++) { // population columns
        for (int yy = 0; yy < latCol; yy++) { // population rows
            //double value1 = x+latCol*y;
            ifng_field_discretized[ifngID][xx + yy*latCol] = ifng_field[int(xx / ifng_lattice_cell_converter)][int(yy / ifng_lattice_cell_converter)]/pow(ifng_lattice_cell_converter, 2);
            ifng_field_discretized[boundifngID][xx + yy*latCol] = bound_ifng[xx][yy][0];
            double value = ifng_field_discretized[ifngID][xx + yy*latCol];
            total_ifng += value/pow(ifng_lattice_cell_converter, 2);
        }
    }
    //  add in updates to propensities of all cells???

}

/* ----------------------------------------------------------------------
   KMC solver on sectors
   can be invoked in serial or parallel
 ------------------------------------------------------------------------- */

void AppLattice::iterate_kmc_sector(double stoptime)
{
  int i,isite,done;
  double dt,timesector;
  double pmax,pmaxall;

  // save ptr to system solver

  Solve *hold_solve = solve;

  int alldone = 0;
  while (!alldone) {
    if (Ladapt) pmax = 0.0;

    for (int iset = 0; iset < nset; iset++) {
      timer->stamp();

      if (nprocs > 1) {
	comm->sector(iset);
	timer->stamp(TIME_COMM);
      }

      solve = set[iset].solve;
      
      propensity = set[iset].propensity;
      i2site = set[iset].i2site;
      int *site2i = set[iset].site2i;

      // update propensities for sites which neighbor a site outside sector
      // necessary since outside sites may have changed
      // attribute this chunk of time to comm, b/c due to decomposition

      int *bsites = set[iset].bsites;
      int *border = set[iset].border;
      int nborder = set[iset].nborder;

      int nsites = 0;
      for (int m = 0; m < nborder; m++) {
	i = border[m];
	isite = i2site[i];
	bsites[nsites++] = isite;
	propensity[isite] = site_propensity(i);
      }
      
      solve->update(nsites,bsites,propensity);
      timer->stamp(TIME_COMM);
      
      // pmax = maximum sector propensity per site

      if (Ladapt) {
	int ntmp = solve->get_num_active();
	if (ntmp > 0) {
	  double ptmp = solve->get_total_propensity();
	  ptmp /= ntmp;
	  pmax = MAX(ptmp,pmax);
	}
      }
      
      // execute events until sector time threshhold reached
      
      done = 0;
      timesector = 0.0;
      while (!done) {
	timer->stamp();
	isite = solve->event(&dt);
	timer->stamp(TIME_SOLVE);
	
	if (isite < 0) done = 1;
	else {
	  timesector += dt;
	  if (timesector >= dt_kmc) done = 1;
	  else {
	    site_event(site2i[isite],ranapp);
	    naccept++;
	  }
	  timer->stamp(TIME_APP);
	}
      }
      
      if (nprocs > 1) {
	comm->reverse_sector(iset);
	timer->stamp(TIME_COMM);
      }
    }

    if (allow_app_update) app_update(dt_kmc);

    // keep looping until overall time threshhold reached

    nsweeps++;

    time += dt_kmc;
    if (time >= stoptime) alldone = 1;
    if (alldone || time >= nextoutput)
      nextoutput = output->compute(time,alldone);
    timer->stamp(TIME_OUTPUT);

    // recompute dt_kmc if adaptive, based on pmax across all sectors

    if (Ladapt) {
      MPI_Allreduce(&pmax,&pmaxall,1,MPI_DOUBLE,MPI_MAX,world);
      if (pmaxall > 0.0) dt_kmc = nstop/pmaxall;
      else dt_kmc = stoptime-time;
      dt_kmc = MIN(dt_kmc,stoptime-time);
    }
  }

  // restore system solver

  solve = hold_solve;
}

/* ----------------------------------------------------------------------
   rejection KMC solver
 ------------------------------------------------------------------------- */

void AppLattice::iterate_rejection(double stoptime)
{
  int i,icolor,nselect,nrange,jset;
  int *site2i;

  // set loop is over:
  // sectors if there are sectors and no colors
  // colors if there are colors and no sectors
  // first nsector sets if there are both sectors and colors

  int nset_loop = nset;
  if (bothflag) nset_loop = nsector;

  int done = 0;
  while (!done) {
    for (int iset = 0; iset < nset_loop; iset++) {
      if (nprocs > 1) {
	timer->stamp();
	if (sectorflag) comm->sector(iset);
	else comm->all();
	timer->stamp(TIME_COMM);
      }

      if (Lmask) boundary_clear_mask(iset);

      timer->stamp();

      // sectors but no colors (could also be no sectors)
      // random selection of sites in iset

      if (sweepflag == RANDOM) {
	site2i = set[iset].site2i;
	nrange = set[iset].nlocal;
	nselect = set[iset].nselect;
	for (i = 0; i < nselect; i++) 
	  sitelist[i] = site2i[ranapp->irandom(nrange) - 1];
	(this->*sweep)(nselect,sitelist);
	nattempt += nselect;

      // sectors but no colors, or colors but no sectors
      // ordered sweep over all sites in iset

      } else if (bothflag == 0) {
	for (i = 0; i < set[iset].nloop; i++)
	  (this->*sweep)(set[iset].nlocal,set[iset].site2i);
	nattempt += set[iset].nselect;

      // sectors and colors
      // icolor loop is over all colors in a sector
      // jset = set that contains sites of one color in one sector
      // ordered sweep over all sites in jset

      } else {
	for (icolor = 0; icolor < ncolors; icolor++) {
	  jset = nsector + iset*ncolors + icolor;
	  for (i = 0; i < set[jset].nloop; i++)
	    (this->*sweep)(set[jset].nlocal,set[jset].site2i);
	  nattempt += set[jset].nselect;
	}
      }

      timer->stamp(TIME_SOLVE);

      if (nprocs > 1) {
	if (sectorflag) comm->reverse_sector(iset);
	else comm->all_reverse();
	timer->stamp(TIME_COMM);
      }
    }

    if (allow_app_update) app_update(dt_rkmc);

    nsweeps++;
    time += dt_rkmc;
    if (time >= stoptime) done = 1;
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }
}

/* ----------------------------------------------------------------------
   iterate the app_update routine only
   app is responsible for doing communciation in app_update()
------------------------------------------------------------------------- */

void AppLattice::iterate_app_update_only(double stoptime,double dt)
{
  int done = 0;
  while (!done) {
    if (allow_app_update) app_update(dt);
    
    time += dt;
    if (time >= stoptime) done = 1;
    if (done || time >= nextoutput) nextoutput = output->compute(time,done);
    timer->stamp(TIME_OUTPUT);
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::sweep_nomask_nostrict(int n, int *list)
{
  for (int m = 0; m < n; m++)
    site_event_rejection(list[m],ranapp);
}

/* ---------------------------------------------------------------------- */

void AppLattice::sweep_mask_nostrict(int n, int *list)
{
  int i;
  for (int m = 0; m < n; m++) {
    i = list[m];
    if (mask[i]) continue;
    site_event_rejection(i,ranapp);
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::sweep_nomask_strict(int n, int *list)
{
  int i;
  for (int m = 0; m < n; m++) {
    i = list[m];
    ranstrict->seed = siteseeds[i];
    site_event_rejection(i,ranstrict);
    siteseeds[i] = ranstrict->seed;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::sweep_mask_strict(int n, int *list)
{
  int i;
  for (int m = 0; m < n; m++) {
    i = list[m];
    if (mask[i]) continue;
    ranstrict->seed = siteseeds[i];
    site_event_rejection(i,ranstrict);
    siteseeds[i] = ranstrict->seed;
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::input_app(char *command, int narg, char **arg)
{
  error->all(FLERR,"Unrecognized command");
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_sector(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal sector command");

  nsector_user = 0;
  if (strcmp(arg[0],"yes") == 0) sectorflag = 1;
  else if (strcmp(arg[0],"no") == 0) sectorflag = 0;
  else {
    sectorflag = 1;
    nsector_user = atoi(arg[0]);
    if (nsector_user != 2 && nsector_user != 4 && nsector_user != 8)
      error->all(FLERR,"Illegal sector command");
  }

  nstop = 1.0;
  tstop = 0.0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nstop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal sector command");
      nstop = atof(arg[iarg+1]);
      if (nstop <= 0.0) error->all(FLERR,"Illegal sector command");
      tstop = 0.0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"tstop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal sector command");
      tstop = atof(arg[iarg+1]);
      if (tstop <= 0.0) error->all(FLERR,"Illegal sector command");
      nstop = 0.0;
      iarg += 2;
    } else error->all(FLERR,"Illegal sector command");
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_sweep(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal sweep command");
  if (strcmp(arg[0],"random") == 0) sweepflag = RANDOM;
  else if (strcmp(arg[0],"raster") == 0) sweepflag = RASTER;
  else if (strcmp(arg[0],"color") == 0) sweepflag = COLOR;
  else if (strcmp(arg[0],"color/strict") == 0) sweepflag = COLOR_STRICT;
  else if (strcmp(arg[0],"none") == 0) sweepflag = NOSWEEP;
  else error->all(FLERR,"Illegal sweep command");

  Lmask = false;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mask") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal sweep command");
      if (strcmp(arg[iarg+1],"no") == 0) Lmask = false;
      else if (strcmp(arg[iarg+1],"yes") == 0) Lmask = true;
      else error->all(FLERR,"Illegal sweep command");
      iarg += 2;
    } else error->all(FLERR,"Illegal sweep command");
  }
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_temperature(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal temperature command");
  temperature = atof(arg[0]);
  if (temperature != 0.0) t_inverse = 1.0/temperature;
}

/* ---------------------------------------------------------------------- */

void AppLattice::set_app_update_only(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal app_update_only command");
  if (strcmp(arg[0],"yes") == 0) app_update_only = 1;    
  else if (strcmp(arg[0],"no") == 0) app_update_only = 0;
  else error->all(FLERR,"Illegal app_update_only command");

  if (app_update_only && !allow_app_update)
    error->all(FLERR,"App does not permit app_update yes");
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void AppLattice::stats(char *strtmp)
{
  char big[8],format[64];
  strcpy(big,BIGINT_FORMAT);

  bigint naccept_all;
  MPI_Allreduce(&naccept,&naccept_all,1,MPI_SPK_BIGINT,MPI_SUM,world);

  if (solve) {
    sprintf(format,"%%10g %%10%s %%10d %%10d",&big[1]);
    sprintf(strtmp,format,time,naccept_all,0,nsweeps);
  } else {
    bigint nattempt_all;
    MPI_Allreduce(&nattempt,&nattempt_all,1,MPI_SPK_BIGINT,MPI_SUM,world);
    sprintf(format,"%%10g %%10%s %%10%s %%10d",&big[1],&big[1]);
    sprintf(strtmp,format,time,naccept_all,nattempt_all-naccept_all,nsweeps);
  }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void AppLattice::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s %10s %10s %10s","Time","Naccept","Nreject","Nsweeps");
}

/* ----------------------------------------------------------------------
   create a subset of owned sites
   insure all ptrs in Set data struct are allocated or NULL
   isector = 0 = all sites (no sector)
   isector > 1 = sites within a sector
   icolor = 0 = all sites (no color)
   icolor > 1 = sites of a certain color
 ------------------------------------------------------------------------- */

void AppLattice::create_set(int iset, int isector, int icolor, Solve *oldsolve)
{
  // sector boundaries

  double xmid = 0.5 * (domain->subxlo + domain->subxhi);
  double ymid = 0.5 * (domain->subylo + domain->subyhi);
  double zmid = 0.5 * (domain->subzlo + domain->subzhi);

  // count sites in subset

  int flag,iwhich,jwhich,kwhich,msector,mcolor;

  int delcolor = delevent + delpropensity;

  int n = 0;
  for (int i = 0; i < nlocal; i++) {
    flag = 1;

    if (isector > 0) {
      if (xyz[i][0] < xmid) iwhich = 0;
      else iwhich = 1;
      if (xyz[i][1] < ymid) jwhich = 0;
      else jwhich = 1;
      if (xyz[i][2] < zmid) kwhich = 0;
      else kwhich = 1;

      if (nsector == 2) msector = iwhich + 1;
      else if (nsector == 4) msector = 2*jwhich + iwhich + 1;
      else msector = 4*kwhich + 2*jwhich + iwhich + 1;

      if (isector != msector) flag = 0;
    }
    
    if (icolor > 0) {
      mcolor = domain->lattice->id2color(id[i],delcolor);
      if (icolor != mcolor) flag = 0;
    }

    if (flag) n++;
  }

  set[iset].nlocal = n;

  // setup site2i for sites in set

  memory->create(set[iset].site2i,n,"app:site2i");

  n = 0;
  for (int i = 0; i < nlocal; i++) {
    flag = 1;

    if (isector > 0) {
      if (xyz[i][0] < xmid) iwhich = 0;
      else iwhich = 1;
      if (xyz[i][1] < ymid) jwhich = 0;
      else jwhich = 1;
      if (xyz[i][2] < zmid) kwhich = 0;
      else kwhich = 1;

      if (nsector == 2) msector = iwhich + 1;
      else if (nsector == 4) msector = 2*jwhich + iwhich + 1;
      else msector = 4*kwhich + 2*jwhich + iwhich + 1;

      if (isector != msector) flag = 0;
    }

    if (icolor > 0) {
      mcolor = domain->lattice->id2color(id[i],delcolor);
      if (icolor != mcolor) flag = 0;
    }

    if (flag) set[iset].site2i[n++] = i;
  }

  // setup i2site for sites in set, only for KMC solver
  // i2site = 0 to nsite-1 for owned points in set, else -1

  if (solve) {
    (int *) memory->create(set[iset].i2site,nlocal+nghost,"app:i2site");
    for (int i = 0; i < nlocal+nghost; i++) set[iset].i2site[i] = -1;
    for (int i = 0; i < set[iset].nlocal; i++) 
      set[iset].i2site[set[iset].site2i[i]] = i;
  } else set[iset].i2site = NULL;

  // allocate propensity array for set

  memory->create(set[iset].propensity,n,"app:propensity");

  // allocate KMC solver for set
  // reuse old solver if provided, else delete it

  if (solve) {
    if (oldsolve) set[iset].solve = oldsolve;
    else set[iset].solve = solve->clone();
  } else {
    set[iset].solve = NULL;
    delete oldsolve;
  }

  // setup border arrays for set
  // nborder = # of sites in sector influenced by site outside sector
  // border = list of these sites, stored as lattice indices
  // bsites = scratch array for use by KMC solver
  // border is only used by KMC solver in sectors and masking
  // bsites is only used by KMC solver in sectors

  if ((solve && sectorflag) || Lmask) {
    set[iset].nborder = find_border_sites(iset);
    if (solve && sectorflag)
      memory->create(set[iset].bsites,set[iset].nborder,"app:bsites");
    else set[iset].bsites = NULL;
  } else {
    set[iset].nborder = 0;
    set[iset].border = NULL;
    set[iset].bsites = NULL;
  }
}

/* ----------------------------------------------------------------------
   free memory inside a set
   return Solver so caller can reuse it if desired
 ------------------------------------------------------------------------- */

Solve *AppLattice::free_set(int iset)
{
  memory->destroy(set[iset].border);
  memory->destroy(set[iset].bsites);
  memory->destroy(set[iset].propensity);
  memory->destroy(set[iset].site2i);
  memory->destroy(set[iset].i2site);
  return set[iset].solve;
}

/* ----------------------------------------------------------------------
   create list of border sites for a set
   border site = site in set with a 1 to Nlayer neighbor outside the set
   neighbor can be another owned site (outside set) or a ghost
   border = lattice index of the sites
 ------------------------------------------------------------------------- */

int AppLattice::find_border_sites(int isector)
{
  int i,j,m;

  int nlayer = delpropensity+delevent;
  int ntotal = nlocal + nghost;
  int nsites = set[isector].nlocal;
  int *site2i = set[isector].site2i;

  // flag sites with -1 that are not in sector
  // flag sites with 0 that are in sector

  int *flag;
  memory->create(flag,ntotal,"app:flag");
  for (i = 0; i < ntotal; i++) flag[i] = -1;
  for (m = 0; m < nsites; m++) flag[site2i[m]] = 0;

  // flag sector sites with -1 that have non-sector neighbor up to nlayer away

  for (int ilayer = 0; ilayer < nlayer; ilayer++) {
    for (m = 0; m < nsites; m++) {
      i = site2i[m];
      if (flag[i]) continue;
      for (j = 0; j < numneigh[i]; j++)
	if (flag[neighbor[i][j]] < 0) break;
      if (j < numneigh[i]) flag[i] = 1;
    }
    for (m = 0; m < nsites; m++) {
      i = site2i[m];
      if (flag[i] > 0) flag[i] = -1;
    }
  }

  // nborder = # of border sites
  // allocate border and fill with site indices

  int nborder = 0;
  for (m = 0; m < nsites; m++) {
    i = site2i[m];
    if (flag[i] < 0) nborder++;
  }

  int *border;
  memory->create(border,nborder,"app:border");

  nborder = 0;
  for (m = 0; m < nsites; m++) {
    i = site2i[m];
    if (flag[i] < 0) border[nborder++] = i;
  }

  memory->destroy(flag);

  set[isector].border = border;
  return nborder;
}
  
/* ----------------------------------------------------------------------
   unset all mask values of owned sites in iset whose propensity
     could change due to events on sites one neighbor outside the set
   border list stores indices of these sites
   their mask value may be out of date, due to state change in other sets
 ------------------------------------------------------------------------- */

void AppLattice::boundary_clear_mask(int iset)
{
  int *border = set[iset].border;
  int nborder = set[iset].nborder;

  for (int m = 0; m < nborder; m++) mask[border[m]] = 0;
}

/* ----------------------------------------------------------------------
   push new site onto stack and assign new id
 ------------------------------------------------------------------------- */

void AppLattice::push_new_site(int i, int* cluster_ids, int id,
					  std::stack<int>* cluststack)
{
  // This can be used to screen out unwanted spin values
  // int isite = iarray[0][i];

  cluststack->push(i);
  cluster_ids[i] = id;
}

/* ----------------------------------------------------------------------
   push connected neighbors of this site onto stack
     and assign current id
   ghost neighbors are masked by id = -1
   previously burned sites are masked by id > 0
 ------------------------------------------------------------------------- */

void AppLattice::push_connected_neighbors(int i, int* cluster_ids, int id,
					  std::stack<int>* cluststack)
{
  int ii;
  int isite = iarray[0][i];

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (iarray[0][ii] == isite && cluster_ids[ii] == 0) {
      cluststack->push(ii);
      cluster_ids[ii] = id;
      domain->pbcshift(xyz[i],xyz[ii]);
    }
  }
}

/* ----------------------------------------------------------------------
   add cluster id of connected ghost sites to neighbor list of cluster
 ------------------------------------------------------------------------- */

void AppLattice::connected_ghosts(int i, int* cluster_ids, 
				  Cluster* clustlist, int idoffset)
{
  int iclust;
  int ii;
  int isite = iarray[0][i];
  int pbcflags[3];

  // check if this was a site that was ignored

  if (cluster_ids[i] == 0) return;

  iclust = cluster_ids[i]-idoffset;
  
  // add ghost cluster to neighbors of local cluster

  for (int j = 0; j < numneigh[i]; j++) {
    ii = neighbor[i][j];
    if (iarray[0][ii] == isite && ii >= nlocal) {
      clustlist[iclust].add_neigh(cluster_ids[ii]);
      domain->set_pbcflags(xyz[i],xyz[ii],pbcflags);
      clustlist[iclust].add_pbcflags(cluster_ids[ii],pbcflags);
    }
  }
}

/* ----------------------------------------------------------------------
   grow per-site arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n 
------------------------------------------------------------------------- */

void AppLattice::grow(int n) // gdag: called by add_site()
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  memory->grow(id,nmax,"app:id");
  memory->grow(xyz,nmax,3,"app:xyz");
  memory->grow(owner,nmax,"app:owner");
  memory->grow(index,nmax,"app:index");

  memory->grow(numneigh,nmax,"app:numneigh");
  if (maxneigh) memory->grow(neighbor,nmax,maxneigh,"app:neighbor");

  for (int i = 0; i < ninteger; i++)
    memory->grow(iarray[i],nmax,"app:iarray");
  for (int i = 0; i < ndouble; i++)
    memory->grow(darray[i],nmax,"app:darray");

  grow_app();
}

/* ----------------------------------------------------------------------
   add an owned site
   called from create_sites command
   grow arrays if necessary
 ------------------------------------------------------------------------- */

void AppLattice::add_site(tagint n, double x, double y, double z) // gdag: called by create_sites.cpp
{
  if (nlocal == nmax) grow(0);

  id[nlocal] = n;
  xyz[nlocal][0] = x;
  xyz[nlocal][1] = y;
  xyz[nlocal][2] = z;

  owner[nlocal] = me;
  index[nlocal] = nlocal;

  for (int i = 0; i < ninteger; i++) iarray[i][nlocal] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][nlocal] = 0.0;

  nlocal++;
}

/* ----------------------------------------------------------------------
   add a ghost site
   called from create_sites command
   grow arrays if necessary
 ------------------------------------------------------------------------- */

void AppLattice::add_ghost(tagint n, double x, double y, double z,
			   int proc, int index_owner)
{
  if (nlocal+nghost == nmax) grow(0);

  int m = nlocal + nghost;

  id[m] = n;
  xyz[m][0] = x;
  xyz[m][1] = y;
  xyz[m][2] = z;

  owner[m] = proc;
  index[m] = index_owner;

  for (int i = 0; i < ninteger; i++) iarray[i][m] = 0;
  for (int i = 0; i < ndouble; i++) darray[i][m] = 0.0;

  nghost++;
}

/* ----------------------------------------------------------------------
   set neighbor connectivity for owned site I
   nvalues = # of neighbors
   called from read_sites command
 ------------------------------------------------------------------------- */

void AppLattice::add_neighbors(int i, int nvalues, char **values)
{
  numneigh[i] = nvalues;
  for (int m = 0; m < nvalues; m++)
    neighbor[i][m] = atoi(values[m]);
}

/* ----------------------------------------------------------------------
   set all values for owned site I
   called from read_sites command
 ------------------------------------------------------------------------- */

void AppLattice::add_values(int i, char **values)
{
  for (int m = 0; m < ninteger; m++) iarray[m][i] = atoi(values[m]);
  for (int m = 0; m < ndouble; m++) darray[m][i] = atof(values[m+ninteger]);
}

/* ----------------------------------------------------------------------
   set single value for owned site I
   type = 0/1 for INT/DOUBLE
   index = 1 to Ninteger or 1 to Ndouble
   called from read_sites command
 ------------------------------------------------------------------------- */

void AppLattice::add_value(int i, int type, int index, char *value)
{
  if (type == 0) iarray[index-1][i] = atoi(value);
  else if (type == 1) darray[index-1][i] = atof(value);
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   lo,hi = inclusive bounds
   5 possibilities:
     (1) i = i to i, (2) * = lo to hi,
     (3) i* = i to hi, (4) *j = lo to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void AppLattice::bounds(char *str, int lo, int hi, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = MAX(atoi(str),lo);
    nhi = MIN(atoi(str),hi);
  } else if (strlen(str) == 1) {
    nlo = lo;
    nhi = hi;
  } else if (ptr == str) {
    nlo = lo;
    nhi = MIN(atoi(ptr+1),hi);
  } else if (strlen(ptr+1) == 0) {
    nlo = MAX(atoi(str),lo);
    nhi = hi;
  } else {
    nlo = MAX(atoi(str),lo);
    nhi = MIN(atoi(ptr+1),hi);
  }
}

/* ----------------------------------------------------------------------
   print connectivity stats
 ------------------------------------------------------------------------- */

void AppLattice::print_connectivity()
{
  int i;

  tagint min = maxneigh;
  tagint max = 0;

  for (i = 0; i < nlocal; i++) {
    min = MIN(min,numneigh[i]);
    max = MAX(max,numneigh[i]);
  }

  tagint minall,maxall;
  MPI_Allreduce(&min,&minall,1,MPI_SPK_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&max,&maxall,1,MPI_SPK_TAGINT,MPI_MAX,world);

  tagint *count = new tagint[maxall+1];
  tagint *countall = new tagint[maxall+1];

  for (i = 0; i <= maxall; i++) count[i] = 0;

  for (i = 0; i < nlocal; i++) count[numneigh[i]]++;
  MPI_Allreduce(count,countall,maxall+1,MPI_SPK_TAGINT,MPI_SUM,world);

  if (me == 0)
    for (i = minall; i <= maxall; i++) {
      if (screen)
	fprintf(screen,"  " TAGINT_FORMAT " sites have %d neighbors\n",\
		countall[i],i);
      if (logfile)
	fprintf(logfile,"  " TAGINT_FORMAT " sites have %d neighbors\n",
		countall[i],i);
    }

  delete [] count;
  delete [] countall;
}

/* ----------------------------------------------------------------------
   sum and print memory usage
   result is only memory on proc 0, not averaged across procs
------------------------------------------------------------------------- */

bigint AppLattice::memory_usage()
{
  bigint bytes = 0;

  bytes += nmax * sizeof(tagint);           // id
  bytes += nmax*3 * sizeof(double);         // xyz
  bytes += ninteger*nmax * sizeof(int);     // iarray
  bytes += ndouble*nmax * sizeof(double);   // darray
  
  bytes += nmax * sizeof(int);              // numneigh
  bytes += nmax*maxneigh * sizeof(int);     // neighbor

  return bytes;
}
