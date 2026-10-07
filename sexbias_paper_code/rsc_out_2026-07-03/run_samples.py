#!/usr/bin/env python3

import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
import sys
import math as m
from scipy.stats import linregress
from datetime import date
from datetime import datetime
import random
import os

# import SPPARKS_init_cond_generator as initWriter

now = datetime.now()

today = date.today()
todayDate = today.strftime("%d-%m-%Y")
current_time = now.strftime("%H,%M,%S")

####################################################   Take input from user #######################################
def checkInputBool(string, print_out=None):
    string = string.strip(" ").lower()
    if string == "yes" or string == "yeah" or string == "true" or string == "ok" or string == "sure" or string == "y" or string == "yah" or string == "t":
        booly = True
    elif string == "no" or string == "nope" or string == "nah" or string == "uhuh" or string == "false" or string == "non" or string == "n" or string == "f" or string == "":
        booly = False
    else:
        booly = False
        print("Input not taken from: " + string + ", input assumed false")

    if print_out is not None:
        print(print_out[booly])

    return booly

dataFileName = 0
printOut = False

ask_questions = False
if ask_questions:
    num_Samples = int(input("How many samples to run? ").strip(" "))
    print("Running "+str(num_Samples)+" samples with different seeds...")

    directory_title = ''

    create_dump_files = checkInputBool(input("Create dump.2D_model.x files?... "), print_out=["Not creating dump files...", "Creating dump files..."])
    dump_filename = ['_', '_dump_'][create_dump_files]

    regional = checkInputBool(input("Do regional simulation?..."), print_out=["Not using regional initialization file... ", "Using regional initialization file..."])
    if regional:
        subregional = False
        # subregional = checkInputBool(input("Do subregional simulation?..."), print_out=["not using subregional initialization file... ", "Using subregional initialization file..."])

    secluded_AR = checkInputBool(input("Seclude ARs?..."), print_out=["Not secluding ARs... ", "Using secluding ARs..."])

    # print(["Using male input file...", "Using female input file..."][regional])
    use_fixedT_ifng_init = False

    # Choose one of the two initializations
    female_or_male = checkInputBool(input("Female (y) or Male (n) initial condition? "), print_out=["Using male initial condition...", "Using female initial condition..."])
    # slide_filename = ['d14f4r1mM', 'd14f4r1'][female_or_male]

    if regional:
        real_or_made_female = False
        fix_T_cells = checkInputBool(input("Should CD8+ T cells be fixed? "),
                                             print_out=["Not CD8+ fixing T cells",
                                                        "Fixing all CD8+ T cells..."])
        if not fix_T_cells:
            use_fixedT_ifng_init = checkInputBool(input("Use 1 week fixed T IFNG initialization? "),
                                             print_out=["Not CD8+ fixing T cells",
                                                        "Fixing all CD8+ T cells..."])  # this will be to use ifng from fixedt distrib from 1 week
        # if female_or_male:
        #     real_or_made_female = checkInputBool(input("Real (y) or Male Made Female (n) initial condition? "),
        #                                     print_out=["Using male made female initial condition...", "Using female initial condition..."])
        #
        # if real_or_made_female:
        #     slide_filename = ['d14m4r4', 'd14f4r4'][female_or_male]
        # else:
        #     slide_filename = ['d14m4r4', 'd14m4r4mF'][female_or_male]

        if subregional:
            slide_filename = ['d14m4r3', 'd14m4r3mF'][female_or_male]
        else:
            slide_filename = ['d14m4r4', 'd14m4r4mF'][female_or_male]
    else:
        real_or_made_female = False
        fix_T_cells = checkInputBool(input("Should CD8+ T cells be fixed? "),
                                             print_out=["Not CD8+ fixing T cells",
                                                        "Fixing all CD8+ T cells..."])
        if not fix_T_cells:
            use_fixedT_ifng_init = checkInputBool(input("Use 1 week fixed T IFNG initialization? "),
                                             print_out=["Not using 1 week IFNG initialization",
                                                             "Using 1 week IFNG initialization..."])  # this will be to use ifng from fixedt distrib from 1 week

        if female_or_male:
            real_or_made_female = checkInputBool(input("Real (y) or Male Made Female (n) initial condition? "),
                                            print_out=["Using male made female initial condition...", "Using female initial condition..."])
        if real_or_made_female:
            slide_filename = ['d14m4', 'd14f4'][female_or_male]
        else:
            slide_filename = ['d14m4', 'd14m4mF'][female_or_male]

    sex_rates = False
    if not fix_T_cells:
        sex_rates = checkInputBool(input("Use female (y) or male (n) input file (rates)?... "), print_out=["Using male input file...", "Using female input file..."])
    # print(["Using male input file...", "Using female input file..."][sex_rates])
    sex_rates_filename = ['mInp_', 'fInp_'][sex_rates]
    if fix_T_cells:
        sex_rates_filename = 'fxdT_'

    directory_title = slide_filename + dump_filename + sex_rates_filename


    print("Running samples with different seeds...")
    # print("Running "+["nonresponder 01RD", "responder 04RD"][res_nonres] + " samples with different seeds...")

    usr_direct = input("Please input description tag of mother directory that will contain all sample directories: ").strip(" ")
    if usr_direct == "":
        usr_direct = "Run_"
else:
    num_Samples = 3

    directory_title = ''

    create_dump_files = True
    dump_filename = ['_', '_dump_'][create_dump_files]

    regional = False
    subregional = False
        # subregional = checkInputBool(input("Do subregional simulation?..."), print_out=["not using subregional initialization file... ", "Using subregional initialization file..."])

    secluded_AR = False

    use_fixedT_ifng_init = False

    # Choose one of the two initializations
    female_or_male = True

    if regional:
        real_or_made_female = False
        fix_T_cells = False
        if not fix_T_cells:
            use_fixedT_ifng_init = False

        if subregional:
            slide_filename = ['d14m4r3', 'd14m4r3mF'][female_or_male]
        else:
            slide_filename = ['d14m4r4', 'd14m4r4mF'][female_or_male]
    else:
        real_or_made_female = False
        fix_T_cells = False
        if not fix_T_cells:
            use_fixedT_ifng_init = False

        if female_or_male:
            real_or_made_female = False # this is where you set either sham or non-sham female
        if real_or_made_female:
            slide_filename = ['d14m4', 'd14f4'][female_or_male]
        else:
            slide_filename = ['d14m4', 'd14m4mF'][female_or_male]

    sex_rates = False
    if not fix_T_cells:
        sex_rates = female_or_male # female rates
    # print(["Using male input file...", "Using female input file..."][sex_rates])
    sex_rates_filename = ['mInp_', 'fInp_'][sex_rates]
    if fix_T_cells:
        sex_rates_filename = 'fxdT_'

    directory_title = slide_filename + dump_filename + sex_rates_filename

    print("Running samples with different seeds...")
    # print("Running "+["nonresponder 01RD", "responder 04RD"][res_nonres] + " samples with different seeds...")

    usr_direct = input(
        "Please input description tag of mother directory that will contain all sample directories: ").strip(" ")
    if usr_direct == "":
        usr_direct = "Run_"

directory_title += usr_direct
newDirName = directory_title + "_"+str(todayDate)+"_"+str(num_Samples)+"_samples"
print("Directory will be titled: "+newDirName)

initial_Spatial_Dist_RandBool = False
# initial_Spatial_Dist_RandBool = checkInputBool(input("Vary initial spatial distribution? "))
# print(["Using a constant initial spatial distribution across all samples...", "Varying initial distribution from sample to sample..."][int(initial_Spatial_Dist_RandBool)])

# if initial_Spatial_Dist_RandBool:
#     init_Canc = int(input("How many initial cancer cells? ").strip(" "))
#     print("Starting with "+str(init_Canc)+" cancer cells...")
#
#     init_Rec_Sites = int(input("How many recruitment sites? ").strip(" "))
#     print("Starting with " + str(init_Rec_Sites) + " recruitment sites...")

init_Canc = 1000
init_Rec_Sites = 100
init_Tc = 0
init_Tn = 0
init_Te = 0
init_mac = 0

# del_dump = True
# del_dump = checkInputBool(input("Delete SPPARKS dump files? "))
# print(["Not deleting dump files...", "Deleting ALL dump files..."][int(del_dump)])
#
# avPlots = checkInputBool(input("Plot average population vs time plots: "))
# print(["Not plotting average population vs time plots...", "Plotting average population vs time plots..."][int(avPlots)])
#
# popCloudPlots = checkInputBool(input("Plot all population vs time plots: "))
# print(["Not plotting all population vs time plots...", "Plotting all population vs time plots..."][int(popCloudPlots)])

##############################################################
def first_changes_to_input_file(fileName, create_dump_files=False, make_female=0.9, init_pos_file_name=None):
    delimiter = [None, '\t', '\s{1,}']
    # os.system("pwd")
    female_diff_factor = 2  # make the input file rates for female (delete AR+ recruitment rates and increase AR- rates accordingly)
    female_recruitment_factor = 0
    if make_female is not None:
        female_recruitment_factor = 1/make_female
        make_female = True

    with open(fileName, 'r') as inp:
        mem = ""
        it = 0
        for ls in inp:
            new_line = ls
            if init_pos_file_name is not None and ls.find("read_sites") >= 0:
                new_line = "read_sites       " + init_pos_file_name + "\n"

            if create_dump_files and ls.find('dump.2D_model.*') >= 0:
                new_line = ls.replace('#', '')
                # print(ls)
            # make male input female
            if make_female:
                if ls.find('recruitment') >= 0 and ls.find('app_style') < 0:
                    if ls.find('AR+') >= 0: # set AR+ T recruitment to 0
                        value = float(ls.split(delimiter[0])[5])
                        new_line = ls.replace(str(value), str(0))
                    elif ls.find('Tp recruitment') >= 0 or ls.find('Tm recruitment') >= 0:
                        value = ls.split(delimiter[0])[5] # get male recruitment value
                        new_value = female_recruitment_factor*float(value)
                        new_line = ls.replace(str(value), str(new_value))
                if ls.find('Tp differentiation') >= 0 and ls.find('app_style') < 0:
                    if ls.find('AR+ Tm') >= 0:
                        value = float(ls.split(delimiter[0])[5])
                        new_line = ls.replace(str(value), str(0))
                    else: # not the GONE interaction, double the differentiation rate to Tm to match total differentiation rate as in male mouse
                        value = ls.split(delimiter[0])[5] # get male recruitment value
                        new_value = female_diff_factor*float(value)
                        new_line = ls.replace(str(value), str(new_value))

            mem += new_line
            # print(mem)
            it += 1
    with open(fileName, 'w') as inp:
        inp.write(mem)

    return 0.0


# old input file editor
def changeInputSeed(fileName, seed_change):
    delimiter = [None, '\t', '\s{1,}']
    # os.system("pwd")
    with open(fileName, 'r') as inp:
        mem = ""
        it = 0
        for ls in inp:
            # print("yo, "+str(it))
            # print(ls.split(delimiter[0]))
            # print(ls)
            new_line = ls
            if it == 2:
                current_seed = int(ls.split(delimiter[0])[1])
                new_seed = current_seed + seed_change
                new_line = ls.replace(str(current_seed), str(new_seed))
                print("\ncurrent seed: "+str(current_seed))
            mem += new_line
            # print(mem)
            it += 1
    with open(fileName, 'w') as inp:
        inp.write(mem)

    return 0.0

# newer input file editor
def changeInputFile(fileName, seed_change, randomize_seed, init_pos_file_name=None, create_dump_files=False):
    delimiter = [None, '\t', '\s{1,}']
    # os.system("pwd")
    with open(fileName, 'r') as inp:
        mem = ""
        it = 0
        for ls in inp:
            # print("yo, "+str(it))
            # print(ls.split(delimiter[0]))
            # print(ls)
            new_line = ls
            if it == 2:
                current_seed = int(ls.split(delimiter[0])[1])
                new_seed = current_seed + seed_change
                if randomize_seed:
                    new_seed = int(random.random() * 10000)
                new_line = ls.replace(str(current_seed), str(new_seed))
                print("\ncurrent seed: "+str(current_seed))
            if init_pos_file_name is not None and ls.find("read_sites") >= 0:
                new_line = "read_sites       " + init_pos_file_name + "\n"

            if create_dump_files and ls.find('dump.2D_model.*') >= 0:
                new_line = ls.replace('#', '')
                # print(ls)

            mem += new_line
            # print(mem)
            it += 1
    with open(fileName, 'w') as inp:
        inp.write(mem)

    return 0.0
# Load modules, create a new directory to put in sample sub-directories, navigate to the new directory
input_file = 'in.simulation'
bash_file = "auto_run.sh"  # "run_sim_0.sh"

sex_of_input_file = ['m', 'f'][sex_rates]
region_of_input_file = ['full', 'regional'][regional]
# input_file_to_copy = 'in.' + sex_of_input_file + '_' + region_of_input_file
input_file_to_copy = 'in.' + region_of_input_file

if fix_T_cells:
    if not regional:
        input_file_to_copy = 'in.full_fixed'
    else:
        input_file_to_copy = 'in.regional_fixed'


# input_file = "in.const_in_cond_simple"

redsky_file = "spk_redsky"
# redsky_file = "spk_redsky_gdb"

init_folder = "sex_bias_initialization_files"

if female_or_male:  # female
    init_file = 'single_cell_ifng_initializationSPPARKS_init'
    init_file = "d7_F_3_SPPARKS_init"
    # init_file = 'd14_F_4_cd8Reloc_SPPARKS_init'
    init_file = 'd14_F_4_SPPARKS_init'
    init_file = 'd14_M_4_madeF_SPPARKS_init'
    # init_file = 'd14_M_4_madeF_secAR_SPPARKS_init'
    if regional:
        # init_file = 'd14_F_4_r1_SPPARKS_init'
        init_file = 'd14_M_4_r4_madeF_SPPARKS_init'  # 'd14_M_4_r2_madeF_noBckgrnd_SPPARKS_init'  #  # 'd14_M_4_r2_madeF_noBckgrnd_SPPARKS_init'  # 'd14_M_4_r1_madeF_SPPARKS_init'
        # init_file = 'd14_M_4_r4_madeF_SPPARKS_initinitIFNG110pp'
        if secluded_AR:
            init_file = 'd14_M_4_r4_madeF_secAR_SPPARKS_init'
        if subregional:
            init_file = 'd14_M_4_r3_madeF_SPPARKS_init'

        if use_fixedT_ifng_init:
            init_file = 'd14_M_4_r4_madeF_fxdT_ifng_t1008_SPPARKS_init'
    else:
        if real_or_made_female:
            init_file = 'd14_F_4_SPPARKS_init'
            if secluded_AR:
                init_file = 'd14_F_4_secAR_SPPARKS_init'
            if use_fixedT_ifng_init:
                init_file = 'd14_F_4_fxdT_ifng_t1008_SPPARKS_init'
        else:
            init_file = 'd14_M_4_madeF_SPPARKS_init'
            if secluded_AR:
                init_file = 'd14_M_4_madeF_secAR_SPPARKS_init'
            if use_fixedT_ifng_init:
                init_file = 'd14_M_4_madeF_fxdT_ifng_t1008_SPPARKS_init'
            # init_file = 'd14_M_4_madeF_SPPARKS_initinitIFNG110pp'


    # init_file = 'd14_F_4_r2_SPPARKS_init'
else:  # male
    # init_file = 'd7_M_1_SPPARKS_init'
    init_file = 'd14_M_4_SPPARKS_init'
    if secluded_AR:
        init_file = 'd14_M_4_secAR_SPPARKS_init'

    if use_fixedT_ifng_init:
        init_file = 'd14_M_4_fxdT_ifng_t1008_SPPARKS_init'

    if regional:
        init_file = 'd14_M_4_r4_SPPARKS_init'  # 'd14_M_4_r2_noBckgrnd_SPPARKS_init' # 'd14_M_4_r4_SPPARKS_init'  # 'd14_M_4_r2_noBckgrnd_SPPARKS_init'  #'d14_M_4_r1_SPPARKS_init'
        # init_file = 'd14_M_4_r2_SPPARKS_init'
        if secluded_AR:
            init_file = 'd14_M_4_r4_secAR_SPPARKS_init'  # 'd14_M_4_r2_noBckgrnd_SPPARKS_init' # 'd14_M_4_r4_SPPARKS_init'  # 'd14_M_4_r2_noBckgrnd_SPPARKS_init'  #'d14_M_4_r1_SPPARKS_init'

        if subregional:
            init_file = 'd14_M_4_r3_SPPARKS_init'

        if use_fixedT_ifng_init:
            init_file = 'd14_M_4_r4_fxdT_ifng_t1008_SPPARKS_init'

doing_lone_t_sim = False
if doing_lone_t_sim:
    init_file = 'loneTcell_SPPARKS_initifng_diff_test_init'
    input_file_to_copy = 'in.singleT_ifng_check'
    # input_file = 'in.singleT_ifng_check'

doing_ifng_pulse_sim = False
if doing_ifng_pulse_sim:
    init_file = 'd14_M_4_ifng_pulse324000_SPPARKS_init_diff_test'
    input_file_to_copy = 'in.singleT_ifng_check'


coarse_graining = False
if coarse_graining:
    input_file = "in.cg"  #  # "in.macInhb_0" "in.const_in_cond_simple"
    bash_file = "run_sim_cg.sh"
    redsky_file = 'spk_redsky_cg'
    # init_file = 'R06RD-R06RD-R06RD-R06RD_lattConst_20um_Patchwork_SPPARKS_init'
    init_folder = "SPPARKS_IMC_init_files_20um_23-01-2024"
    init_file = '06RD_responder_SPPARKS_init'

print("Slide: " + init_file)

########################## Homogenize macrophage distribution for 04RD and 01RD; comment out for regular run_samples ###
# os.system("mkdir " + newDirName)
# # dirAdd = "./" + newDirName + "/"
# os.system("cp " + init_file + " " + newDirName)
#
# os.chdir(newDirName)
#
# initWriter.homogenizePop(init_file, 5, True)
#
# sys.exit("Done")
########################################################################################################################

# os.system("module load GCC/7.3.0-2.30 OpenMPI/3.1.1")
os.system("cp " + input_file_to_copy + ' ' + input_file)
# dirAdd = "./" + newDirName + "/"
os.system("mkdir " + newDirName)
# print("made directory: "+str(newDirName))
# os.system("cp " + bash_file + " " + newDirName + "; cp " + redsky_file + " " + newDirName + "; cp " + input_file + " " + newDirName + "; cp ./" + init_folder + "/" + init_file + " " + newDirName + ";")
os.system("cp " + bash_file + " " + newDirName + "; cp " + redsky_file + " " + newDirName + "; cp " + input_file_to_copy + " " + newDirName + "/" + input_file + "; cp ./" + init_folder + "/" + init_file + " " + newDirName + ";")

os.chdir(newDirName)

make_female_rates = None
if sex_rates:
    # 0.9 is the rate at which we maintain AR- cells in the male, 0.7 is the fraction decrease of T cell recruitment from females to males
    # we multiply by 1/0.9 to increase the recruitment rate of AR- T cells such that total recruitment between males and females match
    # then we multiply by 1/0.7 to make it so that recruitment in the males is 70% of that in the females
    make_female_rates = 0.9 * 0.7

first_changes_to_input_file(fileName=input_file, create_dump_files=create_dump_files, make_female=make_female_rates, init_pos_file_name=init_file)

for i in range(num_Samples):
    # We choose a new seed and make a new sample directory
    if bool(i):
        # changeInputSeed(input_file, 1)
        changeInputFile(input_file, 1, False)  #, init_pos_file_name=init_file, create_dump_files=create_dump_files)
    else:
        changeInputFile(input_file, 0, True)   #, init_pos_file_name=init_file)# , create_dump_files=create_dump_files)

    print("job #"+str(i+1))
    os.system("mkdir sample"+str(i))

    sample_name = "sample"+str(i)
    os.system("cp " + init_file + " " + sample_name + ";")
    os.system("mv " + bash_file + " " + sample_name + "; cp " + redsky_file + " " + sample_name + "; cp " + input_file + " " + sample_name + ";")# cp " + init_file + " sample"+str(i)+";")

    os.chdir(sample_name)
    # now within sample directory
    # if bool(i):
    #     changeInputSeed(input_file, 1)
    os.system("sbatch "+bash_file)
    # Now that we have started the sample run, we leave the directory and take only the bash file with us (leave copy of in.file)
    os.system("mv " + bash_file + " .. ;")  # cp " + input_file + " .. ;")
    os.chdir("..")
# reset seed number... Not necessary
changeInputSeed(input_file, 1 - num_Samples)
print("done")
