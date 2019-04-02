import re
import os
import subprocess
import string
import time
import numpy as np
import argparse

ENVIRONMENT_PATH = "/home/blakei/env_rat-6.3.6.sh"

#def submit_job(shell_script_name):
#    subprocess.check_call("qsub -l cput=150:59:59 {0}".format(shell_script_name), shell = True)

# Pruning down root files to an ntuple with entries needed for coincidence tagging
def PruneforCoincTag():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
. ~/oxsx/bin/compile_with_rat.sh prune_for_coinctag.cpp
./prune_for_coinctag '/data/snoplus/blakei/antinu/mc/rootfiles/Bruce_LAB_5day_flux1000_s*' /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root""").format(ENVIRONMENT_PATH), shell = True)
#arguments:                           input root                                                                       output pruned ntuple

# produces a fake data set to be used as data input for a fit, a single entry _oxsx.ntuple.root file containing oscillated (or unoscillated) EPrompt, produced by applying coincidence tagging on the pruned ntuple made with above func.
def FakeDataEPrompt():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
. ~/oxsx/bin/compile_with_rat.sh prune_to_KE_EPrompt.cpp
./prune_to_KE_EPrompt /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root /data/snoplus/blakei/antinu/test/Bruce5yr1000flux 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297 0.0215 E1""").format(ENVIRONMENT_PATH), shell = True)
#arguments:                    input pruned ntp                                     output _oxsx fake data ntuple(no .root)   nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 s13 E1orKE

def FakeDataKE():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
. ~/oxsx/bin/compile_with_rat.sh prune_to_KE_EPrompt.cpp
./prune_to_KE_EPrompt /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root /data/snoplus/blakei/antinu/test/Bruce5yr1000flux 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297 0.0215 KE""").format(ENVIRONMENT_PATH), shell = True)
#arguments:                    input pruned ntp                                    output _oxsx fake data ntuple(no .root)   nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 s13 E1orKE

# NOT NEEDED FOR FIT, file produces root file w/ histograms showing results of applying coincidence tagging on a pruned ntuple
def OscEPrompt():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
. ~/oxsx/bin/compile_with_rat.sh OscEPrompt.cpp
./OscEPrompt /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root EPromptOut.root 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297 0.0215 MC""").format(ENVIRONMENT_PATH), shell = True)
#arguments:                    input pruned ntp                  output _oxsx fake data ntuple     nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 s13 MCorData looping method

def PruneFlatTree():
    subprocess.check_call(("""source {0}
cd ~/oxsx/util/
./prune_flat_tree /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root EPromptOut.root 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297 0.0215 MC""").format(ENVIRONMENT_PATH), shell = True)

# Currently input data for LH2D is generated within the example, to input data own data, comment out custom produced data part and uncomment out generic datantuple filling part
def LH2DEVindex():
    subprocess.check_call(("""source {0}
cd ~/oxsx/examples/
. ~/oxsx/bin/compile.sh LH2dEVindexPlot.cpp
./LH2dEVindexPlot /data/snoplus/blakei/antinu/test/PrunedBruce5yr1000flux.root USUALLYPUTDATA 3CAD.ratdb /data/snoplus/blakei/antinu/Test.root 1yrflux13CADConstrain 500 4000 850 1300 1. 8. 1.6 2.2 1e6 7.4e-5 0.297""").format(ENVIRONMENT_PATH), shell = True)
#arguments:        input pruned unosc ntp                                Data to Fit          reactoinfofile temp file for filling outputLH2dfilename(no .root) nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 (for output naming)


def OscFit():
    subprocess.check_call(("""source {0}
cd ~/oxsx/examples/
. ~/oxsx/bin/compile.sh OscillationFitTest.cpp
./OscillationFitTest /data/snoplus/blakei/antinu/test/Bruce5yr1000fluxKE_oxsx.root /data/snoplus/blakei/antinu/test/3CAD_36000evsKEds21_7.4e-05_ss12_0.297_ss13_0.0215_oxsx.root 3CAD.ratdb Test.root""").format(ENVIRONMENT_PATH), shell = True)
#arguments:        input pruned unosc ntp                 Data to Fit                   reactoinfofile temp file for filling outputLH2dfilename  nhit1min nhit1max nhit2min nhit2max E1min E1max E2min E2max deltaT d21 s12 (for output naming)


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Select func to run: (-r) \n prunefortag \n fake data \n fakedatake \n osceprompt \n lh2d \n oscfit")
    parser.add_argument("-r", dest="func", help="need to select which function to run", default=0)
    args = parser.parse_args()

    Func = str(args.func)
    if (Func.find("prunefortag") != -1):
        PruneforCoincTag()
    if (Func.find("fakedata") != -1):
        FakeDataEPrompt()
    if (Func.find("fakedatake") != -1):
        FakeDataKE()
    if (Func.find("osceprompt") != -1):
        OscEPrompt()
    if (Func.find("lh2d") != -1):
        LH2DEVindex()
    if (Func.find("oscfit") != -1):
        OscFit()
    
