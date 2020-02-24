#!/bin/tcsh 

set ENV_FILE    = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/environment/1_setup.csh"
set ENV_FILE2   = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/environment/j2root.csh"
set GROOVY_FILE = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/groovy/monitor-kincorr.groovy"
#set OUTPUT_DIR  = "/volatile/clas12/dmriser/farm_out/elastic_jan24/data/03/"
#set DATA_DIR    = "/work/clas12/rg-a/trains/v16_v2/skim8_ep"
#set DATA_DIR = "/lustre/expphy/volatile/clas12/rg-a/production/recon/pass0/calib/v35/recon/005038/"

#mkdir -p $OUTPUT_DIR
source $ENV_FILE
#source $ENV_FILE2
#run-groovy $GROOVY_FILE $DATA_DIR/* 
#cp *.hipo $OUTPUT_DIR

# One job per node
run-groovy $GROOVY_FILE input.hipo 

