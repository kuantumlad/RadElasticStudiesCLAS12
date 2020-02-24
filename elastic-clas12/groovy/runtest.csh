#!/bin/tcsh 

set ENV_FILE    = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/environment/1_setup.csh"
#set ENV_FILE2   = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/environment/j2root.csh"
set GROOVY_FILE = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/groovy/monitor-kincorr.groovy"


setenv BASE_PATH /w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic
setenv JYPATH ${BASE_PATH}/analysis_code
setenv PYTHONPATH ${HIPO_TOOLS}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/work/clas12/kenjo/j2root/lib
#  Python3 and cython 
# setenv PATH /apps/anaconda/python3/anaconda201812/bin:${PATH}

# ROOT Newer Build 
# setenv ROOTSYS /work/clas12/tylern/software/root
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ROOTSYS}/lib
setenv PATH ${ROOTSYS}/bin:${PATH}
setenv PYTHONPATH ${PYTHONPATH}:${ROOTSYS}/lib

#mkdir -p $OUTPUT_DIR


source $ENV_FILE


#run-groovy $GROOVY_FILE $DATA_DIR/* 
#cp *.hipo $OUTPUT_DIR

# One job per node
run-groovy $GROOVY_FILE input.hipo 

