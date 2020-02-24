setenv JYPATH "/work/clas12/kenjo/j2root/target/*"
if (! $?LD_LIBRARY_PATH) then       
  setenv LD_LIBRARY_PATH /work/clas12/kenjo/j2root/lib
else
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/work/clas12/kenjo/j2root/lib
endif
