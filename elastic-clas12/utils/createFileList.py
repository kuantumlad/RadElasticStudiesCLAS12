import os
import glob
import sys

job_number=sys.argv[1]
out_fname=sys.argv[2]

input_dir="/volatile/clas12/osg/bclary/job_"+job_number+"/output/" #/v/lustre2/expphy/volatile/clas12/osg/bclary/job_"+job_number+"/output/"
#f_out = open(out_fname,'w')

n_file_to_merge=[]
cc=0
for _, dirnames, filenames in os.walk(input_dir):
    for dd in dirnames:                
        print(input_dir+"/"+dd+"/dst.hipo")
        n_file_to_merge.append(input_dir+"/"+dd+"/dst.hipo")
        cc+=1
        #f_out.write(input_dir+"/"+dd+"/dst.hipo" + " \n")
    
    print(cc)
    if cc > 9680: break

print("MERGING FILES FOR JOB ")
print("input path: %s " % input_dir)
print("job number: %s" % job_number)
print("output file name: %s" % out_fname )
print("number of files to merge: %d " % len(n_file_to_merge))

print(">> Merging files now")
out_dir="/lustre19/expphy/volatile/clas12/bclary/esepp/cooked/inbend_pass1/" #"/volatile/clas12/bclary/simulation_out/elastic/beam10/esepp/inbend_pass5/"
#list_of_files=" ".join(n_file_to_merge)

split_list =  [n_file_to_merge[x:x+1000] for x in xrange(0, len(n_file_to_merge), 1000)]

nset=0

for ss in split_list:
    list_of_files=" ".join(ss) 
    print(list_of_files)
    temp_out_fname=out_fname+'_set'+str(nset)+'.hipo'
    print('saving to file ' + out_dir+temp_out_fname)   
    os.system("hipo-utils -merge -o {}{} {}".format(out_dir,temp_out_fname,list_of_files) )
    
    nset+=1

