import os
from os import listdir
from os import walk
import time,os.path

def FetchIncompleteRuns(folder,BadRuns=[]):
    runs = []
    for (dirpath, dirnames, filenames) in walk(folder):     
        for file in filenames :
		if( os.path.getsize("%s/%s"%(folder,file)) > 512 and ("raw" in file)) :
			tmp = int(''.join(x for x in file if x.isdigit()))
			if tmp not in BadRuns : 
				runs.append(tmp)

        break
    return runs


def GetRunNumber(file) : 
    return int(''.join(x for x in file if x.isdigit()))
	

runs = FetchIncompleteRuns("/afs/cern.ch/user/m/mbenoit/eos/user/m/mbenoit/RawData/DESY_TB_DATA_December2013")
print  runs

#runs = [2001]

queue = "8nh"
batch_folder = "/afs/cern.ch/user/m/mbenoit/batch"
os.system("mkdir %s"%batch_folder)
os.system("rm -fr  %s/*"%batch_folder)
os.system("mkdir %s/launch"%batch_folder)
launch_folder = "%s/launch"%batch_folder
log_folder = "/afs/cern.ch/user/m/mbenoit/logs"

batch = []



for run in runs :
	#CreateConfigTelAlone(run)
	#CreateConfig(run)
	
	filename="%s/Batch_Run%i.sh"%(launch_folder,run)
	batch.append(filename)
	f=open(filename,'w')						
	
	f.write("source  /afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/setup_CERN_ROOT6_gcc48_PYTHON2.7.sh \n")
	f.write("source  /afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/build_env.sh \n")	
	
	f.write("source ~/setup_EOS.sh \n")
		
		
	#File structure for castor 
	f.write("mkdir output \n")
	f.write("mkdir output/db \n")	
	f.write("mkdir output/logs \n")	
	f.write("mkdir output/results \n")
	f.write("mkdir output/histo \n")
	f.write("mkdir output/lcio-raw \n")
	f.write("mkdir output/raw \n")
	
	f.write("xrdcp -f root://eospublic.cern.ch//eos/user/m/mbenoit/RawData/DESY_TB_DATA_December2013/run%06i.raw output/raw \n"%(run))
	f.write("ls \n")	
	
	f.write("/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/jobsub.py -c /afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/examples/clic_timepix/configDecember2013_EOS.cfg converter %i \n"%(run))	
	f.write("/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/jobsub.py -c /afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/examples/clic_timepix/configDecember2013_EOS.cfg clusearch %i \n"%(run))	
	f.write("/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/jobsub.py -c /afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/examples/clic_timepix/configDecember2013_EOS.cfg hitmaker %i \n"%(run))	
	f.write("/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/jobsub.py -c /afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/examples/clic_timepix/configDecember2013_EOS.cfg align %i \n"%(run))	
	f.write("/afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/jobsub.py -c /afs/cern.ch/eng/clic/TBData/software/ROOT6_gcc48_python2.7/ILCSoft/v01-17-05/Eutelescope/trunk/jobsub/examples/clic_timepix/configDecember2013_EOS.cfg fitter %i \n"%(run))	
	
	
	
	for folder in ["db","histo","logs","results","lcio-raw"] :	
		f.write("xrdcp -fr output/%s root://eospublic.cern.ch//eos/user/m/mbenoit/RawData/DESY_TB_DATA_December2013_results/%s \n"%(folder,folder))			
	#f.write("xrdcp -fr output/histo/tbtrackrun%06i.root root://eospublic.cern.ch//eos/user/m/mbenoit/RawData/DESY_TB_DATA_December2013_results/tbtrack/tbtrackrun%06i.root \n"%(run,run))
	
		

	f.write("ls \n")		
		
	os.system("chmod u+x %s"%filename)
	f.close()
	print filename

	

for job in batch :
    os.system("cd ~/batch/launch")
    run = GetRunNumber(job)
    log = "%s/Run%06i"%(log_folder,run)
    os.system("mkdir %s"%log)
    os.system("bsub -o %s/STDOUT -q %s %s"%(log,queue,job))
