import os,getpass,filecmp
from CleanUp import *

flag=os.getenv("Flag")

if flag == None:
    flag="SNU"

path_jobpre="/data1/"
if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    path_jobpre="/data2/"

ANALYZER_DIR= str(os.getenv("ANALYZER_DIR"))
ANALYZER_LOG= str(os.getenv("ANALYZER_LOG_PATH"))

if not ANALYZER_DIR == "None" :
	datadir = ANALYZER_DIR + "/data/"	
	if not (os.path.exists(datadir)):
		print "This is the first time running Analyzer in this location"
		print "Making data directory in $ANALYZER_DIR"
		os.system("mkdir " + datadir)
        
	outfiledir= ANALYZER_DIR +"/data/output/"
	lumifiledir= ANALYZER_DIR +"/data/Luminosity/"+ str(os.getenv("yeartag"))
	if not os.path.exists(ANALYZER_DIR +"/data/Luminosity/"):
		os.system("mkdir " +ANALYZER_DIR +"/data/Luminosity/")
	if not os.path.exists(lumifiledir):
		os.system("mkdir " + lumifiledir)
	btagfiledir = ANALYZER_DIR +"/data/BTag/"+ str(os.getenv("yeartag"))
	if not os.path.exists(ANALYZER_DIR +"/data/BTag/"):
		os.system("mkdir " +ANALYZER_DIR +"/data/BTag/")
	if not os.path.exists(btagfiledir):
		os.system("mkdir " + btagfiledir)

	if not (os.path.exists(outfiledir)):
		os.system("mkdir " + outfiledir)
		print "Making data/output directory in $ANALYZER_DIR"


	EightTeVdataOne=path_jobpre+"" + getpass.getuser() + "/"+flag+"_SKTreeOutput/"
	EightTeVdataTwo="/data2/" + getpass.getuser() + "/"+flag+"_SKTreeOutput/"
	 
	if os.path.exists(os.getenv("ANALYZER_DIR")+ "/nohup.out"):
		os.system("rm " +os.getenv("ANALYZER_DIR")+ "/nohup.out")


        CleanUpLogs(path_jobpre+flag+"Analyzer_rootfiles_for_analysis/CATAnalyzerStatistics/"+ getpass.getuser()+ "/")
        CleanUpJobLogs(ANALYZER_LOG)
	if os.getenv("HOSTNAME") == "cms.snu.ac.kr":
            CleanUpLogs(path_jobpre+"CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs("/data2/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs("/data7/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs("/data8/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs("/data8/DATA/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs("/data7/DATA/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
            CleanUpLogs(os.getenv("ANALYZER_BATCHLIB_PATH"))
            CleanUpLogs(EightTeVdataOne)
            CleanUpLogs(EightTeVdataTwo)
        else:
            CleanUpLogs("/data4/CAT_SKTreeOutput/" + getpass.getuser()+ "/")
        localfiledir = os.getenv("ANALYZER_FILE_DIR")
	datasetfiledir = os.getenv("ANALYZER_DATASETFILE_DIR")
	txtfiledir = os.getenv("ANALYZER_DIR")+ "/"+flag+"Run/txt/"
	cltxtfiledir = os.getenv("ANALYZER_DIR")+ "/"+flag+"Run/txt/Cluster/"
	seldir =os.getenv("ANALYZER_DIR")+  "/SNUConfig/SelectionConfig/"
	os.system("cp " + localfiledir + "/Luminosity/triggers_catversion_"+str(os.getenv("CATVERSION"))+"* "  + lumifiledir)
	os.system("cp " + localfiledir + "/Luminosity/lumi_catversion_"+str(os.getenv("CATVERSION"))+".txt "  + lumifiledir)
	os.system("cp " + datasetfiledir + "/list_all_mc_"+str(os.getenv("CATVERSION"))+".sh " + txtfiledir)
        if os.getenv("HOSTNAME") == "cms.snu.ac.kr":
            list_sel= ["muons","electrons","jets","fatjets"]
            for x in list_sel:
                if not filecmp.cmp(localfiledir + "/Selection/"+str(x)+".sel",seldir+"/"+str(x)+".sel"):
                    print "#"*50
                    print "File "+str(x)+".sel is out of date...... Updating "
                    print "Differences are:"

                    os.system("diff " + localfiledir + "/Selection/"+str(x)+".sel " + seldir+"/"+str(x)+".sel")
                    os.system("cp " + localfiledir + "/Selection/"+str(x)+".sel " + seldir)
                    print "#"*50
        else:
            os.system("rm " + seldir  +"/electrons.sel") 
            while not os.path.exists(seldir  +"/electrons.sel"):
                os.system("scp 147.47.242.42:/data1/"+flag+"Analyzer_rootfiles_for_analysis/CATAnalysis2016/Selection/*.sel  " + seldir)
                
	#os.system("cp " + localfiledir + "/*.csv " + btagfiledir)
	#os.system("source " +  os.getenv("LQANALYZER_DIR") + "/bin/IncludePrivateSamples.sh")
else:
	print "Area is not setup. Cannot make directories needed for analysis"

