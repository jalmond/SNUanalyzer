import os,datetime,sys,filecmp
from datetime import timedelta

path_admin = os.getenv("ANALYZER_MOD")+"/config.txt"
flag=os.getenv("Flag")

path_jobpre="/data1/"
if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    path_jobpre="/data2/"

def checkLumiFile(backupdir,backup_datelist):
    

    nowtime = datetime.datetime.now()
    diff = datetime.timedelta(days=2)
    PAST = nowtime - diff
    date_now=nowtime.strftime("%d")
    date_past=PAST.strftime("%d")
    month_now=nowtime.strftime("%m")
    month_past=PAST.strftime("%m")

    if not  int(date_now)  % 2 == 0:
        sys.exit()

    snuversion=os.getenv("SNUVERSION")
    lumidir=os.getenv("ANALYZER_DATASETFILE_DIR")
    currentfile=lumidir+"/datasets_snu_SNU_mc_" + snuversion+".txt"
    
    backupfile = lumidir+"/BackUp/"+backup_datelist[len(backup_datelist)-1] + "/datasets_snu_SNU_mc_" + snuversion+".txt"
    if not os.path.exists( currentfile):
        print currentfile + " does not exist"

    if not os.path.exists( backupfile):
        print backupfile + " does not exist"
        return
    
        
    if not filecmp.cmp(currentfile,backupfile):
        print "There is a difference in " + currentfile + " and " + backupfile
        os.system("diff " + currentfile + " " + backupfile)
    
def makeBackUp(backupdir,copylist,backup_datelist):
    
    nowtime = datetime.datetime.now()
    diff = datetime.timedelta(days=2)
    PAST = nowtime - diff
    date_now=nowtime.strftime("%d")
    date_past=PAST.strftime("%d")
    month_now=nowtime.strftime("%m")
    month_past=PAST.strftime("%m")
    
    if  int(date_now)  % 2 == 0:
        print "Making backup for directory " + backupdir
        lastest_backup=backupdir+ nowtime.strftime("%d-%m-%y")+"/"
        if not os.path.exists(lastest_backup):
            os.system("mkdir " + lastest_backup)
        for xcopy in copylist:
            os.system("cp " + xcopy + lastest_backup)
        print "Backing up statfiles in " + lastest_backup
            
        if len(backup_datelist) > 10:
            remove_backup=backupdir+backup_datelist[0]
            if os.path.exists(remove_backup):
                os.system("rm -r " + remove_backup)
    else:
        print "Odd day: checking if previous day was backed up"
        diff1day = datetime.timedelta(days=1)
        PAST1DAY=nowtime - diff1day
        lastest_backup=backupdir+ PAST1DAY.strftime("%d-%m-%y")+"/"
        if not os.path.exists(lastest_backup):
            print "Making directory " + lastest_backup
            os.system("mkdir " + lastest_backup)
        for xcopy in copylist:
            os.system("cp " + xcopy + lastest_backup)
        print "Backing up statfiles in " + lastest_backup

        if len(backup_datelist) > 10:
            remove_backup=backupdir+backup_datelist[0]
            if os.path.exists(remove_backup):
                os.system("rm -r " + remove_backup)


def editconfig(backup_datelist):
    
    nowtime = datetime.datetime.now()
    date_now=nowtime.strftime("%d")
    if not  int(date_now)  % 2 == 0:
        diff1day = datetime.timedelta(days=1)
        PAST1DAY=nowtime - diff1day
        pastformat=PAST1DAY.strftime("%d-%m-%y")
        nobackup=False
        for x in backup_datelist:
            if x == pastformat:
                nobackup=True
        if nobackup:
            return
        else:
            diff1day = datetime.timedelta(days=1)
            PAST1DAY=nowtime - diff1day
            strpastday= PAST1DAY.strftime("%d-%m-%y")
            if len(backup_datelist) > 10:
                mod_file_admin = open(path_admin, "w")
                for xline in copy_admin_file:
                    if not backup_datelist[0] in xline:
                        mod_file_admin.write(xline)
                mod_file_admin.write("backup " + strpastday+"\n")
                mod_file_admin.close()
            else:
                mod_file_admin = open(path_admin, "w")
                for xline in copy_admin_file:
                    mod_file_admin.write(xline)
                mod_file_admin.write("backup " + strpastday+"\n")
                mod_file_admin.close()


    if len(backup_datelist) > 10:
        mod_file_admin = open(path_admin, "w")
        for xline in copy_admin_file:
            if not backup_datelist[0] in xline:
                mod_file_admin.write(xline)
        mod_file_admin.write("backup " + nowtime.strftime("%d-%m-%y")+"\n")
        mod_file_admin.close()        
    else:
        mod_file_admin = open(path_admin, "w")
        for xline in copy_admin_file:
            mod_file_admin.write(xline)
        mod_file_admin.write("backup " + nowtime.strftime("%d-%m-%y")+"\n")    
        mod_file_admin.close()




dobackup=False
backup_date=[]
file_admin = open(path_admin, "r")
copy_admin_file=[]
for line in file_admin:
    copy_admin_file.append(line)
    if "administrator" in line and os.getenv("USER") in line:
        dobackup=True
    if "backup" in line:
        splitline = line.split()
        if len(splitline) == 2:
            backup_date.append(splitline[1])
file_admin.close()

if dobackup:
    lumidir=os.getenv("ANALYZER_DATASETFILE_DIR")
    checkLumiFile(lumidir+"/BackUp/",backup_date)



nowtime2 = datetime.datetime.now()
for xbackup in backup_date:
    if nowtime2.strftime("%d-%m-%y") in xbackup:
        dobackup=False
        print "ADMIN: Backup already made today"
            
if dobackup:
    print "ADMIN: Making backup"
    
    flag=os.getenv("Flag")

    print flag
    copylist=[]
    copylist.append(path_jobpre+flag+"Analyzer_rootfiles_for_analysis/SNUAnalyzerStatistics/MasterFile_v* ")
    copylist.append(path_jobpre+flag+"Analyzer_rootfiles_for_analysis/SNUAnalyzerStatistics/JobSummary* ")
    makeBackUp(path_jobpre+flag+"Analyzer_rootfiles_for_analysis/SNUAnalyzerStatistics/BackUp/",copylist,backup_date)
    copylist2=[]
    lumidir=os.getenv("ANALYZER_DATASETFILE_DIR")
    copylist2.append(lumidir+"/data* ")

    makeBackUp(lumidir+"/BackUp/",copylist2,backup_date)
    
    
    ### change config file to switch dir names
    editconfig(backup_date)
