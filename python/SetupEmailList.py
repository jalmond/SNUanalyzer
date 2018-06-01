import os,sys

path_emailconfig = os.getenv("ANALYZER_MOD") +"emailconfig.txt"                                                                                                               

path_snuconfig= os.getenv("ANALYZER_DIR") + "/bin/snuconfig"

file_snuconfig= open(path_snuconfig,"r")
email_address=""
for line in file_snuconfig:
    if "email" in line:
        splitline = line.split()
        if len(splitline) == 3:
            email_address=splitline[2]
        else:
            print "ERROR in email in bin/snuconfig. Please Fix"
            sys.exit()
file_snuconfig.close()
            

email_in_list=False
file_emailconfig = open(path_emailconfig,"r")
copy_email_list=[]
for line in file_emailconfig:
    copy_email_list.append(line)
    if email_address:
        if email_address in line:
            email_in_list=True

file_emailconfig.close()

if not email_in_list:
    print "Adding email address: " + email_address + " to SNUAnalyzer emaillist"
    file_emailconfig = open(path_emailconfig,"w")
    for xline in copy_email_list:
        xline = xline.split()
        file_emailconfig.write(xline[0])
        file_emailconfig.write("\n")
    file_emailconfig.write(email_address)
    file_emailconfig.close()
