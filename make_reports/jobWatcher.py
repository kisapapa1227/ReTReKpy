import os,time
import subprocess
from subprocess import PIPE

wk="/var/www/html/ReTReKpy/"

def chkJM(uid):
    comp=subprocess.run(['ps','aux'],capture_output=True,text=True)
    for i in comp.stdout.split("\n"):
        if 'jobMaster.py' in i and uid in i:
            return False
    return True

def restartJM(uid):
    comp=subprocess.Popen(['python3','make_reports/jobMaster.py','-restart','-uid',uid],cwd=wk,stdout=PIPE,stderr=PIPE)
    time.sleep(30)
#
# main
#
if os.path.isfile(wk+"jobMaster.lock"):
    with open(wk+"jobMaster.lock") as fp:
        ls=fp.readlines()
        for l in ls:
            uid=l.split("\n")[0]
#            with open(wk+"now.lock","w") as fpp:
#                fpp.write(uid+"\n");
            if chkJM(uid):
                restartJM(uid)
