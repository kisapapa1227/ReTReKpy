import os,sys,time
from datetime import datetime,timezone, timedelta
import sqlite3
import asyncio,subprocess
from rdkit import Chem
from subprocess import PIPE

wk="/var/www/html/ReTReKpy/"
pb="/var/www/html/public/images/"
db="jobMaster.db"

def logs(mes):
    with open(pb+"report/jobMaster.log","a") as fp:
        fp.write(mes+"\n")

def chkJobMaster(cur,conn):
    i="integer";t="text";d="double"
    cur.execute(f'create table if not exists master(uid {t},max_loop {i},running_job {i},db_type {t},exp {d},abort {i});')
    cur.execute(f'create table if not exists slave(uid {t},id {i},smiles {t},name {t},found {i},route {i},done {t},start {d},end {d});')
    conn.commit()

def register(cur,uid,db_type,max_loop):
    cur.execute(f'delete from master where uid="{uid}"')
    cur.execute(f'delete from slave where uid="{uid}"')
    cur.execute(f'insert into master(uid,max_loop,running_job,db_type) values("{uid}","{max_loop}",0,"{db_type}");')
    for i in range(max_loop):
        smiles,name,done,e=analyze(wk+uid+"_async"+str(i+1)+".sh");
        route=int(e)
        cur.execute(f'insert into slave(uid,id,smiles,name,found,route,done) values ("{uid}","{i+1}","{smiles}","{name}",0,{route},"{done}");')
    # about cron:komai
    # lock file
    fn=wk+"jobMaster.lock"
    if os.path.isfile(fn)==False:
        with open(fn,"w") as fp:
            fp.write(uid+"\n")
    else:
        flag=True
        with open(fn,"r") as fp:
            ls=fp.readlines()
            for l in ls:
                if uid in l:
                    flag=False
        if flag:
            with open(fn,"a") as fp:
                fp.write(uid+"\n")
    # cronetab
    flag=True
    comp=subprocess.run(['crontab','-l'],capture_output=True,text=True)
    print(fn)
    for i in comp.stdout.split("\n"):
        if "jobWatcher.py" in i:
            flag=False
    if flag:
        sh=wk+"addCron.sh"
        with open(sh,"w") as fp:
            fp.write("#!/bin/sh\n#\n#\n")
            fp.write('(crontab -l;echo "*/2 * * * * /usr/bin/python3 /var/www/html/ReTReKpy/make_reports/jobWatcher.py")|crontab -')
        comp=subprocess.run(['sh',sh],capture_output=True,text=True)

def analyze (fn):
    smiles='happy_smiles';done='ready';loop=0
    with open(fn) as fp:
        for l in fp:
            if 'smiles:' in l:
                smiles=l.split("smiles:")[1].split('\n')[0]
            if 'abort:' in l:
                done=l.split("abort:")[1].split('\n')[0]
                name='abort'
            if 'python3' in l:
                es=l.split('\n')[0].split(' ')
                smiles="".join(es[2].split("'"))
                smiles="".join(es[2].split("'"))
#                loop=int("".join(es[3].split("'").split('"')))
                loop=int("".join(es[3].split("'")))
                name="".join(es[13].split("'"))
    return smiles,name,done,loop

def jobStart(id,uid,conn,retry=False):
    cur=conn.cursor()
    where=f' where uid="{uid}" and id={id+1}'
    route=0

    if retry==False:
        ret=cur.execute(f'update master set abort=0 where uid="{uid}"')
        conn.commit()
    mes="uid:"+uid+" id:"+str(id)
    ret=cur.execute(f'select done from slave'+ where)
    stat=ret.fetchall()[0][0]
    if stat!='ready' and stat!='Searching':
        return

    sh=wk+uid+"_async"+str(id+1)+".sh"
    ret=cur.execute(f'update master set running_job={id+1};')
    ut=time.time()
    ret=cur.execute(f'update slave set done="Searching", start={ut} '+where)
    ret=cur.execute(f'select name from slave '+where)
    name=ret.fetchall()[0][0]
# read substance.txt
    readTxtForCountingRoute(uid,name,init=True)
#    print(name)
    conn.commit()
    process=subprocess.Popen(["sh",sh],stdout=PIPE,stderr=PIPE)
#    print(sh+" starting")
    running=True
    while running: 
        time.sleep(5)
        running,x=askProcess(uid,name,ut)
        print("Route:"+str(x))
        if route!=x:
            route=x
            abort,nn=updateSlave(uid,route,where,cur,running)
            print("abort : uppdateSlave:"+str(nn)+":"+str(running)+":"+str(abort))
            conn.commit()
#
    if abort:
        if nn<5:
            print("restart process")
            jobStart(id,uid,conn,retry=True)
        else:
            ret=cur.execute(f'update slave set done="abort" '+where)

def getTid(db):
    comp=subprocess.run(['python3',wk+"make_reports/readDb.py","-db_list","-database",db],capture_output=True,text=True)

#    for n,i in enumerate(comp.stdout.split("###")):
#        logs(str(n+1)+"------------->"+i)
    ds=comp.stdout.split("###")[-2].split('##')[0]
    logs(ds)
    return ds

def addDb(uid,ssh,db_type,name):
    sh=ssh.split('_')[1]

    if db_type=="com":
        db="sList.db"
    else:
        db="sList"+uid+".db"

    comp=subprocess.run(['python3',wk+"make_reports/addDb.py","-u",uid,"-s",sh,"-d",pb,"-database",db],capture_output=True,text=True)
    tid=getTid(db)

    src=pb+uid+"/report/"+name+".pdf"
    if db_type=="com":
        dst=pb+"/report/"+tid+".pdf"
    else:
        dst=pb+uid+"/report/"+tid+".pdf"

    logs("src:"+src)
    logs("dst:"+dst)
    comp=subprocess.run(['mv',src,dst],capture_output=True,text=True)
    for i in comp.stdout.split("\n"):
        logs(i)

def updateSlave(uid,r,where,cur,fin):
    now=time.time()

    ret=cur.execute(f'select start, route from slave'+where)
    pt=ret.fetchall()
    ut=pt[0][0];route=pt[0][1]
    exp=ut+(now-ut)/r*route;
    ret=cur.execute(f'update slave set found={r}'+where)
    ret=cur.execute(f'update master set exp={exp} where uid="{uid}"')
    if r==route:
        ret=cur.execute(f'update slave set done="Completed", end={now}'+where)
# komai need:uid,db_type,sh
#    sh=wk+uid+"_async"+str(id+1)+".sh"
#    where=f' where uid="uid" and id={id+1}'
        x=where;
        ret=cur.execute(f'select db_type from master where uid="{uid}"')
        db_type=ret.fetchall()[0][0]
        sid=x.split("id=")[2]
        sh=wk+uid+"_async"+sid+".sh"
        ret=cur.execute(f'select name from slave '+where)
        name=ret.fetchall()[0][0]
        addDb(uid,sh,db_type,name)
    if fin==False and r!=route:
        ret=cur.execute(f'select abort from master where uid="{uid}"')
        x=ret.fetchall()[0][0]
        ret=cur.execute(f'update slave set done="ready"'+where)
        ret=cur.execute(f'update master set abort={x+1} where uid="{uid}"')
        return True,x
    return False,0

def askProcess(uid,name,ut):
    comp=subprocess.run(['ps','aux'],capture_output=True,text=True)
    n=readTxtForCountingRoute(uid,name)
    flag=False
    for i in comp.stdout.split("\n"):
#        print(i)
        if 'exe.py' in i and uid in i:
            flag=True
    return flag, n

def readTxtForCountingRoute(uid,name,init=False):
    n=0
    fi=pb+uid+"/"+name+".txt"
    if os.path.isfile(fi)==False:
        return n
    if init:
        os.remove(fi)
        return n

    with open(fi) as fp:
        for l in fp:
            if 'Route:' in l:
                n=int(l.split('Route:')[1].split('\n')[0])
#                print(n)
    return n

def askList(cur,uid):
    ret=cur.execute(f'select max_loop from master where uid="{uid}"')
    mes="";
    for id in range(ret.fetchall()[0][0]):
        ret=cur.execute(f'select smiles,done from slave where uid="{uid}" and id="{id+1}"')
        d=ret.fetchall()
        if id!=0:
            mes+="\n"
        mes+=d[0][0]+" "+d[0][1]
    return mes

def askRunningId(cur,uid):
    ret=cur.execute(f'select running_job from master where uid="{uid}"')
    rrr=ret.fetchall()
    return rrr[0][0]

def askProgress(cur,uid):
    now=time.time()
    ret=cur.execute(f'select running_job,max_loop from master where uid="{uid}"')
    for d in ret:
        id,max=d[0],d[1]
    ret=cur.execute(f'select start,found,route from slave where uid="{uid}" and id="{id}"')
    JST=timezone(timedelta(hours=+9),'JST')
    for d in ret:
        ut=d[0];route=d[1];max_route=d[2]
    dt=datetime.fromtimestamp(now,JST)
    mess=str(route)+"/"+str(max_route)+";"+dt.strftime('%Y/%m/%d %H:%M:%S')+"(日本標準時間)"

#    dt=datetime.datetime.fromtimestamp(ut,datetime.timezone.utc)

#    print(dt.strftime('%Y/%m/%d %H:%M:%S'))
#    print("running id:"+str(id)+" route:"+str(route))
    mess+=";"+myFormatTime(now-ut)
    ret=cur.execute(f'select exp from master where uid="{uid}"')
    exp=ret.fetchall()[0][0]
#    dt=datetime.fromtimestamp(ut,JST)
#    print("will finish at "+dt.strftime('%Y/%m/%d %H:%M:%S'))
    if route<2 or route==max_route:
        mess+=";--:--:--"
    else:
        mess+=";"+myFormatTime(exp-now)
    return mess

def myFormatTime(dt):
    tt=int(dt)
    return str(tt//3600).zfill(2)+":"+str((tt//60)%60).zfill(2)+":"+str(tt&60).zfill(2)

def setMaxLoop(uid,cur):
    ret=cur.execute(f'select max_loop from master where uid="{uid}"')
    return ret.fetchall()[0][0]

def closeJob(uid):
   #komai 
    fn=wk+"jobMaster.lock"
    uids=[]
    if os.path.isfile(fn)==False:
        with open(fn,"r") as fp:
            ls=fp.readlines()
            for l in ls:
                if not uid in l and len(l)>2:
                    uids.append(l)
    if len(uids)==0:
        os.remove(fn)
        sh=wk+"addCron.sh"
        with open(sh,"w") as fp:
            fp.write("#!/bin/sh\n#\n#\n")
            fp.write('(crontab -l;grep -v jobWatcher.py")|crontab -')
        comp=subprocess.run(['sh',sh],capture_output=True,text=True)
        return

    with open(fn,"w") as fp:
        for u in uids:
            fp.write(u+"\n")

########################
#main
########################

skip=False
db_type="com"

conn=sqlite3.connect(wk+db)
cur=conn.cursor()

chkJobMaster(cur,conn)

stage=""
pid=-1;
max_loop=-1
for n,op in enumerate(sys.argv):
    if skip:
        skip=False;continue
    match op:
        case "-start":
            stage=sys.argv[n]
        case "-ask_db_type":
            stage=sys.argv[n]
        case "-ask":
            stage=sys.argv[n]
        case "-restart":
            stage=sys.argv[n]
        case "-stat":
            stage=sys.argv[n]
        case "-uid":
            uid=sys.argv[n+1]
            skip=True
        case "-pid":
            pid=int(sys.argv[n+1])
            skip=True
        case "-db_type":
            db_type=sys.argv[n+1]
            skip=True
        case "-max_loop":
            max_loop=int(sys.argv[n+1])
            skip=True

match stage:
    case "-restart":
        max_loop=setMaxLoop(uid,cur)

match stage:
    case "-start":
        register(cur,uid,db_type,max_loop)
        conn.commit()
        for id in range(max_loop):
            jobStart(id,uid,conn)
        closeJob(uid)
    case "-restart":
        for id in range(max_loop):
            jobStart(id,uid,conn)
        closeJob(uid)
    case "-ask_db_type":
        ret=cur.execute(f'select db_type from master where uid="{uid}"')
        print(ret.fetchall()[0][0])
    case "-ask":
        id=askRunningId(cur,uid)
        if id==pid:
            mess="pid:"+str(id)+"###"+askProgress(cur,uid)
            print(mess)
        else:
            mess="pid:"+str(id)+"###"+askList(cur,uid)
            print(mess,end="")
