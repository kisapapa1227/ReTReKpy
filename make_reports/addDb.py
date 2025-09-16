import sys,os
import math
import subprocess
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Draw, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

from reportlab.graphics import renderPDF
from reportlab.pdfgen import canvas
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.lib.pagesizes import A4, portrait, landscape
from svglib.svglib import svg2rlg
from reportlab.lib.units import mm
from pptx import Presentation
from pptx.util import Inches, Pt, Cm
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import MSO_AUTO_SIZE
from pptx.dml.color import RGBColor
from PIL import Image
import shutil
import Levenshtein
import sqlite3
import datetime

_web=True
_web=False

db="sList.db"

if _web:
    dr='/var/www/html/public/images/tmp' # working directory
else:
    dr='tmp' # working directory

_magic_x=650


def getInfo(tg,lines):
    key='Route:'+str(tg)
    on=False
    ret=[]
    for l in lines:
        if key in l:
            on=True
            continue
        if on:
            if 'Route' in l:
                if len(ret)<1:
#                    print(l1)
                    ret.append([l1.split(',')[0]])
                return ret
            each=l.split(',')
            ret.append(each[:-1])
    return ret

def log(com):
    with open("/var/www/html/public/images/report/addDb.log","a") as fp:
        fp.write(com+"\n")

def getSmiles(tg,lines):
    key='Route:'+str(tg)
    on=False
    ret=[]
    for l in lines:
        if key in l:
            on=True
            continue
        if on:
            if 'Route' in l:
                if len(ret)<1:
#                    print(l1)
                    ret.append([l1.split(',')[0]])
                return ret
            if not '>' in l:
                l1=l
                continue
            each=l.split(',')
            ret.append(each[:-1])
    return ret

def Debug(c,head=None,sp=None):
    return
    for n,i in enumerate(c):
        if head:
            print(head,end=":")
        if sp:
            print(n,i[sp],end="<--")
        print(n,i)
#    exit()


def chkConnect(n,connect):
    if n<0: 
        return n
    n=connect[n][3]
    while connect[n][1][0]==-1:# ext step: no smiles
        n=connect[n][3]
    return n

def chkJoin(c,s):
    if len(c[2])<2:
        return False
    c[2].sort()
    p1=c[2][-1]
#    print("in checkJoin",p1,s[p1][0],s[p1][1],"<---")
#    print(c[2])
    if s[p1][0]=='':
        return False
    return True

def addCatal(s,p):
    pp=[]
    for m,ss in enumerate(s):
        if ss=='>':
            return pp
        if not m in p:
            pp.append(m)
    return pp

def addCatal2(s,p):
    pp=[]
    for m,ss in enumerate(s):
        if ss=='>':
            return pp
        if not ss in p:
            pp.append(m)
    return pp

def lenForCat(cat,fc):
    sp=0
    for i in cat:
        sp+=fc[i]
    return sp
    
def setBranchLength(connect,goal,fc):
    branch=[]
    right=0

    for n in goal:
        px=0
        while True:
            if not fc:
                sp=len(connect[n][-1])
                connect[n][4]=px=px-(sp+2)
            else:
                if n<len(fc):
                    sp=lenForCat(connect[n][-1],fc[n])+0.5
                    wd=fc[n][connect[n][1][0]]+1
                else:
                    wd=fc[-1]['cheat']
                    sp=2
#   komai
                # for starting material
                #
                connect[n][4]=px=px-(sp+wd)

            if len(connect[n][2])<1:
                if len(branch)<1:
                    break
                else:
                    ww=branch.pop(-1)
                    n,px=ww[0],ww[1]
            elif len(connect[n][2])>1:
                for i in connect[n][2][1:]:
                    branch.append([i,px])
                n=connect[n][2][0]
            else:
                n=connect[n][2][0]
            if right>sp:
                right=sp
    for c in connect:
        c[4]-=right

def initConnect(smiles,fc=False,mode=1):
    connect=[]
    drop={}
    for n,proc in enumerate(smiles):
        e=False
#        print("proc")
        for m, comp in enumerate(proc):
            if ">" in comp:
                e=True
                continue
            if e:
# connect=[smiles of output product,[#ID,], [connected from], connect to, x, y]
                connect.append([comp,[m],[],-1,-1,0,0])
#
# up to here: make template for each steps (arrow)
#
    Debug(connect,head='A')

    for n, prev in enumerate(smiles):
        for m, proc in enumerate(smiles):
            if n !=m and connect[n][0] in proc:
#                print("n",n,m)
                connect[n][3]=m
                connect[m][2].append(n)
#
# up to here: find connection (next arrow)
#
    Debug(connect,head='B')
    goal=[]
    l=len(smiles)
    #komai
    if len(connect)<1:# case of no route: V2.2
        return [],0
    for i in range(l):
        if connect[i][3]==-1:
            goal.append(i)

    if len(goal)>1:
        print("warning multiple goal")

###################################
    for n,c in enumerate(connect):
        while chkJoin(c,connect):
            p1=c[2].pop(-1)
            drop[n]=connect[p1][0]
            p2=c[3]
            p3=len(connect)
            connect[p1][3]=p3
            connect.append(['',[-1],[p1],p2,-1,0,0])
            if p2<0:
                goal.append(p3)
            else:
                connect[p2][2].append(p3)
#
# up to here: make extra stepping stone
#
# p1 -> [,], p2=To, p3=last(inserted)
# in case of p2=-1 --> trouble
#
    Debug(connect,head='C')
#
# up to here: find head : open edge
#
    start=[]
    branch=[]
    connect[n][4]=connect[n][5]=0;
    depth=1
    for n in goal:
        depth-=1
        while True:
            connect[n][5]=depth
            if len(connect[n][2])<1:
                start.append(n)
                if len(branch)<1:
                    break
                else:
                    m=branch.pop(-1)
                    depth=depth-1
                    #n=connect[m][3]
            else:
                m=connect[n][2][0]
                if len(connect[n][2])>1:
                    branch=branch+connect[n][2][1:]
#        print(n,"->",m)
            n=m
#
# up to here: define tree structure
#
    for n in range(len(connect)):
        connect[n].append([])
    Debug(connect,head='D')

#    Debug(connect,head='A')
#    print("keys",drop)
#    print("keys",list(drop.keys()))
# connect=['smiles',[#ID,], [from], to, x, y, start_chem, [catalyser]]
    for n,s in enumerate(smiles):
        c=connect[n]
        p=[]
        if len(c[2])<1:
            init=whichIsFirstMaterial(s,mode)
            c[6]=init
            p.append(init)# will be modified
#            p.append(0)# will be modified
            pp=addCatal(s,p)
            c[7]+=pp
#            print("upper",n)
        else:#
            d=[]
#            print("lower",n)
            for f in c[2]:#
                d.append(connect[f][0])
            if n in list(drop.keys()):
                d.append(drop[n])
            pp=addCatal2(s,d)
            c[7]+=pp
    Debug(connect,head='E')
#
# up to here: register (catlizers) #id 
#
    if fc:
        if len(connect[i][2])<1:
            fc[-1]['cheat']=0
        else:
            i1=connect[i][2][0]
            i2=connect[i1][1][0]
            fc[-1]['cheat']=fc[i1][i2]

    setBranchLength(connect,goal,fc)
#
# inverse ->> n: upstream, m: downstream (goal to sta]rts)
#
#            connect[m][4]=connect[n][4]-(sp+2)
    Debug(connect,head='F')

    return connect,start
#        smile2sgv(comp,"e"+str(n+1)+"x"+str(m+1)+".svg")

def getAllRoutes(all):
    ret=[]
    key='Route:'
    for f in all:
        if key in f:
#            print("--",f,"<--",f.split(key))
            num=int(f.split('\n')[0].split(key)[1])
            ret.append(num)

    return ret

def whichIsFirstMaterial(s,mode):
    #komai
    r=0;maxWeight=0.0
    result=s[-1]
    for i, smiles in enumerate(s):
        if ">" in smiles:
            break
        if mode==3:# number of atom
            mol=Chem.MolFromSmiles(smiles)
            mw=mol.GetNumAtoms()
        elif mode==2:# weigh of atom
#            mw=1-Levenshtein.ratio(smiles,result)
            mol=Chem.MolFromSmiles(smiles)
            mw=rdMolDescriptors.CalcExactMolWt(mol)
        else:# similarity between source and product
            mw=Levenshtein.jaro_winkler(result,smiles)
#        print(smiles,result,mw)
#        mol=Chem.MolFromSmiles(smiles)
#        mw=rdMolDescriptors.CalcExactMolWt(mol)
        if mw>maxWeight:
            maxWeight=mw
            r=i
    return r

def isIncluded(s,p):
    flag=False

    for i in p:
        if len(i)<len(s):
            continue
        flag=True
        for j in s:
            if not j in i:
                flag=False
                break
        if flag:
            return True
    return False

def getMatrix(rt,all):# remove redundant routes
    ret={}

    tg="smiles"
    ret[tg],ppt=getMatrixAgent(rt,all,tg)
    tg="products"
    ret[tg],_=getMatrixAgent(rt,all,tg)
    return ret,ppt

def isInRef(s,ref):
    if s in ref:
        return ref.index(s),False
    n=len(ref)
    ref.append(s)
    return n,True

def getMatrixAgent(rt,all,typ="smiles"):# remove redundant routes

    _smiles={}
    ref=[]
    sAll={}
    sTotalMatch={}
    sSubMatch={}
    sDrop=[]
    ret={}

    routes=rt.copy()

    for r in rt:
        _smiles[r]=getSmiles(r,all)

# small trick
    l=0
    pt=0
    for r in reversed(rt):
        ss=_smiles[r]
        if l<len(ss):
            l=len(ss)
            pt=r

#    print("len",len(_smiles),l,pt)
    if l<1:
        return False

#    print("len",_smiles)
    for ss in reversed(_smiles[pt]):
        if typ=="smiles":
#            print("1")
            ref.append(ss)
        else:
#            print("2",ss)
            ref.append(ss[-1])
#        print("------------------------------ref\n",ss)

    if typ!="smiles":
        ss=_smiles[pt][0][0]
        ref.append(ss)

#    print("Fist Ref",ref,sim)

    for r in reversed(rt):
#        print(r,end="->")
        m=[]
        if typ=="smiles":
            for ss in reversed(_smiles[r]):
                n,_=isInRef(ss,ref)
                m.insert(0,n)
        else:
            for ss in reversed(_smiles[r]):
                n,f=isInRef(ss[-1],ref)
                m.insert(0,n)
            if len(_smiles[r])<1:
                sAll[r]=[]
                continue

            s=_smiles[r][0][0]
            n,f=isInRef(s,ref)
            m.insert(0,n)
        sAll[r]=m
        
    print("type",typ)
    print("the trick",pt,sAll[pt])
    print(sAll)

    ppt=pt
    for r in rt:
        if not sAll[r] in sTotalMatch.values():
            sTotalMatch[r]=sAll[r]
            if sAll[r]==sAll[pt]:
                ppt=r

    ret["total"]=sTotalMatch

    sSubMatch=sTotalMatch.copy()

    for r in sTotalMatch:
        flag=False
        if not r in sSubMatch:
            continue
        for p in sTotalMatch:
            if p == r:
                continue
            if p not in sSubMatch:
                continue
            if len(sTotalMatch[r])>len(sTotalMatch[p]):
                continue
            flag=True
            for i in sTotalMatch[r]:
                if not flag:
                    break
                if not i in sTotalMatch[p]:
                    flag=False
            if flag:
                sSubMatch.pop(r)
                break

    ret["sub"]=sSubMatch

#        if type!="smiles":
#            print("hot now",r,m)

#    print("Total")
#    for r in sTotalMatch.keys():
#        print(r,sTotalMatch[r])

    if typ=="smiles":
        return ret,ppt

#    print("this is submatch",sSubMatch)
    
    return ret,ppt

def getLeftTab():
    return 2

def setWds(c,fc):

    for n in range(len(c)):
        if n<len(fc):
            if len(c[n][2])<1:# starting point
                tg=-1
                for i in fc[n].keys():
                    if tg>0:
                        continue
                    if not i in c[n][-1]:
#                        print("keys",i,fc[n][i])
                        tg=fc[n][i]
                if tg<0:
                    fc[n]['s']=getLeftTab() # depend on
                else:
                    fc[n]['s']=tg # depend on
            else:
                p=c[n][2][0]#connecting node
                if p<0:# new tab
                    fc[n]['s']=getLeftTab() # depend on
                else:
                    q=c[p][1][0]# this is the target
                    fc[n]['s']=fc[p][q]
        else:
            fc.append({'s':2})

def isInIt(s,d):

    for i in s:
        if not i in d:
            return False
    return True

def registerBranch(p,rt):
    ret=[]
    for r in rt:
        if isInIt(p,rt[r]):
            ret.append(r)
    return ret

def checkIn(pp,rt):
    ret=[]
    if type(pp)==int:
        return checkIn2(pp,rt)

    for r in rt:
        if checkIn3(pp,rt[r]):
            ret.append(r)
#    print("checkIn",pp,ret)
    return ret

def checkIn3(p1,p2):

    ll=len(p1)
#    print("xxx",ll,p1,p2)

    for i in range(len(p2)-ll+1):
        flag=True
        for j in range(ll):
#            print(p2[i+j],"vs",p1[j])
            if p2[i+j]!=p1[j]:
                flag=False
        if flag:
            return flag
    return False

def checkIn2(p,rt):
    ret=[]
    for r in rt:
        if p in rt[r]:
            ret.append(r)
    return ret

def branchCheck(s,d):

    for i,j in enumerate(d):
#        print("<--",i,s,j,j==s)
        if j==s:
#            print("--->",i,j,s,index[i])
#            return index[i]
            return i
    return False

def deepArrange(routes):
    n=0
    p=[]
    branch=[]
    for r in routes.values():
        m=max(r)
        if n<m:
            n=m

    p=[n]
    n-=1
    while n>-1:
        p2=checkIn(p,routes)
        p1=checkIn(p+[n],routes)
        p3=checkIn([n],routes)
        if p1!=p2 or p1!=p3:
            branch.append(p)
            p=[n]
        else:
            p.append(n)
        n-=1
    branch.append(p)
    print(branch)

    ok={}
    for r in routes:
        pp=[]
        rr=[]
        for s in routes[r]:
            pp.append(s)
            key=branchCheck(pp,branch)
#            print(pp,branch,key)
#            print("->",pp,key,"in",routes[r])
            if not key is False:
                rr.append(key)
                pp=[]
        ok[r]=rr

    for i in ok:
        print("Route",i,end=":")
        for j in ok[i]:
            print(branch[j],end="->")
        print()
    return [ok,branch]

def dicToString(ss):
    r=""
    for key, value in ss.items():
        r+=str(key)+":"
        for k in value:
            r+=str(k)+","
        r+=";"
    return r

def listToString(ss):
    r=""
    for s in ss:
        ll=len(s)
        for n,l in enumerate(s):
            if type(l)==list:
                for i in l:
                    r+=str(i)+","
            else:
                r+=str(l)
            r+=";"
        r+="##"
    return r

##############################################
#################### main ####################
##############################################

input_dir="./work"
script='async.sh'
user='kisa'

com='-u user -d input_dir -s script.sh'
ops={'-h':com}
skip=False
outp=False
_product_only=True
_product_only=False
_include_subsets=False
_include_subsets=True
total_loop=0
lock=False
lock_file="/var/www/html/public/images/report/addDb.lock"

for n,op in enumerate(sys.argv):
    if skip:
        skip=False;continue
#    print(n,op)
    match op:
        case '-database':
            db=sys.argv[n+1]
            skip=True
        case '-d':
            input_dir=sys.argv[n+1]
            skip=True
        case '-n':
            total_loop=int(sys.argv[n+1])
            skip=True
        case '-u':
            user=sys.argv[n+1]
            skip=True
        case '-lock':
            lock=True
        case '-s':
            script=sys.argv[n+1]
            skip=True
        case '-done':
            outp=sys.argv[n+1]
            skip=True

#print(_routes)

input_file="."
with open(user+"_"+script,"r") as f:
    for l in f:
        if '#name' in l:
            uname=l.split('\n')[0].split('name=')[1]
        if '#email' in l:
            email=l.split('\n')[0].split('email=')[1]
        if 'python3' in l:
            es=l.split('\n')[0].split(' ')
#            for n,i in enumerate(es):
#                print(n,i)
#print("script:"+user+"_"+script)
smi="".join(es[2].split("'"))
loop=es[3]
factor=es[4]
substance="".join(es[13].split("'"))

chk=db+":"+substance
if os.path.isfile(lock_file):
    with open(lock_file,"r") as fp:
        for line in fp:
            if chk in line:
#                print("locked:skip "+chk)
                exit()

if lock!=False:
    with open(lock_file,"a") as fp:
        fp.write(chk+"\n")

opts=""
for i in es[5:11]:
    if i=='':
        opts+='Non,'
    else:
        opts+=i+','
#opts=es[5:11]

# output #5-10

input_file=substance+".txt"
with open(input_dir+"/"+user+"/"+input_file,"r") as f:
    all=f.readlines()
with open(input_dir+"/"+user+"/"+input_file+"_info","r") as f:
    all_info=f.readlines()

allRoutes=getAllRoutes(all)
routeMatrix,ppt=getMatrix(allRoutes,all)

tg="smiles"
if _product_only:
    tg="products"
parent=routeMatrix[tg]

if _include_subsets:
    show="total"
else:
    show="sub"
routes=list(parent[show].keys()) # v2.4 07152025

print("parent",parent,"<---------")
print("routes",routes)

conn=sqlite3.connect(db)
cur=conn.cursor()

cur.execute('create table if not exists searchList(id integer primary key, user text, uname text, email text, smiles text, cSmiles text, date text, substance text, loop integer, factors text, options text);')
cur.execute('create table if not exists parent(id integer, total text, sub text);')
conn.commit()

#smiles=getSmiles(routes[0],all)[-1][-1]
mol=Chem.MolFromSmiles(smi)
cSmiles=Chem.MolToSmiles(mol)
dd=datetime.datetime.now().strftime('%Y%m%d%H%M%S')

if total_loop!=0:
    loop=total_loop
#print(opts)
sql=f'insert into searchList(user, uname, email, smiles, cSmiles, date, substance, loop, factors, options) values("{user}","{uname}","{email}","{smi}","{cSmiles}",{dd},"{substance}",{loop},"{factor}","{opts}")'
cur.execute(sql)
#print(sql)
conn.commit()

sql=f'select id from "searchList" where date = "{dd}";'
ret=cur.execute(sql)
conn.commit()

total=dicToString(parent['total'])
sub=dicToString(parent['sub'])
id=ret.fetchall()[0][0]

#print("total",total)
#print(" sub" ,sub)
sql=f'insert into parent(id, total, sub) values({id},"{total}","{sub}")'
#print(sql)
cur.execute(sql)
conn.commit()

sTable=f"searchTable{id}"

cur.execute(f'create table if not exists "{sTable}" (route int, connect_list text, smiles_list text, info_list text);')
conn.commit()

#komai
for route in routes:
    info=getInfo(route,all_info)
    smiles=getSmiles(route,all)
    connect,start=initConnect(smiles)
    s1=listToString(connect)
    s2=listToString(smiles)
#    print("Route"+str(route)+";+connect",len(connect),connect)
#    print("Route"+str(route)+";+smiles",len(smiles),smiles)
    s3=listToString(info)
    sql=f'insert into {sTable} values({route},"{s1}","{s2}","{s3}");'
#    print("sql",sql)
    cur.execute(sql)
    conn.commit()

print("end")

print(f"substance:{substance}###")
conn.close()
if outp!=False:
    with open(outp,"a") as fp:
        fp.write("addDb_done\n")

if lock:
    lines=""
    with open(lock_file,"r") as fp:
        for line in fp:
            if not chk in line:
                lines=lines+line
    if lines=="":
        os.remove(lock_file)
    else:
        with open(lock_file,"w") as fp:
            fp.write(lines)
