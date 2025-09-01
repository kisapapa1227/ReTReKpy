import sys,os,time
import math
import subprocess
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.Scaffolds import MurckoScaffold
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
from dbTools import *
import re

def testForSvg(statAll):
    print(statAll[69])
    n=1
    fns=[];dst="/var/www/html/public/images/tmp/"
    for s in statAll[69][4]:
        if s!="":
            print(n,s)
            fn="now"+str(n)+".svg"
            smile2svg(s,fn)
            fns.append(fn)
            n+=1

    com=["montage","-tile",str(n-1)+"x1"]
    com=com+fns
    com.append("output.png")

    print(com)
    subprocess.Popen(com,cwd=dst)

def getPrm(cur,pids,sout):
    prm={};prmT={}
    for z in sout:
        prmT[z]={}
    for pid in pids:
        ret=cur.execute(f"""select * from searchList where id='{pid}';""")
        ret=ret.fetchall()
        if len(ret)<1:
            continue
        prm[pid]=ret[0]

    for k,v in prm.items():
#        print("*",k,v)
#        print(type(v))
        for z,x in zip(sout,v):
#            prmT[z].append(x)
            prmT[z][k]=x
#            print("?",z,x)

    return prm,prmT

def easyGetReaction(pid,cur):
    index={};newInfo={}
    sTable=f"searchTable{pid}"
    sql=f'select * from "{sTable}";'
    ans=cur.execute(sql).fetchall()
    for xx in reversed(ans):
        newInfo[int(xx[0])]=getLink(strToSmiles(xx[2]),index)

    newList,idx=getNewList(newInfo)
#    print("-------------->",newList)
#    print("-------------->",idx,len(idx))
    hint=getReactionSummary(newList,idx)
#    print("\n\n\n");
#    print(hint)
    return hint

def checkIt(ss,ps):

    for s in ss:
        if not s in ps:
            return False
    return True

def inIt(h0,h1):
#    print("----->",h0,h1)
    l=len(h0)

    for n,h in enumerate(h1):
        if l<h[0]:
            if checkIt(h0,h[1]):
                return n
    return -1

def dropSubset(hint):
    ret={};rm=[];h1=[]

    for n,h0 in enumerate(hint[0]):
        h1.append([len(h0[1]),h0[1]])

    for n,h0 in enumerate(h1):
        m=inIt(h0[1],h1)
        if m>-1:
            rm.append(n)
    return rm

def getMax(h1):
    ret=0
    for h in h1:
        mx=max(h[1])
        if ret<mx:
            ret=mx
    return ret 

def makeNextIndex(simple_route,routes):
    nextIndex={}
    for route in simple_route:
        nextIndex[route]=routes.index(route)+1
    return nextIndex

def getOpeningMaterial(smiles):
    s=smiles[0]

    mol=Chem.MolFromSmiles(s)
    core=MurckoScaffold.GetScaffoldForMol(mol)
    coreSmiles=Chem.MolToSmiles(core)

    if coreSmiles=="" and smiles[1]!=">":
        s=smiles[1]
        mol=Chem.MolFromSmiles(s)
        core=MurckoScaffold.GetScaffoldForMol(mol)
        coreSmiles=Chem.MolToSmiles(core)
    return s, mol

def getLevel(s,sims):
    flag=False;lv=0
    for p in sims:
        if flag:
            if s==p:
                return lv
            if p=='(':
                lv=lv+1
            elif p==')':
                lv=lv-1
        if s==p:
            if ss=="[":
                continue
            flag=True
        ss=p
    return 0

def findRing(sims):
    flag=True;r="";lv=0
    for s in sims:
        if s.isdigit():
            if flag:# starting point
                if ss=="[":
                    continue
                p=s;flag=False;l=getLevel(s,sims);r=ss+s;continue
            if s==p:# end point
                print("-->",l,r);
                r=r+s;return r
        if not flag:
            if s=="(":
                lv=lv+1;continue
            elif s==")":
                lv=lv-1;continue
            print(lv,l,s,lv<=l)
            if lv <= l and not s.isdigit():
                print("hit ",s)
                r=r+s
        ss=s
    return r

def prevFindRing(sims):
    r="";flag=True;level=0;p=1
    for s in sims:
        if s==")" or s=="]":
            level=level+1;continue
        if s=="(" or s=="[":
            level=level-1;continue
        if level!=0:
            continue
        if s.isdigit():
            if flag:
                p=s;r=ss+s;flag=False
            else:
                if s==p:
                    r=r+s;return r
                else:
                    print("unacceptable",sims)
                    r=r+s
        else:
            if not flag:
                r=r+s
        ss=s
    return r

def numRing(stars):
    level=0;n=0
#    print("starts",stars)
    for s in stars:
        if s=="[":
            level+=1;continue
        if s=="]":
            level-=1;continue
        if level!=0:
            continue
        if s.isdigit():
            n+=1
    return int(n/2)

def getBinder(sims):
    cnt=0;tmp=sims

    while numRing(tmp)>0:
        r=findRing(tmp)
        print("find",cnt,r,tmp);cnt+=1
        p=Chem.MolFromSmiles(tmp)
        x=Chem.MolFromSmiles(r)
#        tmp=Chem.MolToSmiles(Chem.ReplaceCore(p,x,labelByIndex=True))
        tmp=Chem.MolToSmiles(Chem.ReplaceCore(p,x))
#        print("-->",tmp)
    return tmp

def countRing(sims):

    binder=""

    x=re.findall(r'\d',sims)
#    print(x,sims)
    n=int(len(x)/2)
    if n>1:
        binder=getBinder(sims)
    return int(len(x)/2),binder

def outOfRing(mol):
    num=mol.GetNumAtoms()
    n=0
    for i in range(num):
        if not mol.GetAtomWithIdx(i).IsInRing():
            n=n+1
    return n

def getMolSequence(allData):
    mols=[];smis=[]

    smiles,mol=getOpeningMaterial(allData['smiles'][0])
    mols.append(mol)
    smis.append(smiles)
    for s in allData['connect']:
        mols.append(Chem.MolFromSmiles(s[0]))
        smis.append(s[0])
    return mols,smis

def list2str(ss):
    d=""
    for s in ss:
        d=d+"["+str(s)+"]"
    return d

def analyzeCoreChanges(mols):
    step=0;p1="";p3=[];p4=[];p5=[]

    core=MurckoScaffold.GetScaffoldForMol(mols[0])
    coreSmiles=Chem.MolToSmiles(core)
    n=outOfRing(mols[0])
#    p3.append("o"+str(n))
    p3.append(n)
    m,bind=countRing(coreSmiles);p4.append(m);p5.append(bind)
    for mol in mols[1:]:
        nextCore=MurckoScaffold.GetScaffoldForMol(mol)
        nextCoreSmiles=Chem.MolToSmiles(nextCore)
        n=outOfRing(mol)
        if nextCoreSmiles!=coreSmiles:
            p1=p1+"x"
            coreSmiles=nextCoreSmiles;step=step+1
            p3.append(n)
#            p3.append("x"+str(n))
            m,bind=countRing(coreSmiles);p4.append(m);p5.append(bind)
        else:
            p1=p1+"o"
#            p3.append("o"+str(n))
            p3.append(n)
    return [p1,len(mols)-1,p3,p4,p5,step]

def getFirstIds(cur):
    pid=1
    while True:
        try:
            ret=cur.execute(f"""select * from searchList where id='{pid}';""")
        except:
            print("easy")

        if len(ret.fetchall())>0:
            return pid
        pid=pid+1
        if pid>50:
            print("maybe no records in this database")
            exit()

def xmerge(src):
    hand=False;cc=""
    if len(src)<1:
        return "Non"
    if len(src)==2:
        return src[0]+"-"+src[1]
    p=""
    for s in src:
        match s:
            case '(':
                hand=False;cc=cc+s
            case ')':
                hand=False;cc=cc+s
            case '=':
                hand=False;cc=cc+s
            case '.':
                hand=False;cc=cc+s
            case _:
#                if hand and "R" in s:
                if hand:
                    if 'R' in s:
                        cc=cc+"-"+s
                    elif 'R' in p:
                        cc=cc+"-"+s
                    else:
                        cc=cc+s
                else:
                    cc=cc+s
                hand=True
        p=s
    return cc

def xconv(src):
    cc="";depth=0;r="";rs=[]

    l=len(src)
    for n,s in enumerate(src):
        match s:
            case '[':
                depth=depth+1
                continue
            case ']':
                depth=depth-1
                if depth==0:
                    cc=cc+xmerge(rs)+","
                    rs=[]
                continue
            case '*':
                rs.append('R'+src[n-1])
            case _:
                if s.isdecimal():
                    continue
                rs.append(s)
        p=s
    dst=",".join(cc.split(",")[:-1])
#    print("\nsrc",src)
#    print("cc",dst)
    return dst

#########################################################
###################### main #############################
#########################################################
input_dir="./";db="sList.db"
com='-chem chemical_name, -ppt for power point -n [1,2,3] to specify routes -p proc_per_line -d input_path, -one_size : rendering at the same size, -product_only : categolize routes with respct to reaction product only, -show_subsets : show routes including subsets'
com+=" -summary 0/1/2"
ops={'-h':com}

ppt=False;skip=False;hit_id=False
ids=[];mode="normal";full=True
for n,op in enumerate(sys.argv):
    if hit_id:
        if op[0]!='-':
            ids.append(int(op))
        else:
            hit_id=False
    if skip:
        skip=False;continue
    match op:
        case '-h':
            print(sys.argv[0],ops['-h']);exit()
        case '-v':
            print(version);exit()
        case '-ppt':
            ppt=True
        case '-simple':
            full=False
        case '-ids':
            hit_id=True
        case '-database':
            db=sys.argv[n+1]
            skip=True
        case '-forEval':
            mode="forEval"
        case '-d':
            input_dir=sys.argv[n+1]
            skip=True

conn=sqlite3.connect(db)
cur=conn.cursor()

if len(ids)<1 or ids[0]<1:
    ids[0]=getFirstIds(cur)

print("ids",ids)

ret=cur.execute(f"""select * from searchList;""")
descr=cur.description

page=openReport("test.txt",0.5)

sout=[]
for d in descr:
    sout.append(d[0])

prm,prmT=getPrm(cur,ids,sout)

for z in sout:
    print(z,prmT[z])

hint={}
for id in ids:
    hint[id]=easyGetReaction(id,cur)

cutList={}
for id in ids:
    cutList[id]=dropSubset(hint[id])
#    print("hint\n")
#    print(hint[id],len(hint[id]))
#    print("cutList\n")
#    print(cutList[id])

#for h in hint[id]:
#    print("h",h,"\n",len(h),"\n")
#exit()
result={}
for id in ids:
    ret=cur.execute(f"""select user,date from searchList where id='{id}';""")
    ret=ret.fetchall()

    uid=ret[0][0];date=ret[0][1];
    dst="/var/www/html/public/images/"+uid+"/"+date+"/route"

    result[id]={}
    m=getMax(hint[id][0])
#    print(len(hint[id][0]),m,len(hint[id][0])-len(cutList[id]))
    result[id]['routeNumber']=len(hint[id][0])
    result[id]['chemicalNumber']=m
    result[id]['fullRouteNumber']=len(hint[id][0])-len(cutList[id])

#    print(id,hint[id])

ox=40
##############################################################
for id in ids:
    sTable=f"searchTable{id}"
    ans=cur.execute(f"""select * from {sTable};""").fetchall()
    if len(ans[0][1].split(";"))<2:
        print("no route");continue

    allData={};routes=[];newInfo={};index={}
    for xx in ans:
        s={}
        if xx[1]=="":
            continue
        routes.append(int(xx[0]))
        s['connect']=strToConnect(xx[1])
        s['smiles']=strToSmiles(xx[2])
        s['info']=strToSmiles(xx[3])
        allData[int(xx[0])]=s

    for xx in reversed(ans):
        newInfo[int(xx[0])]=getLink(strToSmiles(xx[2]),index)

    newList,idx=getNewList(newInfo)

    simple_route=[]
    for nn in newList:
        simple_route.append(nn[1][0])

    full_route=[]
    for nn in newList:
        full_route=full_route+nn[1]

    if full:
        theRoute=full_route
    else:
        theRoute=simple_route

    hint=getReactionSummary(newList,idx)

    nextIndex=makeNextIndex(theRoute,hint[3])

    print("core structure analyze")
    print("structure changed","number of Ring","number of atoms outside Rings")
    print("<-- route number (database), printed route number (in pdf)")
    statAll={};stat=[{},{},{},{},{}]
    for route in theRoute:
        mols,smis=getMolSequence(allData[route])
#    [coreChange,coreCNum,ringNum,outOfRings]=analyzeCoreChanges(mols)
        r=analyzeCoreChanges(mols)
#        print(r[0],list2str(r[3]),list2str(r[2]),list2str(r[4]),"<--",route,nextIndex[route],r[1])
        statAll[route]=r
        for i in range(5):
            match i:
                case 0:
                    st=r[i]
                case 1:
                    st=r[i]
                case 2:
                    st=list2str(r[i])
                case 3:
                    st=list2str(r[i])
                case 4:
                    st=list2str(r[i])

            if not st in stat[i]:
                stat[i][st]=[route]
            else:
                stat[i][st].append(route)

#    for i in range(5):
#        print(i,"<----------------------------")
#        for key in sorted(stat[i]):
#            print(key,list2str(stat[i][key]))

#testForSvg(statAll):

    print(f"start-forEval for {id}")
    print("key1")
    i=0
    for key in sorted(stat[i]):
        print("##"+key,end="#")
        for route in stat[i][key]:
            print(str(route),end=",")
    print("\nkey2")
    i=1
    for key in sorted(stat[i]):
        print("##"+str(key),end="#")
        for route in sorted(stat[i][key]):
            print(str(route),end=",")
    print("\nkey3")
    i=2
    for key in sorted(stat[i]):
        print("##"+str(key),end="#")
        for route in sorted(stat[i][key]):
            print(str(route),end=",")
    print("\nkey4")
    i=3
    for key in sorted(stat[i]):
        print("##"+str(key),end="#")
        for route in sorted(stat[i][key]):
            print(str(route),end=",")
    print("\nkey5")
    i=4
    for key in sorted(stat[i]):
        print("##"+str(xconv(key)),end="#")
        for route in sorted(stat[i][key]):
            print(str(route),end=",")
    print("\nend-forEval")

#  for drawImages
    print(f"start-forDraw for {id}")
    for key in sorted(stat[0]):
        for route in stat[0][key]:
            ddd=dst+"/"+str(route)+"/"
            print("route"+str(route)+"\npath:"+ddd.split("/var/www/html/public")[1],end=";")
            xx=ox
            connect=allData[route]['connect']
            smiles=allData[route]['smiles']
            info=allData[route]['info']
            fc,image_size=makeSvg(smiles,ddd)
            objs=makeObject(page,connect,image_size)
            for obj in objs:
                for item in obj:
                    x0=item['pos'][0][0];y0=item['pos'][0][1]
                    if item['type']==1:
                        page.drawImage(ddd+"/"+item['svg'][0],x0,y0)
                    else:
                        x1=item['pos'][1][0];y1=item['pos'][1][1]
                        xx=x0;

                        for n,fn in enumerate(item['svg']):
                            xx=xx+page.mg
                            r=image_size[fn]
                            page.drawImage(ddd+fn,xx,y0+page.mg,cat=True)
                            xx+=r[0]
                        page.drawArrow(arrow(x0,y0,x1,y1,page.scale))
            page.stroke()
    print("end-forDraw")
#
#        print(core)
#        print(lastP)
#        print("info",s['info'])
#komai
