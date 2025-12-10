import sys,os,time
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
from dbTools import *

def getMaxMin(y0,btm,top):

    if y0<btm:
        btm=y0
    if y0>top:
        top=y0

    return btm,top

def getPageRange(btm,top,py,page):
    p=[]
    if top-btm<py:
        p.append([btm,top,top-py]);return p
    py1=top;py0=py1-py;py2=py1-py
    p.append([py0,py1,py2]);
    # py-page.mg
    py1=py0
    while py1>btm:
        py0=py1-page.height/page.scale+page.mg
        if py0<btm:
            ppy=py1-page.height/page.scale+page.mg*3
            p.append([btm,py1,ppy]);py1=py0
        else:
            p.append([py0,py1,py0]);py1=py0
    return p

def getAllData(ans,_routes):
    allData={};routes=[];newInfo={};index={}

    for xx in ans:
        s={}
        if xx[1]=="":
            continue
        if len(_routes)>0:
            if not int(xx[0]) in _routes:
                continue

        routes.append(int(xx[0]))
        s['connect']=strToConnect(xx[1])
        s['smiles']=strToSmiles(xx[2])
        s['info']=strToSmiles(xx[3])
        allData[int(xx[0])]=s

    for xx in reversed(ans):
        if len(_routes)>0:
            if not int(xx[0]) in _routes:
                continue
        newInfo[int(xx[0])]=getLink(strToSmiles(xx[2]),index)

    return allData,routes,newInfo,index

def getPids(cur):
    ret=cur.execute(f"""select id from searchList;""")
    ret=ret.fetchall()
    rs=""
    for r in ret:
        if rs=="":
            rs=str(r[0])
        else:
            rs=rs+","+str(r[0])
    return rs

def addLog(mes):
    with open(debug_log,"a") as fp:
        fp.write(mes+"\n")

debug=True
debug=False

version="1.8 08152025"

_web=False
_web=True

if _web:
    dr='/var/www/html/public/images/tmp/' # working directory
else:
    dr='tmp/' # working directory

##############################################
#################### main ####################
##############################################

outp=False
proc_per_line=10 # how many reactions per line roughly ?
_fc=True
_tite=False
_type=3
_product_only=False 
_include_subsets=True
_include_subsets=False
oper="normal"
forced=False

start=time.time()
ppt=False;chem=None;_routes=[]
input_dir="./"
chem=False
com='-chem name, -ppt for power point -n [1,2,3] to specify routes -p proc_per_line -d input_path, -one_size : rendering at the same size, -product_only : categolize routes with respct to reaction product only, -show_subsets : show routes including subsets'
com+=" -summary 0/1/2"
ops={'-h':com}
skip=False
pid=-1

debug_log="/var/www/html/public/images/report/readDb_web.log";

scale=0.15
scale=0.4
scale=0.2
db="sList.db"
fp=open("readDb.log","w")
noRoute=False
lock=False
lock_file="/var/www/html/public/images/report/readDb.lock"
for n,op in enumerate(sys.argv):
    fp.write(f" {op}")
    if skip:
        skip=False;continue
    match op:
        case '-h':
            print(sys.argv[0],ops['-h']);exit()
        case '-v':
            print(version);exit()
        case '-chem':
            chem=sys.argv[n+1]
            fp.write(f" {chem}")
            skip=True
        case '-lock':
            lock=True
        case '-ppt':
            ppt=True
        case '-heavy':
            _type=2
        case '-atom_number':
            _type=3
        case '-one_size':
            _fc=False
        case '-show_subsets':
            _include_subsets=True
        case '-route':
            _routes=getList(sys.argv[n+1])
            fp.write(f" {_routes}")
            skip=True
        case '-tite':
            _tite=True
        case '-s':
            scale=float(sys.argv[n+1])
            fp.write(f" {scale}")
            skip=True
        case '-id':
            pid=sys.argv[n+1]
            fp.write(f" {pid}")
            skip=True
        case '-database':
            db=sys.argv[n+1]
            skip=True
        case '-d':
            input_dir=sys.argv[n+1]
            fp.write(f" {input_dir}")
            skip=True
        case '-product_only':
            _product_only=False
        case '-p':
            proc_per_line=int(sys.argv[n+1])
            fp.write(f" {proc_per_line}")
            skip=True
        case '-db_list':
            oper='db_list'
        case '-thumbnail':
            oper='thumbnail'
        case '-drop':
            oper='drop'
        case '-db':
            oper='db'
        case '-askTids':
            oper='askTids'
        case '-force':
            forced=True
        case '-done':
            outp=sys.argv[n+1]
            skip=True

    fp.write("\n")
fp.close()
conn=sqlite3.connect(db)
cur=conn.cursor()

#print("options:",pid,oper,db)

if _web:
    tmp_dir="/".join(debug_log.split("/")[:-1])
    os.makedirs(tmp_dir,exist_ok=True)

if _web and debug:
    with open(debug_log,"w") as fp:
        fp.write(f"{oper}\n")

if oper=='askTids':
    print(getPids(cur))
    exit()

if oper=='drop':
    dr="/var/www/html/public/images"
    if pid=="all":
        pid=getPids(cur);
  
#    with open(debug_log,"w") as fp:
#        fp.write(pid+":remove\n")

    for p in pid.split(','):
        if int(p)<0:
            continue;
        ret=cur.execute(f"""select user,date from searchList where id='{p}';""")
        ret=ret.fetchall()

        uid=ret[0][0];
        date=ret[0][1];
        dst="/var/www/html/public/images/"+uid+"/"+date+"/"
        if os.path.isdir(dst):
#            with open(debug_log,"w") as fp:
#                fp.write(dst+":remove\n")
            shutil.rmtree(dst)
        if db=="sList.db":
            pdf="/var/www/html/public/images/report/"+str(p)+".pdf"
        else:
            pdf="/var/www/html/public/images/"+uid+"/report/"+str(p)+".pdf"
        if os.path.isfile(pdf):
            os.remove(pdf)

        cur.execute(f"""drop table searchTable{p}""")
        cur.execute(f"""delete from searchList where id={p}""")
        cur.execute(f"""delete from parent where id={p}""")
    conn.commit()
    conn.close()
    exit()

if oper=='thumbnail':
    ret=cur.execute(f"""select id,cSmiles from searchList;""")

    for d in ret.fetchall():
        smile_file=input_dir+"+"+str(d[0])+".svg"
        if os.path.exists(smile_file):
            continue
        smile2svg(d[1],"pid"+str(d[0])+".svg",dr=input_dir,size=200)
    conn.close()
    exit()

ret=cur.execute(f"""select * from searchList;""")
descr=cur.description

sout=""
flag=True
for d in descr:
    if not flag:
        sout+="##"
        flag=False
    sout+=d[0]+"##"
    flag=True
sout+="#"

for id in ret.fetchall():
    flag=True
    for i in id:
         sout+=str(i)+"##"
    sout+="#"

if _web and debug:
    with open(debug_log,"a") as fp:
        fp.write(f"a{oper}\n")

if oper=='db_list':
    conn.close()
    print(sout,end="")
    exit(0)

ret=cur.execute(f"""select substance from searchList where id='{pid}';""")
ret=ret.fetchall()

if len(ret)<1:
    ret=cur.execute(f"""select * from searchList where id='{pid}';""")
    ret=ret.fetchall()
    print(f"No data for id={pid}:",ret)
    conn.close()
    exit()

ret=cur.execute(f"""select user,date,substance from searchList where id='{pid}';""")
ret=ret.fetchall()

uid=ret[0][0];
date=ret[0][1];
chemical=ret[0][2]

if chem:
    name=chem
else:
    name=str(pid)

if ppt:
    output_file=name+'.pptx'
else:
    output_file=name+'.pdf'

log_name="readDb"+str(pid)+oper+".log"

if _web:
    if not os.path.isdir(input_dir):
        os.makedirs(input_dir,exist_ok=True)
    log_name=input_dir+"/readDb"+str(pid)+"normal.log"
    if os.path.isfile(log_name):
        os.remove(log_name)
    index_file=input_dir+"/readDb"+str(pid)+"index.txt"
    if os.path.isfile(index_file):
        os.remove(index_file)

    output_file=input_dir+output_file
#    print(input_dir)
#    exit()

if os.path.isfile(output_file) and not forced:
    conn.close()
    exit()

chk=db+":"+str(pid)
if lock!=False:
    if os.path.isfile(lock_file):
        with open(lock_file,"r") as fp:
            for line in fp:
                if chk in line:
#                print("this-->"+chk)
                    exit()

if lock!=False:
    with open(lock_file,"a") as fp:
        fp.write(chk+"\n")
head=None

if _include_subsets:
    show='total'
    com="."
else:
    show='sub'
    com=" and subsets eliminated."

span=proc_per_line*2
if span<10:
    span=10

_py=500
py=_py
ox=40

print("output_file--------------------------->",output_file)
if scale>0.5:
    scale=0.5
if scale<0.05:
    scale=0.05

page=openReport(output_file,scale,title=chemical)

allId={}
cur=conn.cursor()
ret=cur.execute(f"""select loop from searchList where id={pid};""").fetchall()

if (len(ret)<1):
    print(f"not availabale id:{pid}")
    exit()

maxLoop=ret[0][0]

sTable=f"searchTable{pid}"

sql=f'select * from "{sTable}";'
ans=cur.execute(sql).fetchall()
if len(ans[0][1].split(";"))<2:
    noRoute=True
    head=[chemical,[[ans[0][2].split(";")[0]]],noRoute]

if oper=='db':
    dump(ans,pid,input_dir)
    with open(log_name,"a") as fp:
        fp.write(f"mission ok\n")
    exit()

allData,routes,newInfo,index=getAllData(ans,_routes)

if noRoute:
    head_page(head,page)
    with open(log_name,"a") as fp:
        fp.write(f"mission ok\n")
    page.close()
    with open(log_name,"a") as fp:
        fp.write(f"mission ok\n")
    exit()

newList,idx=getNewList(newInfo)

ret=cur.execute(f"""select * from parent where id={pid};""")
ans=ret.fetchall()

parent=strToParent(ans[0])
hint=getReactionSummary(newList,idx)

    #print("hint",hint[0],hint[1],hint[2])
easy=str(len(parent['total']))+" variations from "+str(len(hint[0]))+" routes over "+str(maxLoop)+" queries"
hint.append(easy)
#easy=str(len(sim))+" routes over "+str(maxLoop)+" queries"

if parent!=False:
    head=[chemical,allData[int(routes[0])]['smiles'],hint]
else:
    head=[chemical,allData[1]['smiles'],hint]
 
if len(_routes)==0:
    head_page(head,page)
else:
    head_page(head,page,arrange=False)

with open(log_name,"w") as fp:
    fp.write("head_page done\n")

sTable=f"searchTable{pid}"
sql=f'select * from "{sTable}";'
ans=cur.execute(sql)

#oy=py-page.mg*5
dst="/var/www/html/public/images/"+uid+"/"+date+"/route"

if len(_routes)>0:
   theList=_routes
else:
   theList=hint[3]

page.scale=scale
py=page.height/page.scale-page.mg
for name,route in enumerate(theList):
    fp=open("/var/www/html/public/images/report/ok.log","a")
    fp.write("Route"+str(route)+"\n")
    fp.close()
    print("---------->  Route"+str(route))
    now=str(int(time.time()-start))+" sec"
    ddd=dst+"/"+str(route)+"/"
    print("images from",ddd)
    with open(index_file,"a") as fp:
        fp.write(f"{name+1}:{route}\n")
    with open(log_name,"a") as fp:
        fp.write(f"{now}:route {route} start\n")
    connect=allData[route]['connect']
    smiles=allData[route]['smiles']
    info=allData[route]['info']
#    print("smiles",smiles)
    print("connect",connect)
    fc,image_size=makeSvg(smiles,ddd)

    if not _fc:
        fc=_fc

    page.setFont()
    objs=makeObject(page,connect,image_size)

#    print("objs",objs)

    stack=[]
    top=-page.height;btm=page.height*10.0

#ugly shift
    ttt=-page.height
    for obj in objs:
        for item in obj:
            if item['type']==1:
                tt=item['pos'][0][1]+item['size'][1]
                if ttt<tt:
                    ttt=tt
            else:
                if ttt<item['pos'][0][1]:
                    ttt=item['pos'][0][1]

    ttt=0
    for obj in objs:
        for item in obj:
            x0=item['pos'][0][0];y0=item['pos'][0][1]
            if item['type']==1:
#                page.drawImage(ddd+item['svg'][0],x0,y0)
                rr=chkProc(ddd+"/svgFile.info",item['svg'][0])
                if not rr[1]:
                    continue
                stack.append([1,x0,y0,ddd+item['svg'][0]])
#                print("this",x0,y0)
                btm,top=getMaxMin(y0,btm,top)
                btm,top=getMaxMin(y0+rr[1][1],btm,top)
            else:
                x1=item['pos'][1][0];y1=item['pos'][1][1]
                xx=x0;
                for n,fn in enumerate(item['svg']):
                    xx=xx+page.mg
                    r=image_size[fn]
#                    page.drawImage(ddd+fn,xx,y0+page.mg)
                    rr=chkProc(ddd+"/svgFile.info",fn)
                    stack.append([1,xx,y0+page.mg,ddd+fn])
                    btm,top=getMaxMin(y0+page.mg,btm,top)
                    btm,top=getMaxMin(y0+page.mg+rr[1][1],btm,top)
                    xx+=r[0]
#                page.drawArrow(arrow(x0,y0,x1,y1,page.scale))
                stack.append([2,x0,y0,x1,y1,page.scale])
                btm,top=getMaxMin(y0,btm,top)
                btm,top=getMaxMin(y1,btm,top)

    if py<top-btm:
        py=page.height/page.scale-page.mg
        page.stroke()

    p_range=getPageRange(btm,top,py,page)

    title=True
    title=False
    for n,xx in enumerate(p_range):
        for i in stack:
            my=i[2]
            if my>xx[0]-10 and my<xx[1]+10:
                if i[0]==1:
                    page.drawImage(i[3],i[1],my-xx[2]-ttt)

        for i in stack:
            my=i[2]
            if my>xx[0]-10 and my<xx[1]+10:
                if i[0]==2:
                    page.drawArrow(arrow(i[1],my-xx[2]-ttt,i[3],i[4]-xx[2]-ttt,i[5]))

        if not title:
            if len(_routes)>0:
                page.drawStringR(ox,py-page.font_size/page.scale,"Route"+str(route))
            else:
                page.drawStringR(ox,py-page.font_size/page.scale,"Route"+str(name+1))

#            print(str(route),py-page.font_size/page.scale)
            title=True

#        page.drawStringR(ox,py-(xx[1]-xx[0])*page.scale,"up"+str(route))
#        page.drawStringR(ox,py-(xx[0]-xx[0])*page.scale,"dwn"+str(route))
#
        if n<len(p_range)-1:
            py=page.height/page.scale-page.mg
            page.stroke()
        else:
            py=py+(xx[0]-xx[1])-page.mg*15

#    print("---------->  Route"+str(route))
#    oy-=(hy*page.scale+page.mg)
    now=str(int(time.time()-start))+" sec"
    with open(log_name,"a") as fp:
        fp.write(f"{now}:route {route} done\n")

#    if name>max_route and debug==True:
#        break
#    break

conn.close()
page.close()
with open(log_name,"a") as fp:
     fp.write(f"mission ok\n")

if outp!=False:
    with open(outp,"a") as fp:
        fp.write("readDb_done\n")

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
