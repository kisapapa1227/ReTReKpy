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

class openReport:
    def __init__(self,file,scale):
        ext=file.split(".")[1]
        if ext=='pdf':
            print("open pdf:",file)
            self.type=ext
            self.scale=scale
            self.page=canvas.Canvas(file,pagesize=landscape(A4))
            self.font_name="Times-Bold"
            self.height=210*mm
            self.width=297*mm
            self.font_size=4*mm
            self._font_size=4*mm
        else:
            self.type='ppt'
            self.scale=scale
            self.page=Presentation()
            self.page.slide_height=Inches(8.27)
            self.page.slide_width =Inches(11.69)
            self.height=self.page.slide_height
            self.sc=Inches(11.69)/830
#            self.sc=Inches(8.0)/800
            blank_slide_layout = self.page.slide_layouts[6]
            self.slide = self.page.slides.add_slide(blank_slide_layout)
            self.file=file
            self.font_name="Times-Bold"
            self.font_size=14 # 14<-16 
            self._font_size=14 #

    def getSpan(self):
#       scale=_magic_x/(span)
        return _magic_x/self.scale
     
    def setFont(self,font_name=None,font_size=None):# size is relative value to default
        if self.type=='pdf':
            if font_name:
                self.font_name=font_name
            if font_size:
                self.font_size=self._font_size*font_size
            self.page.setFont(self.font_name,self.font_size)
        else:
            if font_size:
                self.font_size=self._font_size*font_size

    def stroke(self):
        if self.type=='pdf':
            self.page.showPage()
        else:
            blank_slide_layout = self.page.slide_layouts[6]
            self.slide = self.page.slides.add_slide(blank_slide_layout)

    def drawImage(self,fn,x,y,scale,bt=False):
        if bt:
            y=y+self.scale*0.55
        else:
            y=y-scale*0.4
#            y=y

        if self.type=='pdf':
            drawing=svg2rlg(fn)
            drawing.renderScale=scale/drawing.width
            renderPDF.draw(drawing, canvas=self.page, x=x, y=y)
        else:
            work=dr+"/komai.svg"

            key="2.0px"
            kkk="4.0px"

            with open(fn,"r") as f:
                all=f.readlines()

            with open(work,"w") as f:
                for line in all:
                    o=line.replace(key,kkk)
                    f.write(o)

            gn=fn.split(".svg")[0]+".jpg"
            if _web:
                subprocess.run(['/usr/bin/convert',work,gn])
            else:
                subprocess.run(['convert',work,gn])
            x0,y0=x*self.sc,self.page.slide_height-(y+scale)*self.sc
        
            self.slide.shapes.add_picture(
                gn,left=x0,top=y0,
#                width=scale*self.sc, height=scale*self.sc
                width=scale*self.sc*0.9, height=scale*self.sc*0.9
    )


    def drawLine(self,pt):
        if self.type=='pdf':
            self.page.line(pt[0],pt[1],pt[2],pt[3])
        else:
            x1,y1=pt[0]*self.sc,self.page.slide_height-pt[1]*self.sc
            x2,y2=pt[2]*self.sc,self.page.slide_height-pt[3]*self.sc
            putLine(x1,y1,x2,y2,self.slide.shapes)

    def drawArrow(self,pt):
        if self.type=='pdf':
            self.page.line(pt[0],pt[1],pt[2],pt[3])
            self.page.line(pt[2],pt[3],pt[4],pt[5])
            self.page.line(pt[2],pt[3],pt[6],pt[7])
        else:
            x0,y0=pt[0]*self.sc,self.page.slide_height-pt[1]*self.sc
            x1,y1=pt[2]*self.sc,self.page.slide_height-pt[3]*self.sc
            x2,y2=pt[4]*self.sc,self.page.slide_height-pt[5]*self.sc
            x3,y3=pt[6]*self.sc,self.page.slide_height-pt[7]*self.sc
            group=[]
            group+=[putLine(x0,y0,x1,y1,self.slide.shapes)]
            group+=[putLine(x1,y1,x2,y2,self.slide.shapes)]
            group+=[putLine(x1,y1,x3,y3,self.slide.shapes)]
            self.slide.shapes.add_group_shape(group)

    def getProperLength(self,mag,string):
        if self.type=='pdf':
            tx=string
            while True:
                t=(self.getSpan()-1)*self.scale-stringWidth(tx,self.font_name,self.font_size)
                if t>0:
                    return len(tx)
                tx=tx[0:len(tx)-3]
        else:
            return 50

    def drawPhrase(self,x,y,string,center=False,adj=True):
        if self.type=='pdf':
            drawString(self,x,y,string,center=center,adj=adj)
        else:
            shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,py)
            # komai

    def drawString(self,x,y,string,center=False,adj=True):
        if self.type=='pdf':
            if center:
                x+=self.getSpan()*self.scale/2-stringWidth(string,self.font_name,self.font_size)/2
            self.page.drawString(x,y,string)
        else:
            xx=x*self.sc
            yy=self.page.slide_height-y*self.sc-self.font_size
            if adj:
#                komai
#                shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,Cm(1))
                pw=Cm(0.1)*self.font_size*(len(string)+2)*0.18
                py=Cm(0.1)*self.font_size*0.5
                shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,py)
            else:
                shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),Cm(1),Cm(10))

            tf=shape.text_frame
            tf.word_wrap=True
#            tf.auto_size=MSO_AUTO_SIZE.SHAPE_TO_FIT_TEXT
            p=tf.paragraphs[0]
            run=p.add_run()
            run.text=string
            run.font.name=self.font_name
            run.font.size=Pt(self.font_size)
            #komai

    def close(self):
        if self.type=='pdf':
            self.page.save()
        else:
            self.page.save(self.file)

def putLine(x0,y0,x1,y1,shapes):
    if x1>x0:
        if y1>y0:
           shape=shapes.add_shape(MSO_SHAPE.LINE_INVERSE,x0,y0,x1-x0,y1-y0)
           shape.element.getchildren()[1].getchildren()[0].set('flipV','1')
        else:
           shape=shapes.add_shape(MSO_SHAPE.LINE_INVERSE,x0,y1,x1-x0,-y1+y0)
    else:
        if y1>y0:
           shape=shapes.add_shape(MSO_SHAPE.LINE_INVERSE,x1,y0,-x1+x0,y1-y0)
        else:
           shape=shapes.add_shape(MSO_SHAPE.LINE_INVERSE,x1,y1,-x1+x0,-y1+y0)
           shape.element.getchildren()[1].getchildren()[0].set('flipV','1')
    shape.fill.solid()
    shape.fill.fore_color.rgb=RGBColor(0,255,0)
    shape.line.color.rgb=RGBColor(0,0,0)
    return shape

def picSize():
    return 1.2

def smile2sgv(smiles,name,dr=dr,size=200):
#    if os.path.exists(dr):
#        shutil.rmtree(dr)
    if not os.path.exists(dr):
        os.makedirs(dr)

    mol=Chem.MolFromSmiles(smiles)
#    mw=rdMolDescriptors.CalcExactMolWt(mol)
    mw=mol.GetNumAtoms()
#    mw=Chem.AddHs(mol).GetNumAtoms()
#    fc=int(math.sqrt(mw/5))
    fc=int(math.sqrt(mw)/1.5)
#    print(smiles,mw,fc)
#    fc=int(math.pow(mw/20,0.333))
#    print(mw,fc,name)
    if fc<1:
        fc=1
    view = rdMolDraw2D.MolDraw2DSVG(size*fc,size*fc)

    option=view.drawOptions()
    option.circleAtoms=False
    option.continuousHightlight=False

    tm = rdMolDraw2D.PrepareMolForDrawing(mol)
    view.DrawMolecule(tm)
    view.FinishDrawing()

    svg=view.GetDrawingText()

    with open(dr+"/"+name,'w') as f:
        f.write(svg)

    return fc

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

def makeSgv(smiles,type=0):
    size=100
    fcs=[]
    for n,proc in enumerate(smiles):
        fc={}
        for m,comp in enumerate(proc):
            if ">" in comp:
                continue
#            print(comp,"e"+str(n+1)+"x"+str(m+1)+".svg") #debuf
            fc[m]=smile2sgv(comp,"e"+str(n+1)+"x"+str(m+1)+".svg",size=size)
        fcs.append(fc)

    return fcs

def arrow(x0,y0,x1,y1,scale,skip=1.0,fc=False,ugly_patch=True):

    if not fc:
        return prevArrow(x0,y0,x1,y1,scale,ugly_patch=ugly_patch)
#
# ugly patch
#
    if ugly_patch and x0>x1:
        x1=x0+scale*3
        y1=y0
# end
# anyway
    l=math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
    rt=(l-scale*skip)/l
    sx=10
    pi=3.1415/12.0 # 30 deg

#from: shift full
    x2=x0+skip*scale
    y2=y0+scale/2
#to: shift full 
    x3=x1
    y3=y1+scale/2

    a=math.atan2(y3-y2,x3-x2)

    xx=x3-sx*math.cos(a+pi)
    yy=y3-sx*math.sin(a+pi)
    ret=[x2,y2,x3,y3,xx,yy]

    xx=x3-sx*math.cos(a-pi)
    yy=y3-sx*math.sin(a-pi)

    return ret+[xx,yy]

def prevArrow(x0,y0,x1,y1,scale,skip=1.0,ugly_patch=True):
#
# ugly patch
#
    if ugly_patch and x0>x1:
        x1=x0+scale*3
        y1=y0
# end
# anyway
    l=math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
    rt=(l-scale*skip)/l
    sx=10
    pi=3.1415/12.0 # 30 deg

#from: shift full
    x2=x1+(x0-x1)*(rt-0.10)
    y2=y1+(y0-y1)*(rt-0.10)+scale/2
#to: shift full 
    x3=x0+(x1-x0)*(1.0-0.10)
    y3=y0+(y1-y0)*(1.0-0.10)+scale/2

    a=math.atan2(y1-y0,x1-x0)

    xx=x3-sx*math.cos(a+pi)
    yy=y3-sx*math.sin(a+pi)
    ret=[x2,y2,x3,y3,xx,yy]

    xx=x3-sx*math.cos(a-pi)
    yy=y3-sx*math.sin(a-pi)

    return ret+[xx,yy]

def getName(a,b,dr=dr):
    fn=dr+"/e"+str(a+1)+"x"+str(b)+".svg"#<-- need modification
    return(fn)

def getGoal(smiles,connect):
    goal=[]
    l=len(smiles)
    for i in range(l):
        if connect[i][3]==-1:
            goal.append(i)
    return goal

def getSharpHead(h,c):
    sharpHead=[]

    for n in h:
        branch=[]
        while True:
            if len(c[n][2])<1:#i.e. starting point
                if len(branch)<1:# and no reserved entry
                    break
                else:
                    m=branch.pop(-1) # pop and restart from reserved entry
                    n=c[m][3]
            elif len(c[n][2])>1:# branching 
                branch=branch+c[n][2][1:] # thus reserve others
                m=c[n][2][0] # go forward
            else:
                m=c[n][2][0] # 
            if m>-1 and c[n][4]<c[m][4]:
                sharpHead.append([n,m])
            n=m
            if n<0:# anyway no
                break
    return sharpHead

def rightEdge(i,fc,cnt,l_fc):
    if i[1][0]<0 or cnt>=l_fc:
        return i[4]+len(i[7])+2
    return i[4]+len(i[7])+2+fc[cnt][i[1][0]]

def prevEvalHeight(connect):
    depth=0
    for i in range(len(connect)):
        if depth > connect[i][5]:
            depth=connect[i][5]
    return depth

def evalHeight(c,fc):
    if not fc:
        return prevEvalHeight(c)

    depth=0
    for n,vv in enumerate(fc):
        for v in vv.values():
            w=c[n][5]-v/2
            if depth > w:
                depth=w
    return depth

def compSmiles(s1,s2):
    for a,b in zip(s1,s2):
        if a!=b:
            return False
    return True

def annealHead(sh,c):
    for p in sh:
        m=len(c)
        keep=c[p[1]][0]
        #komai 0830
#        c[p[1]][3]=-2 # right edge out: no connection
        c[p[1]][3]=-2 # right edge out: no connection
        c[p[0]][2]=[-1] # right edge in

def _ccn(i):
    if i<10:
        return 1
    if i<100:
        return 2
    return 3

def getPropString(m):
    if len(m)<2:
        s=str(m[0])
#        w=_ccn(m[0])+2
        w=_ccn(m[0])+2
        return s,w

#    s=str(m[-1])+"-"+str(m[0])
    s=str(m[0])+"-"+str(m[-1])
#    w=_ccn(m[0])+_ccn(m[-1])+2
    w=_ccn(m[0])+_ccn(m[-1])+2
    return s,w

def preDraw(init,proc,branch):

    ixx=0
    yy=1
    pos={}
    deep={}

    for i in proc[init]:
        pos[i]=[ixx,yy]
        s,w=getPropString(branch[i])
        ixx=ixx+w
    deep[init]=yy

    right=ixx
    theLeft=0
    y1=yy

    for j in reversed(proc):
        if j==init:
            continue
        deep[j]=0
        first=True
        prev=True
        curr=True
        ixx=-1

        for i in reversed(proc[j]):
            if i in pos:
                ixx=pos[i][0]
                yy=pos[i][1]
                curr=True
            else:
#                if ixx<0:
#                    ixx=right
                if ixx<0:
                    ixx=theLeft
                curr=False
                s,w=getPropString(branch[i])
                if first:
                    y1=y1+1
                    first=False
                yy=y1

                pos[i]=[ixx-w,yy]

            ixx=pos[i][0]
            if ixx<theLeft:
                theLeft=ixx
            prev=curr

            if deep[j]<yy:
                deep[j]=yy
    if theLeft<0:
        for i in pos:
            pos[i][0]-=theLeft

    return pos,deep

def similarityOnRoute(init,proc,branch,all):

    pos,deep=preDraw(init,proc,branch)
    name={}
    for j in proc:
        if deep[j] in name:
            name[deep[j]]+=","+str(j)
        else:
            name[deep[j]]=str(j)

    sim=[]
    fps=False;
    for i in name:
        for j in name[i].split(','):
            tg_route=j
            ssMax=0
            ss=getSmiles(j,all)
            sl=len(ss)
            if sl>ssMax:
                ssMax,tg_route=sl,j
        fpMin=10
        if j!=tg_route:
            ss=getSmiles(tg_route,all)
        if not fps:
            fps=FingerprintMols.FingerprintMol(Chem.MolFromSmiles(ss[-1][-1]))
        else:
            tg=whichIsFirstMaterial(ss[0],_type)
            fpp=FingerprintMols.FingerprintMol(Chem.MolFromSmiles(ss[0][tg]))
            sm=DataStructs.FingerprintSimilarity(fps,fpp)
#            sim.append([tg_route,sm])
#    return sim easy way
# question....
            for smile in ss:
                fpp=FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smile[-1]))
                ss=DataStructs.FingerprintSimilarity(fps,fpp)
                if sm>ss:
                    sm=ss
            sim.append([tg_route,sm])
    print("sim",sim)
    return sim

def head_page(head,page,hope):
#head=[chemical,getSmiles(route[0],all),ox,oy,["statistics will be ..."]]
    ox=40
    oy=480
    mag=3

    page.setFont(font_size=mag)
    page.drawString(ox,oy,head[0])
    scale=page.scale
#
    smiles=head[1][-1][-1]
    fn="head.svg"
    x0=ox
    py=oy-scale*10.5
#    print("smiles --->",smiles)
    smile2sgv(smiles,fn,size=100)
    page.drawImage(dr+"/"+fn,x0,py,scale*10)

    pt=0;l=len(smiles)
#    py=py-mag*page.scale
    mg=2
    page.setFont(font_size=mg)
    wd=page.getProperLength(mg,smiles)
#    wd=7
    tl=smiles
    n=int(math.ceil(l/wd))
    py=(n+0.5)*page.font_size
    while pt<l:
        end=min(pt+wd,l)
        tl=smiles[pt:end]
        py-=page.font_size
        page.drawString(ox,py,tl,center=True)
        pt+=wd

#    print("head",type(head[2][0]))
    if type(head[2][0])==type(True):
        page.drawString(250,500,head[2][1])
        page.close()
        exit()

    proc=head[-1][0]
    branch=head[-1][1]
    init=head[-1][2]
    sear=head[-1][3]

    ox=250
    oy=550
    sx=100
    sy=50
    sy=35

    page.drawString(ox,oy,"Drawing summary")
#
# draw initial one
#
    ll=0
    for i in proc[init]:
        _,w=getPropString(branch[i])
        ll+=w

    if ox+(ll+8)*(sx)>_magic_x:
        sx=(_magic_x-ox)/(ll+8)
#komai
    if ll>30:
        page.setFont(font_size=1.5)
    print("length",ll)

    page.drawString(ox+(2)*sx,oy-(1)*sy,"reactions")

#    ox=ox+sx*5
    ox=ox+sx*0
    yy=oy-sy*2

    first=True
    scale=20
    
    pos,deep=preDraw(init,proc,branch)
#
# the first line
#
    for i in pos:
        ix=pos[i][0]
        iy=pos[i][1]
        s,w=getPropString(branch[i])
        page.drawString(ox+(ix+2)*sx,oy-(iy+1)*sy,s)

    ddd={}
    for i in proc:
        first=True
        for j in proc[i]:
            x1=pos[j][0]
            y1=pos[j][1]
            if not first:
                draw=False
                if j0 in ddd:
                    if j in ddd[j0]:
                        ddd[j0].append(j)
                        draw=True
                else:
                    ddd[j0]=[j]
                    draw=True
#                if draw:
                if True:
                    page.drawArrow(arrow(ox+(x0+w-1)*sx,oy-(y0+1)*sy-4,ox+(x1+1.9)*sx,oy-(y1+1)*sy-4,scale,fc=True)) 
            s,w=getPropString(branch[j])
            first=False
            x0=x1
            y0=y1
            j0=j

    name={}
    for j in proc:
        if deep[j] in name:
            name[deep[j]]+=","+str(j)
        else:
            name[deep[j]]=str(j)

    ll=0
    for j in name:
        if ll<len(name[j]):
            ll=len(name[j])

    if ll>23:
        page.setFont(font_size=1)
    elif ll>10:
        page.setFont(font_size=1.5)
    else:
        page.setFont(font_size=2)

    bt=0
    for j in name:
        print(j,name[j])
        if bt<j:
            bt=j
        if len(name[j])>20:
            page.drawString(_magic_x+sx*2,oy-(j+1)*sy,name[j][:22])
            page.drawString(_magic_x+sx*2,oy-(j+1)*sy-sy/2,name[j][22:])
        else:
            page.drawString(_magic_x+sx*2,oy-(j+1)*sy,name[j])

    page.setFont(font_size=2)
    page.drawString(_magic_x+sx*3,oy-(1)*sy,"route #id")

    if oy-(bt+2)*sy-sy<0:
        py=oy-sy
    else:
        py=oy-(bt+2)*sy-sy

    page.drawString(ox,py,sear)
    py-=sy*0.5

    page.setFont(font_size=1)

    add="Routes start at similarity of "
    for i in sorted(hope):
        add+="{:.2f}".format(i[1])+"("+str(i[0])+"),"
        if len(add)>90:
            page.drawString(ox,py,add)
            py-=sy*0.01
            add=""

    if len(add)>0:
       page.drawString(ox,py-sy*0.5,add)

    page.stroke()
    page.setFont(font_size=1)

# komai2
def getList(s):
    ret=[]
    s=s.split('[')[1].split(']')[0].split(',')
    for i in s:
        ret.append(int(i))
    return ret

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

    print("len",len(_smiles),l,pt)
    if l<1:
        return False

    for ss in reversed(_smiles[pt]):
        if typ=="smiles":
            ref.append(ss)
        else:
            ref.append(ss[-1])

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

def getTheBottom(c,f):
    py=evalHeight(c,fc)
    return py-3

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

for n,op in enumerate(sys.argv):
    if skip:
        skip=False;continue
#    print(n,op)
    match op:
        case '-d':
            input_dir=sys.argv[n+1]
            skip=True
        case '-u':
            user=sys.argv[n+1]
            skip=True
        case '-s':
            script=sys.argv[n+1]
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
            for n,i in enumerate(es):
                print(n,i)
smi="".join(es[2].split("'"))
loop=es[3]
factor=es[4]
substance="".join(es[13].split("'"))

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
parent=routeMatrix["smiles"]

routes=list(parent['sub'].keys())

#print("parent",parent,"<---------")

conn=sqlite3.connect(db)
cur=conn.cursor()

cur.execute('create table if not exists searchList(id integer primary key, user text, uname text, email text, smiles text, cSmiles text, date text, substance text, loop integer, factors text, options text);')
cur.execute('create table if not exists parent(id integer, total text, sub text);')
conn.commit()


#smiles=getSmiles(routes[0],all)[-1][-1]
mol=Chem.MolFromSmiles(smi)
cSmiles=Chem.MolToSmiles(mol)
dd=datetime.datetime.now().strftime('%Y%m%d%H%M%S')

print(opts)
cur.execute(f'insert into searchList(user, uname, email, smiles, cSmiles, date, substance, loop, factors, options) values("{user}","{uname}","{email}","{smi}","{cSmiles}",{dd},"{substance}",{loop},"{factor}","{opts}")')
conn.commit()

sql=f'select id from "searchList" where date = "{dd}";'
ret=cur.execute(sql)
conn.commit()

total=dicToString(parent['total'])
sub=dicToString(parent['sub'])
id=ret.fetchall()[0][0]

print("total",total)
print(" sub" ,sub)
sql=f'insert into parent(id, total, sub) values({id},"{total}","{sub}")'
print(sql)
cur.execute(sql)
conn.commit()

print("----------------------->")

sTable=f"searchTable{id}"

cur.execute(f'create table if not exists "{sTable}" (route int, connect_list text, smiles_list text, info_list text);')
conn.commit()

#komai
#print(routes)
for route in routes:
    info=getInfo(route,all_info)
    smiles=getSmiles(route,all)
    connect,start=initConnect(smiles)
    s1=listToString(connect)
    s2=listToString(smiles)
#    print("---->",route)
#    print("connect",connect)
#    print("connect(s1)",s1)
#    print(smiles)
#    print("smiles(s2)",s2)
#    print(info)
    s3=listToString(info)
    sql=f'insert into {sTable} values({route},"{s1}","{s2}","{s3}");'
#    print("sql",sql)
    cur.execute(sql)
    conn.commit()

print("end")

print(f"substance:{substance}###")
conn.close()
