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

max_route=7
debug=False
debug=True

version="1.3 12062024"

_web=False
_web=True

if _web:
    dr='/var/www/html/public/images/tmp/' # working directory
else:
    dr='tmp/' # working directory

_magic_x=650
_magic_x=550

class openReport:
    def __init__(self,file,scale):
        self.left=210*mm*0.05
        self.top=297*mm*0.05
        self.height=210*mm*0.9
        self.width=297*mm*0.9
        self.scale=scale
        self.mg=20

        ext=file.split(".")[1]
        if ext=='pdf':
            print("open pdf:",file)
            self.type=ext
            self.page=canvas.Canvas(file,pagesize=landscape(A4))
            self.font_name="Times-Bold"
            self.font_size=4*mm
            self._font_size=4*mm
            self.P2PRatio=1.0
        else:
            self.type='ppt'
            self.page=Presentation()
            self.page.slide_height=Inches(8.27)
            self.page.slide_width =Inches(11.69)
            blank_slide_layout = self.page.slide_layouts[6]
            self.slide = self.page.slides.add_slide(blank_slide_layout)
            self.file=file
            self.font_name="Times-Bold"
            self.font_size=12 # 14<-16 
            self._font_size=12 #
#            self.P2PRatio=self.width/297*mm
            self.P2PRatio=self.page.slide_width/(297*mm)

    def getSpan(self):
#       scale=_magic_x/(span)
        return _magic_x/self.scale
     
    def setFont(self,font_name=None,font_size=None):# size is relative value to default
        if self.type=='pdf':
            if font_name:
                self.font_name=font_name
            if font_size:
                self.font_size=int(self._font_size*font_size)
                print("in setFont",self.font_size)
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

    def getImageSize(self,fn):
        if _web:
            ret=subprocess.run(['/usr/bin/identify',fn],capture_output=True,text=True).stdout
        else:
            ret=subprocess.run(['identify',fn],capture_output=True,text=True).stdout
        xy=ret.split(" ")[2].split("x")
        return float(xy[0]),float(xy[1])

    def drawImage(self,fn,x,y):
        if 'branch' in fn:
            return 0,0

        ix,iy=self.getImageSize(fn)

        ppx=page.left+x*page.scale;ppy=page.top+y*page.scale
        if self.type=='pdf':
            drawing=svg2rlg(fn)
            drawing.renderScale=self.scale
            renderPDF.draw(drawing, canvas=self.page,x=ppx,y=ppy)
        else:
            key="2.0px"
            kkk="4.0px"

            with open(fn,"r") as f:
                all=f.readlines()

            gn=fn.split(".svg")[0]+".jpg"
            mag=360
            if _web:
                subprocess.run(['/usr/bin/convert','-density',f'{mag}',fn,gn])
            else:
                subprocess.run(['convert','-density',f'{mag}',fn,gn])

            if ix>iy:
                sx=1.0;sy=iy/ix
            else:
                sy=ix/iy;sy=1.0
            ppx=page.left+x*page.scale;ppy=page.top+y*page.scale
            width=self.scale*self.P2PRatio*ix;height=self.scale*self.P2PRatio*iy
            x0,y0=ppx*self.P2PRatio,self.page.slide_height-(ppy)*self.P2PRatio-height
#            x0,y0=x*self.sc,self.page.slide_height-(y+iy)*self.sc
#            width=self.scale*self.P2PRatio*0.9*sx;height=self.scale*self.P2PRatio*0.9*sy
            self.slide.shapes.add_picture(
                gn,left=x0,top=y0,
                width=width, height=height
    )
#            printF(gn,x0,y0,width,height)
#            printF("(x,y):",x,y)
        return ix,iy


    def drawLine(self,pt):
        x0=self.left+pt[0]*self.scale;y0=self.top+pt[1]*self.scale
        x1=self.left+pt[2]*self.scale;y1=self.top+pt[3]*self.scale
        if self.type=='pdf':
            self.page.line(x0,y0,x1,y1)
        else:
            x0,y0=x0*self.P2PRatio,self.page.slide_height-y0*self.P2PRatio
            x1,y1=x1*self.P2PRatio,self.page.slide_height-y1*self.P2PRatio
            putLine(x0,y0,x1,y1,self.slide.shapes)

    def drawArrow(self,pt):
        x0=self.left+pt[0]*page.scale;y0=self.top+pt[1]*page.scale
        x1=self.left+pt[2]*page.scale;y1=self.top+pt[3]*page.scale
        x2=self.left+pt[4]*page.scale;y2=self.top+pt[5]*page.scale
        x3=self.left+pt[6]*page.scale;y3=self.top+pt[7]*page.scale

        if self.type=='pdf':
            self.page.line(x0,y0,x1,y1)
            self.page.line(x1,y1,x2,y2)
            self.page.line(x1,y1,x3,y3)
        else:
            x0,y0=x0*self.P2PRatio,self.page.slide_height-y0*self.P2PRatio
            x1,y1=x1*self.P2PRatio,self.page.slide_height-y1*self.P2PRatio
            x2,y2=x2*self.P2PRatio,self.page.slide_height-y2*self.P2PRatio
            x3,y3=x3*self.P2PRatio,self.page.slide_height-y3*self.P2PRatio
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

    def drawStringR(self,x,y,string,center=False,adj=True):
        xx=self.left+x*self.scale
        if center:
            xx-=stringWidth(string,self.font_name,self.font_size)/2
        yy=self.top+y*self.scale
        if self.type=='pdf':
            self.page.drawString(xx,yy,string)
        else:
            xx=xx*self.P2PRatio
            yy=self.page.slide_height-yy*self.P2PRatio-self.font_size
            if adj:
                pw=Cm(0.1)*self.font_size*(len(string)+2)*0.20
                py=Cm(0.1)*self.font_size*0.5
                shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,py)
            else:
                pw=self.page.slide_width*0.9
                sw=self.font_size*(len(string)+2)*0.18/pw+1
                shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,py)

            tf=shape.text_frame
            tf.word_wrap=True
#            tf.auto_size=MSO_AUTO_SIZE.SHAPE_TO_FIT_TEXT
            p=tf.paragraphs[0]
            run=p.add_run()
            run.text=string
            run.font.name=self.font_name
            run.font.size=Pt(self.font_size)
            sw=self.font_size*(len(string)+2)*0.18/pw+1
            shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,py)

    def drawStringX(self,x,y,string,center=False,adj=True):
        if self.type=='pdf':
            if center:
                x+=self.getSpan()*self.scale/2-stringWidth(string,self.font_name,self.font_size)/2
            self.page.drawString(x,y,string)
        else:
            xx=x*self.P2PRatio
            yy=self.page.slide_height-y*self.P2PRatio-self.font_size
            if adj:
#                shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,Cm(1))
                pw=Cm(0.1)*self.font_size*(len(string)+2)*0.20
                py=Cm(0.1)*self.font_size*0.5
                shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,py)
            else:
                pw=self.page.slide_width*0.9
                sw=self.font_size*(len(string)+2)*0.18/pw+1
                shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,py)

            tf=shape.text_frame
            tf.word_wrap=True
#            tf.auto_size=MSO_AUTO_SIZE.SHAPE_TO_FIT_TEXT
            p=tf.paragraphs[0]
            run=p.add_run()
            run.text=string
            run.font.name=self.font_name
            run.font.size=Pt(self.font_size)

    def close(self):
        if self.type=='pdf':
            self.page.save()
        else:
            self.page.save(self.file)
    def isPdf(self):
        if self.type=='pdf':
            return True
        return False

def fm(x):
    return "{:.1f}".format(x)

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

def cropSvg(fn):
    if not os.path.isfile(fn):
        return [0,0]

    fno=fn+".tmp";fo=open(fno,"w")

    with open(fn,"r") as fp:
        for l in fp:
            fo.write(l)
            if 'viewBox' in l:
                w=l.split("width='")[1].split('px')[0]
                h=l.split("height='")[1].split('px')[0]

    x=[float(w),0]
    y=[float(h),0]
    fo.close()
    inPath=False
    with open(fn,"r") as fp:
        for l in fp:
            if "d='" in l:
                inPath=True
            if inPath:
                if "d=" in l:
                    if '/>' in l:
                        inPath=False
                    l=l.split("d=")[1].split("'")[1].split("\n")[0]
                f,v=extF(l)
                if f:
                    px=[]
                    for p in v:
                        try:
                            px.append(float(p))
                        except:
                            px.append(0.0)
                        if len(px)==2:
                            repv(x,y,px)
                            px=[]
                else:
                    inPath=False
    mg=10

    x0=int(x[0]-mg);x1=int(x[1]+mg)
    y0=int(y[0]-mg);y1=int(y[1]+mg)

    sx=int(x1-x0)
    sy=int(y1-y0)

    ff=open(fn,"w")
    with open(fno,"r") as fp:
        for l in fp:
            if "<rect " in l:
                ff.write(l)
            elif "stroke-width:" in  l:
                ff.write(l.replace("stroke-width:2.0","stroke-width:1.0"))
            elif "viewBox" in  l:
                ff.write(f"width='{sx}px' height='{sy}px' viewBox='{x0} {y0} {sx} {sy}'>\n")
            else:
                ff.write(l)
    ff.close()
    return [sx,sy]

def smile2svg(smiles,name,dr=dr,size=200):
#    if os.path.exists(dr):
#        shutil.rmtree(dr)
    if not os.path.exists(dr):
        os.makedirs(dr)
    
    mol=Chem.MolFromSmiles(smiles)
    mw=mol.GetNumAtoms()
    fc=int(math.sqrt(mw)/1.2)

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

    r=cropSvg(dr+"/"+name)

    return fc,r

#    os.remove(fo)

def Debug(c,head=None,sp=None):
    return
    for n,i in enumerate(c):
        if head:
            print(head,end=":")
        if sp:
            print(n,i[sp],end="<--")
        print(n,i)
#    exit()

def printF(head,*args):
    print(head,end='->')
    for arg in args:
        print(round(arg,2),end=" ")
    print()

def makeSgv(smiles,type=0):
    size=100
    fsize={}
    fcs=[]
    for n,proc in enumerate(smiles):
        fc={}
        for m,comp in enumerate(proc):
            if ">" in comp:
                continue
            fn="e"+str(n+1)+"x"+str(m+1)+".svg"
            fc[m],fsize[fn]=smile2svg(comp,fn,size=size)
        fcs.append(fc)

    return fcs,fsize

def arrowR(x0,y0,x1,y1):

    l=math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
    sx=10
    pi=3.1415/12.0 # 30 deg

    x2=x0;y2=y0
    x3=x1;y3=y1

    a=math.atan2(y1-y0,x1-x0)
    x4=x3-sx*math.cos(a+pi);y4=y3-sx*math.sin(a+pi)
    x5=x3-sx*math.cos(a-pi);y5=y3-sx*math.sin(a-pi)

    return [x2,y2,x3,y3,x4,y4,x5,y5]

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

    xx=x3-sx*math.cos(a+pi);yy=y3-sx*math.sin(a+pi)
    ret=[x2,y2,x3,y3,xx,yy]

    xx=x3-sx*math.cos(a-pi);yy=y3-sx*math.sin(a-pi)

    return ret+[xx,yy]

def prevArrow(x0,y0,x1,y1,scale,skip=0.0,ugly_patch=True):
#
# ugly patch
#
    if ugly_patch and x0>x1:
        x1=x0+scale*3
        y1=y0
# end

    l=math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
    rt=(l-scale*skip)/l
    sx=10
    pi=3.1415/12.0 # 30 deg

    x2=x1+(x0-x1)*(rt+0.10);y2=y1+(y0-y1)*(rt-0.10)
    x3=x0+(x1-x0)*(1.0+0.10);y3=y0+(y1-y0)*(1.0-0.10)

    a=math.atan2(y1-y0,x1-x0)
    x4=x3-sx*math.cos(a+pi);y4=y3-sx*math.sin(a+pi)
    x5=x3-sx*math.cos(a-pi);y5=y3-sx*math.sin(a-pi)

    return [x2,y2,x3,y3,x4,y4,x5,y5]

def getName(a,b,dr=dr):
    fn=dr+"e"+str(a+1)+"x"+str(b)+".svg"#<-- need modification
    return(fn)

def getGoal(smiles,connect):
    goal=[]
    l=len(smiles)
    for i in range(l):
        if connect[i][3]==-1:
            goal.append(i)
    return goal

def _ccn(i):
    if i<10:
        return 1
    if i<100:
        return 2
    return 3

# komai..
def getPropString(page,m):
    if len(m)<2:
        s=str(m[0])
        w=stringWidth(s,page.font_name,page.font_size)*page.scale
        return s,w
    s=str(m[0])+"-"+str(m[-1])
    w=stringWidth(s,page.font_name,page.font_size)*page.scale
    return s,w

def preDraw(page,init,proc,branch):

    ixx=0
    yy=1
    pos={}
    deep={}

    for i in proc[init]:
        pos[i]=[ixx,yy]
        s,w=getPropString(page,branch[i])
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
                s,w=getPropString(page,branch[i])
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

def similarityOnRoute(page,init,proc,branch,allData):

    pos,deep=preDraw(page,init,proc,branch)
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
#            print(j,"-------->",allData.keys())
#            print(allData[int(j)],"<--------")
            ss=allData[int(j)]['smiles']
            sl=len(ss)
            if sl>ssMax:
                ssMax,tg_route=sl,j
        fpMin=10
        ss=allData[int(j)]['smiles']
#        print("<=============================>",ss)

        if not fps:
            fps=FingerprintMols.FingerprintMol(Chem.MolFromSmiles(ss[-1][-1]))

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
    return sim

def easyTrans(page,x,y,half=False,lf=1.2):
    nx,ny=x/page.scale,y/page.scale+(page.height-page.font_size)/page.scale
    if half:
        ny+=page.font_size/page.scale/2.0
    return nx,ny,y-page.font_size*lf

def head_page(head,page,hope):
    ox=40
    oy=480
    mag=3

    nx,ny,_=easyTrans(page,0,0)
    ll=len(head[0])

    page.setFont(font_size=mag)
    if ll>7:
        page.setFont(font_size=mag*7.0/float(ll))
    page.drawStringR(nx,ny,head[0])
    if ll>7:
        page.setFont(font_size=mag)

    scale=page.scale
#
    smiles=head[1][-1][-1]
    fn="head.svg"
    x0=ox
    smile2svg(smiles,fn,size=100)
    scale=page.scale
    page.scale=1.0
    page.drawImage(dr+"/"+fn,x0,page.top)
    page.scale=scale

    pt=0;l=len(smiles)
#    py=py-mag*page.scale
    mg=2
    page.setFont(font_size=mg)
    wd=page.getProperLength(mg,smiles)
#    wd=7
    tl=smiles
    py=-page.height
    while pt<l:
        end=min(pt+wd,l)
        tl=smiles[pt:end]
        nx,ny,_=easyTrans(page,page.width/2,py)
        page.drawStringR(nx,ny,tl,center=True)
        pt+=wd
        py+=page.font_size

    if type(head[2][0])==type(True):
        nx,ny,_=easyTrans(page,page.width/2,-page.height/2)
        page.drawStringR(nx,ny,head[2][1])
        page.close()
        exit()

    proc=head[-1][0]
    branch=head[-1][1]
    init=head[-1][2]
    sear=head[-1][3]

    oy=550
    sy=50
    sy=35

#    print("font",page.font_size,page.scale)
    py=0
    nx,ny,py=easyTrans(page,page.width/2,py)
    page.drawStringR(nx,ny,"Reaction route summary",center=True)
#
# draw initial one
#
    ll=0
    for i in proc[init]:
        _,w=getPropString(page,branch[i])
        ll+=w

    sx=14
    nx,ny,_=easyTrans(page,page.width*0.8,py)
    page.drawStringR(nx,ny,"route #id")
    
    mes="reactions"
    nx,ny,py=easyTrans(page,page.width*0.2,py)
    page.drawStringR(nx,ny,mes)

    pos,deep=preDraw(page,init,proc,branch)

#
# the first line
    name={}
    for j in proc:
        if deep[j] in name:
            name[deep[j]]+=","+str(j)
        else:
            name[deep[j]]=str(j)
#
    y_skip=1.0
    if len(name)>10:
        y_skip=10/float(len(name))

    t_max=0
    for i in pos:
        if t_max<pos[i][0]:
            t_max=pos[i][0]

    x_skip=1
    t_ox=ox=page.width*0.2
    ll=0
    for i in pos:
        ix=pos[i][0]
        iy=pos[i][1]
        s,w=getPropString(page,branch[i])
# dummy 2 lines
#        s=t_s
#        ix=0
        pt=ix+w
        if ll<pt:
            ll=pt

    box=[32,40]
    if ll>box[1]:
        x_skip=box[0]/ll
    elif ll>box[0]:
        t_ox=ox-(ll-box[0])*sx
        x_skip=box[0]/ll

#    if x_skip<y_skip:
#        page.setFont(font_size=x_skip*mg)
#    else:
#        page.setFont(font_size=y_skip*mg)
    page.setFont(font_size=y_skip*mg)

#    page.setFont(font_size=2.5)
    tmp_scale=page.scale
    page.scale=1
#    print("fontsize",page.font_size);

#    x_skip=y_skip=1.0
#    print("pos",pos)
#    print("bos",branch)
#    print("proc",proc)
    for i in pos:
        ix=pos[i][0]
        iy=pos[i][1]
        s,w=getPropString(page,branch[i])
#        nx,ny,_=easyTrans(page,t_ox+(ix*x_skip)*sx,py+(-iy)*page.font_size)
        nx,ny,_=easyTrans(page,t_ox+(ix*x_skip)*sx,py+(1-iy)*page.font_size)
        page.drawStringR(nx,ny,s)
#        page.drawLine([nx,ny,nx+w,ny])

    for i in proc:
        draw=False
        for j in proc[i]:
            s,w=getPropString(page,branch[j])
            x1=pos[j][0]
            y1=pos[j][1]
# 
#        nx,ny,_=easyTrans(page,t_ox+(ix*x_skip)*sx,py+(1-iy)*page.font_size*x_skip)
#        page.drawLine([nx,ny,nx+w,ny])
#
            nx1,ny1,_=easyTrans(page,t_ox+((x1-0.2)*x_skip)*sx,
                    py+(1.25-y1)*page.font_size)
            if draw:
               # page.drawLine([nx0,ny0,nx1,ny1])
                page.drawArrow(arrowR(nx0,ny0,nx1,ny1))
            nx0,ny0,_=easyTrans(page,t_ox+((x1+0.5)*sx)*x_skip,
                    py+(1.25-y1)*page.font_size)
            nx0+=w
            draw=True

    page.scale=tmp_scale
    bt=0
    for j in name:
        if bt<j:
            bt=j
#        if len(name[j])>20:
#            page.drawString(_magic_x+sx*2,oy-(j+1)*sy,name[j][:22])
#            page.drawString(_magic_x+sx*2,oy-(j+1)*sy-sy/2,name[j][22:])
#        else:
        nx,ny,tpy=easyTrans(page,page.width*0.8,py+(1-j)*page.font_size)
        page.drawStringR(nx,ny,name[j])
    tpy=py+(1-bt)*page.font_size
# now
#    if oy-(bt+2)*sy-sy<0:
#        py=oy-sy
#    else:
#        py=oy-(bt+2)*sy*t_skip-sy

    tpy-=page.font_size
    nx,ny,tpy=easyTrans(page,page.width*0.2,tpy-page.font_size)
    page.drawStringR(nx,ny,sear)

    page.setFont(font_size=1)

    add="Routes start at similarity of "
    for i in sorted(hope):
        add+="{:.2f}".format(i[1])+"("+str(i[0])+"),"
        if len(add)>90:
            nx,ny,_=easyTrans(page,page.width*0.2,tpy-page.font_size)
            page.drawStringR(nx,ny,add)
            py-=sy*0.01
            add=""

    if len(add)>0:
       nx,ny,_=easyTrans(page,page.width*0.2,tpy-page.font_size)
       page.drawStringR(nx,ny,add)

    page.stroke()
    page.setFont(font_size=1)

def getList(s):
    ret=[]
    s=s.split('[')[1].split(']')[0].split(',')
    for i in s:
        ret.append(int(i))
    return ret

def whichIsFirstMaterial(s,mode):
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
        if mw>maxWeight:
            maxWeight=mw
            r=i
    return r

def checkIn(pp,rt):
    ret=[]
    if type(pp)==int:
        return checkIn2(pp,rt)

    for r in rt:
        if checkIn3(pp,rt[r]):
            ret.append(r)
    return ret

def checkIn3(p1,p2):

    ll=len(p1)

    for i in range(len(p2)-ll+1):
        flag=True
        for j in range(ll):
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
        if j==s:
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

    ok={}
    for r in routes:
        pp=[]
        rr=[]
        for s in routes[r]:
            pp.append(s)
            key=branchCheck(pp,branch)
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

def strToParent(ss):
    ret={}
    n={1:'total',2:'sub'}

    for i in range(1,3):
        r={}
        for s in ss[i].split(';')[0:-1]:
            key=int(s.split(':')[0])
            values=s.split(':')[1]
            p=[]
            for q in values.split(',')[0:-1]:
                p.append(int(q))
            r[key]=p
        ret[n[i]]=r
    return ret

def strToSmiles(ss):
    ret=[]
    for s in ss.split('##')[0:-1]:
        item=[]
        for p in s.split(';')[0:-1]:
            item.append(p)
        ret.append(item)
    return ret

def dummyParent(p):
    print("total",p['total'])
    for k in p['total'].keys():
        continue
    return k


def strToList(ii):
        w=[]
        for i in ii.split(',')[0:-1]:
            w.append(int(i))
        return w

def strToConnect(ss):
    ret=[]
    for s in ss.split('##')[0:-1]:
        item=[]
        ii=s.split(';')
        item.append(ii[0])
        item.append(strToList(ii[1]))
        item.append(strToList(ii[2]))
        for n in range(3,7):
            item.append(int(ii[n]))
#        item.append('A')
        item.append(strToList(ii[7]))
        ret.append(item)
    return ret

def dump(fall,pid,odir):
    with open(odir+str(pid)+'db.txt','w') as fp:
        for xx in fall:
            fp.write('route '+str(xx[0])+'\n')
            fp.write(xx[2]+'\n')
            fp.write(xx[3]+'\n')

def repv(x,y,p):
    if x[0]>p[0]:
        x[0]=p[0]
    if x[1]<p[0]:
        x[1]=p[0]
    if y[0]>p[1]:
        y[0]=p[1]
    if y[1]<p[1]:
        y[1]=p[1]

def frac(ss,k):
    rr=[]
    for s in ss:
        if k in s:
            tt=s.split(k)
            for t in tt:
                if not t=="":
                    rr.append(t)
        else:
            rr.append(s)
    return rr

def xsplit(v,s):
    vv=[]
    for rx in v:
        for v in rx.split(s):
            if v!='':
                vv.append(v)
    return vv

def extF(l):
    if "class='atom" in l:
        return False,[]

    v=frac([l],"M")
    v=frac(v,"Q")
    v=frac(v,"L")

    vv=xsplit(v,',')
    vvv=xsplit(vv,' ')
    if "/>" in vvv[-1]:
        return False, False
    return True,vvv

def getSubObj(page,n,c,fsize,onset=False):

    ret={}
    if onset:
        fn=getName(n,c[n][6]+1,dr='')
    else:
        fn=getName(n,c[n][1][0]+1,dr='')

    if fn in fsize:
        r=fsize[fn]
    else:
        r=[page.mg,page.mg]

    ret['type']=1;
#    py=page.height-page.mg*5
    py=page.height
#komai...Need..
#    pos=[];pos.append([0,(py-r[1]/2.0-page.font_size*3)/page.scale])
    wy=-r[1]/2.0
    pos=[];pos.append([0,wy])
    ret['pos']=pos
    ret['svg']=[fn];ret['size']=r
#    komai...
    return ret

def calcBranch(obj,tag):
    base=0
    x=obj[tag]['pos'][0]
    y=obj[tag]['pos'][1]+obj[tag]['size'][1]/2.0

    for o in obj:
        tx=o['pos'][0]
        ty=o['pos'][1]
        if tx<x:
            if base>ty:
                base=ty
    return x,base

def getCatObj(page,n,c,fsize):
    l=0
    fns=[]
    height=0
    ret={}

    for j in c[n][-1]:# catalytic materials
        fn=getName(n,j+1,dr='')
        fns.append(fn)
        r=fsize[fn]
# will be check nicely
        h=(page.mg+r[1])
        if height<h:
            height=h
        l+=(page.mg+r[0])
    l+=(page.mg)
    if len(c[n][-1])<1:
        l+=page.mg*5
#Need..
    ret['type']=2;ret['pos']=[[0,0],[l,0]];ret['svg']=fns;ret['size']=[l,height]
    return ret

def setItem(a,b,c,d):
    ret={}
    ret['type']=a;ret['pos']=b;ret['svg']=c;ret['size']=d
    return ret

def shiftForBranch(objs):

    for i in range(1,len(objs)):
        x0=objs[i][-1]['pos'][0][0]
        x1=objs[i][ 0]['pos'][0][0]+objs[i][ 0]['size'][0]
        ty=-_py;by=_py
        for item in objs[i]:
            if item['type']==2:
                yy=item['pos'][0][1]+item['size'][1]
                if ty<yy:
                    ty=yy
        for item in objs[0]:
            if item['type']==1:
                p0=item['pos'][0][0]
                p1=item['pos'][0][0]+item['size'][0]
                if p0<x1 and x0<p1:
                    yy=item['pos'][0][1]-item['size'][1]/2
                    if by>yy:
                        by=yy
        shift=ty-by
        fool=0
        for item in objs[i]:
            if fool==0:
                item['pos'][0][1]-=shift
                fool+=1
                continue
            item['pos'][0][1]-=shift
            if item['type']==2:
                item['pos'][1][1]-=shift

def lineBraker(objs):
    #komai
    for j in range(len(objs[0])-1,-1,-1):
        item=objs[0][j]
        xx=item['pos'][0][0]+item['size'][0]
        if xx>page.width/page.scale*0.9:
            cutIt(objs,item['pos'][0][0])

def cutIt(objs,ok):
    theWidth=page.width/page.scale*0.9
    by=_py
    ty=-py
    for obj in objs:
        for item in obj:
            xx=item['pos'][0][0]+item['size'][0]
            if xx<theWidth:
                if item['type']==1:
                    yy=item['pos'][0][1]-item['size'][1]/2.0
                else:
                    yy=item['pos'][0][1]
                if by>yy:
                    by=yy
            else:
                if item['type']==2:
                    yy=item['pos'][0][1]+item['size'][1]
                else:
                    yy=item['pos'][0][1]+item['size'][1]/2.0
                if ty<yy:
                    ty=yy

    shift=ty-by+page.mg
#    print("shift:",shift,ty,by,page.mg)
    for obj in objs:
        for item in obj:
            xx=item['pos'][0][0]+item['size'][0]
            if xx>theWidth:
                item['pos'][0][0]-=(ok-page.width*0.1)
                item['pos'][0][1]-=shift
                if item['type']==2:
                    item['pos'][1][0]-=(ok-page.width*0.1)
                    item['pos'][1][1]-=shift

def easyGetBottom(objs):
    hy=_py
    top=-_py
    for obj in objs:
        for item in obj:
            if item['type']==1:
                yy=item['pos'][0][1]-item['size'][1]/2.0
                if hy>yy:
                    hy=yy
            else:
                yy=item['pos'][0][1]+item['size'][1]
                if top<yy:
                    top=yy
    return top-hy,top,hy

def makeObject(page,connect,fsize):

    line=0
    objs=[];hint=[];tmp_goal=[];xx=0
#
# find goal
#
    for i in range(len(connect)):
        if connect[i][3]<0:
            goal=n=i
            break
#
# temporal base y=0
#
    tag=0;obj=[]
#    fn=getName(n,connect[n][6]+1,dr='') # goal
    hint.append([]) # to which the line connect (line,tag)

    l=len(connect[n][2])

    if l>1:
        tmp_goal.append([connect[n][2][1],line,tag])

    while n>-1:
# substance object
        xx-=page.mg
        item=getSubObj(page,n,connect,fsize)
        # set postion
        xx-=item['size'][0]
        item['pos'][0][0]=xx
        obj.append(item);tag+=1

# catalysis object
        l=len(connect[n][2])
        if l>1:
            tmp_goal.append([connect[n][2][1],line,tag])

        xx-=page.mg
#        xx-=page.mg
        item=getCatObj(page,n,connect,fsize)
        # set postion
        item['pos'][1][0]=xx
        xx-=item['size'][0]
        item['pos'][0][0]=xx
        obj.append(item);tag+=1

        if l<1:
            xx-=page.mg
            item=getSubObj(page,n,connect,fsize,onset=True)
            # set postion
            xx-=item['size'][0]
            item['pos'][0][0]=xx
            obj.append(item);tag+=1
            objs.append(obj);line+=1
            break
        n=connect[n][2][0]

    while len(tmp_goal)>0:
        r=tmp_goal.pop()
        n=r[0];obj=[];tag=0
#        hint.append([r[1],r[2]]) # to which the line connect (line,tag)

        item=objs[r[1]][r[2]-1] # connecting position
#        x0=item['pos'][0][0]-page.mg
        x0=item['pos'][0][0]-page.mg
        y0=item['pos'][0][1]+item['size'][1]/2
#        y0=item['pos'][0][1]+item['size'][1]/2
        item=objs[r[1]][r[2]]
        ln=item['size'][0]
        xx=x0-ln
        item=setItem(2,[[xx,0],[x0,y0]],[],[ln,0]);
        obj.append(item);tag+=1

        n=connect[n][2][0]
        xx-=page.mg
        item=getSubObj(page,n,connect,fsize)
        xx-=item['size'][0]
        item['pos'][0][0]=xx
        obj.append(item);tag+=1

        if len(connect[n][2])<1:
            n=-1
        else:
            n=connect[n][2][0]

        while n>-1:
# catalysis object
            l=len(connect[n][2])
            if l>1:
                tmp_goal.append([connect[n][2][1],line,tag])

            xx-=page.mg

            item=getCatObj(page,n,connect,fsize)
            # set postion
            item['pos'][1][0]=xx
            xx-=item['size'][0]
            item['pos'][0][0]=xx
            obj.append(item);tag+=1
# substance object
            xx-=page.mg
            item=getSubObj(page,n,connect,fsize)
            xx-=item['size'][0]
            item['pos'][0][0]=xx
            obj.append(item);tag+=1
# up to now
            if l<1:
                objs.append(obj);line+=1
                break
            n=connect[n][2][0]

    xx=0
    for obj in objs:
        for item in obj:
            if xx>item['pos'][0][0]:
                xx=item['pos'][0][0]

    for obj in objs:
        for item in obj:
            for i in range(len(item['pos'])):
                item['pos'][i][0]-=xx

    if len(objs)>1:
        shiftForBranch(objs)

    lineBraker(objs)

    return objs

##############################################
#################### main ####################
##############################################

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
com='-chem chemical_name, -ppt for power point -n [1,2,3] to specify routes -p proc_per_line -d input_path, -one_size : rendering at the same size, -product_only : categolize routes with respct to reaction product only, -show_subsets : show routes including subsets'
com+=" -summary 0/1/2"
ops={'-h':com}
skip=False
pid=-1

debug_log="/var/www/html/public/images/report/readDb_web.log";

scale=0.15
scale=0.4
scale=0.2
fp=open("readDb.log","w")
for n,op in enumerate(sys.argv):
    fp.write(f" {op}")
    if skip:
        skip=False;continue
    match op:
        case '-h':
            print(sys.argv[0],ops['-h']);exit()
        case '-v':
            print(version);exit()
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
        case '-n':
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
        case '-force':
            forced=True

    fp.write("\n")
fp.close()
db="sList.db"
conn=sqlite3.connect(db)
cur=conn.cursor()

#print("options:",pid,oper)

if _web and debug:
    with open(debug_log,"w") as fp:
        fp.write(f"{oper}\n")

if oper=='drop':
    if int(pid)<0:
        exit()
    cur.execute(f"""drop table searchTable{pid}""")
    cur.execute(f"""delete from searchList where id={pid}""")
    cur.execute(f"""delete from parent where id={pid}""")
    conn.commit()
    conn.close()
    exit()

#komai
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

chemical=ret[0][0]

if ppt:
    output_file=str(pid)+'.pptx'
else:
    output_file=str(pid)+'.pdf'

log_name="readDb"+str(pid)+oper+".log"

if _web:
    if os.path.isdir(input_dir):
        os.makedirs(input_dir,exist_ok=True)
    log_name=input_dir+"/readDb"+str(pid)+"normal.log"
    output_file=input_dir+output_file

print(output_file)
if os.path.isfile(output_file) and not forced:
    conn.close()
    exit()

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

page=openReport(output_file,scale)

allId={}
cur=conn.cursor()
ret=cur.execute(f"""select loop from searchList where id={pid};""").fetchall()

print("ret..",ret)
if (len(ret)<1):
    print(f"not availabale id:{pid}")
    exit()

maxLoop=ret[0][0]

sTable=f"searchTable{pid}"
sql=f'select * from "{sTable}";'
ans=cur.execute(sql)

if oper=='db':
    dump(ans.fetchall(),pid,input_dir)
    with open(log_name,"a") as fp:
        fp.write(f"mission ok\n")
    exit()

allData={}
routes=[]
for xx in ans.fetchall():
    s={}
    routes.append(int(xx[0]))
    s['connect']=strToConnect(xx[1])
    s['smiles']=strToSmiles(xx[2])
    s['info']=strToSmiles(xx[3])
    allData[xx[0]]=s

ret=cur.execute(f"""select * from parent where id={pid};""")
parent=strToParent(ret.fetchall()[0])
ppt=dummyParent(parent)

hint=deepArrange(parent[show])
hint.append(ppt)

#print(head)

sim=similarityOnRoute(page,hint[2],hint[0],hint[1],allData) 
print(sim,"sim")
#easy=str(len(parent['total']))+" variations from "+str(len(sim))+" routes over "+str(maxLoop)+" queries"
easy=str(len(sim))+" routes over "+str(maxLoop)+" queries"

hint.append(easy)

if parent!=False:
    head=[chemical,allData[int(routes[0])]['smiles'],hint]
else:
    head=[chemical,allData[int(routes[0])]['smiles'],hint]
head_page(head,page,sim)

with open(log_name,"w") as fp:
    fp.write("head_page done\n")

sTable=f"searchTable{pid}"
sql=f'select * from "{sTable}";'
ans=cur.execute(sql)

#oy=py-page.mg*5
oy=page.height
for route in routes:
    now=str(int(time.time()-start))+" sec"
    with open(log_name,"a") as fp:
        fp.write(f"{now}:route {route} start\n")
    connect=allData[route]['connect']
    smiles=allData[route]['smiles']
    info=allData[route]['info']

    fc,image_size=makeSgv(smiles)

    if not _fc:
        fc=_fc

    page.setFont()
    xx=ox

#    p.append([xx,bt,up])
    objs=makeObject(page,connect,image_size)

    hy,top,btm=easyGetBottom(objs)
    top+=page.font_size/page.scale+page.mg

    if oy-hy*page.scale-page.mg<0:
        page.stroke()
        oy=page.height

    oyy=oy/page.scale
    for obj in objs:
        for item in obj:
            x0=item['pos'][0][0];y0=item['pos'][0][1]+oyy-top
            if item['type']==1:
                page.drawImage(dr+"/"+item['svg'][0],x0,y0)
            else:
                x1=item['pos'][1][0];y1=item['pos'][1][1]+oyy-top
                xx=x0;

                for n,fn in enumerate(item['svg']):
                    xx=xx+page.mg
                    r=image_size[fn]
                    page.drawImage(dr+"/"+fn,xx,y0+page.mg)
                    xx+=r[0]
                page.drawArrow(arrow(x0,y0,x1,y1,page.scale))
## komai111
    page.drawStringR(ox,(oy-page.font_size)/page.scale,"Route"+str(route))
#    page.drawLine([0,oyy,ox,oyy])
    print("obj",obj)
    print("---------->  Route"+str(route))
    oy-=(hy*page.scale+page.mg)
    now=str(int(time.time()-start))+" sec"
    with open(log_name,"a") as fp:
        fp.write(f"{now}:route {route} done\n")
    if route>max_route and debug==True:
        break

with open(log_name,"a") as fp:
    fp.write(f"mission ok\n")
if page.isPdf():
    conn.close()
    page.close()
    exit()

conn.close()
page.close()
