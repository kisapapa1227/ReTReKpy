import os
import subprocess


var="/var/www/html/public/images"
limit=2000000 # 2GB
limit=200000 # 200MB
limit=20000000 # 20GB

def getInfo(dr):
    fn=dr+"/svgFile.info"
    d=False;s=0
    with open(fn) as f:
        lines=f.readlines()
    if 'date' in lines[-2]:
        d=lines[-2].split('\n')[0].split(':')[1]
    if 'size' in lines[-1]:
        s=lines[-1].split('\n')[0].split(':')[1]
    return d,s

def getNextDir(dr):
    ret=[]
    for d in os.listdir(dr):
        if os.path.isdir(dr+"/"+d):
            ret.append(dr+"/"+d)
    return ret

l2=[]
ret=getNextDir(var)
for r in ret:#depth 1
    l2=l2+getNextDir(r)

l1=l2;l2=[]
for r in l1:#depth 1
    l2=l2+getNextDir(r)

l1=l2;l2=[]
for r in l1:#depth 1
    l2=l2+getNextDir(r)

#print("l2",l2)
tg=[]
for l in l2:
    tt="/".join(l.split("/")[0:8])
    if not tt in tg:
        tg.append(tt)

#print(tg)

l1=l2;l2=[];fs={};d=''
for t in tg:
    ss=0
    for r in l1:
        if t in r:
            d,s=getInfo(r)
            ss=ss+int(s)
    if d!='':
        fs[t]=[d,ss]

sorted=sorted(fs.items(),key=lambda x:x[1][0])
#print(sorted)
xs=[]
for f in sorted:
    xs.append([f[0],f[1][0],f[1][1]])

total=0
for x in xs:
    total=total+x[2]

print('total',total)
if total<limit:
    exit()

for x in xs:
    total=total-x[2]
    print("remove ",x[0],str(x[2]),"kbyte")
    ret=subprocess.run(['rm','-rf',x[0]],capture_output=True,text=True).stdout
    if total<limit:
        print("total",total)
        exit()
    
#for f in fs:
#    print(f,fs[f])
