24a25,26
> import sqlite3
> import datetime
26d27
< _web=False
27a29,31
> _web=False
> 
> db="sList.db"
31,32d34
<     if not os.path.exists(dr):
<         os.makedirs(dr)
162a165,170
>     def drawPhrase(self,x,y,string,center=False,adj=True):
>         if self.type=='pdf':
>             drawString(self,x,y,string,center=center,adj=adj)
>         else:
>             shape=self.slide.shapes.add_textbox(xx,yy-Pt(self.font_size),pw,py)
> 
164d171
< #        komai
173d179
< #                komai
181a188
>             tf.word_wrap=True
250a258,275
> def getInfo(tg,lines):
>     key='Route:'+str(tg)
>     on=False
>     ret=[]
>     for l in lines:
>         if key in l:
>             on=True
>             continue
>         if on:
>             if 'Route' in l:
>                 if len(ret)<1:
> #                    print(l1)
>                     ret.append([l1.split(',')[0]])
>                 return ret
>             each=l.split(',')
>             ret.append(each[:-1])
>     return ret
> 
295a321,322
> #    print("in checkJoin",p1,s[p1][0],s[p1][1],"<---")
> #    print(c[2])
367a395
> #        print("proc")
382a411
> #                print("n",n,m)
390a420,422
>     #komai
>     if len(connect)<1:# case of no route: V2.2
>         return [],0
440a473
> #        print(n,"->",m)
462a496
> #            print("upper",n)
464a499
> #            print("lower",n)
500a536
> #            print(comp,"e"+str(n+1)+"x"+str(m+1)+".svg") #debuf
614,686d649
< def align(c,smiles,s,fc,tite):
<     connect=c
<     m=0
<     l_fc=len(fc)
<     annealed=[]
<     for i in c:
<         if m>i[4]:
<             m=i[4]
<     for i in c:
<         i[4]-=m;
< 
<     loop=True
<     depth=0# cleaned floor depth
<     while loop:
<         d=0
<         y0=100
<         y1=-y0
<         x0=s*2
<         loop=False
<         ccc=0
<         for cnt,i in enumerate(c):
<             xx0=rightEdge(i,fc,cnt,l_fc)
<             if xx0 > s:# arrow starting point!
<                 loop=True
<                 if x0 > i[4]:
<                     x0=i[4]
<                     ccc=cnt
<             elif i[5]<=depth:
<                 if y1 < i[5]:
<                     y1=i[5] # smaller 
<                 if y0 > i[5]:
<                     y0=i[5] # bigger
<         if loop:
<             for cnt,i in enumerate(c):
<                 xx0=rightEdge(i,fc,cnt,l_fc)
<                 if xx0 > s:
<                     i[4] -= x0
<                     i[5] += (y0-y1-1)
<         depth+=(y0-y1-1)
< 
<         sharpHead=getSharpHead(getGoal(smiles,c),connect)
<         annealHead(sharpHead,connect)
< 
<         if tite:
<             span=6
<         else:
<             span=4
<         deep=[0]
<         for n,vv in enumerate(fc):
<             for v in vv.values():
<                 i=-c[n][5]
<                 # ugly
<                 if i<0 or i>5:
<                     i=5
<                     c[n][5]=-i
<                 while i>=len(deep):
< #                if i>=len(deep):
<                     deep.append(0)
<                 if deep[i]<v:
<                     deep[i]=v
<         dx=[0]
<         for n in range(len(deep)-1):
<             dx.append((deep[n]+deep[n+1])/span+dx[n])
< 
<         for i in c:
<             try:
<                 if -i[5]<len(dx):
<                     i[5]=-int(dx[-i[5]])
<                 else:
<                     i[5]=i[4]
<             except:
<                     i[5]=i[4]
<         
837c800
< #    print("sim",sim)
---
>     print("sim",sim)
905c868
< #    print("length",ll)
---
>     print("length",ll)
971c934
< #        print(j,name[j])
---
>         print(j,name[j])
1105c1068
< #    print("len",len(_smiles),l,pt)
---
>     print("len",len(_smiles),l,pt)
1131a1095,1097
>             if len(_smiles[r])<1:
>                 sAll[r]=[]
>                 continue
1138,1140c1104,1106
< #    print("type",typ)
< #    print("the trick",pt,sAll[pt])
< #    print(sAll)
---
>     print("type",typ)
>     print("the trick",pt,sAll[pt])
>     print(sAll)
1302c1268
< #    print(branch)
---
>     print(branch)
1318,1322c1284,1288
< #    for i in ok:
< #        print("Route",i,end=":")
< #        for j in ok[i]:
< #            print(branch[j],end="->")
< #        print()
---
>     for i in ok:
>         print("Route",i,end=":")
>         for j in ok[i]:
>             print(branch[j],end="->")
>         print()
1324a1291,1313
> def dicToString(ss):
>     r=""
>     for key, value in ss.items():
>         r+=str(key)+":"
>         for k in value:
>             r+=str(k)+","
>         r+=";"
>     return r
> 
> def listToString(ss):
>     r=""
>     for s in ss:
>         ll=len(s)
>         for n,l in enumerate(s):
>             if type(l)==list:
>                 for i in l:
>                     r+=str(i)+","
>             else:
>                 r+=str(l)
>             r+=";"
>         r+="##"
>     return r
> 
1329,1342c1318,1322
< chemical="losartan"
< proc_per_line=10 # how many reactions per line roughly ?
< _fc=True
< _tite=False
< _type=3
< _id=False
< _product_only=False 
< _include_subsets=True
< 
< #print(sys.argv)
< ppt=False;chem=None;_routes=[]
< input_dir="./"
< com='-chem chemical_name, -ppt for power point -n [1,2,3] to specify routes -p proc_per_line -d input_path, -one_size : rendering at the same size, -product_only : categolize routes with respct to reaction product only, -show_subsets : show routes including subsets'
< com+=" -summary 0/1/2"
---
> input_dir="./work"
> script='async.sh'
> user='kisa'
> 
> com='-u user -d input_dir -s script.sh'
1344a1325
> 
1350,1368c1331,1332
<         case '-h':
<             print(sys.argv[0],ops['-h']);exit()
<         case '-ppt':
<             ppt=True
<         case '-heavy':
<             _type=2
<         case '-atom_number':
<             _type=3
<         case '-one_size':
<             _fc=False
<         case '-show_subsets':
<             _include_subsets=True
<         case '-n':
<             _routes=getList(sys.argv[n+1])
<             skip=True
<         case '-tite':
<             _tite=True
<         case '-sub_id':
<             _id=sys.argv[n+1]
---
>         case '-database':
>             db=sys.argv[n+1]
1373,1376c1337,1338
<         case '-product_only':
<             _product_only=True
<         case '-p':
<             proc_per_line=int(sys.argv[n+1])
---
>         case '-u':
>             user=sys.argv[n+1]
1378,1379c1340,1341
<         case '-chem':
<             chemical=sys.argv[n+1]
---
>         case '-s':
>             script=sys.argv[n+1]
1382,1398d1343
< input_file=chemical+'.txt'
< if _id:
<     output_file=chemical+str(_id)
< else:
<     output_file=chemical
< 
< if ppt:
<     output_file=output_file+'.pptx'
< else:
<     output_file=output_file+'.pdf'
< 
< if _web:
<     if not os.path.exists(input_dir+'/report'):
<         os.makedirs(input_dir+'/report')
<     output_file=input_dir+'/report/'+output_file
< print(input_file)
< print(output_file)
1400d1344
< print(proc_per_line)
1402c1346,1368
< head=None
---
> input_file="."
> with open(user+"_"+script,"r") as f:
>     for l in f:
>         if '#name' in l:
>             uname=l.split('\n')[0].split('name=')[1]
>         if '#email' in l:
>             email=l.split('\n')[0].split('email=')[1]
>         if 'python3' in l:
>             es=l.split('\n')[0].split(' ')
>             for n,i in enumerate(es):
>                 print(n,i)
> smi="".join(es[2].split("'"))
> loop=es[3]
> factor=es[4]
> substance="".join(es[13].split("'"))
> 
> opts=""
> for i in es[5:11]:
>     if i=='':
>         opts+='Non,'
>     else:
>         opts+=i+','
> #opts=es[5:11]
1404c1370
< print(input_dir+"/"+input_file)
---
> # output #5-10
1406c1372,1373
< with open(input_dir+"/"+input_file,"r") as f:
---
> input_file=substance+".txt"
> with open(input_dir+"/"+user+"/"+input_file,"r") as f:
1407a1375,1376
> with open(input_dir+"/"+user+"/"+input_file+"_info","r") as f:
>     all_info=f.readlines()
1412,1429c1381,1383
< tg="smiles"
< if _product_only:
<     tg="products"
< parent=routeMatrix[tg]
< 
< if _include_subsets:
<     show='total'
<     com="."
< else:
<     show='sub'
<     com=" and subsets eliminated."
< 
< print("parent",parent['total'])
< routes=list(parent[show].keys())
< 
< #easy="All output:"+str(len(allRoutes))+", Net output:"+str(len(parent['total']))+", Prime output:"+str(len(parent['sub']))
< 
< print(f"Categolized in respect to {tg}"+com)
---
> #print("=========================================")
> #print(routeMatrix,ppt)
> #exit()
1431,1432c1385,1386
< for r in parent[show]:
<     print(r,parent[show][r])
---
> tg="smiles"
> parent=routeMatrix["smiles"]
1434,1439c1388
< print("--->",parent[show][r])
< if len(parent[show][r])<1:
<     hint=[False,"No search results"]
< else:
<     hint=deepArrange(parent[show])
<     hint.append(ppt)
---
> routes=list(parent['sub'].keys())
1441,1463c1390,1391
< span=proc_per_line*2
< if span<10:
<     span=10
< scale=_magic_x/(span)
< 
< _py=500
< py=_py
< ox=40
< 
< if len(_routes)>0:
<     routes=_routes
< 
< page=openReport(output_file,scale)
< 
< #print(head)
< sim=similarityOnRoute(hint[2],hint[0],hint[1],all) 
< easy=str(len(parent['total']))+" variations from "+str(len(sim))+" routes over "+str(len(allRoutes))+" queries"
< hint.append(easy)
< if parent!=False:
<     head=[chemical,getSmiles(routes[0],all),hint]
< else:
<     head=[chemical,getSmiles(1,all),hint]
< head_page(head,page,sim)
---
> print("parent",parent,"<---------")
> print("routes",routes)
1465,1467c1393,1394
< for route in routes:
<     print(f"Route:{route}",end="..")
<     oy=py-scale
---
> conn=sqlite3.connect(db)
> cur=conn.cursor()
1469c1396,1398
<     smiles=getSmiles(route,all)
---
> cur.execute('create table if not exists searchList(id integer primary key, user text, uname text, email text, smiles text, cSmiles text, date text, substance text, loop integer, factors text, options text);')
> cur.execute('create table if not exists parent(id integer, total text, sub text);')
> conn.commit()
1471c1400,1403
<     fc=makeSgv(smiles)
---
> #smiles=getSmiles(routes[0],all)[-1][-1]
> mol=Chem.MolFromSmiles(smi)
> cSmiles=Chem.MolToSmiles(mol)
> dd=datetime.datetime.now().strftime('%Y%m%d%H%M%S')
1473,1474c1405,1407
<     if not _fc:
<         fc=_fc
---
> print(opts)
> cur.execute(f'insert into searchList(user, uname, email, smiles, cSmiles, date, substance, loop, factors, options) values("{user}","{uname}","{email}","{smi}","{cSmiles}",{dd},"{substance}",{loop},"{factor}","{opts}")')
> conn.commit()
1476,1481c1409,1411
<     connect,start=initConnect(smiles,fc=fc,mode=_type)
<     Debug(connect,head='G',sp=5)
< #
<     align(connect,smiles,span-3,fc,_tite)
<     Debug(connect,head='H',sp=5)
<     hy=evalHeight(connect,fc)
---
> sql=f'select id from "searchList" where date = "{dd}";'
> ret=cur.execute(sql)
> conn.commit()
1483,1484c1413,1415
<     if fc:
<         setWds(connect,fc)
---
> total=dicToString(parent['total'])
> sub=dicToString(parent['sub'])
> id=ret.fetchall()[0][0]
1486,1490c1417,1422
<     if oy+hy*scale*2<0:
<         page.stroke()
<         py=_py
<         oy=py-scale
<     page.setFont()
---
> print("total",total)
> print(" sub" ,sub)
> sql=f'insert into parent(id, total, sub) values({id},"{total}","{sub}")'
> print(sql)
> cur.execute(sql)
> conn.commit()
1492,1495c1424
<     tmp_goal=[]
<     for i in range(len(connect)):
<         x0=connect[i][4]*scale+ox
<         y0=connect[i][5]*scale*2+oy
---
> sTable=f"searchTable{id}"
1497,1517c1426,1427
<         if fc:
<             wds=fc[i]['s']
<         if len(connect[i][2])<1:
<             fn=getName(i,connect[i][6]+1)# need
< #            print("starting",i,fn,x0,y0)
<             if not fc:
<                 page.drawImage(fn,x0,y0,scale) # for starting material
<             else:
< #                wds=fc[i][connect[i][6]]
<                 page.drawImage(fn,x0,y0,scale*wds) # for starting material
<         n=connect[i][3]
< #    if n<0:
< #        continue
<         if n<0:# for goal
<             if not fc:
<                 xx= x0+(len(connect[i][7])+2)*scale
<             else:
<                 if i<len(fc):
<                     sp=len(connect[i][7])+2
<                 else:
<                     sp=len(connect[i][7])+2
---
> cur.execute(f'create table if not exists "{sTable}" (route int, connect_list text, smiles_list text, info_list text);')
> conn.commit()
1519,1537c1429,1446
< #                    xx=x0+(fc[-1]['cheat']+fc[i][connect[i][6]])*scale
<                 xx=x0+(wds+sp)*scale
<             yy= y0
<             if len(tmp_goal)<1 and n==-1:
<                 tmp_goal=[xx,yy]
<         else:
<             xx=connect[n][4]*scale+ox
<             yy=connect[n][5]*scale*2+oy
< 
<         if connect[i][1][0]>-1:
< #        if connect[i][3]>-1:
<             fn=getName(i,connect[i][1][0]+1)
< #            print("result",i,fn,xx,yy)
<             if not fc:
<                 page.drawImage(fn,xx,yy,scale) # for resultant material
<                 xp=x0+scale*1.5
<             else:
<                 page.drawImage(fn,xx,yy,scale*fc[i][connect[i][1][0]]) # for resultant material
<                 xp=x0+scale*(0.5+wds)
---
> #komai
> for route in routes:
>     info=getInfo(route,all_info)
>     smiles=getSmiles(route,all)
>     connect,start=initConnect(smiles)
>     s1=listToString(connect)
>     s2=listToString(smiles)
> #    print("---->",route)
> #    print("connect",connect)
> #    print("connect(s1)",s1)
> #    print(smiles)
> #    print("smiles(s2)",s2)
> #    print(info)
>     s3=listToString(info)
>     sql=f'insert into {sTable} values({route},"{s1}","{s2}","{s3}");'
> #    print("sql",sql)
>     cur.execute(sql)
>     conn.commit()
1539,1574d1447
<         for j in connect[i][-1]:
<             fn=getName(i,j+1)
< #            page.drawImage(fn,xp,y0+step*0.55,scale) # for catalytic material
<             if not fc:
<                 page.drawImage(fn,xp,y0,scale,bt=True) # for catalytic material
<                 xp=xp+scale
<             else:
< #                page.drawImage(fn,xp,y0,scale*fc[i][j],bt=True) # for catalytic material
< #                fx=fc[i][j] if fc[i][j]<2 else 1
<                 fx=1
<                 page.drawImage(fn,xp,y0,scale*fx,bt=True) # for catalytic material
<                 xp=xp+fc[i][j]*scale
< #        page.drawArrow(arrow(x0,y0,xx,yy,scale)) 
<         if n<0 and len(tmp_goal)>0:# for goal
<             xx=tmp_goal[0]
<             yy=tmp_goal[1]
< ## komai
< #
<         if not fc:
<             page.drawArrow(arrow(x0,y0,xx,yy,scale)) 
<         else:
<             skip=wds
< #            print("arrow:",n,wds,x0,xx,xx-x0,skip)
<             page.drawArrow(arrow(x0,y0,xx,yy,scale,skip=skip,fc=fc))
< 
<         if fc:
<             py=oy+getTheBottom(connect,fc)*scale
<         else:
<             if py>y0-scale:
<                 py=y0-scale
< 
<     page.drawString(ox,oy+scale*1.5,"Route"+str(route))
< 
< #for i in range(10):
< #    page.drawLine([0,i,10,i])
< page.close()
1575a1449,1451
> 
> print(f"substance:{substance}###")
> conn.close()
