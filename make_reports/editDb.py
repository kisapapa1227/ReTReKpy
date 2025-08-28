import sqlite3
import sys


def myForm(x,startWith=1):
    ret='"'+str(x[startWith])+'"'
    for n in range(startWith+1,len(x)):
        ret=ret+',"'+str(x[n])+'"'
    return ret

def getMaxPid(cur):
    pid_max=0
    ret=cur.execute(f"""select * from searchList""")

    for d in ret.fetchall():
#        print("pid",d[0])
        if pid_max<d[0]:
            pid_max=d[0]
    return pid_max

def rmRecords(src,ids):
    conn=sqlite3.connect(src);cur=conn.cursor()
    for p in ids:
        sTable=f"searchTable{p}"
        print(sTable)
        cur.execute(f"""drop table {sTable}""")
        cur.execute(f"""delete from searchList where id={p}""")
        cur.execute(f"""delete from parent where id={p}""")
    conn.commit()
    conn.close()

def cpRecords(src,dst,ids):
    conn=sqlite3.connect(src);cur=conn.cursor()
    oconn=sqlite3.connect(dst);ocur=oconn.cursor()

    ocur.execute('create table if not exists searchList(id integer primary key, user text, uname text, email text, smiles text, cSmiles text, date text, substance text, loop integer, factors text, options text);')
    ocur.execute('create table if not exists parent(id integer, total text, sub text);')
    oconn.commit()

    pid=getMaxPid(ocur)+1

    subs={}

    for p in ids:
        ret=cur.execute(f"""select * from searchList where id={p}""")
        sTable=f"searchTable{p}"

        d=ret.fetchall()[0]
#
# about searchTable
#
        subs[p]=d[7]
        fields='"'+str(pid)+'",'+myForm(d)
        osTable=f"searchTable{pid}"
        ocur.execute(f'create table if not exists "{osTable}" (route int, connect_list text, smiles_list text, info_list text);')

        sql=f'insert into searchList(id, user, uname, email, smiles, cSmiles, date, substance, loop, factors, options) values({fields})'
        ocur.execute(sql)
#
# about parents
#
        ret=cur.execute(f"""select * from parent where id={p}""")
        d=ret.fetchall()[0];total=d[1];sub=d[2]
        sql=f'insert into parent(id, total, sub) values({pid},"{total}","{sub}")'
        ocur.execute(sql)
#
# about sTable 
#
        ret=cur.execute(f"""select * from {sTable}""")
        for d in ret.fetchall():
            fields=myForm(d,startWith=0)
            sql=f'insert into {osTable} values({fields});'
            ocur.execute(sql)
        oconn.commit()
        pid=pid+1
    oconn.close()
    return subs

########################################################
################# main #################################
########################################################
operation="cp";dst=False;src=False;hit_id=False;skip=False
ids=[]

for n,op in enumerate(sys.argv):
    if hit_id:
        if op[0]!='-':
            ids.append(int(op))
        else:
            hit_id=False
    if skip:
        skip=False;continue
    match op:
        case '-rm':
            operation="rm"
        case '-mv':
            operation="mv"
        case '-src':
            src=sys.argv[n+1]
            skip=True
        case '-dst':
            dst=sys.argv[n+1]
            skip=True
        case '-ids':
            hit_id=True

print("operation",operation)
print("ids",ids)
print("dst",dst)
print("src",src)

if operation=="cp" or operation=="mv":
    subs=cpRecords(src,dst,ids)
    for key in subs:
        print("substance"+str(key)+":"+subs[key])

if operation=="mv" or operation=="rm":
    rmRecords(src,ids)
