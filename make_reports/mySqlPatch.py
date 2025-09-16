import mysql.connector

conn=mysql.connector.connect(
        host="localhost",
        user="root",
        password="password"
        )

cursor=conn.cursor()

sql="select host,user from mysql.user;"
cursor.execute(sql)

cursor.close()
conn.close()
