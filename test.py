#!/usr/bin/env python
import pandas as pd
import sqlalchemy.engine
from lsst_transients.utils import database_io



regfile = '/home/jrblaum/lsst_transients/test/supernova_test.reg'

with open(regfile) as f:
    i=0
    strs=[]
    for line in f:
        if i==0:
            pass
        else:
           trimmed =line[0:len(line)-1]
           strs.append(trimmed)
        i+=1
    indices =range(1,i) 
    series = pd.Series(strs, index=indices)
    reg_dataframe = pd.DataFrame.from_dict({'ds9_info': series})

dbname = 'test1.db'
engine = sqlalchemy.engine.create_engine('sqlite:///%s' % dbname)
connection = engine.connect()
reg_dataframe.to_sql("reg_dataframe", connection, if_exists='replace', index_label='regID')
connection.close()

dbname = 'test2.db'
db = database_io.SqliteDatabase(dbname)
db.connect()

db.insert_dataframe(reg_dataframe, 'reg_dataframe')

for r in range(5):
    flux_dataframe = pd.DataFrame(columns=['flux', 'err'], dtype=float)
    db.insert_dataframe(flux_dataframe, 'flux_table_%i' % (r + 1))
    
cond_dataframe = pd.DataFrame(columns=['date (modified Julian)','duration (s)','seeing (")'], dtype=float)
db.insert_dataframe(cond_dataframe, 'cond_table') 

db.disconnect()
