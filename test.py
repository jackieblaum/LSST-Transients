#!/usr/bin/env python
import pandas as pd
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
    
dbname = 'test.db'
db = database_io.SqliteDatabase(dbname)
db.connect()

db.insert_dataframe(reg_dataframe, 'reg_dataframe')

db.disconnect()
