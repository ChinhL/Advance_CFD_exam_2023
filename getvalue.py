import pandas as pd
import math
import os
import numpy as np
def getData(filePath):
    if not os.path.isfile(filePath):
        raise IOError('No file with path %s' % filePath)
    header = []
    values = []
    with open(filePath, 'r') as f:
        i = 1;
        for line in f:
            if line[0] == '#':
                header = line.split()
            else:    
                split_line = [float(x) for x in line.split()]
                values.append(np.asarray(split_line))
                i = i+1
    return {'header': header,
            'values': np.asarray(values)}
def getpressure(A):
    pinlet=getData(A+'/postProcessing/swakExpression_p_inlet/0/p_inlet')
    poutlet=getData(A+'/postProcessing/swakExpression_p_outlet/0/p_outlet')
    outcome=pd.DataFrame(index=range(len(pinlet['values'])))
    outcome['time']=np.nan
    outcome['inlet']=np.nan
    outcome['outlet']=np.nan
    outcome['drop']=np.nan
    for i in range (0,len(outcome)):
        outcome['inlet'][i]=pinlet['values'][i][1]
        outcome['outlet'][i]=poutlet['values'][i][1]
        outcome['drop'][i]=pinlet['values'][i][1]-poutlet['values'][i][1] 
        outcome['time'][i]=pinlet['values'][i][0]
    return outcome
def bulkT(A):
    Tubulk=getData(A+'/postProcessing/swakExpression_Tbulk_up/0/Tbulk_up')
    ubulk=getData(A+'/postProcessing/swakExpression_velocity_bulk/0/velocity_bulk')
    outcome=pd.DataFrame(index=range(len(Tubulk['values'])))
    outcome['time']=np.nan
    outcome['Tuintergral']=np.nan
    outcome['uintergral']=np.nan
    outcome['Tbulk']=np.nan
    for i in range (0,len(outcome)):
        if i==0:
            outcome['Tuintergral'][i]=Tubulk['values'][i][1]
            outcome['uintergral'][i]=ubulk['values'][i][1]
            outcome['Tbulk'][i]=0
            outcome['time'][i]=Tubulk['values'][i][0]
        else:
            outcome['Tuintergral'][i]=Tubulk['values'][i][1]
            outcome['uintergral'][i]=ubulk['values'][i][1]
            outcome['Tbulk'][i]=Tubulk['values'][i][1]/ubulk['values'][i][1] 
            outcome['time'][i]=Tubulk['values'][i][0]
    return outcome


