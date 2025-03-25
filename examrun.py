import getvalue as gv
import matplotlib.pyplot as plt
pressuredrop=gv.getpressure('/home/chinhle/OpenFOAM/chinhle-8/run/examtask')
pressuredrop.to_csv('/home/chinhle/OpenFOAM/chinhle-8/run/pressuredrop.csv',index=False)
bulktemperature=gv.bulkT('/home/chinhle/OpenFOAM/chinhle-8/run/examtask')
bulktemperature.to_csv('/home/chinhle/OpenFOAM/chinhle-8/run/bulktemperature.csv',index=False)
plt.plot(bulktemperature['time'],bulktemperature['Tbulk'])
plt.xlabel('time')
plt.ylabel('bulk temperature')
plt.show()