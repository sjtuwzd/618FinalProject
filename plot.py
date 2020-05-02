import numpy as np
import matplotlib.pyplot as plt
import sys

with open("data/results.txt") as results_file:
    results_data = results_file.read()

results_data = results_data.split('\n')
del results_data[-1]

results_data = [float(x) for x in results_data]

blue_hex='#042247'
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Portfolio Return (%)', color=blue_hex)
ax.set_ylabel('Frequency', color=blue_hex)
#plt.axis((2,50,0,1200));

title = "No Title Provided"
if (len(sys.argv) > 1):
    title = sys.argv[1];
fig.suptitle(title, fontsize=20, fontweight='bold', color=blue_hex)

y,binEdges,_=plt.hist(results_data, bins=30, histtype='stepfilled', color=blue_hex);
bincenters = 0.5*(binEdges[1:]+binEdges[:-1]);
plt.plot(bincenters,y,'-', color='#885E23', linewidth=4);

plt.show()
