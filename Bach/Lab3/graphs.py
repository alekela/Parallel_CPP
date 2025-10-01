import os
from time import sleep
import matplotlib.pyplot as plt
import numpy as np


sizes = [1001, 100001]
Ns = list(range(1, 9))
times = {}

legend = ["ideal"] + [f"h={10 / (sizes[i] - 1)}" for i in range(len(sizes))]
fig = plt.figure()
ax = fig.add_subplot()
fig2 = plt.figure()
ax2 = fig2.add_subplot()

ax.set_xlabel("Threads num")
ax.set_ylabel("T0/T")
ax.plot(range(1, max(Ns) + 1), range(1, max(Ns) + 1))
ax.set_title("Speed up")
ax.grid()
ax.set_xlim(0, max(Ns))
ax.set_ylim(0, max(Ns))

ax2.set_xlabel("Threads num")
ax2.set_ylabel("Efficiency")
ax2.plot(range(1, max(Ns) + 1), [1 for _ in range(max(Ns))])
ax2.set_title("Efficiency")
ax2.grid()
ax2.set_xlim(0, max(Ns))
# ax2.set_ylim(0, max(Ns))


option = 'threads'

for size in sizes:
	times[size] = []
	for n in Ns:
		with open(f"res/out_threads_{n}_size_{size}.out") as f:
			for line in f:
				if "Time of working" in line:
					times[size].append(float(line.split(":")[1]))
	t0 = times[size][0]
	for i in range(len(times[size])):
		times[size][i] = t0 / times[size][i]
	ax.plot(Ns, times[size])
	ax2.plot(Ns, np.array(times[size]) / np.array(Ns))

with open(f"Results_{option}.txt", 'w') as f:
	for size in sizes:
		f.write(f"N={size}\n")
		f.write('\t'.join(map(str, Ns)) + '\n')
		f.write('\t'.join(map(str, times[size])) + '\n\n')

ax.legend(legend)
ax2.legend(legend)
fig.savefig(f"Speed_up_{option}.png")
fig2.savefig(f"Efficiency_{option}.png")

	