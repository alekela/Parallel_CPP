import os
from time import sleep
import matplotlib.pyplot as plt

sizes = [1001, 100001]
Ns = list(range(1, 9))
times = {}

plt.figure()
legend = ["ideal"] + [f"real N={size}" for size in sizes]
plt.xlabel("Threads num")
plt.ylabel("T0/T")
plt.plot(range(1, max(Ns) + 1), range(1, max(Ns) + 1))
plt.title("Speed up")
plt.grid()

option = ''

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
	plt.plot(Ns, times[size])

with open(f"Results_{option}.txt", 'w') as f:
	for size in sizes:
		f.write(f"N={size}\n")
		f.write('\t'.join(map(str, Ns)) + '\n')
		f.write('\t'.join(map(str, times[size])) + '\n\n')
plt.legend(legend)
plt.xlim(0, max(Ns))
plt.ylim(0, max(Ns))
plt.savefig(f"out_{option}.png")
	