import os
from time import sleep
import matplotlib.pyplot as plt

sizes = [2000, 10000, 50000]
Ns = list(range(1, 9))
times = {}

plt.figure()
legend = ["ideal"] + [f"real N={size}" for size in sizes]
plt.xlabel("N_proc")
plt.ylabel("T/T0")
plt.plot(range(1, max(Ns) + 1), range(1, max(Ns) + 1))
plt.title("Acceleration")

for size in sizes:
	times[size] = []
	for n in Ns:
		with open("run") as f:
			data = f.readlines()
		new_data = []
		for line in data:
			'''if "#SBATCH --tasks-per-node" in line:
				line = line.split("=")[0] + "=" + str(n) + '\n'
			if "#SBATCH --output" in line:
				line = line.split("=")[0] + "=" + f"out_n_{n}_size_{size}.out" + '\n'
			'''
			if "mpirun" in line:
				line = line.split()
				line[2] = str(n)
				line[4] = str(size)
				line[6] = f"out_n_{n}_size_{size}.out"
				line = " ".join(line)
			new_data.append(line)

		with open("run", 'w') as f:
			for i in new_data:
				f.write(i)
		os.system("bash run")
		sleep(15 * int(size / sizes[0]))
		print(f"out_n_{n}_size_{size}.out done")
		with open(f"out_n_{n}_size_{size}.out") as f:
			for line in f:
				if "Time of working" in line:
					times[size].append(float(line.split(":")[1]))
	t0 = times[size][0]
	for i in range(len(times[size])):
		times[size][i] = t0 / times[size][i]
	plt.plot(Ns, times[size])

with open("Results_1.txt", 'w') as f:
	for size in sizes:
		f.write(f"N={size}\n")
		f.write('\t'.join(map(str, Ns)) + '\n')
		f.write('\t'.join(map(str, times[size])) + '\n\n')
#plt.legend(legend)
plt.xlim(0, max(Ns))
plt.ylim(0, max(Ns))
plt.savefig("out_1.png")
	