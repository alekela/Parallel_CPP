import os
from time import sleep
import matplotlib.pyplot as plt

Ns = [128, 512, 1024]
ps = list(range(1, 13))
times = {}

plt.figure()
legend = ["ideal"] + [f"real N={N}" for N in Ns]
plt.xlabel("procs")
plt.ylabel("T0/T")
plt.plot(range(1, max(ps) + 1), range(1, max(ps) + 1))
plt.title("Speed up")

for N in Ns:
	times[N] = []
	for p in ps:
		with open("srun") as f:
			data = f.readlines()
		new_data = []
		for line in data:
			if "#SBATCH --tasks-per-node" in line:
				line = line.split("=")[0] + "=" + str(p) + '\n'
			if "#SBATCH --output" in line:
				line = line.split("=")[0] + "=" + f"out_p_{p}_N_{N}.out" + '\n'
			if "mpirun $NAME" in line:
				line = " ".join(line.split()[:2]) + " " + str(N)
			new_data.append(line)

		with open("srun", 'w') as f:
			for i in new_data:
				f.write(i)
		os.system("sbatch srun")

		flag = True
		sleep(5)
		while flag:
			with open(f"out_p_{p}_N_{N}.out") as f:
				data = f.readlines()
			for line in data:
				if "Time of working" in line:
					flag = False
			sleep(2)
		print(f"starting analyze out_p_{p}_N_{N}.out")
		with open(f"out_p_{p}_N_{N}.out") as f:
			for line in f:
				if "Time of working" in line:
					times[N].append(float(line.split(":")[1]))
	t0 = times[N][0]
	for i in range(len(times[N])):
		times[N][i] = t0 / times[N][i]
	plt.plot(ps, times[N])

with open("Results.txt", 'w') as f:
	for N in Ns:
		f.write(f"N={N}\n")
		f.write('\t'.join(map(str, ps)) + '\n')
		f.write('\t'.join(map(str, times[N])) + '\n\n')
plt.legend(legend)
plt.xlim(0, max(ps))
plt.ylim(0, max(ps))
plt.savefig("Speed_up_graph.png")
	