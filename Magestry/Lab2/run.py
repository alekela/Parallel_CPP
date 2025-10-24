import os
from time import sleep
import matplotlib.pyplot as plt

Ns = [128, 512, 1024]
ps = list(range(1, 13))
times = {}

flag_write = input("Rewrite results? (y/n): ")
while flag_write not in ['y', 'n']:
	flag_write = input("Rewrite results? (y/n): ")

if flag_write == 'y':
	with open("Results.csv", 'w') as f:
		f.write("N,p,time\n")
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
					line = line.split("=")[0] + "=" + f"Outs/out_p_{p}_N_{N}.out" + '\n'
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
				with open(f"Outs/out_p_{p}_N_{N}.out") as f:
					data = f.readlines()
				for line in data:
					if "Time of working" in line:
						flag = False
				sleep(2)

			with open(f"Outs/out_p_{p}_N_{N}.out") as f:
				for line in f:
					if "Time of working" in line:
						times[N].append(float(line.split(":")[1]))
		with open("Results.csv", 'a') as f:
			for i in range(len(ps)):
				f.write(f"{N},{ps[i]},{times[N][i]}\n")

with open("Results.csv") as f:
	title = f.readline()
	data = f.readlines()
data = list(map(lambda x: x.strip().split(','), data))
ddata = {f"{line[0]},{line[1]}": float(line[2]) for line in data}
times = {}
for N in Ns:
	times[N] = []
	for p in ps:
		times[N].append(ddata[f"{N},{p}"])

plt.figure()
legend = ["ideal"] + [f"real N={N}" for N in Ns]
plt.xlabel("procs")
plt.ylabel("T0/T")
plt.plot(range(1, max(ps) + 1), range(1, max(ps) + 1))
plt.title("Speed up")

for N in Ns:
	t0 = times[N][0]
	for i in range(len(times[N])):
		times[N][i] = t0 / times[N][i]
	plt.plot(ps, times[N])

plt.legend(legend)
plt.grid()
plt.xlim(0, max(ps))
plt.ylim(0, max(ps))
plt.savefig("Speed_up_graph.png")

plt.figure()
legend = ["ideal"] + [f"real N={N}" for N in Ns]
plt.xlabel("procs")
plt.ylabel("Efficiency")
plt.plot(range(1, max(ps) + 1), [1 for _ in range(max(ps))])
plt.title("Efficiency")

for N in Ns:
	for i in range(len(times[N])):
		times[N][i] = times[N][i] / ps[i]
	plt.plot(ps, times[N])

plt.legend(legend)
plt.grid()
plt.xlim(0, max(ps))
plt.ylim(0, 1.2)
plt.savefig("Efficiency_graph.png")
	