import os
from time import sleep
import matplotlib.pyplot as plt

Ns = [128, 512, 1024]
cases = [(i, 1) for i in range(1, 13)] + [(1, i) for i in range(1, 13)] + [(2, 6), (3, 4), (4, 3), (6, 2)]
cases = [(i, 1) for i in range(1, 13)]
kind = 'procs' # 'threads'

flag_write = input("Rewrite results? (y/n): ")
while flag_write not in ['y', 'n']:
	flag_write = input("Rewrite results? (y/n): ")

if flag_write == 'y':
	with open("Results.csv", 'w') as f:
		f.write("N,p,threads,time\n")
	for N in Ns:
		for case in cases:
			p = case[0]
			thread = case[1]

			with open("run_test") as f:
				data = f.readlines()
			new_data = []
			for line in data:
				if "#SBATCH --ntasks-per-node" in line:
					line = line.split("=")[0] + "=" + str(p) + '\n'
				if "#SBATCH --output" in line:
					line = line.split("=")[0] + "=" + f"Outs/out_p_{p}_t_{thread}_N_{N}.out" + '\n'
				if "cpus-per-task" in line:
					line = line.split("=")[0] + "=" + str(thread) + '\n'
				if "srun --mpi pmix" in line:
					line = line.strip().split()
					line[-1] = str(N)
					line = " ".join(line) + '\n'
				new_data.append(line)
					
			with open("run_test", 'w') as f:
				for i in new_data:
					f.write(i)
			os.system("sbatch run_test")

			sleep(5)
			while True:
				with open(f"Outs/out_p_{p}_t_{thread}_N_{N}.out") as f:
					data = f.readlines()
				data = "\n".join(data)
				if "Time of working" in data:
					time_var = float(data.split(":")[1].split("\n")[0])
					break
				sleep(2)

			with open("Results.csv", 'a') as f:
				f.write(f"{N},{p},{thread},{time_var}\n")


with open("Results.csv") as f:
	title = f.readline()
	data = f.readlines()
data = list(map(lambda x: x.strip().split(','), data))
times = {}
for N in Ns:
	times[N] = []
	for case in cases:
		for line in data:
			if line[:3] == [str(N), str(case[0]), str(case[1])]:
				times[N].append(float(line[3]))

print(times)
index = 0 if kind == 'procs' else 1
ps = list(map(lambda x: x[index], cases))
plt.figure()
legend = ["ideal"] + [f"real N={N}" for N in Ns]
plt.xlabel(kind)
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
plt.savefig(f"Speed_up_graph_{kind}.png")

plt.figure()
legend = ["ideal"] + [f"real N={N}" for N in Ns]
plt.xlabel(kind)
plt.ylabel("Efficiency")
plt.plot(range(1, max(ps) + 1), [1 for _ in range(max(ps))])
plt.title("Efficiency")

for N in Ns:
	for i in range(len(times[N])):
		times[N][i] = times[N][i] / ps[i]
	plt.plot(ps, times[N])

plt.legend(legend)
plt.grid()
plt.xlim(0, 12)
plt.ylim(0, 1.2)
plt.savefig(f"Efficiency_graph_{kind}.png")
	