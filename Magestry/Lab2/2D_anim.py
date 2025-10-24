import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import imageio
from imageio.v2 import imread


names = ["Out_csvs_N_100_p_1", "Out_csvs_N_100_p_9"]

for filename in names:
	times = []
	sep = ','
	if f"Pics_{filename}" in os.listdir():
		shutil.rmtree(f"Pics_{filename}")
	os.mkdir(f"Pics_{filename}")
	os.mkdir(os.path.join(f"Pics_{filename}", "U_exp"))
	os.mkdir(os.path.join(f"Pics_{filename}", "V_exp"))
	os.mkdir(os.path.join(f"Pics_{filename}", "U_theory"))
	os.mkdir(os.path.join(f"Pics_{filename}", "V_theory"))

	for file in os.listdir(filename):
		if file[:5] == "Time=":
			# Получение времени из имени файла
			time = file.split("=")[1][:-4]
			times.append(time)

			with open(os.path.join(filename, file)) as f :
				# time = float(f.readline().split(":")[1])
				# time = float(f.readline())
				title = f.readline().split(sep)
				data = f.readlines()
			data = list(map(lambda x : list(map(float, x.split(sep))), data))
			try:
				data.sort(key=lambda x: (x[0], x[1]))
			except Exception as e:
				print(e)
				s = []
				for i in range(len(data)):
					s.append(len(data[i]))
					if len(data[i]) != 6:
						print(data[i], i)
				print(set(s))
				print(file)
			
			x = np.array(list(map(lambda x : x[0], data)))
			y = np.array(list(map(lambda x : x[1], data)))

			u_exp = np.array(list(map(lambda x : x[2], data)))
			v_exp = np.array(list(map(lambda x : x[3], data)))
			u_theory = np.array(list(map(lambda x : x[4], data)))
			v_theory = np.array(list(map(lambda x: x[5], data)))
			min_u_exp = min(u_exp)
			max_u_exp = max(u_exp)
			min_v_exp = min(v_exp)
			max_v_exp = max(v_exp)
			min_u_theory = min(u_theory)
			max_u_theory = max(u_theory)
			min_v_theory = min(v_theory)
			max_v_theory = max(v_theory)

			Lx = 1
			Ly = len(y)
			for i in range(1, len(x)) :
				if y[i] == y[0] :
					Ly = i
					break
			
			Lx = len(x) // Ly
			x = x.reshape(Ly, Lx)
			y = y.reshape(Ly, Lx)
			u_exp = u_exp.reshape(Ly, Lx)
			v_exp = v_exp.reshape(Ly, Lx)
			u_theory = u_theory.reshape(Ly, Lx)
			v_theory = v_theory.reshape(Ly, Lx)
			# print(x)
			# print(y)
	
			fig, ax = plt.subplots(1, 1)
			c = ax.pcolormesh(x, y, u_exp, cmap = 'jet')
			c.set_clim(min_u_exp, max_u_exp)

			fig.colorbar(c, ax = ax)
			ax.set_title(f't = {time} s', loc='right')
			plt.savefig(os.path.join(f"Pics_{filename}", "U_exp", f"Time={time}.png"))

			fig, ax = plt.subplots(1, 1)
			c = ax.pcolormesh(x, y, v_exp, cmap = 'jet')
			c.set_clim(min_v_exp, max_v_exp)

			fig.colorbar(c, ax = ax)
			ax.set_title(f't = {time} s', loc='right')
			plt.savefig(os.path.join(f"Pics_{filename}", "V_exp", f"Time={time}.png"))

			fig, ax = plt.subplots(1, 1)
			c = ax.pcolormesh(x, y, u_theory, cmap = 'jet')
			c.set_clim(min_u_theory, max_u_theory)

			fig.colorbar(c, ax = ax)
			ax.set_title(f't = {time} s', loc='right')
			plt.savefig(os.path.join(f"Pics_{filename}", "U_theory", f"Time={time}.png"))

			fig, ax = plt.subplots(1, 1)
			c = ax.pcolormesh(x, y, v_theory, cmap = 'jet')
			c.set_clim(min_v_theory, max_v_theory)
			fig.colorbar(c, ax = ax)
			ax.set_title(f't = {time} s', loc='right')
			plt.savefig(os.path.join(f"Pics_{filename}", "V_theory", f"Time={time}.png"))
			plt.close("all")



	
	species = ["U_exp", "V_exp", "U_theory", "V_theory"]
	for s in species:
		images = []
		times.sort(key=lambda x: float(x))
		for time in times:
			images.append(imread(os.path.join(f"Pics_{filename}", s, f"Time={time}.png")))
		imageio.mimsave(f"Res_{filename}_{s}.gif", images)
		print(f"Anim saved in Res_{filename}_{s}.gif.")

