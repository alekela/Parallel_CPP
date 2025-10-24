import numpy as np
import matplotlib.pyplot as plt


with open("Out_csvs/Time=3.000000.csv") as f :
	title = f.readline().split(',')
	data = f.readlines()
data = list(map(lambda x : list(map(float, x.split(','))), data))
data.sort(key=lambda x: (x[0], x[1]))

data_x_0_25_1 = np.array(list(filter(lambda x: x[0] == 0.245, data)))
data_x_0_25_2 = np.array(list(filter(lambda x: x[0] == 0.255, data)))

data_y_0_25_1 = np.array(list(filter(lambda x: x[1] == 0.245, data)))
data_y_0_25_2 = np.array(list(filter(lambda x: x[1] == 0.255, data)))

data_x_0_75_1 = np.array(list(filter(lambda x: x[0] == 0.745, data)))
data_x_0_75_2 = np.array(list(filter(lambda x: x[0] == 0.755, data)))

data_y_0_75_1 = np.array(list(filter(lambda x: x[1] == 0.745, data)))
data_y_0_75_2 = np.array(list(filter(lambda x: x[1] == 0.755, data)))

data_x_0_25 = (data_x_0_25_1 + data_x_0_25_2) / 2
data_y_0_25 = (data_y_0_25_1 + data_y_0_25_2) / 2
data_x_0_75 = (data_x_0_75_1 + data_x_0_75_2) / 2
data_y_0_75 = (data_y_0_75_1 + data_y_0_75_2) / 2

figs = [plt.figure() for _ in range(4)]
axes = [figs[i].add_subplot() for i in range(4)]

axes[0].set_xlabel("y")
axes[0].plot(data_x_0_25[:, 1], data_x_0_25[:, 2])
axes[0].plot(data_x_0_25[:, 1], data_x_0_25[:, 4])
axes[0].plot(data_x_0_25[:, 1], data_x_0_25[:, 3])
axes[0].plot(data_x_0_25[:, 1], data_x_0_25[:, 5])
axes[0].legend(["U_exp", "U_theory", "V_exp", "V_theory"])

axes[1].set_xlabel("y")
axes[1].plot(data_x_0_75[:, 1], data_x_0_75[:, 2])
axes[1].plot(data_x_0_75[:, 1], data_x_0_75[:, 4])
axes[1].plot(data_x_0_75[:, 1], data_x_0_75[:, 3])
axes[1].plot(data_x_0_75[:, 1], data_x_0_75[:, 5])
axes[1].legend(["U_exp", "U_theory", "V_exp", "V_theory"])

axes[2].set_xlabel("x")
axes[2].plot(data_y_0_25[:, 0], data_y_0_25[:, 2])
axes[2].plot(data_y_0_25[:, 0], data_y_0_25[:, 4])
axes[2].plot(data_y_0_25[:, 0], data_y_0_25[:, 3])
axes[2].plot(data_y_0_25[:, 0], data_y_0_25[:, 5])
axes[2].legend(["U_exp", "U_theory", "V_exp", "V_theory"])

axes[3].set_xlabel("x")
axes[3].plot(data_y_0_75[:, 0], data_y_0_75[:, 2])
axes[3].plot(data_y_0_75[:, 0], data_y_0_75[:, 4])
axes[3].plot(data_y_0_75[:, 0], data_y_0_75[:, 3])
axes[3].plot(data_y_0_75[:, 0], data_y_0_75[:, 5])
axes[3].legend(["U_exp", "U_theory", "V_exp", "V_theory"])

figs[0].savefig("x=0_25.png")
figs[1].savefig("x=0_75.png")
figs[2].savefig("y=0_25.png")
figs[3].savefig("y=0_75.png")

