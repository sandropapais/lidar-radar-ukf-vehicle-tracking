# import csv
#
# with open('fusion/est_out.csv', newline='') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
#     for row in spamreader:
#         print(', '.join(row))

# # first get all lines from file
# with open('radar/nis_out.csv', 'r') as f:
#     lines = f.readlines()
#
# # remove spaces
# lines = [line.replace(' ', '') for line in lines]
#
# # finally, write lines in the file
# with open('radar/nis_out.csv', 'w') as f:
#     f.writelines(lines)

# import matplotlib.pyplot as plt
# import matplotlib.cbook as cbook
#
# import numpy as np
# import pandas as pd
#
# msft = pd.read_csv('fusion/est_out.csv')
# msft.plot("x_est")
#
# fname2 = cbook.get_sample_data('fusion/est_out.csv', asfileobj=False)
# with cbook.get_sample_data('fusion/est_out.csv') as file:
#     array = np.loadtxt(file)
#
# fig, axs = plt.subplots(2, sharex=True)
# axs[0].plot(array[:, 0], array[:, 1])
# axs[0].set(ylabel='$f(x)=x^2$')
# axs[1].plot(array[:, 0], array[:, 2])
# axs[1].set(xlabel='$x$', ylabel='$f(x)=x^3$')

import math
import matplotlib.pyplot as plt
import csv
import numpy as np

t = []
x_est = []
y_est = []
y_est = []
v_est = []
yaw_est = []
yawdot_est = []
x_var = []
y_var = []
v_var = []
yaw_var = []
yawdot_var = []
x_tru = []
y_tru = []
vx_tru = []
vy_tru = []

with open('fusion/est_out.csv', 'r') as csvfile:
    out_reader = csv.reader(csvfile, delimiter=',')

    for row in out_reader:
        if row[0] == '0':
            t.append(float(row[1]))
            x_est.append(float(row[2]))
            y_est.append(float(row[3]))
            v_est.append(float(row[4]))
            yaw_est.append(float(row[5]))
            yawdot_est.append(float(row[6]))
            x_var.append(float(row[7]))
            y_var.append(float(row[8]))
            v_var.append(float(row[9]))
            yaw_var.append(float(row[10]))
            yawdot_var.append(float(row[11]))
            x_tru.append(float(row[12]))
            y_tru.append(float(row[13]))
            vx_tru.append(float(row[14]))
            vy_tru.append(float(row[15]))

x_est = np.array(x_est)
y_est = np.array(y_est)
v_est = np.array(v_est)
yaw_est = np.array(yaw_est)
yawdot_est = np.array(yawdot_est)
x_var = np.array(x_var)
y_var = np.array(y_var)
v_var = np.array(v_var)
yaw_var = np.array(yaw_var)
yawdot_var = np.array(yawdot_var)
x_tru = np.array(x_tru)
y_tru = np.array(y_tru)
vx_tru = np.array(vx_tru)
vy_tru = np.array(vy_tru)
v_tru = np.sqrt(np.multiply(vx_tru, vx_tru) + np.multiply(vy_tru, vy_tru))

fig, axs = plt.subplots(3, 2)
axs[0, 0].plot(t, x_est-x_tru)
axs[0, 0].plot(t, x_est-x_tru+3*np.sqrt(x_var), '--r')
axs[0, 0].plot(t, x_est-x_tru-3*np.sqrt(x_var), '--r')
axs[0, 0].set_ylabel('X Position (m)')
axs[0, 1].plot(t, y_est-y_tru)
axs[0, 1].plot(t, y_est-y_tru+3*np.sqrt(y_var), '--r')
axs[0, 1].plot(t, y_est-y_tru-3*np.sqrt(y_var), '--r')
axs[0, 1].set_ylabel('Y Position (m)')
axs[1, 0].plot(t, v_tru)
axs[1, 0].plot(t, v_est, '--r')
#axs[1, 0].plot(t, v_est-v_tru-3*np.sqrt(v_var), '--r')
axs[1, 0].set_ylabel('Speed (m/s)')
axs[2, 0].set_xlabel('Times (s)')
axs[2, 1].set_xlabel('Times (s)')
fig.suptitle(r'State Estimation Error and $3\sigma$')
plt.show()

# plt.plot(t, x_est, label="x")
# plt.plot(t, y_est, label="y")
# plt.xlabel('Times (s)')
# plt.ylabel('Estimated Position (m)')
# plt.title('Vehicle Estimated States')
# plt.legend()
# plt.show()