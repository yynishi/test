import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anm
import csv

fig = plt.figure()

ob = []
with open('csv/obstacle.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        ob.append(row)
ob = [[float(v) for v in row] for row in ob]

target_course = []
with open('csv/target_course.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        target_course.append(row)
target_course = [[float(v) for v in row] for row in target_course]

low_speed_car = []
with open('csv/low_speed_car.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        low_speed_car.append(row)
low_speed_car = [[float(v) for v in row] for row in low_speed_car]

path_x = []
path_y = []
flag = 0
with open('csv/frenet_path.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        row = [s for s in row if s != '']
        if flag == 0:
            path_x.append(row)
            flag = 1
        else:
            path_y.append(row)
            flag = 0
path_x = [[float(v) for v in row] for row in path_x]
path_y = [[float(v) for v in row] for row in path_y]

car = []
with open('csv/car.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        car.append(row)
car = [[float(v) for v in row] for row in car]

target_fp = []
with open('csv/target_fp.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        target_fp.append(row)
target_fp = [[float(v) for v in row] for row in target_fp]


fig = plt.figure(figsize = (10, 10))

area = 4.0
def update(i, ):
    if i != 0:
        plt.cla()                      # 現在描写されているグラフを消去

    plt.plot(ob[:][0], ob[:][1], "xk", markersize = "8")
    plt.plot(target_course[0][:], target_course[1][:])
    plt.plot(path_x[int(i / 8)], path_y[int(i / 8)], "-or")
    plt.plot(target_fp[0][i], target_fp[1][i], "v", color = "b")
    plt.plot(low_speed_car[0][i], low_speed_car[1][i], "s", markersize = "10")
    plt.plot(car[0][:i], car[1][:i], linewidth = "5")
    plt.xlim(car[0][i] - area, car[0][i] + area)
    plt.ylim(car[1][i] - area, car[1][i] + area)
    plt.plot()


ani = anm.FuncAnimation(fig, update, interval = 50, frames = len(low_speed_car[0]))

# ani.save('anim.mp4', writer="ffmpeg")
plt.show()