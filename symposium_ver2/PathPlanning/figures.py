# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import copy
import math
import sys
import os
import time
import csv
import pprint
from matplotlib.font_manager import FontProperties

font = {'family' : 'DejaVu Sans'}
matplotlib.rc('font', **font)

def main():
    wx = [0.0, 70.0]
    wy = [0.0, 0.0]
    ob = []
    states = []

    with open('csv/obstacle.csv') as f:
        reader = csv.reader(f)
        for row in reader:
            ob.append(row)
    ob = [[float(v) for v in row] for row in ob]

    with open('csv/experiment4.csv') as f:
        reader = csv.reader(f)
        for row in reader:
            states.append(row)
    states = [[float(v) for v in row] for row in states]

    plt.cla()
    # for stopping simulation with the esc key.
    plt.gcf().canvas.mpl_connect(
        'key_release_event',
        lambda event: [exit(0) if event.key == 'escape' else None])
    plt.plot(wx, wy, linewidth="2", label=u"目標パス", linestyle="dashed")
    plt.plot(ob[0][:], ob[1][:], "xk", label=u"障害物")
    plt.plot(states[0][:], states[1][:], linewidth="2", label=u"走行経路")
    fp = FontProperties(fname='/usr/share/fonts/truetype/takao-gothic/TakaoGothic.ttf', size=20)
    plt.xlabel("x[m]", fontsize=13)
    plt.xticks(fontsize=13)
    plt.ylabel("y[m]", fontsize=13)
    plt.yticks(fontsize=13)
    plt.legend(prop=fp, bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0)
    plt.xlim(-1, 71)
    plt.ylim(-20, 20)
    plt.grid(True)
    plt.pause(0.0001)
    # plt.show()
    plt.savefig("figures/experiment4.pdf")


if __name__ == '__main__':
    main()