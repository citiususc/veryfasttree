#!/usr/bin/env python3

import os.path

import numpy
import numpy as np
import sys
try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Module matplotlib required. Please install with pip install matplotlib", file=sys.stderr)
    exit(-1)

def parse(path, step):
    table = list()
    units = {"m": 1, "g": 2, "t": 3}
    with open(path) as file:
        for line in file:
            fields = line.strip().replace(',', '.').split(" ")
            for i in range(len(fields)):
                if fields[i][-1] in units.keys():
                    fields[i] = float(fields[i][:-1]) * 1000 ** units[fields[i][-1]]
                else:
                    fields[i] = float(fields[i])

            table.append(np.array(fields))

    if step > 1:
        table2 = list()
        table2.append(table[0])
        for i in range(1, len(table), step):
            nrow = table[i]
            for j in range(1, step):
                if i + j == len(table):
                    break
                nrow += table[i + j]
            nrow /= j + 1
            table2.append(nrow)
        table = table2

    t = np.arange(0, len(table), 1)
    men = (np.array([f[1] for f in table]) + np.array([f[2] for f in table])) / 1000000
    cpu = np.array([f[-1] for f in table]) / 100

    fig, ax1 = plt.subplots()
    title = f"{os.path.splitext(os.path.basename(path))[0]}\n" \
            f"Avg CPU: {round(numpy.average(cpu), 2)}      Max CPU: {round(max(cpu), 2)}\n" \
            f"Max Memory: {round(max(men), 2)} (GB)"
    plt.title(title)

    color = 'tab:blue'
    ax1.set_xlabel('')
    ax1.set_ylabel('CPU', color=color)
    ax1.plot(t, cpu, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:red'
    ax2.set_ylabel('Memory (GB)', color=color)
    ax2.plot(t, men, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()
    plt.savefig(path + "." + str(step) + ".png")


def main():
    if len(sys.argv) < 2:
        print("use <file|folder> [step]")
        exit(-1)

    path = sys.argv[1]
    if len(sys.argv) > 2:
        step = int(sys.argv[2])
    else:
        step = 1

    if os.path.isfile(path):
        parse(path, step)
    else:
        for file in os.listdir(path):
            if file.endswith(".prof"):
                try:
                    parse(os.path.join(path, file), step)
                except:
                    pass


if __name__ == '__main__':
    main()
