#!/usr/bin/env python3

import sys
import os


def show(name, secs, out):
    d = int(secs / 3600 / 24)
    h = int(secs / 3600) % 24
    m = int(secs / 60) % 60
    s = round(secs % 60, 2)

    row = name + ": "
    if d > 0:
        row += str(d) + " days, "
    if h > 0:
        row += str(h) + " hours, "
    if m > 0:
        row += str(m) + " minutes, "
    row += str(s) + " seconds"
    print(row, file=out)


def parse(file, short):
    total = float(0)
    table = []
    row_check = False
    relative = False
    last_secs = 0

    with open(file) as finput, open(file + '.sum', 'w') as out:
        for line in finput:
            if "seconds" in line and line.strip()[0].isdigit():
                secs_s, name = line.split("seconds")[0:2]
                name = name[1:].strip()
                secs = float(secs_s)
                if not relative:
                    secs, last_secs = secs - last_secs, secs
                total += secs
                if len(table) == 0:
                    row_check = True
                    table.append((name, secs))
                else:
                    rname = table[-1][0]
                    rsecs = table[-1][1]

                    if short and "Tree partitioned" in name:
                        table[-1] = (rname, rsecs + secs)
                        row_check = True
                        continue

                    if row_check:
                        for i in range(min(len(rname), len(name))):
                            if rname[i] != name[i]:
                                break

                        if i > 0 and "," in name and not short:
                            i = name.index(',') + 1

                        if i > 2 and name[i - 2] == ' ' and short:
                            i -= 1

                        if i > 4 and (name[i - 1] == ' ' or name[i - 1] == ','):
                            table[-1] = (rname[:i].strip(), rsecs + secs)
                            row_check = short
                            continue

                    if len(rname) <= len(name) and rname == name[:len(rname)]:
                        table[-1] = (rname, rsecs + secs)
                    else:
                        row_check = True
                        table.append((name.split(",")[0], secs))
            elif len(table) == 0:
                relative = relative or "-relative-progress" in line
                print(line.strip(), file=out)

        tab = 0
        for name, secs in table:
            if len(name) > tab:
                tab = len(name)

        for name, secs in table:
            show(" " * (tab - len(name)) + name, secs, out)
        show(" " * (tab - len("Total")) + "Total", total, out)


def main():
    if len(sys.argv) < 2:
        print("use <file|folder> [-s]\n  -s: shorter summary")
        exit(-1)

    path = sys.argv[1]
    short = len(sys.argv) > 2 and "-s" == sys.argv[2]

    if os.path.isfile(path):
        parse(path, short)
    else:
        for file in os.listdir(path):
            if file.endswith(".out"):
                try:
                    parse(os.path.join(path, file), short)
                except:
                    pass


if __name__ == '__main__':
    main()
