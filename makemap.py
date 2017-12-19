#!/usr/bin/env python3

# Create mapfiles for island speciation model
# Ludwig Leidinger 2017 <l.leidinger@gmx.net>

def print_map(xlen, ylen, landtype):
    print("#", landtype, ":")
    id = 0
    mintemp = 273
    maxtemp = 303
    tempstep = round((maxtemp - mintemp) / ylen)
    temp = mintemp
    for x in range(xlen):
        temp += tempstep
        for y in range(ylen):
            id += 1
            print(id, x, y, temp, landtype)
    print()

def add_ocean():
    print("# dummy island:")
    id = 0
    temp = 373
    x = 1000
    y = 0
    print(id, x, y, "island")
    print()

print_map(9, 9, "continent")
add_ocean()
