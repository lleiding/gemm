#!/usr/bin/env python3

# Create mapfiles for island speciation model
# Ludwig Leidinger 2017 <l.leidinger@gmx.net>

import argparse
import math

def print_map(xlen, ylen, landtype, xpos, ypos, ident):
    print("#", landtype, ":")
    mintemp = 273
    maxtemp = 303
    temp = mintemp
    if landtype == "continent":
        tempstep = round((maxtemp - mintemp) / ylen)
        for x in range(xpos, xpos + xlen):
            temp += tempstep
            for y in range(ypos, ypos + ylen):
                ident += 1
                print(ident, x, y, temp, landtype)
    else: # perhaps just calculate distance to shortedge?
        shortedge = min([xlen, ylen])
        steps = math.ceil(shortedge / 2)
        temp = 298
        step = 0
        for x in range(xpos, xpos + xlen):
            step += 1
            for y in range(ypos, ypos + ylen):
                ident += 1
                print(ident, x, y, temp, landtype)
    print()

def add_ocean():
    print("# dummy island:")
    ident = 9999
    temp = 1273
    x = 1000
    y = 0
    print(ident, x, y, temp, "island")
    print()

def print_all():
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--width", type = int, default = 9,
                        help = "number of cells in x direction")
    parser.add_argument("-b", "--height", type = int, default = 9,
                        help = "number of cells in y direction")
    parser.add_argument("-l", "--land", type = str, default = "continent",
                        choices = ["continent", "island"], help = "type of landmass")
    parser.add_argument("-t", "--time", type = int, default = 1000,
                        help = "number of timesteps for which definition is valid")
    parser.add_argument("-x", "--longitute", type = int, default = 0,
                        help = "x position")
    parser.add_argument("-y", "--latitude", type = int, default = 0,
                        help = "y position")
    parser.add_argument("-i", "--identifier", type = int, default = 0,
                        help = "start of identifying number")
    args = parser.parse_args()
    print("# timesteps:")
    print(args.time)
    print()
    print_map(args.width, args.height, args.land, args.longitute, args.latitude, args.identifier)
    args.land == "continent" and add_ocean()

print_all()
