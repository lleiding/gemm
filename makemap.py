#!/usr/bin/env python3

# Create mapfiles for island speciation model
# Ludwig Leidinger 2017 <l.leidinger@gmx.net>

import argparse
import math
import random

def print_map(xlen, ylen, landtype, xpos, ypos, ident, isol, precon):
    print("#", landtype, ":")
    mintemp = 273
    maxtemp = 303
    minprec = 0
    maxprec = 1000 if precon else minprec
    precstep = (maxprec - minprec) / (ylen - 1)
    temp = mintemp
    if landtype == "continent":
        prec = 0
        tempstep = round((maxtemp - mintemp) / ylen)
        for x in range(xpos, xpos + xlen):
            prec = 0
            temp += tempstep
            for y in range(ypos, ypos + ylen):
                ident += 1
                print(ident, x, y, temp, landtype, "no", prec)
                prec += precstep
    else:
        temp = 298
        tempstep = 3 # CAVE!
        for x in range(xpos, xpos + xlen):
            prec = 0
            for y in range(ypos, ypos + ylen):
                ident += 1
                mindist = min([abs(x - xpos), abs(y - ypos), abs(xpos + xlen - x - 1), abs(ypos + ylen - y - 1)])
                localtemp = temp - mindist * tempstep
                isolated = "isolated" if random.random() < isol else "no"
                print(ident, x, y, localtemp, landtype, isolated, prec)
                prec += precstep
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
                        help = "number of timesteps for which definition is valid. A value of zero turns off output")
    parser.add_argument("-x", "--longitute", type = int, default = 0,
                        help = "x position")
    parser.add_argument("-y", "--latitude", type = int, default = 0,
                        help = "y position")
    parser.add_argument("-n", "--identifier", type = int, default = 0,
                        help = "start of identifying number")
    parser.add_argument("-i", "--isolation", type = float, default = 0.0,
                        help = "proportion of randomly distributed, isolated patches on island")
    parser.add_argument("-p", "--precipitation", action = "store_true", default = False,
                        help = "turns on additional environment niche")
    parser.add_argument("--nodummy", action = "store_true", default = False,
                        help = "turns off creation of dummy island")
    args = parser.parse_args()
    if args.time > 0:
        print("# timesteps:")
        print(args.time)
        print()
    print_map(args.width, args.height, args.land, args.longitute, args.latitude, args.identifier, args.isolation, args.precipitation)
    args.land == "continent" and not args.nodummy and add_ocean()

print_all()
