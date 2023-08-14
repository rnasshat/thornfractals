from math import pi, sin, cos
import random
import numpy as np
from PIL import Image

def thorn(resolution = (512, 512), iterations = 256, c = (0.01, -0.01), plane = (-pi, pi, -pi, pi), bailout = 10000):
    """ Implementation of the "Thorn fractal" or "Secant Sea" as described on http://paulbourke.net/fractals/thorn/ """
    a = np.zeros((resolution[1], resolution[0], 3)) #Array indices are backwards relative to image coordinates

    for y in range(resolution[1]):
        stepy = plane[2] + y*(plane[3]-plane[2]) / resolution[1]
        for x in range(resolution[0]):
            stepx = plane[0] + x*(plane[1]-plane[0]) / resolution[0]

            (outx, outy) = (stepx, stepy) #TODO break up into several functions for later multithreading / extension
            numiter = 0
            for i in range(iterations):
                (prevx, prevy) = (outx, outy)
                numiter = i
                
                try:
                    outx = prevx / cos(prevy) + c[0]
                    outy = prevy / sin(prevx) + c[1]
                except ZeroDivisionError:
                    break

                if outx*outx + outy*outy > bailout:
                    break

            out = (float(numiter)/(iterations-1)) #TODO parameterize colorization
            a[y, x, 0] = (sin(2*pi * out * 5 - 0.25*pi)+1)/2
            a[y, x, 1] = (sin(2*pi * out * 7 - 0.25*pi)+1)/2
            a[y, x, 2] = (sin(2*pi * out * 13 - 0.25*pi)+1)/2


    im = Image.fromarray((a * 255).astype(np.uint8))

    fname = "thorn{0}x{1}px{2}i c{3} {4}.png".format(resolution[0], resolution[1], iterations, c, plane)
    im.save(fname)

#control inputs for reference against images on the linked site
thorn(resolution = (1024, 768), c = (0.102, -0.04), plane = (-pi, pi, -0.5*pi, 0.5*pi))

def unitcircle(theta = 0.0):
    return (cos(theta), sin(theta))

def unitcirclefuzzy(theta = 0.0):
    return (cos(theta)+random.uniform(-0.1, 0.1), sin(theta)+random.uniform(-0.1, 0.1))

def tests():
    thorn((2048, 2048), c=(-0.1, 0.15), plane=(0, pi, -0.5*pi, 0.5*pi))
    thorn((2048, 2048), c=(-0.1, 0.15), plane=(0.75, 1.2566370614359172, -0.7853981633974483, -0.1))
    thorn((2048, 2048), c=(-0.1, 0.15), plane=(0.333, 1.2566370614359172, -1.0995574287564276, 0))
    thorn((2048, 2048), c=(-0.1, 0.15), plane=(-pi, 0, -0.5*pi, 0.5*pi))
    thorn((4096, 4096), c=(0,0), plane=(-pi, pi, -0.5*pi, 1.5*pi))
    thorn((4096, 4096), c=((0.1*pi, -0.05*pi)), plane=(-pi, pi, -0.5*pi, 1.5*pi))
