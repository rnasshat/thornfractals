from math import pi, sin, cos
import random
import numpy as np
from PIL import Image

def thorn(x = 0.0, y = 0.0, c = (0.01, -0.01), iterations = 256, bailout = 10000):
    numiter = 0
    for i in range(iterations):
        (prevx, prevy) = (x, y)
        numiter = i
        
        try:
            x = prevx / cos(prevy) + c[0]
            y = prevy / sin(prevx) + c[1]
        except ZeroDivisionError:
            break

        if x**2 + y**2 > bailout:
            break
        
    return numiter

def colorfunc(val = 0.0, frequency = (5, 7, 13), phaseshift = (-0.25, -0.25, -0.25)):
    return np.array((
        (sin(2*pi * val * frequency[0] + phaseshift[0]*pi)+1)/2,
        (sin(2*pi * val * frequency[1] + phaseshift[1]*pi)+1)/2,
        (sin(2*pi * val * frequency[2] + phaseshift[2]*pi)+1)/2
    ))

def drawfractal(
    resolution = (512, 512), c = (0.01, -0.01), plane = (-pi, pi, -pi, pi), iterations = 256, bailout = 10000,
    rgbfreq = (5, 7, 13), rgbphase = (-0.25, -0.25, -0.25), supersample = 1
):
    """Implementation of the "Thorn fractal" or "Secant Sea" as described on http://paulbourke.net/fractals/thorn/
    plane format is (xmin, xmax, ymin, ymax)
    coordinates start in the upper left corner, positive y is down
    """
    a = np.zeros((resolution[1], resolution[0], 3)) #Array indices are backwards relative to image coordinates

    stepy = (plane[3]-plane[2]) / resolution[1]
    stepx = (plane[1]-plane[0]) / resolution[0]
    for y in range(resolution[1]):
        cury = plane[2] + y*stepy
        for x in range(resolution[0]):
            curx = plane[0] + x*stepx

            out = np.zeros(3)
            for ysamples in range(supersample):
                for xsamples in range(supersample):
                    out += colorfunc(
                        thorn(curx + xsamples*(stepx/supersample), cury + ysamples*(stepy/supersample), c, iterations, bailout)/(iterations-1),
                        rgbfreq,
                        rgbphase
                        )
            out /= (supersample**2)
                    
            a[y, x, 0] = out[0]
            a[y, x, 1] = out[1]
            a[y, x, 2] = out[2]

        print(y+1, "/", resolution[1], " rows complete", end="\r")

    print("\n")
    
    im = Image.fromarray((a * 255).astype(np.uint8))

    fname = "thorn{0}x{1}px{2}i c{3} {4}.png".format(resolution[0], resolution[1], iterations, c, plane)
    im.save(fname)

def unitcircle(theta = 0.0):
    return (cos(theta), sin(theta))

def unitcirclefuzzy(theta = 0.0):
    return (cos(theta)+random.uniform(-0.1, 0.1), sin(theta)+random.uniform(-0.1, 0.1))

def tests():
    drawfractal((1024, 768), (0.102, -0.04), (-pi, pi, -0.5*pi, 0.5*pi))    #control inputs for reference against images on the linked site
    
    drawfractal((1600, 900), (-0.1, 0.15), (0, pi, -0.5*pi, 0.5*pi), supersample=4)
    drawfractal((2048, 2048), (-0.1, 0.15), (0.75, 1.2566370614359172, -0.7853981633974483, -0.1), rgbphase=(1.85, 0.333, 0.666), supersample=2)
    drawfractal((2048, 2048), (-0.1, 0.15), (0.333, 1.2566370614359172, -1.0995574287564276, 0), rgbphase=(-0.069, -0.5, 0.555), supersample=2)
    drawfractal((2048, 2048), (-0.1, 0.15), (-pi, 0, -0.5*pi, 0.5*pi), rgbfreq=(2, 5, 23))
    drawfractal((4096, 4096), (0,0), (-pi, pi, -0.5*pi, 1.5*pi), rgbphase=(1.46, 1.46, 1.46))
    drawfractal((4096, 4096), (0.1*pi, -0.05*pi), (-pi, pi, -0.5*pi, 1.5*pi), rgbphase=(1.46, 1.46, 1.46))
