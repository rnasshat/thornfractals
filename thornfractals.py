from math import pi, sin, cos
import random
import numpy as np
from PIL import Image

#FIXME code style is all over the place

def safe_div(a, b):
    return a / b if b else 0

def thorn(x = 0.0, y = 0.0, c = (0.01, -0.01), iterations = 256, bailout = 10000):
    numiter = 0
    for i in range(iterations):
        (prevx, prevy) = (x, y)
        numiter = i
        
        x = safe_div(prevx, cos(prevy)) + c[0]
        y = safe_div(prevy, sin(prevx)) + c[1]

        if x**2 + y**2 > bailout:
            break
        
    return numiter

def thorn_alt(x = 0.0, y = 0.0, c = (0.0, 0.0), iterations = 256, bailout = 10000):
    numiter = 0
    for i in range(iterations):
        prev_c = c
        numiter = i
        
        c = (safe_div(prev_c[0], cos(prev_c[1])) + x, safe_div(prev_c[1], sin(prev_c[0])) + y)

        if c[0]**2 + c[1]**2 > bailout:
            break
        
    return numiter

def rgb_sines(val = 0.0, frequency = (5, 7, 13), phaseshift = (0, 0, 0)):
    return np.array((
        (sin(val * frequency[0]*2*pi + phaseshift[0]*2*pi)+1)/2,
        (sin(val * frequency[1]*2*pi + phaseshift[1]*2*pi)+1)/2,
        (sin(val * frequency[2]*2*pi + phaseshift[2]*2*pi)+1)/2
    ))

def drawfractal(
    fractal_func = thorn, color_func = rgb_sines,
    resolution = (512, 512), c = (0.0, 0.0), plane = (-pi, pi, -pi, pi), iterations = 256, bailout = 10000,
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
                    out += color_func(
                        fractal_func(curx + xsamples*(stepx/supersample), cury + ysamples*(stepy/supersample), c, iterations, bailout)/(iterations-1),
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

def tests(): #TODO clean this up
    drawfractal(resolution=(1024, 768), c=(0.102, -0.04), plane=(-pi, pi, -0.5*pi, 0.5*pi))
    #control inputs for reference against images on the linked site
    
    drawfractal(resolution=(1600, 900), c=(-0.1, 0.15), plane=(0, pi, -0.5*pi, 0.5*pi), supersample=4)
    drawfractal(resolution=(2048, 2048), c=(-0.1, 0.15), plane=(0.75, 1.2566370614359172, -0.7853981633974483, -0.1), rgbphase=(1.85/2, 0.333/2, 0.666/2), supersample=2)
    drawfractal(resolution=(2048, 2048), c=(-0.1, 0.15), plane=(0.333, 1.2566370614359172, -1.0995574287564276, 0), rgbphase=(-0.069/2, -0.5/2, 0.555/2), supersample=2)
    drawfractal(resolution=(2048, 2048), c=(-0.1, 0.15), plane=(-pi, 0, -0.5*pi, 0.5*pi), rgbfreq=(2, 5, 23))
    drawfractal(resolution=(4096, 4096), c=(0,0), plane=(-pi, pi, -0.5*pi, 1.5*pi), rgbphase=(1.46/2, 1.46/2, 1.46/2))
    drawfractal(resolution=(4096, 4096), c=(0.1*pi, -0.05*pi), plane=(-pi, pi, -0.5*pi, 1.5*pi), rgbphase=(1.46/2, 1.46/2, 1.46/2))
