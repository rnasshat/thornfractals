from math import pi, tau, sin, cos
import random
import numpy as np
from PIL import Image

#FIXME code style is all over the place

def safe_div(a, b):
    return a / b if b else 0

def linear_to_srgb(x):
    if x >= 0.0031308:
        return min(max(1.055 * (x**(1.0/2.4)) - 0.055, 0.0), 1.0)
    else:
        return min(max(12.92 * x, 0.0), 1.0)

def srgb_to_linear(x):
    x = min(max(x, 0.0), 1.0)
    if x >= 0.04045:
        return ((x + 0.055)/(1.055))**2.4
    else:
        return x / 12.92

def oklch_to_rgb(color = (1, 0, 0)): #Modified from code on https://bottosson.github.io/posts/oklab/
    a = color[1] * cos(np.deg2rad(color[2]))
    b = color[1] * sin(np.deg2rad(color[2]))

    l_ = color[0] + 0.3963377774 * a + 0.2158037573 * b
    m_ = color[0] - 0.1055613458 * a - 0.0638541728 * b
    s_ = color[0] - 0.0894841775 * a - 1.2914855480 * b
    
    l = l_*l_*l_
    m = m_*m_*m_
    s = s_*s_*s_
    
    return (
        +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
		-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
		-0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s)

def thorn(x = 0.0, y = 0.0, c = (0.01, -0.01), iterations = 256, bailout = 10000):
    """Implementation of the "Thorn fractal" or "Secant Sea" as described on http://paulbourke.net/fractals/thorn/
    plane format is (xmin, xmax, ymin, ymax)
    coordinates start in the upper left corner, positive y is down
    """
    numiter = 0
    for i in range(iterations):
        (prevx, prevy) = (x, y)
        numiter = i
        
        x = safe_div(prevx, cos(prevy)) + c[0]
        y = safe_div(prevy, sin(prevx)) + c[1]

        if x**2 + y**2 > bailout:
            break
        
    return numiter/(iterations-1)

def thorn_alt(x = 0.0, y = 0.0, c = (0.0, 0.0), iterations = 256, bailout = 10000):
    """Variant on the thorn fractal ("Mandelbrot" type)"""
    numiter = 0
    for i in range(iterations):
        prev_c = c
        numiter = i
        
        c = (safe_div(prev_c[0], cos(prev_c[1])) + x, safe_div(prev_c[1], sin(prev_c[0])) + y)

        if c[0]**2 + c[1]**2 > bailout:
            break
        
    return numiter/(iterations-1)

def rgb_sines(val = 0.0, frequency = (5, 7, 13), phase = (-0.25, -0.25, -0.25)):
    """coloring function"""
    return np.array((
        (sin(val * frequency[0]*tau + phase[0]*tau)+1)/2,
        (sin(val * frequency[1]*tau + phase[1]*tau)+1)/2,
        (sin(val * frequency[2]*tau + phase[2]*tau)+1)/2))
    
def oklch_cycle(val = 0.0, l = 1, c = 0.34, frequency = 360, offset = 0):
    """coloring function"""
    return np.array(oklch_to_rgb((l, c, val * frequency + offset)))

def drawfractal(
    resolution = (512, 512), plane = (-pi, pi, -pi, pi), supersample = 1, *,
    fractal_func = thorn, color_func = rgb_sines,
    fractal_func_args = None, color_func_args = None):
    
    if color_func_args is None:
        color_func_args = {}
    if fractal_func_args is None:
        fractal_func_args = {}
    a = np.zeros((resolution[1], resolution[0], 3)) #Array indices are backwards relative to image coordinates

    stepy = (plane[3]-plane[2]) / resolution[1]
    stepx = (plane[1]-plane[0]) / resolution[0]
    for y in range(resolution[1]):
        cury = plane[2] + y*stepy
        for x in range(resolution[0]):
            curx = plane[0] + x*stepx

            out = color_func( #color_func expects a single float followed by arguments
                            fractal_func( #fractal_func expects two floats x, y followed by arguments
                                        curx,
                                        cury,
                                        **fractal_func_args),
                            **color_func_args)
            if supersample > 1:
                for ysamples in range(supersample + 1):
                    for xsamples in range(supersample + 1):
                        out += color_func(
                                        fractal_func(
                                                    curx - stepx/2.0 + xsamples*stepx/supersample,
                                                    cury - stepy/2.0 + ysamples*stepy/supersample,
                                                    **fractal_func_args),
                                        **color_func_args)
                out /= (supersample + 1)**2 + 1
            
            a[y, x, 0] = linear_to_srgb(out[0])
            a[y, x, 1] = linear_to_srgb(out[1])
            a[y, x, 2] = linear_to_srgb(out[2])
        print(y+1, "/", resolution[1], " rows complete", end="\r")
    print("\n")
    
    im = Image.fromarray((a * 255).astype(np.uint8))

    fname = "{0}{1}x{2}px {3} plane={4}.png".format(
                                                    fractal_func.__name__,
                                                    resolution[0],
                                                    resolution[1],
                                                    fractal_func_args["c"] if "c" in fractal_func_args else "",
                                                    plane)
    im.save(fname)

def unitcircle(theta = 0.0):
    return (cos(theta), sin(theta))

def unitcirclefuzzy(theta = 0.0):
    return (cos(theta)+random.uniform(-0.1, 0.1), sin(theta)+random.uniform(-0.1, 0.1))

def tests():
    drawfractal( #control inputs
                resolution=(1024, 768),
                plane=(-pi, pi, -0.5*pi, 0.5*pi),
                fractal_func_args={"c": (0.102, -0.04)})
    drawfractal(
                resolution=(1600, 900),
                plane=(0, pi, -0.5*pi, 0.5*pi),
                supersample=4,
                fractal_func_args={"c": (-0.1, 0.15)})
    drawfractal(
                resolution=(2048, 2048),
                plane=(0.75, 1.2566370614359172, -0.7853981633974483, -0.1),
                supersample=2,
                fractal_func_args={"c": (-0.1, 0.15)},
                color_func_args={"phase": (1.85/2, 0.333/2, 0.666/2)})
    drawfractal(
                resolution=(2048, 2048),
                plane=(0.333, 1.2566370614359172, -1.0995574287564276, 0),
                supersample=2,
                fractal_func_args={"c": (-0.1, 0.15)},
                color_func_args={"phase": (-0.069/2, -0.5/2, 0.555/2)})
    drawfractal(
                resolution=(2048, 2048),
                plane=(-pi, 0, -0.5*pi, 0.5*pi),
                fractal_func_args={"c": (-0.1, 0.15)},
                color_func_args={"phase": (-0.125, -0.125, -0.125), "frequency": (2, 5, 23)})
    drawfractal(
                resolution=(4096, 4096),
                plane=(-pi, pi, -0.5*pi, 1.5*pi),
                fractal_func_args={"c": (0,0)},
                color_func_args={"phase": (1.46/2, 1.46/2, 1.46/2)})
    drawfractal(
                resolution=(4096, 4096),
                plane=(-pi, pi, -0.5*pi, 1.5*pi),
                fractal_func_args={"c": (0.1*pi, -0.05*pi)},
                color_func_args={"phase": (1.46/2, 1.46/2, 1.46/2)})
    drawfractal(
                resolution=(1024, 1024),
                supersample=3,
                fractal_func=thorn_alt)
    drawfractal(
                resolution=(1024, 1024),
                supersample=3,
                fractal_func=thorn_alt,
                fractal_func_args={"c": (1, 0.1)},
                color_func_args={"frequency": (7, 5, 7)})
