import ctypes
import numpy as np
import timeit

lib = ctypes.CDLL('./libfhistjl.so')

lib.hist1d.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_long, # input and length
        ctypes.POINTER(ctypes.c_double), ctypes.c_long, # bincoutns and length
        ctypes.c_double, ctypes.c_double, ctypes.c_double # binedges
        ]
lib.hist1d.restype = ctypes.c_int

def jlhist(a, bins, range):
    bincounts = np.zeros(bins)
    step = (range[1] - range[0]) / bins

    res = lib.hist1d(
            input_data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.c_long(len(input_data)), 
            bincounts.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.c_long(len(bincounts)), 
            ctypes.c_double(range[0]), 
            ctypes.c_double(step),
            ctypes.c_double(range[1])
            )
    return bincounts


input_data = np.random.randn(10**4)

np_timer = timeit.Timer(lambda: np.histogram(input_data, bins=10, range=(0.0, 1.0)))
print(np_timer.timeit(1000)/1000)

jl_timer = timeit.Timer(lambda: jlhist(input_data, bins=10, range=(0.0, 1.0)))
print(jl_timer.timeit(1000)/1000)
