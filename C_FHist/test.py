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


for N_input in [10**3, 10**4, 10**5, 10**6]:
    input_data = np.random.randn(N_input)
    print("=====================================")
    print(f"Input size: {N_input}")
    print("All close:", np.allclose(np.histogram(input_data, bins=10, range=(0.0, 1.0))[0], jlhist(input_data, bins=10, range=(0.0,
                                                                                                                 1.0))))
    np_timer = timeit.Timer(lambda: np.histogram(input_data, bins=10, range=(0.0, 1.0)))
    jl_timer = timeit.Timer(lambda: jlhist(input_data, bins=10, range=(0.0, 1.0)))
    print(f"Numpy    time (μs): {np.min(np_timer.repeat(number=2, repeat=500)) * 1000}")
    print(f"FHist.jl time (μs): {np.min(jl_timer.repeat(number=2, repeat=500)) * 1000}")
