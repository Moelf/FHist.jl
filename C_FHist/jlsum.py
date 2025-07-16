import ctypes

lib = ctypes.CDLL('./libjlsum.so')
lib.jlsum.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_long,
        ]
lib.jlsum.restype = ctypes.c_double
def jlsum(a):
    res = lib.jlsum(
            a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.c_long(len(a)), 
            )
    return res


import numpy as np
import timeit
for N_input in [10**4, 10**5, 10**6]:
    input_data = np.random.rand(N_input)
    print("=====================================")
    print(f"Input size: {N_input}")
    print("np.isclose:", np.isclose(np.sum(input_data), jlsum(input_data)))
    np_timer = timeit.Timer(lambda: np.sum(input_data))
    jl_timer = timeit.Timer(lambda: jlsum(input_data))
    print(f"Numpy    time (μs): {np.min(np_timer.repeat(number=5, repeat=500)) / 5 * 1000}")
    print(f"FHist.jl time (μs): {np.min(jl_timer.repeat(number=5, repeat=500)) / 5 * 1000}")
