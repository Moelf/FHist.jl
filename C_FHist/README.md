## Setup Julia nightly:
Follow instruction at https://julialang.org/install/
```
curl -fsSL https://install.julialang.org | sh
juliaup add nightly
```

## How to test

```bash
julia +nightly --project=.. ~/.julia/juliaup/julia-nightly/share/julia/juliac/juliac.jl --output-lib libfhistjl.so --compile-ccallable --experimental --trim C_FHist.jl

> pixi run python test.py
=====================================
Input size: 1000
All close: True
Numpy    time (μs): 0.04149973392486572
FHist.jl time (μs): 0.007666647434234619
=====================================
Input size: 10000
All close: True
Numpy    time (μs): 0.10570790618658066
FHist.jl time (μs): 0.02658367156982422
=====================================
Input size: 100000
All close: True
Numpy    time (μs): 0.8217496797442436
FHist.jl time (μs): 0.4670834168791771
=====================================
Input size: 1000000
All close: True
Numpy    time (μs): 10.170875117182732
FHist.jl time (μs): 7.914000190794468
```
