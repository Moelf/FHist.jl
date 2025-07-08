## Setup Julia nightly:
Follow instruction at https://julialang.org/install/
```
curl -fsSL https://install.julialang.org | sh
juliaup add 1.12
```

## How to test

```bash
>julia +1.12 --startup-file=no --project=.. -e "using Pkg; Pkg.update()"
# the the -beta4+0 is the latest version at the time of writing, change it to the latest version
# the x64.linux.gnu need to be changed to your platform
>julia +1.12 --project=.. ~/.julia/juliaup/julia-1.12.0-beta4+0.x64.linux.gnu/share/julia/juliac.jl --output-lib libfhistjl.so --compile-ccallable --experimental --trim C_FHist.jl

>pixi ls -x
Package  Version  Build            Size     Kind   Source
numpy    2.3.1    py313h41a2e72_0  6.3 MiB  conda  numpy

> pixi run python test.py
```

## Results on macOS ARM (Mac Mini M4)
```
> julia +1.12 -e "using InteractiveUtils; versioninfo()"
Julia Version 1.12.0-beta4
Commit 600ac61d3d2 (2025-06-05 07:03 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: macOS (arm64-apple-darwin24.0.0)
  CPU: 10 × Apple M4
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, apple-m1)
  GC: Built with stock GC

> pixi run python test.py
=====================================
Input size: 1000
All close: True
Numpy    time (μs): 0.021641701459884644
FHist.jl time (μs): 0.003958120942115784
=====================================
Input size: 10000
All close: True
Numpy    time (μs): 0.05824156105518341
FHist.jl time (μs): 0.011725164949893951
=====================================
Input size: 100000
All close: True
Numpy    time (μs): 0.438716821372509
FHist.jl time (μs): 0.0859750434756279
=====================================
Input size: 1000000
All close: True
Numpy    time (μs): 3.94724179059267
FHist.jl time (μs): 0.829133577644825
```

## Results on Linux x86_64 (Ryzen 9 3900X)

```
> julia +1.12 -e "using InteractiveUtils; versioninfo()"
Julia Version 1.12.0-beta4
Commit 600ac61d3d2 (2025-06-05 07:03 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 24 × AMD Ryzen 9 3900X 12-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver2)
  GC: Built with stock GC

> pixi run python test.py
=====================================
Input size: 1000
All close: True
Numpy    time (μs): 0.05090880440548062
FHist.jl time (μs): 0.007097609341144562
=====================================
Input size: 10000
All close: True
Numpy    time (μs): 0.19860140746459365
FHist.jl time (μs): 0.020549201872199774
=====================================
Input size: 100000
All close: True
Numpy    time (μs): 1.3177349930629134
FHist.jl time (μs): 0.19442759221419692
=====================================
Input size: 1000000
All close: True
Numpy    time (μs): 6.48591200588271
FHist.jl time (μs): 1.8390380078926682
```
