## Setup Julia nightly:
Follow instruction at https://julialang.org/install/
```
curl -fsSL https://install.julialang.org | sh
juliaup add nightly
```

## How to test

```bash
>julia +nightly --project=.. ~/.julia/juliaup/julia-nightly/share/julia/juliac/juliac.jl --output-lib libfhistjl.so --compile-ccallable --experimental --trim C_FHist.jl

>pixi ls -x
Package  Version  Build            Size     Kind   Source
numpy    2.3.1    py313h41a2e72_0  6.3 MiB  conda  numpy

> pixi run python test.py
```


## Results on macOS ARM (Mac Mini M4)

```
> julia +nightly -e "using InteractiveUtils; versioninfo()"
Julia Version 1.13.0-DEV.817
Commit c0cc1e1022b (2025-07-04 12:09 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: macOS (arm64-apple-darwin24.0.0)
  CPU: 10 × Apple M4
  LLVM: libLLVM-20.1.2 (ORCJIT, apple-m4)

> pixi run python test.py
=====================================
Input size: 1000
All close: True
Numpy    time (μs): 0.01876652240753174
FHist.jl time (μs): 0.0034248456358909607
=====================================
Input size: 10000
All close: True
Numpy    time (μs): 0.05549192428588867
FHist.jl time (μs): 0.010091438889503479
=====================================
Input size: 100000
All close: True
Numpy    time (μs): 0.4449000582098961
FHist.jl time (μs): 0.08234996348619461
=====================================
Input size: 1000000
All close: True
Numpy    time (μs): 3.905266709625721
FHist.jl time (μs): 0.7887914776802063
```


## Results on Linux x86_64 (Ryzen 9 3900X)

```
> julia +nightly -e "using InteractiveUtils; versioninfo()"
Julia Version 1.13.0-DEV.819
Commit 4846c3d938b (2025-07-04 15:43 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 24 × AMD Ryzen 9 3900X 12-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-20.1.2 (ORCJIT, znver2)
  GC: Built with stock GC

> pixi run python test.py
=====================================
Input size: 1000
All close: True
Numpy    time (μs): 0.04856220039073378
FHist.jl time (μs): 0.00662259990349412
=====================================
Input size: 10000
All close: True
Numpy    time (μs): 0.2016977989114821
FHist.jl time (μs): 0.04962040111422539
=====================================
Input size: 100000
All close: True
Numpy    time (μs): 1.333067798987031
FHist.jl time (μs): 0.12432239891495554
=====================================
Input size: 1000000
All close: True
Numpy    time (μs): 6.55430739861913
FHist.jl time (μs): 1.6349460027413443
```
