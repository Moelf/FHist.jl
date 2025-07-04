## Setup Julia nightly:
Follow instruction at https://julialang.org/install/
```
curl -fsSL https://install.julialang.org | sh
juliaup add nightly
```

## How to test

```bash
julia +nightly --project=.. ~/.julia/juliaup/julia-nightly/share/julia/juliac/juliac.jl --output-lib libfhistjl.so --compile-ccallable --experimental --trim C_FHist.jl

# these are numpy.histogram time and julia time respectively
python test.py
5.117958411574364e-05
3.2396082766354086e-05
```
