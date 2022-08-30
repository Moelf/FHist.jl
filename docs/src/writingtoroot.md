### Write out a histogram
Checkout configuration for PythonCall.jl: https://cjdoris.github.io/PythonCall.jl/stable/pythoncall/#pythoncall-config

Most importantly, you probably want to set:
```julia
ENV["JULIA_PYTHONCALL_EXE"] = readchomp(`which python`)
```
before the `using PythonCall` line. Especially if you're using LCG or Athena or CMSSW environment.

```julia
julia> using PythonCall, FHist

julia> np = pyimport("numpy")

julia> up = pyimport("uproot")

julia> h = Hist1D(rand(100))
              ┌                              ┐ 
   [0.0, 0.2) ┤██████████████████████▊ 22      
   [0.2, 0.4) ┤██████████████████▋ 18          
   [0.4, 0.6) ┤█████████████████▋ 17           
   [0.6, 0.8) ┤██████████████████▋ 18          
   [0.8, 1.0) ┤██████████████████████████  25  
              └                              ┘ 
edges: 0.0:0.2:1.0
bin counts: [22, 18, 17, 18, 25]
total count: 100

julia> pywith(up.recreate("./example.root")) do file
           file["myhist"] = np.array(bincounts(h)), np.array(binedges(h))
       end;
(<py array([22, 18, 17, 18, 25])>, <py array([0. , 0.2, 0.4, 0.6, 0.8, 1. ])>)
```
