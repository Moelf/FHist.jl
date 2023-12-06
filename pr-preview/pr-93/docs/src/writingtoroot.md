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

### Writing out a histogram with errors (sumW2)
[ref](https://github.com/scikit-hep/uproot5/issues/696#issuecomment-1235918878)

Unfourtunately `uproot` does not support this yet, so we have to use pyROOT:
```julia
using PythonCall
ROOT = pyimport("ROOT")

# h is a FHist.jl histogram but you just need two arrays
bc = bincounts(h)
be = binerrors(h)

file = ROOT.TFile("/tmp/example.root", "recreate")
# 100 bins from 0 to 1
th1d = ROOT.TH1D("blah", "blah blah", 100, 0, 1)
for i in eachindex(bc)
    # ROOT has under/overflow bins outside of normal range
    # so we're skipping `0`-th bin by adhering Julia's 1-based index
    th1d.SetBinContent(i, bc[i])
    # similarly, we're skipping overflow bin
    th1d.SetBinError(i, be[i])
end
th1d.Write()
file.Close()
```

we can veryfy the round trip gives us back the original number:
```julia
julia> binerrors(h)
100-element Vector{Float64}:
 0.0
 0.000339284442870306
 0.0003276982770686565
 0.0020320752377659085
...

julia> round_h = UnROOT.parseTH(ROOTFile("/tmp/example.root")["blah"]);

# the third array we read back is the sumW2, so we take sqrt() to recover error
julia> sqrt.(round_h[3])
100-element Vector{Float64}:
 0.0
 0.000339284442870306
 0.0003276982770686565
 0.0020320752377659085
 0.0021085552210994537
...
```
