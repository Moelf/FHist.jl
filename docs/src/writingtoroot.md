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
We need https://github.com/scikit-hep/hist for this:

```julia
julia> using PythonCall, FHist

julia> function to_pyhist(h::Hist1D)
    pyhist = pyimport("hist")
    bes, bcs, sumw2s = binedges(h), bincounts(h), sumw2(h)
    h1 = pyhist.Hist.new.Variable(collect(bes)).Weight()
    for idx in eachindex(bcs, sumw2s)
        h1.view()[idx-1] = (bcs[idx], sumw2s[idx])
    end
    h1
end

julia> pywith(up.recreate("./example.root")) do file
           file["myhist"] = to_pyhist(h)
       end;
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
