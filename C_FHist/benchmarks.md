## On a AMD Zen2 system:

### JULIA_CPU_TARGET='znver2;generic' julia +1.12
```
Input size: 1000000
All close: True
Numpy    time (μs): 6.558584200683981
FHist.jl time (μs): 1.7967560008401051
```

### JULIA_CPU_TARGET='generic;znver2' julia +1.12
```
Input size: 1000000
All close: True
Numpy    time (μs): 6.538515799911693
FHist.jl time (μs): 2.6081674004672095
```

### JULIA_CPU_TARGET='znver2;generic' julia +nightly
```
Input size: 1000000
All close: True
Numpy    time (μs): 6.641940600820817
FHist.jl time (μs): 1.5721895993920043
```

### JULIA_CPU_TARGET='generic;znver2' julia +nightly
```
Input size: 1000000
All close: True
Numpy    time (μs): 6.617425798322074
FHist.jl time (μs): 2.40858779870905
```
