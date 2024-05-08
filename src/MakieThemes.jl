"""

## Example
```
with_theme(ATLASTHEME) do
    h1 = Hist1D(randn(10^4))
    hist(h1; label="atlas style histogram")
end
```
"""
const ATLASTHEME = MakieCore.Attributes(
      Axis = (
              xtickalign=true, ytickalign=true,
              xticksmirrored=true, yticksmirrored=true,
              xminortickalign=1, yminortickalign=1,
              xticksize=10, yticksize=10,
              xminorticksize=6, yminorticksize=6,
              xgridvisible = false, ygridvisible = false,
              xminorticksvisible = true, yminorticksvisible = true,
             ),
      Colorbar = (
                  colormap = :haline,
                  highclip = :red,
                  lowclip = :black
                 )
     )
