"""

## Example
```
with_theme(ATLASTHEME) do
    h1 = Hist1D(randn(10^4))
    hist(h1; label="atlas style histogram")
end
```
"""
ATLASTHEME = nothing
