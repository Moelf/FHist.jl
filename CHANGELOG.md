## 0.12.0
- Exact bin assignment for values on uniform bin edges, `NaN`/`Inf` handling unified across uniform and non-uniform binnings.
- `Hist2D`/`Hist3D` `nentries` no longer counts discarded entries; `empty!` resets `nentries`.
- `Hist3D` auto binning fixed (was one bin per axis).
- Arithmetic keeps uniform bin edges; `*` allows negative bin contents.
- New 3D methods: `rebin`, `restrict`, `append!`, `mean`/`std`/`median`; edge based `rebin` for 2D/3D; `integral`/`normalize` with `width` for 2D/3D.
- GPU histogramming extension (`Hist1D(gpu_array; ...)`, `gpu_bincounts`).

## 0.9.00
- `normalize(; width=true)` becomes the default instead of `width=false`, see: https://github.com/Moelf/FHist.jl/issues/78
