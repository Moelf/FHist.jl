### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ be0b45f2-86dc-11ee-2dce-2994fb540eec
let
    docs_dir = dirname(dirname(@__DIR__))
    pkg_dir = dirname(docs_dir)

    using Pkg: Pkg
    Pkg.activate(docs_dir)
    Pkg.develop(; path=pkg_dir)
    Pkg.instantiate()
end;

# ╔═╡ 5a884f80-7c0a-4345-b2d3-1785b1b5a46a
using Statistics

# ╔═╡ 58b588b6-4b8f-447d-b205-e31b52e99e4e
using FHist

# ╔═╡ cb7cecc9-e7ec-4555-83ef-5f70a84b4e15
using Plots

# ╔═╡ a354ab2a-e809-4b0f-829d-1f6b33ed5634

md"""
# Documents some of the functions of the FHist package.
"""

# ╔═╡ 2d857d08-f9c2-4a60-86dd-2afbd3b599f6
begin
	a = randn(100000);
	h1 = Hist1D(a, -4:0.1:4)
	Plots.plot(h1)
end

# ╔═╡ 75fbe0da-30f0-42e7-aba6-c37bd8da71db
md"""
    sample(h::Hist1D, n::Int=1)

Sample a histogram's with weights equal to bin count, `n` times.
The sampled values are the bins' lower edges.
"""

# ╔═╡ 64eae0d9-63d7-49f2-9855-9a35246f8c85
ns = sample(h1, n=100)

# ╔═╡ 56fde2f7-2f30-4e9a-ae94-6b665b99d64f
hs = Hist1D(ns, -4:0.5:4)

# ╔═╡ 4f96c103-5284-4ffb-8927-562d88dffd40
md"""
    bincounts(h::Hist1D)

Get the bin counts of a histogram.
"""

# ╔═╡ 6ac9e44c-f72c-4d27-adb4-c51faa92debd
bincounts(h1)

# ╔═╡ ae5c165f-31fa-43f9-8118-32f478b13e49
md"""
    binedges(h::Hist1D)

Get the bin edges of a histogram.
"""

# ╔═╡ 6e171319-93a2-4923-acd3-d09002548163
binedges(h1)

# ╔═╡ dd1d3023-65b4-48cc-b8de-2abc2ebb06d3
md"""
    bincenters(h::Hist1D)

Get the bin centers of a histogram.
"""

# ╔═╡ 45631c9d-a19b-43ee-84c6-97985fd5e8af
bincenters(h1)

# ╔═╡ 9f4a2537-59f8-4ab8-9822-85efa1b1668e
md"""
    binerrors(f::T, h::Hist1D) where T<:Function = f.(h.sumw2)
    binerrors(h::Hist1D) = binerrors(sqrt, h)

Calculate the bin errors from `sumw2` with a Gaussian default.
"""

# ╔═╡ d56f6346-cc84-469f-9921-c0de5f523819
binerrors(h1)

# ╔═╡ 75d72814-e864-4a6e-bf4b-e45882f10518
md"""
    nbins(h::Hist1D)

Get the number of bins of a histogram.
"""

# ╔═╡ 7caea1bb-ee35-406e-85bc-7cd234c6ad82
nbins(h1)

# ╔═╡ 3a164df3-7685-4299-b5ec-3e0787e10e69
md"""
    integral(h::Hist1D; width=false)

Get the integral a histogram; `width` means multiply each bincount
by their bin width when calculating the integral.

!!! warning
    Be aware of the approximation you make
    when using `width=true` with histogram with overflow bins, the overflow
    bins (i.e. the left/right most bins) width will be taken "as is".
"""

# ╔═╡ 6d545c67-e6d9-4464-bb79-a813a0ed1d7f
integral(h1, width=false)

# ╔═╡ aeabcde5-d2cb-4173-af04-824e21f224cd
md"""
    Statistics.mean(h)
    Statistics.std(h)
    Statistics.median(h)
    Statistics.quantile(h::Hist1D, p)

Compute statistical quantities based on the bin centers weighted
by the bin counts.

When the histogram is `Hist2D`, return tuple instead, e.g `(mean(project(h, :x)), mean(project(h, :y)))` etc.
"""

# ╔═╡ 7a7f21cb-54b0-474b-a77d-e3feb17120e9
mean(h1) 


# ╔═╡ 17aafd33-756d-443c-abd6-ab2926d71cb3
std(h1)

# ╔═╡ fa898cda-8504-43f4-b8db-391a497cad6e
median(h1)

# ╔═╡ d2f82cac-d710-4dcd-9710-e816372c06f6
quantile(h1,1)

# ╔═╡ a5a832b0-3ba5-4d13-8f2a-95fdbe4b925b
md"""
    lookup(h::Hist1D, x)

For given x-axis value `x`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""

# ╔═╡ d9a00887-b8a4-419b-806f-f150f8bbe813
lookup(h1, .2)

# ╔═╡ fa1b47db-94bb-4956-9857-7f3481069008
lookup(h1, 2)

# ╔═╡ 314ca5ce-b1d3-48e5-8f45-b2a7a05b8248
md"""
    normalize(h::Hist1D; width=true)

Create a normalized histogram via division by `integral(h)`, when `width==true`, the
resultant histogram has area under the curve equals 1.

!!! warning
    Implicit approximation is made when using `width=true` with histograms
    that have overflow bins: the overflow data lives inthe left/right most bins
    and the bin width is taken "as is".
"""

# ╔═╡ 79d4ab1f-cd3f-498c-9e01-fe944df7b94d
hn1 = normalize(h1);

# ╔═╡ 77b95061-0e0c-4dcc-97b9-ff5764564b96
typeof(hn1)

# ╔═╡ d926a1d2-fbc4-4d49-be39-7f2bbd484097
Plots.plot(hn1)

# ╔═╡ 3a9d52a4-8028-4bb9-a72b-a3d24e0db589
md"""
    cumulative(h::Hist1D; forward=true)

Create a cumulative histogram. If `forward`, start
summing from the left.
"""

# ╔═╡ 611690ff-3821-49ad-b6d8-1f0c053a6c39
cumulative(h1)

# ╔═╡ 69446fed-b9c3-4c1a-8c56-5b8331a73734
md"""
    rebin(h::Hist1D, n::Int=1)
    rebin(n::Int) = h::Hist1D -> rebin(h, n)

Merges `n` consecutive bins into one.
The returned histogram will have `nbins(h)/n` bins.
"""

# ╔═╡ d07dbae1-5917-42b4-8da2-10903d8dcd27
rebin(h1, 2)

# ╔═╡ e5cf45cf-798a-4c44-8536-fc3974582f34
md"""
    bayes_rebin_edges(h::Hist1D; prior=BayesHistogram.Geometric(0.995))

Find optimal bin edges for a histogram using Bayesian rebinning algorithm.
This function only find edges, it doesn't return a new histogram.

For possible priors, see [`BayesHistogram.jl`](https://github.com/francescoalemanno/BayesHistogram.jl/blob/main/src/BayesHistogram.jl).

"""

# ╔═╡ 4cb8a35c-727c-4ead-94d1-ebb890f40b14
function bayes_rebin_edges(h::Hist1D; prior=BayesHistogram.Geometric(0.995))
    old_edges = binedges(h)
    length(old_edges) < 4 && error("too little bins to rebin")
    fake_xs = [first(old_edges); bincenters(h); last(old_edges)]
    weights = [0; bincounts(h); 0]
    res = BayesHistogram.bayesian_blocks(fake_xs; weights=weights, prior=prior)
    return res.edges
end

# ╔═╡ 05476b2a-fbb1-4572-ae23-c1116dbbb386
FHist.bayes_rebin_edges(h1)

# ╔═╡ Cell order:
# ╠═be0b45f2-86dc-11ee-2dce-2994fb540eec
# ╠═5a884f80-7c0a-4345-b2d3-1785b1b5a46a
# ╠═58b588b6-4b8f-447d-b205-e31b52e99e4e
# ╠═cb7cecc9-e7ec-4555-83ef-5f70a84b4e15
# ╠═a354ab2a-e809-4b0f-829d-1f6b33ed5634
# ╠═2d857d08-f9c2-4a60-86dd-2afbd3b599f6
# ╠═75fbe0da-30f0-42e7-aba6-c37bd8da71db
# ╠═64eae0d9-63d7-49f2-9855-9a35246f8c85
# ╠═56fde2f7-2f30-4e9a-ae94-6b665b99d64f
# ╠═4f96c103-5284-4ffb-8927-562d88dffd40
# ╠═6ac9e44c-f72c-4d27-adb4-c51faa92debd
# ╠═ae5c165f-31fa-43f9-8118-32f478b13e49
# ╠═6e171319-93a2-4923-acd3-d09002548163
# ╠═dd1d3023-65b4-48cc-b8de-2abc2ebb06d3
# ╠═45631c9d-a19b-43ee-84c6-97985fd5e8af
# ╠═9f4a2537-59f8-4ab8-9822-85efa1b1668e
# ╠═d56f6346-cc84-469f-9921-c0de5f523819
# ╠═75d72814-e864-4a6e-bf4b-e45882f10518
# ╠═7caea1bb-ee35-406e-85bc-7cd234c6ad82
# ╠═3a164df3-7685-4299-b5ec-3e0787e10e69
# ╠═6d545c67-e6d9-4464-bb79-a813a0ed1d7f
# ╠═aeabcde5-d2cb-4173-af04-824e21f224cd
# ╠═7a7f21cb-54b0-474b-a77d-e3feb17120e9
# ╠═17aafd33-756d-443c-abd6-ab2926d71cb3
# ╠═fa898cda-8504-43f4-b8db-391a497cad6e
# ╠═d2f82cac-d710-4dcd-9710-e816372c06f6
# ╠═a5a832b0-3ba5-4d13-8f2a-95fdbe4b925b
# ╠═d9a00887-b8a4-419b-806f-f150f8bbe813
# ╠═fa1b47db-94bb-4956-9857-7f3481069008
# ╠═314ca5ce-b1d3-48e5-8f45-b2a7a05b8248
# ╠═79d4ab1f-cd3f-498c-9e01-fe944df7b94d
# ╠═77b95061-0e0c-4dcc-97b9-ff5764564b96
# ╠═d926a1d2-fbc4-4d49-be39-7f2bbd484097
# ╠═3a9d52a4-8028-4bb9-a72b-a3d24e0db589
# ╠═611690ff-3821-49ad-b6d8-1f0c053a6c39
# ╠═69446fed-b9c3-4c1a-8c56-5b8331a73734
# ╠═d07dbae1-5917-42b4-8da2-10903d8dcd27
# ╠═e5cf45cf-798a-4c44-8536-fc3974582f34
# ╠═4cb8a35c-727c-4ead-94d1-ebb890f40b14
# ╠═05476b2a-fbb1-4572-ae23-c1116dbbb386
