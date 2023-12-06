### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ be0b45f2-86dc-11ee-2dce-2994fb540eec
# hideall
let
    docs_dir = dirname(dirname(@__DIR__))
    pkg_dir = dirname(docs_dir)

    using Pkg: Pkg
    Pkg.activate(docs_dir)
    Pkg.develop(; path=pkg_dir)
    Pkg.instantiate()
end;

# ╔═╡ 5a884f80-7c0a-4345-b2d3-1785b1b5a46a
using FHist, Statistics, Plots

# ╔═╡ d604a07e-2179-4570-b6bb-df3b8f6cfc7e
md"## Hist1D"

# ╔═╡ 6ebafb79-448b-425a-890a-39ae454816be
md"""
Let's generate three dummy histograms sampled from three different distributions:
1. unit Gaussian with 1000 points
2. unit Gaussian with 10000 points
3. unit Gaussians with 1000 points and a mean of 0.5
"""

# ╔═╡ 316167ab-4d59-4978-a3b0-539e23c1eaaa
begin
	h1 = Hist1D(randn(10^3), -2:0.3:2)
	h2 = Hist1D(randn(10^4), -2:0.3:2)
	h3 = Hist1D(randn(10^3) .+ 0.5, -2:0.3:2)
end

# ╔═╡ 2d857d08-f9c2-4a60-86dd-2afbd3b599f6
begin
	Plots.plot(h1; size=(600, 500), legend=:topleft)
	Plots.plot!(h3)
end

# ╔═╡ 68d705dd-d0d9-4394-8d3c-f8b2113e6cef
md"## Hist2D"

# ╔═╡ 53eacafc-71ac-4e70-a238-581335ac4729
begin
		h2d = Hist2D((randn(10000), randn(10000)))
		p2d = plot(h2d)
end

# ╔═╡ Cell order:
# ╠═be0b45f2-86dc-11ee-2dce-2994fb540eec
# ╠═5a884f80-7c0a-4345-b2d3-1785b1b5a46a
# ╟─d604a07e-2179-4570-b6bb-df3b8f6cfc7e
# ╟─6ebafb79-448b-425a-890a-39ae454816be
# ╠═316167ab-4d59-4978-a3b0-539e23c1eaaa
# ╠═2d857d08-f9c2-4a60-86dd-2afbd3b599f6
# ╟─68d705dd-d0d9-4394-8d3c-f8b2113e6cef
# ╠═53eacafc-71ac-4e70-a238-581335ac4729
