### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 4cf8f502-eba2-11ec-275f-3f1858e4c2e6
# hideall
let
    docs_dir = dirname(dirname(@__DIR__))
    pkg_dir = dirname(docs_dir)

    using Pkg: Pkg
    Pkg.activate(docs_dir)
    Pkg.develop(; path=pkg_dir)
    Pkg.instantiate()
end;

# ╔═╡ 82a6fb18-e9c9-450d-9c93-e05bd6cf5859
using FHist, CairoMakie

# ╔═╡ 180e65b9-13db-4b0c-8b16-369b6978cec0
md"""
Let's generate three dummy histograms sampled from three different distributions:
1. unit Gaussian with 1000 points
2. unit Gaussian with 10000 points
3. unit Gaussians with 1000 points and a mean of 0.5
"""

# ╔═╡ 8491f5c3-93b6-4670-8f6e-0d08d7afbf75
begin
	h1 = Hist1D(randn(10^3); binedges = -2:0.3:2)
	h2 = Hist1D(randn(10^4); binedges = -2:0.3:2)
	h3 = Hist1D(randn(10^3) .+ 0.5; binedges = -2:0.3:2)
end

# ╔═╡ da61f742-b530-4aa4-9253-a23f6731f892
md"""
To plot a single histogram, you can do any of the following:
"""

# ╔═╡ 4c09d3a5-a034-49fb-aba3-1891c6f90b5a
begin
	fig = Figure()
	p1, ax1= plot(fig[1, 1],h1; label = "plot(h1)")
	p2, ax2 = stephist(fig[1, 2], h1, label = "stairs(h1)")
	p3, ax3 = hist(fig[2, 1], h1, label = "hist(h1)")
	p4, ax4 = errorbars(fig[2, 2], h1, color=:black, whiskerwidth=6, label = "errorbars(h1)")
	axislegend.([p1, p2, p3, p4])
	fig
end

# ╔═╡ 91a8b01a-54e7-4955-bb56-1803d3f2576a
begin
	h1_log = Hist1D(rand(10^5);binedges = [0.001, 0.01, 0.1, 1])
	ax_log = plot(h1_log; axis=(;xscale=log10, yscale=log10), label = "log-log scale")
	xlims!(0.0001, nothing)
	ylims!(1, nothing)
	axislegend()
	ax_log
end

# ╔═╡ 3b4d5b69-dbbb-49dd-b7ad-9ef1605f76eb
md"""
When you're studying a single histogram, it's hepful to have statbox (CERN ROOT enjoyer?)

We provide a `statbox!()` function to add such box:
"""

# ╔═╡ 99f8bcde-b823-4c5b-8820-22c403c21d6b
begin
	f_s = stephist(h1)
	statbox!(f_s, h1)
	f_s
end

# ╔═╡ c774440d-c7c1-4b75-b115-64038992174d
md"""
`statbox` also works for 2D histogram as one'd expect:
"""

# ╔═╡ f5630dae-870e-45e8-995a-db946687031b
begin
	fig2d = Figure()
	h2d = Hist2D((randn(10000), randn(10000)))
	plot(fig2d[1,2], h2d)
	statbox!(fig2d, h2d; position=(1,1))
	Colorbar(fig2d[1,3])
	fig2d
end

# ╔═╡ 254c736a-79e3-4cf5-a662-77471c5aa465
begin
	even_edges = 0:0.1:1
	uneven_edges = [0, 0.2, 0.5, 1]
	data = (randn(10000), randn(10000))
	h2d_even =   Hist2D(data; binedges = (even_edges, even_edges))
	h2d_uneven = Hist2D(data; binedges = (uneven_edges, even_edges))
	
	f2_uneven = Figure()
	plot(f2_uneven[2,1], h2d_even; axis=(title="even x-binning", ))
	plot(f2_uneven[1,1], h2d_uneven; axis=(title="uneven x-binning", ))
	Colorbar(f2_uneven[:,2])
	f2_uneven
end

# ╔═╡ d549a735-12a4-443a-98c7-be5cce4e6789
md"""
You can freely combine these by using the bang `!` version to plot things on top of each other, for example this is some plot you might see:
"""

# ╔═╡ 12b44438-2af4-4ca5-8569-57de2cadc607
begin
    pos_h2 = restrict(h2, 0, Inf)
    hist(pos_h2;
         label = "fake bkg", 
         axis=(
		   limits=(nothing, nothing, 0.1, nothing),
               yminorticks=IntervalsBetween(5),
               yminorticksvisible = true,
               yscale=log10
              )
        )

    errorbars!(pos_h2; whiskerwidth=7)
    stephist!(h3; color=:red, label = "fake sig")
    axislegend()
    current_figure()
end

# ╔═╡ 386b8f0e-ad21-4af3-b2b6-b3412f44315c
md"""
But it's also pretty common to have multiple background components (processes), in such case, we often want a stacked histogram with error bar representing the sum of each component added up correctly.

We also have a recipe for that:
"""

# ╔═╡ 70a5ba34-53c5-4e12-93b9-49eace436524
with_theme(ATLASTHEME) do
	f, a, p = stackedhist([h2, h2, h1]; )
	labels = ["proc1", "blah", "third one"]
	elements = [PolyElement(polycolor = p.attributes.color[][i]) for i in 1:length(labels)]
	title = "Legend title"
	Legend(f[1,2], elements, labels, title)
	f
end

# ╔═╡ 089fb6a1-3f2a-45f4-b4bf-4013dee3a6da
md"""
Notice the `with_theme(ATLASTHEME)`, you can apply a bunch of appearance tweaks for the entire plot.
"""

# ╔═╡ a3dd8f72-d292-4f85-be76-2647b9433b42
md"""
Finally, we have a recipe for ratio plot which is argubaly useful for the bottom pannel in many plots:
"""

# ╔═╡ b1188446-d16b-4e0c-93be-c5e3d0de379c
begin
    f_ratio, a, p = stephist(h1; label="h1")
	errorbars!(a, h1)
    stephist!(a, h3; label="h3")
	errorbars!(a, h3)
	
    axislegend()
    ratioax = Axis(f_ratio[2, 1], aspect = 5.5, xlabel = "my x-axis", ylabel="h1 / h3")
    linkxaxes!(ratioax, a)
    hidexdecorations!(a)
    ratiohist!(ratioax, h1/h3; color=:black, errors=true)
    rowsize!(f_ratio.layout, 1, Aspect(1, 0.5))
    f_ratio
end

# ╔═╡ 2f97d05b-9103-451e-8fce-bb7458acf2ba
md"""
# Shading/Hatching errorbar band

"""

# ╔═╡ b114c8a7-5f01-44fe-9660-b5021d359399
let
    f, a, p = stackedhist([h1, h1 ]; error_color=(:black, 0.5))
	labels = ["proc1", "blah"]
	elements = [PolyElement(polycolor = p.attributes.color[][i]) for i in 1:length(labels)]
	title = "Legend title"
	Legend(f[1,2], elements, labels, title)
	f
end

# ╔═╡ fae56f0f-bac4-4dea-b33c-d45e6f3b51fc
let
    f, a, p = stackedhist([h1, h1 ]; error_color=Pattern('/'))
	f
end

# ╔═╡ Cell order:
# ╠═4cf8f502-eba2-11ec-275f-3f1858e4c2e6
# ╠═82a6fb18-e9c9-450d-9c93-e05bd6cf5859
# ╠═180e65b9-13db-4b0c-8b16-369b6978cec0
# ╠═8491f5c3-93b6-4670-8f6e-0d08d7afbf75
# ╠═da61f742-b530-4aa4-9253-a23f6731f892
# ╠═4c09d3a5-a034-49fb-aba3-1891c6f90b5a
# ╠═91a8b01a-54e7-4955-bb56-1803d3f2576a
# ╠═3b4d5b69-dbbb-49dd-b7ad-9ef1605f76eb
# ╠═99f8bcde-b823-4c5b-8820-22c403c21d6b
# ╠═c774440d-c7c1-4b75-b115-64038992174d
# ╠═f5630dae-870e-45e8-995a-db946687031b
# ╠═254c736a-79e3-4cf5-a662-77471c5aa465
# ╠═d549a735-12a4-443a-98c7-be5cce4e6789
# ╠═12b44438-2af4-4ca5-8569-57de2cadc607
# ╟─386b8f0e-ad21-4af3-b2b6-b3412f44315c
# ╠═70a5ba34-53c5-4e12-93b9-49eace436524
# ╠═089fb6a1-3f2a-45f4-b4bf-4013dee3a6da
# ╠═a3dd8f72-d292-4f85-be76-2647b9433b42
# ╠═b1188446-d16b-4e0c-93be-c5e3d0de379c
# ╟─2f97d05b-9103-451e-8fce-bb7458acf2ba
# ╠═b114c8a7-5f01-44fe-9660-b5021d359399
# ╠═fae56f0f-bac4-4dea-b33c-d45e6f3b51fc
