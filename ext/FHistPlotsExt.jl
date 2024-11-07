module FHistPlotsExt
using FHist, RecipesBase, Statistics

isdefined(Base, :get_extension) ? (using Plots) : (using ..Plots)

#
@recipe function f(h::Hist1D)
    seriestype --> :barbins
    N = nentries(h)
    M = round(mean(h); sigdigits= 2)
    S = round(std(h); sigdigits= 2)
    label --> "Entries = $N\nMean = $M\nStd Dev = $S\nOverflow = $(h.overflow)"
    x:= binedges(h)
    y:= bincounts(h)
    ()
end

@recipe function f(h::Hist2D)
    seriestype --> :bins2d
    x := binedges(h)[1]
    y := binedges(h)[2]
    z := (surf = bincounts(h), )
    ()
end


"""
    WithData(factor, support)
    WithData(bins, N::Int = 1)
    WithData(h::Hist1D)

A struct that helps plotting a model curve on top of the binned data with correct normalization.
If the model is normalized to 1, the scaling factor accounts for the bin width and number of entries in the histogram.
For a model curve normalized to the total number of events, default option `N=1` is provided.

# Examples:
```julia
julia> data = randn(1000);
julia> h = Hist1D(data; binedges=range(-5,5, 100));
julia> model_fun(x) = length(data) * exp(-x^2 / 2) / sqrt(2π);
julia> #
julia> let
    plot(h, seriestype=:stepbins)
    plot!(model_fun, WithData(h.binedges[1]), lw=2)
end  # displays the plot
```
"""
struct WithData
    factor::Float64
    support::Tuple{Float64, Float64}
end
WithData(bins, N::Int = 1) = WithData(N * diff(bins)[1], (bins[1], bins[end]))
WithData(h::Hist1D) = WithData(h.binedges[1], Int(sum(h.bincounts)))
#
@recipe function f(model_fun::Function, scale::WithData; n_points = 100)
    normalized_function(x) = model_fun(x) * scale.factor
    xv = range(scale.support..., n_points)
    return (normalized_function, xv)
end


"""
    curvedfitwithpulls(h, best_model;, data_scale_curve = true)

Plots the histogram with Poisson uncertainties, the model by a curve, and a pull distribution.

# Example
```julia
julia> h = Hist1D(data; binedges=1.1:0.1:2.5))
julia> curvedfitwithpulls(h0, best_model, xlab = "X-axis", ylab = "Y-axis")
````

The recipe does something along the lines of the following code:
```julia
WD = data_scale_curve ? WithData(h) : WithData(h.binedges[1], 1)

# pulls
mismatches = h0.bincounts - best_model.(bincenters(h0)) .* WD.factor
hpull = Hist1D(; binedges, bincounts = mismatches ./ binerrors(h0))
dx = (step(binedges)) / 2

# plotting
plot(layout = grid(2, 1, heights = (0.8, 0.2)), link = :x)
plot!(sp = 1, best_model, WD, lw = 2, lc = :cornflowerblue)
scatter!(sp = 1, bincenters(h0), h0.bincounts, yerror = binerrors(h0), xerror = dx)
scatter!(sp = 2, bincenters(hpull), hpull.bincounts, yerror = 1, xerror = dx)
```
"""
@userplot CurvedFitWithPulls

@recipe function f(h::CurvedFitWithPulls; data_scale_curve = true)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: Hist1D) ||
       !(typeof(h.args[2]) <: Function)
        error("Marginal CurvedFitWithPulls should be given a histogram and a function.  Got: $(typeof(h.args))")
    end
    h, model = h.args
    binedges = h.binedges[1]
    WD = data_scale_curve ? WithData(h) : WithData(binedges, 1)

    mismatches = h.bincounts - best_model.(bincenters(h)) .* WD.factor
    hpull = Hist1D(; binedges, bincounts = mismatches ./ binerrors(h))
    dx = (binedges[2] - binedges[1]) / 2

    # set up the subplots
    legend := false
    link := :x
    framestyle := :box
    grid := false
    layout := grid(2, 1, heights = (0.8, 0.2))

    xlimits := (WD.support)

    # model curve
    @series begin
        subplot := 1
        linewidth := 1.5
        best_model, WD
    end
    markersize := 1.5

    # data
    @series begin
        xticks := false
        xguide := ""
        ylimits := (0, :auto)
        bottom_margin := -4.5mm
        markercolor --> :black
        seriestype := :scatter
        subplot := 1
        yerror := binerrors(h0)
        xerror := dx
        bincenters(h), h.bincounts
    end

    # pull
    @series begin
        linecolor := :gray
        subplot := 2
        seriestype := :hline
        [0.0]
    end
    @series begin
        bottom_margin := 0mm
        markercolor --> :black
        yguide := "Pull"
        guidefontvalign := :center
        seriestype := :scatter
        ylims := (-7, 7)
        subplot := 2
        yerror := 1.0
        xerror := dx
        bincenters(hpull), hpull.bincounts
    end
end

end

