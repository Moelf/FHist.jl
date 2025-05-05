using FHist, BenchmarkTools

const uniform_xs = rand(10^6)
const gauss_xs = randn(10^6)
const weights = rand(10^6)

function hist1d_loop(xs; binedges)
    hist = Hist1D(; binedges)
    for x in xs
        push!(hist, x)
    end
    hist
end

function hist1d_loop_atomic(xs; binedges)
    hist = Hist1D(; binedges)
    for x in xs
        atomic_push!(hist, x)
    end
    hist
end

suite = BenchmarkGroup()

suite["Uniform Binning"] = BenchmarkGroup()
suite["Non-uniform Binning"] = BenchmarkGroup()

suite["Uniform Binning"]["One-shot Uniform input"] = @benchmarkable Hist1D(uniform_xs; binedges = 0:0.1:1)
suite["Uniform Binning"]["One-shot Gaussian input"] = @benchmarkable Hist1D(gauss_xs; binedges = 0:0.1:1)
suite["Uniform Binning"]["One-shot Uniform input with weights"] = @benchmarkable Hist1D(uniform_xs; weights, binedges = 0:0.1:1)
suite["Uniform Binning"]["One-shot Gaussian input with weights"] = @benchmarkable Hist1D(gauss_xs; weights, binedges = 0:0.1:1)
suite["Uniform Binning"]["push!-loop Uniform input"] = @benchmarkable hist1d_loop(uniform_xs; binedges = 0:0.1:1)
suite["Uniform Binning"]["atomic_push!-loop Gaussian input"] = @benchmarkable hist1d_loop_atomic(gauss_xs; binedges = 0:0.1:1)

suite["Non-Uniform Binning"]["One-shot Uniform input"] = @benchmarkable Hist1D(uniform_xs; binedges = [0, 0.1, 0.3, 0.5, 1.0])
suite["Non-Uniform Binning"]["One-shot Gaussian input"] = @benchmarkable Hist1D(gauss_xs; binedges = [0, 0.1, 0.3, 0.5, 1.0])
suite["Non-Uniform Binning"]["One-shot Uniform input with weights"] = @benchmarkable Hist1D(uniform_xs; weights, binedges = [0, 0.1, 0.3, 0.5, 1.0])
suite["Non-Uniform Binning"]["One-shot Gaussian input with weights"] = @benchmarkable Hist1D(gauss_xs; weights, binedges = [0, 0.1, 0.3, 0.5, 1.0])

tune!(suite)

run(suite; verbose=true, seconds = 2)
