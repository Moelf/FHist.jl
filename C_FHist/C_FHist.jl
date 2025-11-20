module C_FHist

using FHist: Hist1D, _fast_bincounts_1d!

Base.@ccallable function hist1d(input::Ptr{Cdouble}, Ninput::Clong, bincounts::Ptr{Cdouble}, Nbincounts::Clong, start::Cdouble, step::Cdouble, stop::Cdouble)::Cvoid
    np_input = unsafe_wrap(Array{Float64}, input, Ninput, own=false)
    np_bincounts = unsafe_wrap(Array{Float64}, bincounts, Nbincounts, own=false)
    binedges = start:step:stop
    h = Hist1D(; bincounts=np_bincounts, binedges)
    _fast_bincounts_1d!(h, np_input, binedges)
    return nothing
end

end
