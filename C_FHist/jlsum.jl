Base.@ccallable function jlsum(p::Ptr{Cdouble}, N::Clong)::Cdouble
    np_input = unsafe_wrap(Array{Float64}, p, N, own=false)
    return sum(np_input)
end
