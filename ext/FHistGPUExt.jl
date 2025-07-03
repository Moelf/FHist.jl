module FHistGPUExt

using FHist
import FHist: Hist1D, gpu_bincounts

using KernelAbstractions as KA
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace, @Const

@kernel unsafe_indices = true function histogram_naive_kernel!(histogram_output, @Const(output_L), @Const(input_raw), @Const(weights), @Const(firstr), @Const(invstep))
    gid = @index(Group, Linear)
    lid = @index(Local, Linear)

    gs = Int32(prod(@groupsize()))
    tid = (gid - Int32(1)) * gs + lid
    if tid <= length(input_raw)
        x = input_raw[tid]
        cursor = floor(Int32, (x - firstr) * invstep)
        bin = cursor + Int32(1)
        @atomic histogram_output[bin] += isnothing(weights) ? one(Float32) : weights[tid]
    end
end

@kernel unsafe_indices = true function histogram_sharemem_v2_kernel!(histogram_output, ::Val{N},
    @Const(input_raw), @Const(weights), @Const(firstr), @Const(invstep)
) where {N}
    gid = @index(Group, Linear)
    lid = @index(Local, Linear)

    @uniform gs = Int32(prod(@groupsize()))
    tid = (gid - Int32(1)) * gs + lid
    max_oid = Int32(length(histogram_output))

    shared_histogram = @localmem eltype(histogram_output) (N)

    # Setting shared_histogram to 0
    min_element = Int32(1)
    while min_element < N
        oid = min_element + lid - Int32(1)
        if oid <= max_oid
            @inbounds shared_histogram[oid] = Int32(0)
        end
        min_element += gs
    end
    @synchronize()

    # Defining bin on shared memory and writing to it if possible
    if tid <= length(input_raw)
        x = input_raw[tid]
        cursor = floor(Int32, (x - firstr) * invstep)
        bin = cursor + Int32(1)
        if bin >= Int32(1) && bin <= N
            @atomic shared_histogram[bin] += isnothing(weights) ? one(Float32) : weights[tid]
        end
    end
    @synchronize()

    min_element = Int32(1)
    while min_element < N
        oid = min_element + lid - Int32(1)
        if oid <= max_oid
            @atomic histogram_output[oid] += shared_histogram[oid]
        end

        min_element += gs
    end
end

# `data` already on gpu
function gpu_bincounts(data; weights=nothing, sync=false, binedges::AbstractRange, blocksize=512, backend)
    kernel! = histogram_sharemem_v2_kernel!(backend, (blocksize,))

    if !isnothing(weights)
        @assert get_backend(weights) == backend "Weights must be on the same backend as histogram_output"
    end

    cu_bincounts = KA.zeros(backend, Float32, length(binedges) - 1)
    firstr = Float32(first(binedges))
    invstep = Float32(inv(step(binedges)))

    kernel!(cu_bincounts, Val(length(cu_bincounts)), data, weights, firstr, invstep, ndrange=size(data))
    if sync
        synchronize(backend)
    end
    return cu_bincounts
end

end
