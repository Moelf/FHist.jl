module SafeLog

const SAFE_LOG = Ref(false)

"""
    set_safe_log(x::Bool)

Turn on/off SAFE_LOG, to plot histograms without/with negative values. By default is `false`.
"""
function set_safe_log(x::Bool)
    SAFE_LOG[]=x
end


"""
    _clip_for_log(c_vec)

Internal fuction to clip counts of histogram. 
"""
function _clip_counts!(c_vec)
    min_positive = eps()
    @. c_vec = max(c_vec, min_positive)
end


"""
    _clip_counts!(c_vec,el_vec, eh_vec)

Internal fuction to clip counts and errors of histogram. 
"""
function _clip_counts!(c_vec, el_vec, eh_vec)
    
    # Set the clipping, and make copy of starting counts
    min_positive = eps()
    c_vec_def = copy(c_vec)

    # Clip bin counts
    @. c_vec = max(c_vec, min_positive)

    # clip lower errors
    mask =  c_vec - el_vec .< min_positive
    el_vec[mask] = c_vec[mask] .- min_positive

    # clip higher errors
    @. eh_vec = max(eh_vec - (c_vec - c_vec_def), min_positive)
end
 
end
