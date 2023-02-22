function _svg(h::Hist1D)
    paddingx, paddingy = 0.05, 0.10
    framewidth, frameheight = 250, 200
    fill = "#ffffff00" # 0 alpha
    strokecolor, strokewidth = "black", 1
    strokewidth = 1
    _c, _e = copy(bincounts(h)), copy(binedges(h))

    hasneg = minimum(_c) < 0
    zerolineyfrac = 0.0
    if hasneg
        zerolineyfrac = -minimum(_c)/(maximum(_c)-minimum(_c))
        _c .-= minimum(_c)
    end
    ys = frameheight * ((2*paddingy-1)/maximum(_c) .* _c .+ (1-paddingy))
    xs = framewidth * (
        (1 - 2 * paddingx) / (maximum(_e) - minimum(_e))
        .* (_e .- minimum(_e))
        .+ paddingx
    )
    points = [(paddingx*framewidth, (1-paddingy)*frameheight)] # bottom left
    for i in 1:(length(xs)-1)
        push!(points, (xs[i],ys[i])) # left bin edge
        push!(points, (xs[i+1],ys[i])) # right bin edge
    end
    push!(points, ((1-paddingx)*framewidth,(1-paddingy)*frameheight))
    if !hasneg
        push!(points, points[1]) # close path
    end
    pathstr = join(["$(x),$(y)" for (x,y) in points],",")
    zeroliney = hasneg ? frameheight*(1-zerolineyfrac) : frameheight*(1-paddingy)
    xlow, xhigh = round.(extrema(_e), sigdigits=10)
    return """
    <svg width="$(framewidth)" height="$(frameheight)" version="1.1" xmlns="http://www.w3.org/2000/svg">
        <polyline points="$(pathstr)" stroke="$(strokecolor)" fill="$(fill)" stroke-width="$(strokewidth)"/>
        <polyline points="$(framewidth*paddingx),$(zeroliney),$(framewidth*(1-paddingx)),$(zeroliney)" stroke="black" stroke-width="1"/>
        <text x="$(framewidth*paddingx)" y="$(frameheight*(1-0.5*paddingy))" dominant-baseline="middle" text-anchor="start" fill="black">$(xlow)</text>
        <text x="$(framewidth*(1-paddingx))" y="$(frameheight*(1-0.5*paddingy))" dominant-baseline="middle" text-anchor="end" fill="black">$(xhigh)</text>
    </svg>
    """
end

function _svg(h::Hist2D)
    paddingx1, paddingy1 = 0.15, 0.15 # left, bottom
    paddingx2, paddingy2 = 0.05, 0.05 # right, top
    framewidth, frameheight = 250, 200

    (xlow, xhigh), (ylow, yhigh) = extrema.(binedges(h))
    function transform(x::Real, y::Real)
        xfrac = (x - xlow)/(xhigh - xlow)
        xdraw = xfrac * framewidth*(1-paddingx1-paddingx2) + framewidth*paddingx1
        yfrac = (y - ylow)/(yhigh - ylow)
        ydraw = frameheight - (yfrac * frameheight*(1-paddingy1-paddingy2) + frameheight*paddingy1)
        return xdraw, ydraw
    end

    function colorscale(x::Real)
        x = isnan(x) ? 0. : x
        # x is fraction between [0,1]
        # 6th deg polynomial fit to viridis color palette
        # See https://gist.github.com/aminnj/dba84777718613d3a37291a43659feff
        # to generate other color palettes. Black and white is trivial:
        #   return 255 .* (1-x,1-x,1-x)
        r = @evalpoly(x, 0.2731, 0.1270, -0.3617, -4.7456, 6.7092, 4.2491, -5.2678)
        g = @evalpoly(x, 0.0039, 1.3810, 0.3969, -6.4246, 15.3241, -14.7345, 4.9587)
        b = @evalpoly(x, 0.3305, 1.3726, 0.3948, -20.7431, 59.7287, -68.3565, 27.4051)
        return round.(Int, 255 .* (r,g,b))
    end

    counts = bincounts(h)
    mincount, maxcount = extrema(counts)
    counts = clamp.(counts, 0, maxcount)
    ex, ey = binedges(h)
    sx, sy = size(counts)
    rectlines = String[]
    textstyle = """ dominant-baseline="middle" text-anchor="middle" font-size="85%" font-family="sans-serif" """
    for i in 1:sx, j in 1:sy
        tx1, ty1 = ceil.(Int, transform(ex[i], ey[j+1]))
        tx2, ty2 = ceil.(Int, transform(ex[i+1], ey[j]))
        tw, th = tx2-tx1, ty2-ty1
        c = counts[i,j]
        (c == 0) && continue
        r,g,b = colorscale((c-mincount)/(maxcount-mincount))
        rcolor = "rgb($(r),$(g),$(b))"
        line = """<rect x="$(tx1)" y="$(ty1)" width="$(tw)" height="$(th)" fill="$(rcolor)" stroke="none" />"""
        if (sx <= 15) && (sy <= 15)
            tcolor = ((0.299*r + 0.587*g + 0.114*b) < 0.60*255) ? "#fff" : "#000"
            line = """<g>$(line)<text class="svgplotlabels" x="$(tx1+tw/2)" y="$(ty1+th/2)" fill="$(tcolor)" $(textstyle)>$(c)</text></g>"""
        end
        push!(rectlines, line * "\n")
    end

    return """
    <svg width="$(framewidth)" height="$(frameheight)" version="1.1" xmlns="http://www.w3.org/2000/svg" class="svgplot">
        <style>
          .svgplotlabels { display: none; }
          .svgplot g:hover text { display: block; cursor: default; }
        </style>
        $(join(rectlines))
        <rect x="$(round(Int,framewidth*paddingx1))" y="$(round(Int,frameheight*paddingy2))" width="$(round(Int,framewidth*(1-paddingx1-paddingx2)))" height="$(frameheight*(1-paddingy1-paddingy2))" fill="none" stroke="#000" stroke-width="1" />
        <text x="$(framewidth*paddingx1)" y="$(frameheight*(1-0.5*paddingy1))" fill="black" $(textstyle)>$(minimum(ex))</text>
        <text x="$(framewidth*(1-paddingx2))" y="$(frameheight*(1-0.5*paddingy1))" fill="black" $(textstyle)>$(maximum(ex))</text>
        <text x="$(framewidth*0.5*paddingx1)" y="$(frameheight*(1-paddingy1))" fill="black" $(textstyle)>$(minimum(ey))</text>
        <text x="$(framewidth*0.5*paddingx1)" y="$(frameheight*paddingy2)" fill="black" $(textstyle)>$(maximum(ey))</text>
    </svg>
    """
end


function _svg(h::Hist3D)
    framewidth = 250
    frameheight = 200

    x1, x2 = 60, 180
    y1, y2 = 30, 130
    z1, z2 = 10, 120

    # cabinet projection
    α = pi/4
    P = [1 0 1/2*cos(α) ; 0 1 1/2*sin(α); 0 0 0]
    projection(points) = map(p->begin
                                 x,y,z = (P*p)
                                 [round(x),round(frameheight-y)]
                             end, points)
    flatstring(points) = join([points...;], " ")

    frontface = projection([[x1,y1,z1],[x2,y1,z1],[x2,y2,z1],[x1,y2,z1],[x1,y1,z1]])
    rightface = projection([[x2,y1,z2],[x2,y2,z2],[x2,y2,z1],[x2,y1,z1],[x2,y1,z2]])
    topface = projection([[x1,y2,z2],[x1,y2,z1],[x2,y2,z1],[x2,y2,z2],[x1,y2,z2]])
    backface = projection([[x1,y1,z2],[x2,y1,z2],[x2,y2,z2],[x1,y2,z2],[x1,y1,z2]])
    bottomface = projection([[x1,y1,z2],[x1,y1,z1],[x2,y1,z1],[x2,y1,z2],[x1,y1,z2]])
    xaxisline = projection([[x1,y1,z1],[x2,y1,z1]])
    yaxisline = projection([[x1,y1,z1],[x1,y2,z1]])
    zaxisline = projection([[x1,y2,z1],[x1,y2,z2]])

    (xlow, xhigh), (ylow, yhigh), (zlow, zhigh) = extrema.(binedges(h))
    function data_to_pixels(x, y, z)
        xt = (x - xlow)/(xhigh - xlow)*(x2-x1)+x1
        yt = (y - ylow)/(yhigh - ylow)*(y2-y1)+y1
        zt = (z - zlow)/(zhigh - zlow)*(z2-z1)+z1
        return xt, yt, zt
    end

    boxes = ""
    if prod(nbins(h)) < 5000
        counts = bincounts(h)
        maxcount = maximum(counts)
        counts = clamp.(counts, 0, maxcount)
        ex, ey, ez = binedges(h)
        sx, sy, sz = size(counts)
        lines = String[]
        for i in 1:sx, j in 1:sy, k in 1:sz
            c = counts[i,j,k]
            (c == 0) && continue
            x1t,y1t,z1t = data_to_pixels(ex[i], ey[j], ez[k])
            x2t,y2t,z2t = data_to_pixels(ex[i+1], ey[j+1], ez[k+1])
            sf = c/maxcount # normalize box size by count relative to maximum
            zf = 1 - (0.5*(ez[k]+ez[k+1])-zlow) / (zhigh-zlow)
            opacity = 0.3 + 0.3*zf # more transparent as you go toward the rear
            dx, dy, dz = x2t - x1t, y2t - y1t, z2t - z1t
            x1t += dx*(1-sf)/2
            x2t -= dx*(1-sf)/2
            y1t += dy*(1-sf)/2
            y2t -= dy*(1-sf)/2
            z1t += dz*(1-sf)/2
            z2t -= dz*(1-sf)/2
            face1 = projection([[x1t,y1t,z1t],[x2t,y1t,z1t],[x2t,y2t,z1t],[x1t,y2t,z1t],[x1t,y1t,z1t]])
            face2 = projection([[x2t,y1t,z2t],[x2t,y2t,z2t],[x2t,y2t,z1t],[x2t,y1t,z1t],[x2t,y1t,z2t]])
            face3 = projection([[x1t,y2t,z2t],[x1t,y2t,z1t],[x2t,y2t,z1t],[x2t,y2t,z2t],[x1t,y2t,z2t]])
            push!(lines, """<polyline points="$(flatstring(face1))" stroke="none" fill="#888" opacity="$(opacity)"/>\n""")
            push!(lines, """<polyline points="$(flatstring(face2))" stroke="none" fill="#ccc" opacity="$(opacity)"/>\n""")
            push!(lines, """<polyline points="$(flatstring(face3))" stroke="none" fill="#aaa" opacity="$(opacity)"/>\n""")
        end
        boxes = join(lines)
    end

    c1, c2, c3 = "#1f77b4", "#d62728", "#2ca02c"
    (x1text, x2text), (y1text, y2text), (z1text, z2text) = map(x->(first(x),last(x)), binedges(h))
    textstyle = """ font-size="85%" font-family="sans-serif" """
    return """
    <svg width="$(framewidth)" height="$(frameheight)" version="1.1" xmlns="http://www.w3.org/2000/svg">
        <polyline points="$(flatstring(backface))" stroke="gray" stroke-width="1" stroke-dasharray="4" fill="none"/>
        <polyline points="$(flatstring(bottomface))" stroke="gray" stroke-width="1" stroke-dasharray="4" fill="none"/>
        $(boxes)
        <polyline points="$(flatstring(frontface))" stroke="black" stroke-width="1" fill="none"/>
        <polyline points="$(flatstring(rightface))" stroke="black" stroke-width="1" fill="none"/>
        <polyline points="$(flatstring(topface))" stroke="black" stroke-width="1" fill="none"/>

        <polyline points="$(flatstring(xaxisline))" stroke="$(c1)" stroke-width="2" fill="none"/>
        <polyline points="$(flatstring(yaxisline))" stroke="$(c2)" stroke-width="2" fill="none"/>
        <polyline points="$(flatstring(zaxisline))" stroke="$(c3)" stroke-width="2" fill="none"/>

        <text x="$(mean(xaxisline[1])+10)" y="$(xaxisline[1][2]+5)" fill="$(c1)" dominant-baseline="hanging" text-anchor="middle" $(textstyle)>x</text>
        <text x="$(xaxisline[1][1])" y="$(xaxisline[1][2]+5)" fill="$(c1)" dominant-baseline="hanging" text-anchor="left" $(textstyle)>$(x1text)</text>
        <text x="$(xaxisline[2][1])" y="$(xaxisline[2][2]+5)" fill="$(c1)" dominant-baseline="hanging" text-anchor="end" $(textstyle)>$(x2text)</text>

        <text x="$(yaxisline[1][1]-5)" y="$(mean(yaxisline[1]))" fill="$(c2)" dominant-baseline="text-bottom" text-anchor="end" $(textstyle)>y</text>
        <text x="$(yaxisline[1][1]-5)" y="$(yaxisline[1][2])" fill="$(c2)" dominant-baseline="text-bottom" text-anchor="end" $(textstyle)>$(y1text)</text>
        <text x="$(yaxisline[2][1]-5)" y="$(yaxisline[2][2])" fill="$(c2)" dominant-baseline="hanging" text-anchor="end" $(textstyle)>$(y2text)</text>

        <text x="$((zaxisline[1][1]+zaxisline[2][1])/2-10)" y="$((zaxisline[1][2]+zaxisline[2][2])/2)" fill="$(c3)" dominant-baseline="text-bottom" text-anchor="end" $(textstyle)>z</text>
        <text x="$(zaxisline[1][1]+15)" y="$(zaxisline[1][2]-4)" fill="$(c3)" dominant-baseline="text-bottom" text-anchor="start" $(textstyle)>$(z1text)</text>
        <text x="$(zaxisline[2][1]+5)" y="$(zaxisline[2][2]+3)" fill="$(c3)" dominant-baseline="hanging" text-anchor="start" $(textstyle)>$(z2text)</text>
    </svg>
         """
end

function Base.show(io::IO, h::Hist1D)
    compact = get(io, :compact, false)
    if compact
        print(io, "Hist1D{$(eltype(bincounts(h)))}, ")
        print(io, "edges=$(repr(binedges(h), context=:limit => true)), ")
        print(io, "integral=$(integral(h))")
    else
        println(io, "edges: ", binedges(h))
        println(io, "bin counts: ", bincounts(h))
        print(io, "total count: ", integral(h))
    end
end

function Base.show(io::IO, m::MIME"text/html", h::Hist1D)
    println(io, """
    <div style="display: flex;">
        <div style="float:left; margin:5px">$(_svg(h))</div>
        <div style="float:left; margin:5px; max-width: 50%; display:flex; justify-content:center; align-items:center;">
            <ul>
                <li>edges: $(repr(binedges(h), context=:limit => true))</li>
                <li>bin counts: $(repr(bincounts(h), context=:limit => true))</li>
                <li>maximum count: $(maximum(bincounts(h)))</li>
                <li>total count: $(integral(h))</li>
            </ul>
        </div>
    </div>
    """)
end

function Base.show(io::IO, h::Hist2D)
    compact = get(io, :compact, false)
    if compact
        print(io, "Hist2D{$(eltype(bincounts(h)))}, ")
        print(io, "edges=$(repr(binedges(h), context=:limit => true)), ")
        print(io, "integral=$(integral(h))")
    else
        println(io, "edges: ", binedges(h))
        println(io, "bin counts: ", bincounts(h))
        print(io, "total count: ", integral(h))
    end
end

function Base.show(io::IO, m::MIME"text/html", h::Hist2D)
    println(io, """
    <div style="display: flex;">
        <div style="float:left; margin:5px">$(_svg(h))</div>
        <div style="float:left; margin:5px; max-width: 50%; display:flex; justify-content:center; align-items:center;">
            <ul>
                <li>edges: $(repr(binedges(h), context=:limit => true))</li>
                <li>bin counts: $(repr(bincounts(h), context=:limit => true))</li>
                <li>maximum count: $(maximum(bincounts(h)))</li>
                <li>total count: $(integral(h))</li>
            </ul>
        </div>
    </div>
    """)
end

function Base.show(io::IO, h::Hist3D)
    print(io, "Hist3D{$(eltype(bincounts(h)))}, ")
    print(io, "edges=$(repr(binedges(h), context=:limit => true)), ")
    print(io, "integral=$(integral(h))")
end

function Base.show(io::IO, m::MIME"text/html", h::Hist3D)
    println(io, """
    <div style="display: flex;">
        <div style="float:left; margin:5px">$(_svg(h))</div>
        <div style="float:left; margin:5px; max-width: 50%; display:flex; justify-content:center; align-items:center;">
            <ul>
                <li>edges: $(repr(binedges(h), context=:limit => true))</li>
                <li>bin counts: $(repr(bincounts(h), context=:limit => true))</li>
                <li>maximum count: $(maximum(bincounts(h)))</li>
                <li>total count: $(integral(h))</li>
            </ul>
        </div>
    </div>
    """)
end
