module PlottingCairoMakieExt
#using .CairoMakie

isdefined(Base, :get_extension) ? (using CairoMakie) : (using ..CairoMakie)

#FIXME remove piracy
function CairoMakie.Cairo.CairoPattern(color::Makie.AbstractPattern)
    # the Cairo y-coordinate are fliped
    color_matrix = replace!(Makie.to_image(color), Makie.RGBA{Float32}(1.0f0,1.0f0,1.0f0,0.0f0) => Makie.RGBA{Float32}(0.0f0,0.0f0,0.0f0,0.0f0))
    bitmappattern = reverse!(Makie.ARGB32.(color_matrix); dims=2)
    cairoimage = CairoMakie.Cairo.CairoImageSurface(bitmappattern)
    cairopattern = CairoMakie.Cairo.CairoPattern(cairoimage)
    return cairopattern
end
end

