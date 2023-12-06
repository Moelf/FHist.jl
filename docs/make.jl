using Documenter, FHist
using Pluto: Configuration.CompilerOptions
using PlutoStaticHTML

notebooks = [
    "Makie Plotting",
]

include("build.jl")

build()
md_files = markdown_files()
T = [t => f for (t, f) in zip(notebooks, md_files)]

makedocs(;
    modules=[FHist],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = 5*10^6
        # assets = ["assets/logo.ico"],
    ),
    pages=[
        "Introduction" => "index.md",
        "APIs" => "api.md",
        "Writing to `.root`" => "writingtoroot.md",
        "Tutorials" => T,
    ],
    sitename="FHist.jl",
    authors="Jerry Ling"
)

deploydocs(;
    repo="github.com/Moelf/FHist.jl",
    push_preview = true
)
