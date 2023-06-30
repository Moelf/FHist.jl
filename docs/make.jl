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
    ),
    pages=[
        "Introduction" => "index.md",
        "APIs" => "api.md",
        "Writing to `.root`" => "writingtoroot.md",
        "Tutorials" => T,
    ],
    repo="https://github.com/Moelf/FHist.jl/blob/{commit}{path}#L{line}",
    sitename="FHist.jl",
    authors="Jerry Ling",
    assets=String[],
)

deploydocs(;
    repo="github.com/Moelf/FHist.jl",
    push_preview = true
)
