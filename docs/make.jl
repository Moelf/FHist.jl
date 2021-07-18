using Documenter, FHist

makedocs(;
    modules=[FHist],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        # assets = ["assets/logo.ico"],
    ),
    pages=[
        "Introduction" => "index.md",
    ],
    repo="https://github.com/Moelf/FHist.jl/blob/{commit}{path}#L{line}",
    sitename="FHist.jl",
    authors="Jerry Ling",
    assets=String[],
)

deploydocs(;
    repo="github.com/Moelf/FHist.jl",
)
