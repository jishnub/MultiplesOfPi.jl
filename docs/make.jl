using Documenter
using MultiplesOfPi

DocMeta.setdocmeta!(MultiplesOfPi, :DocTestSetup, :(using MultiplesOfPi); recursive=true)

makedocs(;
    modules=[MultiplesOfPi],
    authors="Jishnu Bhattacharya",
    repo="https://github.com/jishnub/MultiplesOfPi.jl/blob/{commit}{path}#L{line}",
    sitename="MultiplesOfPi.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jishnub.github.io/MultiplesOfPi.jl",
        assets=String[],
    ),
    pages=[
        "Reference" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jishnub/MultiplesOfPi.jl",
)
