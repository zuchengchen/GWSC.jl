using Documenter, GWSC

makedocs(;
    modules=[GWSC],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/bingining/GWSC.jl/blob/{commit}{path}#L{line}",
    sitename="GWSC.jl",
    authors="Zu-Cheng Chen",
    assets=String[],
)
