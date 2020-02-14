using Documenter, GWSC

makedocs(
    doctest = false,
    sitename = "GWSC.jl",
    modules = [GWSC],
    pages = ["Home" => "index.md",],
    repo = "https://github.com/bingining/GWSC.jl/blob/{commit}{path}#L{line}",
    authors = "Zu-Cheng Chen",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
)

deploydocs(repo = "github.com/bingining/GWSC.jl.git", push_preview = true)
