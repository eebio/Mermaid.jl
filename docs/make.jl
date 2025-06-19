using Documenter, Mermaid, CommonSolve
using DocumenterInterLinks

links = InterLinks(
    "CommonSolve" => "https://docs.sciml.ai/CommonSolve/dev/",
    "DiffEq" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    "ModelingToolkit" => "https://docs.sciml.ai/ModelingToolkit/stable/",
    "Symbolics" => "https://docs.sciml.ai/Symbolics/stable/",
)

PAGES = [
    "Introduction" => "index.md",
    "Tutorial" => "tutorial.md",
    "Examples" => ["examples/mtk.md"],
    "Mermaid Interface" => "interface.md",
    "Is Mermaid right for me?" => "is_mermaid_right_for_me.md",
    "API" => "API.md",
]

makedocs(sitename="Mermaid.jl", repo=Remotes.GitHub("mjowen", "Mermaid.jl"), modules=[CommonSolve, Mermaid], checkdocs = :exports,
    pages = PAGES, plugins = [links])

deploydocs(
    repo = "github.com/mjowen/Mermaid.jl",
)
