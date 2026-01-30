using Documenter, Mermaid, CommonSolve
using DocumenterInterLinks
# Load extensions to be documented
using Agents, OrdinaryDiffEq, MethodOfLines, Surrogates, Flux

links = InterLinks(
    "CommonSolve" => "https://docs.sciml.ai/CommonSolve/dev/",
    "DiffEq" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    "ModelingToolkit" => "https://docs.sciml.ai/ModelingToolkit/stable/",
    "Symbolics" => "https://docs.sciml.ai/Symbolics/stable/",
    "MethodOfLines" => "https://docs.sciml.ai/MethodOfLines/dev/",
    "PythonCall" => "https://juliapy.github.io/PythonCall.jl/stable/"
)

PAGES = [
    "Introduction" => "index.md",
    "Tutorial" => "tutorial.md",
    "Examples" => [
        "examples/mtk.md",
        "examples/duplicated_components.md",
        "examples/out_of_sync.md",
        "examples/spatial_maps.md",
        "examples/surrogates.md",
        "examples/external_components.md",
    ],
    "Mermaid Interface" => "interface.md",
    "Is Mermaid right for me?" => "is_mermaid_right_for_me.md",
    "API" => "API.md"
]

modules = [CommonSolve, Mermaid,
    Base.get_extension(Mermaid, :AgentsExt),
    Base.get_extension(Mermaid, :DiffEqExt),
    Base.get_extension(Mermaid, :MethodOfLinesExt),
    Base.get_extension(Mermaid, :SurrogatesExt),
]

format = Documenter.HTML(assets = ["assets/logo.ico"])

makedocs(sitename = "Mermaid.jl", format = format,
    repo = Remotes.GitHub("mjowen", "Mermaid.jl"), modules = modules, checkdocs = :exports,
    pages = PAGES, plugins = [links])

deploydocs(
    repo = "github.com/mjowen/Mermaid.jl",
)
