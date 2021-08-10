import Literate

corefiles = "main.jl"

vignetteconfig = Dict(
    "repo_root_url" => "https://github.com/TanyaS08/project_wombat",
    "codefence" => Pair("````julia", "````"),
    "flavor" => Literate.FranklinFlavor(),
    "credit" => false
)

Literate.markdown("main.jl", "vignettes"; config=vignetteconfig)
