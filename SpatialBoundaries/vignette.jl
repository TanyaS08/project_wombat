import Literate

corefiles = "main.jl"

vignetteconfig = Dict(
    "repo_root_url" => "https://github.com/TanyaS08/project_wombat",
    "codefence" => Pair("````julia", "````"),
    "flavor" => Literate.FranklinFlavor(),
    "credit" => false
)

Literate.markdown("main.jl", "vignettes"; config = vignetteconfig)

# write README

README = readlines("_README.md")
append!(README, readlines(joinpath("vignettes", "main.md")))

open("README.md", "w") do readme
    for line in README
        println(readme, line)
    end 
end