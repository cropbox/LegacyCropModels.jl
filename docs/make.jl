using Documenter
using LegacyCropModels

makedocs(
    sitename = "LegacyCropModels.jl",
    remotes = nothing,
    format = Documenter.HTML(
        edit_link = "main",
        repolink = "https://github.com/cropbox/LegacyCropModels.jl",
    ),
    pages = [
        "Home" => "index.md",
    ],
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/cropbox/LegacyCropModels.jl.git",
        devbranch = "main",
    )
end
