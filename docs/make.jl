using Documenter
using SDtoolbox

makedocs(
    sitename = "SDtoolbox",
    format = Documenter.HTML(prettyurls = false),
    modules = [SDtoolbox],
    pages = [
        "Home" => "index.md",
        "Postshock" => "postshock.md",
        "ZND" => "znd.md",
        "CV" => "cv.md",
        "Cell size calculation" => "cell_size.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/daharn/SDtoolbox.jl"
)
