using SeededAlignment
using Documenter

DocMeta.setdocmeta!(
    SeededAlignment,
    :DocTestSetup,
    :(using SeededAlignment);
    recursive=true,
)

makedocs(;
    modules = [SeededAlignment],
    authors = "Willó Corry and contributors",
    sitename = "SeededAlignment.jl",
    format = Documenter.HTML(;
        canonical = "https://MurrellGroup.github.io/SeededAlignment.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        #"Tutorial" => "tutorial.md"
        #"Customizing Alignments" => "customizing alignments.md",
        "Types" => "types.md",
        "API Reference" => "api.md",
    ],
    doctest = false,
)

deploydocs(;
    repo="github.com/MurrellGroup/SeededAlignment.jl",
    devbranch="main",
)