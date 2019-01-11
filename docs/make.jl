push!(LOAD_PATH, "../src/")

using Documenter, GeMM

makedocs(sitename="Island Model",
         modules = [GeMM],
         pages = ["index.md",
                  "io.md",
                  "aux.md",
                  "model.md"],
         Documenter.HTML(prettyurls=false))
