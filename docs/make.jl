push!(LOAD_PATH, "../src/")

using Documenter, GeMM

makedocs(sitename="Island Speciation Model",
         modules = [GeMM],
         pages = ["index.md",
                  "model.md",
                  "io.md",
                  "aux.md"])
