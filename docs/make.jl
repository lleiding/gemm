push!(LOAD_PATH, "../src/")

using Documenter, GeMM

makedocs(sitename="Island Speciation Model",
         modules = [GeMM],
         pages = ["index.md",
                  "io.md",
                  "aux.md",
                  "model.md"])
