push!(LOAD_PATH, "../src/")

using Documenter, GeMM

makedocs(sitename="Island Speciation Model",
         modules = [GeMM],
         pages = ["index.md",
                  "framework.md",
                  "io.md",
                  "processes.md",
                  "extensions.md"])
