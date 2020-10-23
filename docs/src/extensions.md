# Extensions

Originally, GeMM was designed to model island plant communities, but it has since
been adapted to various other study questions. These files include functions that
are not needed by the "core" model, but were added for these other studies.

## habitatchange.jl

Investigate the effect of temporal habitat change on plant communities.

```@autodocs
Modules = [GeMM]
Pages = ["habitatchange.jl"]
```

## invasion.jl

Carry out species invasions experiments in an island setting.

<!--FIXME Documenter.jl doesn't recognise function definitions that are not
	on the top level - like the `invade!()` functions. -->

```@autodocs
Modules = [GeMM]
Pages = ["invasion.jl"]
```

## zosterops.jl

Adapt GeMM to investigate eco-evolutionary population dynamics of East-African
Zosterops populations in the Taita Hills, Kenya.

```@autodocs
Modules = [GeMM]
Pages = ["zosterops.jl"]
```
