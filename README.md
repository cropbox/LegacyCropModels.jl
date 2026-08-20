<a href="https://github.com/cropbox/Cropbox.jl"><img src="https://github.com/cropbox/Cropbox.jl/raw/main/docs/src/assets/logo.svg" alt="Cropbox" width="150"></a>

# LegacyCropModels.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://cropbox.github.io/LegacyCropModels.jl/stable/)
[![Latest Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://cropbox.github.io/LegacyCropModels.jl/dev/)

[LegacyCropModels.jl](https://github.com/cropbox/LegacyCropModels.jl) provides
declarative implementations of selected LINTUL models using the
[Cropbox](https://github.com/cropbox/Cropbox.jl) crop modeling framework. The
package includes LINTUL 1 for potential production, LINTUL 2 for water limited
production, and a spring wheat implementation of LINTUL 3 with water and
nitrogen dynamics.

The model systems retain the original LINTUL variable identifiers and add
descriptive aliases, physical units, and variable documentation. A shared
calendar based weather system supplies the daily inputs used by all three
models. It keeps cyclic calendar day of year separate from the monotonically
increasing FST simulation day used by source timing conditions.
`Lintul3SavedStateModel` is also provided for comparison with the saved
development state of the selected historical source program.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/cropbox/LegacyCropModels.jl")
```

## Getting Started

See the [documentation](https://cropbox.github.io/LegacyCropModels.jl/dev/) for
weather data preparation, model configuration, and simulation examples.

Variable declarations and their documentation can also be inspected directly:

```julia
using Cropbox
using LegacyCropModels

look(Lintul1System, :GLAE)
look(Lintul3System, :NNI)
```

## License

LegacyCropModels.jl is released under the MIT License. Historical FST source
programs and generated Fortran files are not included.
