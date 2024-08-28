# DES-SAFTVRMie
# Julia packages version
name = "Clapeyron"
uuid = "7c7805af-46cc-48c9-995b-ed0ed2dc909a"
authors = ["Hon Wa Yew <yewhonwa@gmail.com>", "Pierre Walker <pjwalker@caltech.edu>", "Andrés Riedemann <andres.riedemann@gmail.com>"]
version = "0.6.2"

[deps]
BlackBoxOptim = "a134a8b2-14d6-55f6-9291-3336d3ab0209"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
OrderedCollections = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
PackedVectorsOfVectors = "7713531c-48ef-4bdd-9821-9ff7a8736089"
PositiveFactorizations = "85a6dd25-e78a-55b7-8502-1745935b8125"
PrecompileTools = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
Preferences = "21216c6a-2e73-6563-6e65-726566657250"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
Scratch = "6c6a2e73-6563-6170-7368-637461726353"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
StableTasks = "91464d47-22a1-43fe-8b7f-2d57ee82463f"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
UUIDs = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[weakdeps]
CoolProp = "e084ae63-2819-5025-826e-f8e611a84251"
EoSSuperancillaries = "c1bf003f-4e47-49d9-bdfd-5a4051db3c04"
JutulDarcy = "82210473-ab04-4dce-b31b-11573c4f8e0a"
Metaheuristics = "bcdb8e00-2c21-11e9-3065-2b553b22f898"
MultiComponentFlash = "35e5bd01-9722-4017-9deb-64a5d32478ff"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[extensions]
ClapeyronCoolPropExt = "CoolProp"
ClapeyronJutulDarcyExt = "JutulDarcy"
ClapeyronMetaheuristicsExt = "Metaheuristics"
ClapeyronMultiComponentFlashExt = "MultiComponentFlash"
ClapeyronSuperancillaries = "EoSSuperancillaries"
ClapeyronSymbolicsExt = "Symbolics"
ClapeyronUnitfulExt = "Unitful"

[compat]
BlackBoxOptim = "^0.6.2"
CSV = "0.10"
CoolProp = "0.2"
DiffResults = "^1.1"
Downloads = "1"
EoSSuperancillaries = "1.2.3"
FillArrays = "^0.12, 0.13, 1"
ForwardDiff = "^0.10"
JSON3 = "1"
JutulDarcy = "0.2"
LinearAlgebra = "1"
LogExpFunctions = "^0.3.24"
Metaheuristics = "2, 3"
MultiComponentFlash = "1.1.13"
NLSolvers = "0.5"
OrderedCollections = "^1.5"
PackedVectorsOfVectors = "^0.1.2"
PositiveFactorizations = "^0.2"
PrecompileTools = "1"
Preferences = "1"
Roots = "2.1"
Scratch = "^1.1"
SparseArrays = "1"
StableTasks = "0.1"
StaticArrays = "^1.5.9"
Symbolics = "5"
Tables = "^1.8"
UUIDs = "1"
Unitful = "^1.12"
julia = "1.6"

[extras]
CoolProp = "e084ae63-2819-5025-826e-f8e611a84251"
EoSSuperancillaries = "c1bf003f-4e47-49d9-bdfd-5a4051db3c04"
JutulDarcy = "82210473-ab04-4dce-b31b-11573c4f8e0a"
Metaheuristics = "bcdb8e00-2c21-11e9-3065-2b553b22f898"
MultiComponentFlash = "35e5bd01-9722-4017-9deb-64a5d32478ff"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[targets]
test = ["CoolProp", "EoSSuperancillaries", "MultiComponentFlash", "Test", "Unitful"]
