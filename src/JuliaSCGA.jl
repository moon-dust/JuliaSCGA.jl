module JuliaSCGA

include("UnitCell.jl")
export UnitCell, addInteraction!, setInteractionOnsite!, setField!, addBasisSite!, getDist
export getFourier_iso, solveLambda_iso, getCorr_iso
export getFourier_aniso, solveLambda_aniso, getCorr_aniso
export getChiT_iso

end # module
