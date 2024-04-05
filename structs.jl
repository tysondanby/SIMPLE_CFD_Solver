abstract type CFDProblem end

mutable struct SIMPLEProblem <: CFDProblem
    geometry
    boundaryconditions
    discretizationsettings
    mesh
    valuefunctions
    sourcefunctions
    valuesatprincipalnodes
    valuesatstaggerednodes
end

abstract type GeometryDefinition end

mutable struct OneDimensionalChannel <: GeometryDefinition
    length
    areafunction
end

abstract type BoundaryCondition end


    