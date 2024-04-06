abstract type CFDProblem end
abstract type GeometryDefinition end
abstract type BoundaryCondition end
abstract type Settings end
abstract type Mesh end
abstract type Node end

mutable struct SIMPLEProblem{T1} <: CFDProblem
    geometry::GeometryDefinition
    boundaryconditions::Vector{BoundaryCondition}
    discretizationsettings::Settings
    mesh::Mesh
    initialpressurefunction::Function
    constantfunctions::Vector{Function} #Things like rho, nu,k
    initialufunction::Function
    initialvfunction::Function
    sourcefunctions::Vector{Function}
    pressurerelax::T1
    urelax::T1
    vrelax::T1
end

mutable struct OneDimensionalChannel{T1} <: GeometryDefinition
    length::T1
    areafunction::Function
end

mutable struct ConstantBoundaryCondition{T1} <: BoundaryCondition
    value::String #"P" or "u" or "v"
    position::Vector{T1}
    value1::T1 # value1 = value2 * dT/dx + value3 * T
end

mutable struct PlanarBoundaryCondition{T1} <: BoundaryCondition
    value::String #"P" or "u" or "v"
    position::Vector{T1}
    normal::Vector{T1} #Direction of positive dT/dx as shown below.
    value1::T1 # value1 = value2 * dT/dx + value3 * T
    value2::T1
    value3::T1
end

mutable struct OneDimensionalMesher <: Settings
    n::Int #Number of principal nodes
    spacingfunction::Function #A one to one function, takes in evenly spaced values on 0:1 and makes them into nondimentional node positions.
end

mutable struct EmptyMesh <: Mesh end

mutable struct PressureNode{T1} <: Node
    position::Vector{T1}
    neighborN::Int
    areaN::T1
    neighborS::Int
    areaS::T1
    neighborE::Int
    areaE::T1
    neighborW::Int
    areaW::T1
    volume::T1
    boundarycondition::Int  #index of SIMPLEProblem.boundaryconditions corresponding to the BC. if no BC, left zero.
    value::T1
end

mutable struct StaggeredNode{T1} <: Node
    position::Vector{T1}
    forwardneighbor::Int #Indicies of the other type of node. If this node is principal, the indicies are for staggered nodes and vice-versa
    forwardarea::T1 #CS areas associated with each neighbor
    backwardneighbor::Int #Indicies of the other type of node. If this node is principal, the indicies are for staggered nodes and vice-versa
    backwardarea::T1 #CS areas associated with each neighbor
    volume::T1
    boundarycondition::Int #index of SIMPLEProblem.boundaryconditions corresponding to the BC. if no BC, left zero.
    value::T1
end

mutable struct SIMPLEMesh <: Mesh
    pressurenodes::Vector{PressureNode}
    unodes::Vector{StaggeredNode}
    vnodes::Vector{StaggeredNode}
end

