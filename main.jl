using LinearAlgebra, Plots
include("structs.jl")
include("functions.jl")
include("mesh.jl")
include("solve.jl")
include("postprocessing.jl")

function homework_7()
    n = 51 #number of pressure nodes
    α_P =0.01
    α_u = 0.01
    mdot = 1.0 #A guess
    function ρ(x)
        return 1.0
    end
    function A(x)
        return 0.5 - .2*x[1]
    end
    function P(x) #A guess
        return 10.0 - 5.0*x[1]
    end
    function u(x)
        return mdot/(ρ(x)*A(x))
    end
    dummyBC = ConstantBoundaryCondition("P",[0.0,0.0,0.0],10.0)

    geometry = OneDimensionalChannel(2.0,A)
    boundaryconditions = fill(dummyBC,0)
    settings = OneDimensionalMesher(n,linear)
    mesh = EmptyMesh()
    sources = fill(u,0)
    
    
    problem = SIMPLEProblem1D(geometry,boundaryconditions,settings,mesh,P,[ρ],u,sources,α_P,α_u)

    mesh!(problem)
    steadysolve_1D!(problem)
    xp,ps = extractpressures(problem)
    xu,us = extractvelocities(problem)
    p1 = plot(xp,ps,ylims=(0,10))
    mdot = ρ(xu[end])*us[end]*A(xu[end])

    function exact(pos)
        x = [pos,0.0,0.0]
        return 10 - (0.5*ρ(x)*0.44721^2)/(ρ(x)*A(x))^2
    end

    plot!(xp,@. exact(xp))

    return [p1],mdot,problem
end

plots,mdot, problem = homework_7()
mdot