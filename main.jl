using LinearAlgebra, Plots, FLOWMath, SparseArrays
include("structs.jl")
include("functions.jl")
include("visualization.jl")
include("mesh.jl")
include("solve.jl")
include("postprocessing.jl")

function homework_7()
    n = 201 #number of pressure nodes
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
    mdot = ρ(xu[end])*us[end]*A(xu[end])

    function exact(pos)
        x = [pos,0.0,0.0]
        return 10 - (0.5*ρ(x)*0.44721^2)/(ρ(x)*A(x))^2
    end

    p1 = plot(xp,[ps (@. exact(xp))],ylims=(0,10),xlabel = "Position (m)", ylabel = "Pressure (Pa)",label=["SIMPLE" "Theoretical"])

    return [p1],mdot,problem
end

#plots,mdot, problem = homework_7()
#mdot
L = .05
H = 0.01
uin = 0.001
P0 = 0.0
entrance = LabelBoundaryCondition("Left",uin)
top = LabelBoundaryCondition("Top",0.0)
bottom = LabelBoundaryCondition("Bottom",0.0)
outlet = LabelBoundaryCondition("Right",P0)
grid = 30
n = 5*grid
m = 4*grid
α = 0.5
α_P =α
α_u = α
α_v = α

function ρ(x)
    return 1000.0
end
function mu(x)
    return 0.001  
end 
function P(x) #A guess
    return .007 - x[1]*(.007/L)
end
function u(x)
    return .001
end
function v(x)
    return 0.0#.0001
end

geometry = TwoDimensionalChannel(L,H)
boundaryconditions = [entrance, top, bottom, outlet]
settings = TwoDimensionalMesher(n,m,linear)
mesh = EmptyMesh()
sources = fill(u,0)

problem = SIMPLEProblem2D(geometry,boundaryconditions,settings,mesh,P,[ρ,mu],u,v,sources,α_P,α_u,α_v)

mesh!(problem)
steadysolve_2D!(problem)
p5 = heatmap(testus[:,2:end]*1000,aspect_ratio = 5/20,clims = (0,1.5),cbar_title = " \nHorizontal Velocity (mm/s)",right_margin=4Plots.mm,xticks = [],yticks=[], c = :jet1)
p6 = heatmap(testvs[:,1:end]*1000,aspect_ratio = 5/20,clims = (-.35,.35),cbar_title = " \nVertical Velocity (mm/s)",right_margin=4Plots.mm,xticks = [],yticks=[], c = :jet1)
p7 = heatmap(testPs[:,1:end]*1000,aspect_ratio = 5/20,clims = (0,8),cbar_title = " \nRelative Pressure (mPa)",right_margin=4Plots.mm,xticks = [],yticks=[], c = :jet1)