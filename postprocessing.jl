function extractpressures(p::SIMPLEProblem1D)
    n = p.discretizationsettings.n
    pressures = zeros(n)
    xs = zeros(n)
    for i = 1:1:n
        pressures[i] = p.mesh.pressurenodes[i].value
        xs[i] = p.mesh.pressurenodes[i].position[1]
    end
    return xs,pressures
end

function extractvelocities(p::SIMPLEProblem1D)
    n = p.discretizationsettings.n
    velocities = zeros(n-1)
    xs = zeros(n-1)
    for i = 1:1:n-1
        velocities[i] = p.mesh.unodes[i].value
        xs[i] = p.mesh.unodes[i].position[1]
    end
    return xs,velocities
end

function extractvelocities(p::SIMPLEProblem2D)
    n = p.discretizationsettings.n
    m = p.discretizationsettings.m
    us = zeros(m,n+1)
    vs = zeros(m+1,n)
    for i = m:-1:1
        for j = 1:1:n+1
            us[i,j] = p.mesh.unodes[(j-1)*m+i].value
        end
    end
    for i = m+1:-1:1
        for j = 1:1:n
            vs[i,j] = p.mesh.vnodes[(j-1)*(m+1)+i].value
        end
    end
    return us,vs
end

function extractpressures(p::SIMPLEProblem2D)
    n = p.discretizationsettings.n
    m = p.discretizationsettings.m
    Ps = zeros(m,n)
    for i = m:-1:1
        for j = 1:1:n
            Ps[i,j] = p.mesh.pressurenodes[(j-1)*m+i].value
        end
    end
    return Ps
end

function wallshear(p::SIMPLEProblem2D)
    n = p.discretizationsettings.n
    m = p.discretizationsettings.m
    mu = p.constantfunctions[2]([0.0,0.0,0.0])
    dx = norm(p.mesh.unodes[1].position - p.mesh.unodes[m+1].position)
    dy = norm(p.mesh.unodes[1].position - p.mesh.unodes[2].position)
    tau = 0.0
    #tau = tau + mu*.5*dx*p.mesh.unodes[1].value/(dy/2)
    #tau = tau + mu*.5*dx*p.mesh.unodes[m].value/(dy/2)
    for i = 2:1:(n-1)
        tau = tau + mu*dx*p.mesh.unodes[(i-1)*m+1].value/(dy/2)
        tau = tau + mu*dx*p.mesh.unodes[i*m].value/(dy/2)
    end
    return tau
end

function pressureplot(p::SIMPLEProblem2D)
    n = p.discretizationsettings.n
    m = p.discretizationsettings.m
    mu = mu = p.constantfunctions[2]([0.0,0.0,0.0])
    dx = norm(p.mesh.unodes[1].position - p.mesh.unodes[m+1].position)
    dy = norm(p.mesh.unodes[1].position - p.mesh.unodes[2].position)
    Ps = []
    xs = []
    j = Int(round(m/2))
    for i = 1:1:n
        push!(xs,p.mesh.pressurenodes[(i-1)*m+j].position[1])
        push!(Ps,p.mesh.pressurenodes[(i-1)*m+j].value)
    end
    return plot(xs*1000,Ps*1000,xlabel = "Axial Position (mm)",ylabel = "Pressure (mPa)")
end