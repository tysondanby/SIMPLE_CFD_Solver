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