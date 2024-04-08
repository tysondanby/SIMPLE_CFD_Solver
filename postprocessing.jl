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