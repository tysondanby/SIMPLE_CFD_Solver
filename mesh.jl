
function mesh!(problem::SIMPLEProblem1D)
    A = problem.geometry.areafunction
    P = problem.initialpressurefunction
    n = problem.discretizationsettings.n
    u = problem.initialufunction
    L = problem.geometry.length
    spacing = problem.discretizationsettings.spacingfunction
    pts = collect(0:1/(n-1):1)
    pressurepoints = []
    velocitypoints = []
    push!(pressurepoints,[spacing(pts[1])*L,0.0,0.0])
    for i = 1:1:(n-1)
        push!(pressurepoints,[spacing(pts[i+1])*L,0.0,0.0])
        push!(velocitypoints,(pressurepoints[i] + pressurepoints[i+1])/2)
    end
    pressurenodes = fill(PressureNode(zeros(3),0,0.0,0,0.0,0,0.0,0,0.0,0.0,0,0.0),n)
    velocitynodes = fill(StaggeredNode(zeros(3),0,0.0,0,0.0,0.0,0,0.0),n-1)
    #Compute positions and neighbors and values
    pressurenodes[1] = PressureNode(pressurepoints[1],0,0.0,0,0.0,1,0.0,0,0.0,0.0,0,P(pressurepoints[1]))
    velocitynodes[1] = StaggeredNode(velocitypoints[1],2,0.0,1,0.0,0.0,0,u(velocitypoints[1]))
    for i = 2:1:n-1
        #pressurenodes[i].position = pressurepoints[i]
        #pressurenodes[i].neighborE = deepcopy(i)
        #pressurenodes[i].neighborW = deepcopy(i-1)
        #pressurenodes[i].value = P(pressurepoints[i])
        pressurenodes[i] = PressureNode(pressurepoints[i],0,0.0,0,0.0,deepcopy(i),0.0,deepcopy(i-1),0.0,0.0,0,P(pressurepoints[i]))
        #velocitynodes[i].position = velocitypoints[i]
        #velocitynodes[i].forwardneighbor = deepcopy(i+1)
        #velocitynodes[i].backwardneighbor = deepcopy(i)
        #velocitynodes[i].value = u(velocitypoints[i])
        velocitynodes[i] = StaggeredNode(velocitypoints[i],deepcopy(i+1),0.0,deepcopy(i),0.0,0.0,0,u(velocitypoints[i]))
    end
    pressurenodes[n] = PressureNode(pressurepoints[n],0,0.0,0,0.0,0,0.0,n-1,0.0,0.0,0,P(pressurepoints[n]))
    #Compute areas and volumes
    for node in pressurenodes
        posW = node.position
        posE = node.position
        if node.neighborW != 0
            posW = velocitynodes[node.neighborW].position
        end
        if node.neighborE != 0
            posE = velocitynodes[node.neighborE].position
        end
        AW = A(posW)
        AE = A(posE)

        node.volume = norm(posE-posW)*(AE+AW)/2
        node.areaW = AW
        node.areaE = AE
    end
    for node in velocitynodes
        posW = pressurenodes[node.backwardneighbor].position
        posE = pressurenodes[node.forwardneighbor].position
        AW = A(posW)
        AE = A(posE)

        node.volume = norm(posE-posW)*(AE+AW)/2
        node.backwardarea = AW
        node.forwardarea = AE
    end
    problem.mesh=SIMPLEMesh1D(pressurenodes,velocitynodes)
    global testmesh = SIMPLEMesh1D(pressurenodes,velocitynodes)
end