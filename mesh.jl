
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


function mesh!(problem::SIMPLEProblem2D)
    P = problem.initialpressurefunction
    n = problem.discretizationsettings.n
    m = problem.discretizationsettings.m
    u = problem.initialufunction
    v = problem.initialvfunction
    L = problem.geometry.length
    H = problem.geometry.height
    spacing = problem.discretizationsettings.spacingfunction #TODO: ONLY works for linear()
    ptsn = collect(0:1/(n-1):1)#collect(1/(2*n-1):2/(2*n-1):1)
    ptsm = collect(0:1/(m-1):1)#collect((1/(2*m)):1/m:(2*m - 1)/(2*m))
    pressurepoints = []
    upoints = []
    vpoints = []
    #Define pressure locations
    for i = 1:1:n
        for j = 1:1:m
            push!(pressurepoints,[spacing(ptsn[i])*L,spacing(ptsm[j])*H,0.0])
        end
    end
    horizontalspacing = norm(pressurepoints[1+m]-pressurepoints[1])
    verticalspacing = norm(pressurepoints[2]-pressurepoints[1])
    #Define u locations
    for j = 1:1:m
        push!(upoints,pressurepoints[j] - [horizontalspacing/2.0, 0.0, 0.0])
    end
    for i = 1:1:n
        for j = 1:1:m
            push!(upoints,pressurepoints[(i-1)*m+j] + [horizontalspacing/2.0, 0.0, 0.0])
        end
    end
    #Define v locations
    for i = 1:1:n
        push!(vpoints,pressurepoints[1+(i-1)*m] - [0.0, verticalspacing/2.0, 0.0])
        for j = 1:1:m
            push!(vpoints,pressurepoints[j+(i-1)*m] + [0.0, verticalspacing/2.0, 0.0])
        end
    end
    pressurenodes = fill(PressureNode(zeros(3),0,0.0,0,0.0,0,0.0,0,0.0,0.0,0,0.0),n*m)
    unodes = fill(StaggeredNode(zeros(3),0,0.0,0,0.0,0.0,0,0.0),n*m+m)
    vnodes = fill(StaggeredNode(zeros(3),0,0.0,0,0.0,0.0,0,0.0),n*m+n)
    #Compute positions and neighbors and values and areas and volumes
    #pressure first
    for i = 1:1:n
        for j = 1:1:m
            location = pressurepoints[(i-1)*m+j]
            neighborN = (i-1)*(m+1)+j+1
            neighborS = (i-1)*(m+1)+j
            neighborE = (i)*m+j
            neighborW = (i-1)*m+j
            areaX = horizontalspacing
            areaY = verticalspacing
            volume = horizontalspacing*verticalspacing
            bc = 0
            if i == n
                bc = 4#TODO: only works for the specific scenario lined out
            end
            pressurenodes[(i-1)*m+j] = PressureNode(location,neighborN,areaY,neighborS,areaY,neighborE,areaX,neighborW,areaX,volume,bc,P(location))
        end
    end
    #U velocities
    for j = 1:1:m
        location = upoints[j]
        forwardnode = j
        backwardnode = 0
        area = verticalspacing
        volume = horizontalspacing*verticalspacing
        bc = 1
        unodes[j] = StaggeredNode(location,forwardnode,area,backwardnode,area,volume,bc,u(location))
    end
    for i = 2:1:n
        for j = 1:1:m
            location = upoints[(i-1)*m+j]
            forwardnode = (i-1)*m+j
            backwardnode = (i-2)*m+j
            area = verticalspacing
            volume = horizontalspacing*verticalspacing
            bc = 0
            if j == 1
                bc = 3
            elseif j == m
                bc = 2
            end
            unodes[(i-1)*m+j] = StaggeredNode(location,forwardnode,area,backwardnode,area,volume,bc,u(location))
        end
    end
    for j = 1:1:m
        location = upoints[n*m+j]
        forwardnode = 0
        backwardnode = (n-1)*m+j
        area = verticalspacing
        volume = horizontalspacing*verticalspacing
        bc = 4
        unodes[n*m+j] = StaggeredNode(location,forwardnode,area,backwardnode,area,volume,bc,u(location))
    end
    #V velocities
    for i = 1:1:n
        location = vpoints[1+(i-1)*(m+1)]
        forwardnode = (i-1)*m +1
        backwardnode = 0
        area = horizontalspacing
        volume = horizontalspacing*verticalspacing/2
        bc = 3
        vnodes[1+(i-1)*(m+1)] = StaggeredNode(location,forwardnode,area,backwardnode,area,volume,bc,v(location))
        for j = 1:1:(m-1)
            location = vpoints[(i-1)*(m+1)+j+1]
            forwardnode = (i-1)*m + j + 1
            backwardnode = (i-1)*m + j
            area = horizontalspacing
            volume = horizontalspacing*verticalspacing
            bc = 0
            if i == 1
                bc = 1
            elseif i == n
                bc = 4
            end
            vnodes[(i-1)*(m+1)+j+1] = StaggeredNode(location,forwardnode,area,backwardnode,area,volume,bc,v(location))
        end
        location = vpoints[i*(m+1)]
        forwardnode = 0
        backwardnode = i*m
        area = horizontalspacing
        volume = horizontalspacing*verticalspacing/2
        bc = 2
        vnodes[i*(m+1)] = StaggeredNode(location,forwardnode,area,backwardnode,area,volume,bc,u(location))
    end

    problem.mesh=SIMPLEMesh(pressurenodes,unodes,vnodes)
    global testmesh = SIMPLEMesh(pressurenodes,unodes,vnodes)
end