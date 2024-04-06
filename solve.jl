function steadysolve_1D(p::SIMPLEProblem)
    pressurelocations = []
    for node in p.mesh.pressurenodes
        push!(pressurelocations,node.position)
    end
    ulocations = []
    for node in p.mesh.unodes
        push!(ulocations,node.position)
    end
    #Initial guesses of P* and u* and v* and....
        	#u v first
            for i = 1:1:length(ulocations)
                p.mesh.unodes[i].value=p.initialufunction(pos)
            end
            for i = 1:1:length(pressurelocations)
                p.mesh.pressurenodes[i].value = p.initialpressurefunction(pressurelocations[i])
            end

    rho = p.constantfunction[1]
    A = p.geometry.areafunction

    convergencecriteria = 1.0

    while convergencecriteria > 1.0E-6 #while not converged TODO
        #Find constants for momentum equations meanwhile finding parameter d at each (staggered) node.
        #Calculate velocities
        n = length(ulocations)
        Au = zeros(n,n)
        b = zeros(n,1)
        d = zeros(n)
        Su = zeros(n)
        nodeP = p.mesh.unodes[1]
        P = nodeP.position
        nodew = p.mesh.pressurenodes[nodeP.backwardneighbor]
        nodee = p.mesh.pressurenodes[nodeP.forwardneighbor]
        nodeE = p.mesh.unodes[nodee.neighborE]
        w = nodew.position
        e = nodee.position
        uA = nodeP.value*A(P)/A(w)
        Fw = rho(w)*uA*A(w)
        Fe = 0.5*(nodeP.value + nodeE.value)*rho(e)*A(e)
        ae = maximum([0,-Fe])
        ap = Fe + Fw*0.5*(A(P)/A(w))^2
        Su[1] = (nodew.value - nodee.value)*A(P) + Fw*(A(P)/A(w))*nodeP.value
        Au[1,1] = ap
        b[1] = Su[1]
        d[1] = A(P)/ap
        for i = 2:1:(n-1)
            nodeP = p.mesh.unodes[i]
            P = nodeP.position
            nodew = p.mesh.pressurenodes[nodeP.backwardneighbor]
            nodee = p.mesh.pressurenodes[nodeP.forwardneighbor]
            nodeW = p.mesh.unodes[nodew.neighborW]
            nodeE = p.mesh.unodes[nodee.neighborE]
            w = nodew.position
            e = nodee.position
            Fw = 0.5*(nodeP.value + nodeW.value)*rho(w)*A(w)
            Fe = 0.5*(nodeP.value + nodeE.value)*rho(e)*A(e)
            aw = maximum([0,Fw]) #upwinding
            ae = maximum([0,-Fe])
            ap = aw + ae + Fe - Fw
            Su[i] = (nodew.value - nodee.value)*A(P)
            Au[i,i-1] = -aw
            Au[i,i] = ap
            Au[i,i+1] = -ae
            b[i] = Su[i]
            d[i] = A(P)/ap
        end
        nodeP = p.mesh.unodes[end]
        P = nodeP.position
        nodew = p.mesh.pressurenodes[nodeP.backwardneighbor]
        nodee = p.mesh.pressurenodes[nodeP.forwardneighbor]
        nodeW = p.mesh.unodes[nodew.neighborW]
        nodeE = p.mesh.unodes[nodee.neighborE]
        w = nodew.position
        e = nodee.position
        Fw = 0.5*(nodeP.value + nodeW.value)*rho(w)*A(w)
        Fe = nodeP.value*rho(P)*A(P)
        aw = maximum([0,Fw]) #upwinding
        ae = 0
        ap = aw + ae + Fe - Fw
        Su[end] = (nodew.value - nodee.value)*A(P)
        Au[end,i-1] = -aw
        Au[end,end] = ap
        b[end] = Su[end]
        d[end] = A(P)/ap
        ustar = Au \ b

        #Find constants for pressure correction equations (at staggered nodes)
        n = length(pressurelocations)
        A = zeros(n,n)
        b = zeros(n,1)
        for i = 2:1:(n-1)
            nodeP = p.mesh.pressurenodes[i]
            P = nodeP.position
            nodew = p.mesh.unodes[nodeP.neighborW]
            nodee = p.mesh.unodes[nodeP.neighborE]
            w = nodew.position
            e = nodee.position
            aw = rho(w)*d[nodeP.neighborW]*A(w)
            ae = rho(e)*d[nodeP.neighborE]*A(e)
            Fwstar = rho(w)*ustar[nodeP.neighborW]*A(w)
            Festar = rho(e)*ustar[nodeP.neighborE]*A(e)
            ap = aw + ae
            A[i,i-1] = -aw
            A[i,i] = ap
            A[i,i+1] = -ae
            b[i] = Fwstar - Festar
        end
        A[1,1] = 1.0
        A[end,end] = 1.0
        pprime = A \ b

        #Correct pressures (at pressure nodes)
        for i = 1:1:length(p.mesh.pressurenodes)
            p.mesh.pressurenodes[i] = p.mesh.pressurenodes[i] + p.pressurerelax*pprime[i]
        end
        #Correct velocities (at staggered nodes)
        for i = 1:1:length(p.mesh.unodes)
            pprimewest = pprime[p.mesh.unodes[i].backwardneighbor]
            pprimeeast = pprime[p.mesh.unodes[i].forwardneighbor]
            ucalculated = ustar[i] + d[i]*(pprimewest-pprimeeast)
            p.mesh.unodes[i] = (1.0-p.urelax)*p.mesh.unodes[i] + p.urelax*ucalculated
        end

        convergencecriteria = maximum(Au*ustar-Su)
    end
end

function steadysolve_2D(p::SIMPLEProblem)
    pressurelocations = []
    for node in p.mesh.pressurenodes
        push!(pressurelocations,node.position)
    end
    ulocations = []
    for node in p.mesh.unodes
        push!(ulocations,node.position)
    end
    vlocations = []
    for node in p.mesh.vnodes
        push!(vlocations,node.position)
    end
    #Initial guesses of P* and u* and v* and....
        	#u v first
    ufield = zeros(length(ulocations))
    for i = 1:1:length(ulocations)
        ufield[i]=p.initialufunction(pos)
    end
    vfield = zeros(length(vlocations))
    for i = 1:1:length(vlocations)
        vfield[i]=p.initialufunction(pos)
    end
    pressurefield = zeros(length(pressurelocations))
    for i = 1:1:length(pressurelocations)
        pressurefield[i] = p.initialpressurefunction(pressurelocations[i])
    end

    rho = p.constantfunction[1]
    nu = p.constantfunction[2]

    while true #while not converged TODO
        #Find constants for momentum equations meanwhile finding parameter d at each (staggered) node.
        #Calculate velocities
        n = length(ulocations)
        A = zeros(n,n)
        b = zeros(n,1)
        for i = 1:1:n
            nodeP = p.mesh.unodes[i]
            P = nodeP.position
            if nodeP.boundarycondition == 0
                #TODO: handle differently if forwardneighbor or backwardneighbor are 0
                nodew = p.mesh.pressurenodes[nodeP.backwardneighbor]
                nodee = p.mesh.pressurenodes[nodeP.forwardneighbor]
                nodeW = p.mesh.unodes[nodew.neighborW]
                nodeE = p.mesh.unodes[nodee.neighborE]
                nodenw = p.mesh.vnodes[nodew.neighborN]
                nodesw = p.mesh.vnodes[nodew.neighborS]
                nodene = p.mesh.vnodes[nodee.neighborN]
                nodese = p.mesh.vnodes[nodee.neighborS]
                w = nodew.position
                e = nodee.position
                W = nodeW.position
                E = nodeE.position
                Fw = 0.5*(rho(P)*nodeP.value + rho(W)*nodeW.value)
                Fe = 0.5*(rho(P)*nodeP.value + rho(E)*nodeE.value)
                A[i,i] = stg#TODO
                b[i] = stg#TODO
            else #TODO only BCs of first kind
                A[i,i] = 1.0
                b[i] = p.boundaryconditions[nodeP.boundarycondition].value1
            end
        end
        
        #Find constants for pressure correction equations (at staggered nodes)
        #Correct pressures (at pressure nodes)
        #Correct velocities (at staggered nodes)

    end
end