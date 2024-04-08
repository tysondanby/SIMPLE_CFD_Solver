function steadysolve_1D!(p::SIMPLEProblem)
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
                p.mesh.unodes[i].value=p.initialufunction(p.mesh.unodes[i].position)
            end
            for i = 1:1:length(pressurelocations)
                p.mesh.pressurenodes[i].value = p.initialpressurefunction(pressurelocations[i])
            end

    rho = p.constantfunctions[1]
    A = p.geometry.areafunction
    uold1 = p.mesh.unodes[1].value
    convergencecriteria = 1.0
    itters = 0#itters < 3#
    while convergencecriteria > 1.0E-6 #while not converged TODO
        itters = itters +1
        #println()
        #println("Begin itteration $itters")
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
        Su[1] = (10.0 - nodee.value)*A(P) + Fw*(A(P)/A(w))*nodeP.value #TODO: here, a BC is applied that is unique to HW7
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
        w = nodew.position
        e = nodee.position
        Fw = 0.5*(nodeP.value + nodeW.value)*rho(w)*A(w)
        Fe = nodeP.value*rho(P)*A(P)
        aw = maximum([0,Fw]) #upwinding
        ae = 0
        ap = aw + ae + Fe - Fw
        Su[end] = (nodew.value - nodee.value)*A(P)
        Au[end,end-1] = -aw
        Au[end,end] = ap
        b[end] = Su[end]
        d[end] = A(P)/ap
        ustar = Au \ b
        
        #Find constants for pressure correction equations (at staggered nodes)
        n = length(pressurelocations)
        Ap = zeros(n,n)
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
            Ap[i,i-1] = -aw
            Ap[i,i] = ap
            Ap[i,i+1] = -ae
            b[i] = Fwstar - Festar
        end
        Ap[1,1] = 1.0
        Ap[end,end] = 1.0
        pprime = Ap \ b
        #println("Coefficient matrix for p'")
        #println(Ap)
        #println("b matrix")
        #println(b)
        #println("p'")
        #println(pprime)
        #Correct pressures (at pressure nodes)
        for i = 1:1:length(p.mesh.pressurenodes)
            p.mesh.pressurenodes[i].value = p.mesh.pressurenodes[i].value + p.pressurerelax*pprime[i]
        end
        
        #Correct velocities (at staggered nodes)
        ucalculated = similar(ustar)
        for i = 1:1:length(p.mesh.unodes)
            pprimewest = pprime[p.mesh.unodes[i].backwardneighbor]
            pprimeeast = pprime[p.mesh.unodes[i].forwardneighbor]
            ucalculated[i] = ustar[i] + d[i]*(pprimewest-pprimeeast)
            p.mesh.unodes[i].value = (1.0-p.urelax)*p.mesh.unodes[i].value + p.urelax*ucalculated[i]
        end
        #update first pressure node
        posa = p.mesh.pressurenodes[1].position
        node1 = p.mesh.unodes[p.mesh.pressurenodes[1].neighborE]
        pos1 = node1.position
        u1 = ucalculated[1]
        pcalc = 10.0 - 0.5*rho(posa)*(u1^2)*(A(pos1)/A(posa))^2#TODO: This is specific to the homework_7
        p.mesh.pressurenodes[1].value = p.mesh.pressurenodes[1].value*(1-p.pressurerelax)+p.pressurerelax*pcalc
        #p.mesh.pressurenodes[1].value = pcalc
        #println(u1)
        #for i = 1:1:n
        #    println(p.mesh.pressurenodes[i].value)
        #    #pos = p.mesh.unodes[i].position
        #    #println(rho(pos)*p.mesh.unodes[i].value*A(pos))
        #end
        #println("velocity:")
        #println(ustar)
        #for i = 1:1:n-1
        #    println(p.mesh.unodes[i].value)
        #end

        convergencecriteria = maximum(Au*ucalculated-Su)
        #println(convergencecriteria)
    end
end

function steadysolve_2D!(p::SIMPLEProblem)
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