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
    n = problem.discretizationsettings.n
    m = problem.discretizationsettings.m
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
            for i = 1:1:length(ulocations)
                p.mesh.unodes[i].value=p.initialufunction(p.mesh.unodes[i].position)
            end
            for i = 1:1:length(vlocations)
                p.mesh.vnodes[i].value=p.initialvfunction(p.mesh.vnodes[i].position)
            end
            for i = 1:1:length(pressurelocations)
                p.mesh.pressurenodes[i].value = p.initialpressurefunction(pressurelocations[i])
            end
        us, vs = extractvelocities(problem)
        global p3 = heatmap(us,clims = (0,0.005))
        global p4 = heatmap(vs,clims = (0,0.005))
    rho = p.constantfunctions[1]([0.0,0.0,0.0])
    mu = p.constantfunctions[2]([0.0,0.0,0.0])
    deltay = norm(p.mesh.unodes[1].position - p.mesh.unodes[2].position)#TODO assumes even grid spacing and uniform in x and y
    deltax = norm(p.mesh.unodes[1].position - p.mesh.unodes[1+m].position)
    convergencecriteria = 1.0
    itters = 0#itters < 3#
    nu = length(ulocations)
    Au = zeros(nu,nu)
    bu = zeros(nu,1)
    du = zeros(nu)
    Su = zeros(nu)
    
    nv = length(vlocations)
    Av = zeros(nv,nv)
    bv = zeros(nv,1)
    dv = zeros(nv)
    Sv = zeros(nv)

    np = length(pressurelocations)
    Ap = zeros(np,np)
    bp = zeros(np,1)
    
    ustar = zeros(nu)
    vstar = zeros(nv)
    pprime = zeros(np)
    vcalculated = similar(vstar)
    deltaus = similar(ustar)
    ucalculated = similar(ustar)


    while itters < 1000#convergencecriteria > 1.0E-6 #while not converged TODO
        
        itters = itters +1
        #Find constants for momentum equations meanwhile finding parameter d at each (staggered) node.
        #Calculate velocities
        
        for i = 1:1:nu
            nodeP = p.mesh.unodes[i]
            P = nodeP.position
            if nodeP.boundarycondition == 1
                Au[i,i] = 1.0
                bu[i] = p.boundaryconditions[1].value
                du[i] = 0.0 #Prescribed velocity means no pressure gradient will change it.
            elseif nodeP.boundarycondition == 4
                #TODO: check if this works. I set the outer node equal to the inner one.
                nodew = p.mesh.pressurenodes[nodeP.backwardneighbor]
                Au[i,i] = 1.0
                Au[i,nodew.neighborW] = -1.0
                du[i] = 0.0#du[i-m]#
            elseif nodeP.boundarycondition == 2 #Top
                nodew = p.mesh.pressurenodes[nodeP.backwardneighbor]
                nodee = p.mesh.pressurenodes[nodeP.forwardneighbor]
                nodeW = p.mesh.unodes[nodew.neighborW]
                nodeE = p.mesh.unodes[nodee.neighborE]
                nodene = p.mesh.vnodes[nodee.neighborN]
                nodenw = p.mesh.vnodes[nodew.neighborN]
                nodese = p.mesh.vnodes[nodee.neighborS]
                nodesw = p.mesh.vnodes[nodew.neighborS]
                nodeS = p.mesh.unodes[p.mesh.pressurenodes[nodese.backwardneighbor].neighborW]
                Fw = 0.5*(nodeP.value + nodeW.value)*rho*deltay
                Fe = 0.5*(nodeP.value + nodeE.value)*rho*deltay
                Fn = 0.5*(nodene.value + nodenw.value)*rho*deltax
                Fs = 0.5*(nodese.value + nodesw.value)*rho*deltax
                Dx = mu / deltax
                Dy = mu / deltay
                aw = (Dx*deltay + maximum([0,Fw])) #upwinding
                ae = (Dx*deltay + maximum([0,-Fe]))
                as = (Dy*deltax + maximum([0,Fs])) #upwinding
                ap = aw + ae + as + Fe - Fw + Fn - Fs + (Dy*deltax + maximum([0,-Fn]))
                Su[i] = (nodew.value - nodee.value)*deltay
                Au[i,nodee.neighborE] = -ae
                Au[i,nodew.neighborW] = -aw
                Au[i,p.mesh.pressurenodes[nodese.backwardneighbor].neighborW] = -as
                Au[i,i] = ap
                bu[i] = Su[i]
                du[i] = deltay/ap
            elseif nodeP.boundarycondition == 3 #Bottom
                nodew = p.mesh.pressurenodes[nodeP.backwardneighbor]
                nodee = p.mesh.pressurenodes[nodeP.forwardneighbor]
                nodeW = p.mesh.unodes[nodew.neighborW]
                nodeE = p.mesh.unodes[nodee.neighborE]
                nodene = p.mesh.vnodes[nodee.neighborN]
                nodenw = p.mesh.vnodes[nodew.neighborN]
                nodese = p.mesh.vnodes[nodee.neighborS]
                nodesw = p.mesh.vnodes[nodew.neighborS]
                nodeN = p.mesh.unodes[p.mesh.pressurenodes[nodene.forwardneighbor].neighborW]
                Fw = 0.5*(nodeP.value + nodeW.value)*rho*deltay
                Fe = 0.5*(nodeP.value + nodeE.value)*rho*deltay
                Fn = 0.5*(nodene.value + nodenw.value)*rho*deltax
                Fs = 0.5*(nodese.value + nodesw.value)*rho*deltax
                Dx = mu / deltax
                Dy = mu / deltay
                aw = (Dx*deltay + maximum([0,Fw])) #upwinding
                ae = (Dx*deltay + maximum([0,-Fe]))
                an = (Dy*deltax + maximum([0,-Fn]))
                ap = aw + ae + an + Fe - Fw + Fn - Fs + (Dy*deltax + maximum([0,Fs]))
                Su[i] = (nodew.value - nodee.value)*deltay
                Au[i,nodee.neighborE] = -ae
                Au[i,nodew.neighborW] = -aw
                Au[i,p.mesh.pressurenodes[nodene.forwardneighbor].neighborW] = -an
                Au[i,i] = ap
                bu[i] = Su[i]
                du[i] = deltay/ap
            else
                nodew = p.mesh.pressurenodes[nodeP.backwardneighbor]
                nodee = p.mesh.pressurenodes[nodeP.forwardneighbor]
                nodeW = p.mesh.unodes[nodew.neighborW]
                nodeE = p.mesh.unodes[nodee.neighborE]
                nodene = p.mesh.vnodes[nodee.neighborN]
                nodenw = p.mesh.vnodes[nodew.neighborN]
                nodese = p.mesh.vnodes[nodee.neighborS]
                nodesw = p.mesh.vnodes[nodew.neighborS]
                nodeN = p.mesh.unodes[p.mesh.pressurenodes[nodene.forwardneighbor].neighborW]
                nodeS = p.mesh.unodes[p.mesh.pressurenodes[nodese.backwardneighbor].neighborW]
                Fw = 0.5*(nodeP.value + nodeW.value)*rho*deltay
                Fe = 0.5*(nodeP.value + nodeE.value)*rho*deltay
                Fn = 0.5*(nodene.value + nodenw.value)*rho*deltax
                Fs = 0.5*(nodese.value + nodesw.value)*rho*deltax
                Dx = mu / deltax
                Dy = mu / deltay
                aw = (Dx*deltay + maximum([0,Fw])) #upwinding
                ae = (Dx*deltay + maximum([0,-Fe]))
                as = (Dy*deltax + maximum([0,Fs])) #upwinding
                an = (Dy*deltax + maximum([0,-Fn]))
                ap = aw + ae + an + as + Fe - Fw + Fn - Fs
                Su[i] = (nodew.value - nodee.value)*deltay
                Au[i,nodee.neighborE] = -ae
                Au[i,nodew.neighborW] = -aw
                Au[i,p.mesh.pressurenodes[nodene.forwardneighbor].neighborW] = -an
                Au[i,p.mesh.pressurenodes[nodese.backwardneighbor].neighborW] = -as
                Au[i,i] = ap
                bu[i] = Su[i]
                du[i] = deltay/ap
            end
        end
        ustar = Au \ bu

        
        for i = 1:1:nv
            nodeP = p.mesh.vnodes[i]
            P = nodeP.position
            if (nodeP.boundarycondition == 1) || (nodeP.boundarycondition == 2) || (nodeP.boundarycondition == 3) || (nodeP.boundarycondition == 4)
                Av[i,i] = 1.0
                bv[i] = 0.0
                dv[i] = 0.0 #Prescribed velocity means no pressure gradient will change it.
            else
                nodes = p.mesh.pressurenodes[nodeP.backwardneighbor]
                noden = p.mesh.pressurenodes[nodeP.forwardneighbor]
                nodeN = p.mesh.vnodes[noden.neighborN]
                nodeS = p.mesh.vnodes[nodes.neighborS]
                nodene = p.mesh.unodes[noden.neighborE]
                nodenw = p.mesh.unodes[noden.neighborW]
                nodese = p.mesh.unodes[nodes.neighborE]
                nodesw = p.mesh.unodes[nodes.neighborW]
                nodeE = p.mesh.vnodes[p.mesh.pressurenodes[nodene.forwardneighbor].neighborS]
                nodeW = p.mesh.vnodes[p.mesh.pressurenodes[nodenw.backwardneighbor].neighborS]
                Fw = 0.5*(nodesw.value + nodenw.value)*rho*deltay
                Fe = 0.5*(nodese.value + nodene.value)*rho*deltay
                Fn = 0.5*(nodeP.value + nodeN.value)*rho*deltax
                Fs = 0.5*(nodeP.value + nodeS.value)*rho*deltax
                Dx = mu / deltax
                Dy = mu / deltay
                aw = (Dx*deltay + maximum([0,Fw])) #upwinding
                ae = (Dx*deltay + maximum([0,-Fe]))
                as = (Dy*deltax + maximum([0,Fs])) #upwinding
                an = (Dy*deltax + maximum([0,-Fn]))
                ap = aw + ae + an + as + Fe - Fw + Fn - Fs
                Sv[i] = (noden.value - nodes.value)*deltax
                Av[i,p.mesh.pressurenodes[nodene.forwardneighbor].neighborS] = -ae
                Av[i,p.mesh.pressurenodes[nodenw.backwardneighbor].neighborS] = -aw
                Av[i,noden.neighborN] = -an
                Av[i,nodes.neighborS] = -as
                Av[i,i] = ap
                bv[i] = Sv[i]
                dv[i] = deltax/ap
            end
        end
        vstar = Av \ bv
        
        #Find constants for pressure correction equations (at staggered nodes)
        for i = 1:1:np
            nodeP = p.mesh.pressurenodes[i]
            P = nodeP.position
            if nodeP.boundarycondition == 4
                Ap[i,i] = 1.0
                bp[i] = 0.0#p.boundaryconditions[4].value
            else
                nodeWindex = p.mesh.unodes[nodeP.neighborW].backwardneighbor
                nodeEindex = p.mesh.unodes[nodeP.neighborE].forwardneighbor
                nodeSindex = p.mesh.vnodes[nodeP.neighborS].backwardneighbor
                nodeNindex = p.mesh.vnodes[nodeP.neighborN].forwardneighbor
                aw = rho*du[nodeP.neighborW]*deltay
                ae = rho*du[nodeP.neighborE]*deltay
                Fwstar = rho*ustar[nodeP.neighborW]*deltay
                Festar = rho*ustar[nodeP.neighborE]*deltay
                as = rho*dv[nodeP.neighborS]*deltax
                an = rho*dv[nodeP.neighborN]*deltax
                Fsstar = rho*vstar[nodeP.neighborS]*deltax
                Fnstar = rho*vstar[nodeP.neighborN]*deltax
                ap = aw + ae + an + as
                if nodeWindex == 0 #TODO: verify that these if statements do what you think they do
                    ap = ap -2*aw
                    ae = ae + aw
                    Ap[i,nodeEindex] = -ae
                elseif nodeEindex == 0
                    ap = ap -2*ae
                    aw = aw + ae
                    Ap[i,nodeWindex] = -aw
                else
                    Ap[i,nodeEindex] = -ae
                    Ap[i,nodeWindex] = -aw
                end
                if nodeSindex == 0 #TODO: verify that these if statements do what you think they do
                    ap = ap -2*as
                    an = an + as
                    Ap[i,nodeNindex] = -an
                elseif nodeNindex == 0
                    ap = ap -2*an
                    as = as + an
                    Ap[i,nodeSindex] = -as
                else
                    Ap[i,nodeSindex] = -as
                    Ap[i,nodeNindex] = -an
                end
                Ap[i,i] = ap
                bp[i] = Fwstar - Festar  + Fsstar - Fnstar
            end
        end
        pprime = Ap \ bp


        for i = 1:1:length(p.mesh.pressurenodes)
            p.mesh.pressurenodes[i].value = p.mesh.pressurenodes[i].value + p.pressurerelax*pprime[i]
        end
        
        #Correct velocities (at staggered nodes)
        for i = 1:1:length(p.mesh.unodes)
            if p.mesh.unodes[i].backwardneighbor == 0 || p.mesh.unodes[i].forwardneighbor == 0
                ucalculated[i] = ustar[i]
                deltaus[i] = 0.0
            else
                pprimewest = pprime[p.mesh.unodes[i].backwardneighbor]
                pprimeeast = pprime[p.mesh.unodes[i].forwardneighbor]
                ucalculated[i] = ustar[i] - du[i]*(pprimewest-pprimeeast)
                deltaus[i] = abs(ucalculated[i]-ustar[i])
            end
            p.mesh.unodes[i].value = (1.0-p.urelax)*p.mesh.unodes[i].value + p.urelax*ucalculated[i]
        end
        for i = 1:1:length(p.mesh.vnodes)
            if p.mesh.vnodes[i].backwardneighbor == 0 || p.mesh.vnodes[i].forwardneighbor == 0
                vcalculated[i] = vstar[i]
            else
                pprimewest = pprime[p.mesh.vnodes[i].backwardneighbor]
                pprimeeast = pprime[p.mesh.vnodes[i].forwardneighbor]
                vcalculated[i] = vstar[i] + dv[i]*(pprimewest-pprimeeast)
            end
            p.mesh.vnodes[i].value = (1.0-p.vrelax)*p.mesh.vnodes[i].value + p.vrelax*vcalculated[i]
        end

        convergencecriteria = maximum(deltaus)#maximum(Au*ucalculated-Su)
        if itters/10 == round(itters/10)
            println("Itteration $itters"*" with criteria: $convergencecriteria")
            global testus, testvs = extractvelocities(problem)
            global testPs = extractpressures(problem)
            global p1 = heatmap(us)#,clims = (0,0.005))
            global p2 = heatmap(vs)#,clims = (0,0.005))
        end
    end
end
