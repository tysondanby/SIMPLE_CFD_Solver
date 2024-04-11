

function uplot(p)
    L = p.geometry.length
    H = p.geometry.height
    n=p.discretizationsettings.n
    m=p.discretizationsettings.m
    x = 0.0:L*(2*n-1)/(2*n)
    y = 1/(2*m)*H:(2*m-1)/(2*m)*H
    us, ~ = extractvelocities(p)

    xpt = range(0, L, length=100)
    ypt = range(0, H, length=100)

    function z(xval,yval) 
        return interp2d(akima, x, y, us, xval, yval)
    end

    return heatmap(xpt,ypt,z)
end

function vplot(p)
    L = p.geometry.length
    H = p.geometry.height
    n=p.discretizationsettings.n
    m=p.discretizationsettings.m
    x = 1/(2*m-1)*L:L
    y = 0:H
    us, vs = extractvelocities(p)

    xpt = collect(1/(2*m-1)*L:L/100:L)
    ypt = collect(0:H/100:H)

    zpt =interp2d(akima, x, y, vs, xpt, ypt)
    #function z(xval,yval) 
    #    return interp2d(akima, x, y, vs, xval, yval)
    #end

    return heatmap(xpt,ypt,zpt)#, clims = (0.0,0.01))
end