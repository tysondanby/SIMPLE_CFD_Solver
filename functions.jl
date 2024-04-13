function contains(v::Vector{T1},value::T1) where T1 <: Any
    for element in v
        if value == element
            return true
        end
    end
    return false
end

function linear(x)
    return x
end

function checkcoefficients(list)
    for element in list
        if element < 0.0
            #@warn("Negative coefficient")
        end
    end
end