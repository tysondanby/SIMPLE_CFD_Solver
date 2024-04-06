function contains(v::vector{T1},value::{T1}) where T1 <: Any
    for element in v
        if value == element
            return true
        end
    end
    return false
end