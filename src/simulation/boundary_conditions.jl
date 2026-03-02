# Functions for applying boundary conditions to particle positions
function periodic_BC_array!(xy::AbstractArray{Float64,2}, L::Real, Np::Int)
    @inbounds for i in 1:Np
        x, y = xy[i, 1], xy[i, 2]

        if abs(x) > L/2
            xy[i, 1] = x - sign(x) * L
        end

        if abs(y) > L/2
            xy[i, 2] = y - sign(y) * L
        end
    end
    return nothing
end