𝑀[90,50:55]

sort(𝑀; dims = 1)

filter(isnumeric, 𝑀)

function Base.sortperm(A::AbstractMatrix, dim::Integer)
    P = mapslices(sortperm, A, dim)
    if dim == 1
        for j = 1:size(P,2)
            offset = (j-1) * size(P,1)
            for i = 1:size(P,1)
                P[i,j] += offset
            end
        end
    else # if dim == 2
        for j = 1:size(P,2)
            for i = 1:size(P,1)
                P[i,j] = (P[i,j] - 1) * size(P,1) + i
            end
        end
    end
    return P
end