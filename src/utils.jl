
macro varname(variable)
    string(variable)
end

" Expand S-expression for macro expression "
macro msexpr(value)
    return quote
        Meta.show_sexpr(@macroexpand $(value))
    end
end


" Checks array's size is equal to expected, raise ArgumentError otherwise "
macro checksize(expected, array)
    return quote
        actual = size($(esc(array)))
        expected = $(esc(expected))
        if expected != actual
            throw(ArgumentError("""Incorrect dimension for array:
                expected $expected != actual $actual"""))
        end
    end
end

" cdf of univariate unit normal distribution "
function normcdf(x)
    cdf(Normal(), x)
end

" Inverse of univariate gassian distribution's cumulative distribution function "
function invnormcdf(p; μ=0, σ=1)
    return quantile(Normal(μ, σ), p)
end


" Find difference operator of matrix or vector `B`. Results stored in `A`
    ⚈ no allocation
    ⚈ no array bound checcking
"
function diff!(A::AbstractMatrix, B::AbstractMatrix; dims::Integer)
    if dims == 1
        for i=1:size(A,1)-1
            for j=1:size(A,2)
                @inbounds A[i,j] = B[i+1,j] - B[i,j]
            end
        end
    elseif dims == 2
        for i=1:size(A,1)
            for j=1:size(A,2)-1
                @inbounds A[i,j] = B[i,j+1] - B[i,j]
            end
        end
    else
        throw(ArgumentError("dimension must be 1 or 2, got $dims"))
    end
end
