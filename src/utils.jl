
import Statistics: quantile
import Distributions: Normal


" Checks array's size is equal to expected, raise ArgumentError otherwise "
macro checksize(expected, array)
    return quote
        name = summary($(esc(array)))
        actual = size($(esc(array)))
        expected = $(esc(expected))
        if expected != actual
            throw(ArgumentError("""Incorrect dimension for $name:
                expected $expected != actual $actual"""))
        end
    end
end

" Inverse of univariate gassian distribution's cumulative distribution function "
function invnormcdf(p; μ=0, σ=1)
    return quantile(Normal(μ, σ), p)
end
