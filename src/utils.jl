
import Statistics: quantile
import Distributions: Normal

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

" Inverse of univariate gassian distribution's cumulative distribution function "
function invnormcdf(p; μ=0, σ=1)
    return quantile(Normal(μ, σ), p)
end

" Read 1 float from a line at a time from a file and populates `estimates` "
function get_estimates!(estimates::Vector{Float64}, filename::String)
    open(filename, "r") do io
        for line in eachline(io)
            append!(estimates, parse(Float64, line))
        end
    end
end
