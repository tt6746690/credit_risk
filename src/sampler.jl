include("utils.jl")

"""
Slice sampling

    References
    Iain Murray's matlab impl for multivariate slice sampling:
        https://homepages.inf.ed.ac.uk/imurray2/teaching/09mlss/slice_sample.m
    David MacKay's book:
        http://www.inference.org.uk/itprnn/book.pdf
    Radford Neal's R impl for univariate slice sampling:
        https://www.cs.toronto.edu/~radford/ftp/slice-R-prog
    Alternative impl that is not well supported, thorough:
        https://github.com/brian-j-smith/Mamba.jl/blob/release-0.11/src/samplers/slice.jl

x0: Initial point
f: The density function to sample from
N: The number of simulate samples

step_limit (m): Limit on steps (default Inf)
width (w): The size of step for creating interval [L, R] (defaults to 1)
burn: The number of samples to discard before saving (defaults to 0)
thin: Discards `thin-1` samples and return the next to avoid correlation (defaults to 1)
logpdf: `true` if g is the log of probability density trying to sample from (defaults to `false`)

"""
function slicesample(
    x0::Array{Float64, 1},
    f,
    N::Int;
    step_limit::Int = typemax(Int),
    width::Union{Array{Float64, 1}, Float64} = 1.,
    burn::Int = 0,
    thin::Int = 1,
    logpdf::Bool = false)

    D = size(x0, 1)
    if isa(width, Float64)
        width = fill(width, D)
    end

    if logpdf == false
        g(x) = log(f(x))
    else
        g = f
    end

    @checksize (D,) x0
    @checksize (D,) width

    L = zeros(D)
    R = zeros(D)
    x = zeros(D)
    x′= zeros(D)
    samples = zeros(D, N)

    x[:] = x0

    # Log density at initial point
    gx = g(x)

    for i in 1:(N*thin+burn)

        # Determine the slice level, in log terms.
        logy = log(rand()) + gx

        # Sample from each dimension separately
        for d in 1:D
            x′[:] = x
            L[:] = x
            R[:] = x

            # Find interval x ∈ [L, R] to sample from
            u = width[d] * rand()
            L[d] = x[d] - u
            R[d] = x[d] + (width[d] - u)

            # Step out
            if step_limit == Inf
                while (g(L) > logy)
                    L[d] -= width[d]
                end

                while (g(R) > logy)
                    R[d] += width[d]
                end
            elseif step_limit > 1
                J = floor(step_limit * rand())
                K = (step_limit - 1) - J

                while (J > 0) && (g(L) > logy)
                    L[d] -= width[d]
                    J -= 1
                end

                while (K > 0) && (g(R) > logy)
                    R[d] += width[d]
                    K -= 1
                end
            end

            # Sample from [L, R], shrinking it on each rejection
            #     until we have found good sample x' in the slice
            while true
                x′[d] = L[d] + rand() * (R[d] - L[d])
                gx = g(x′)

                gx >= logy && break

                if x′[d] > x[d]
                    R[d] = x′[d]
                else
                    L[d] = x′[d]
                end
            end

            x[d] = x′[d]
        end

        if i > burn && ((i-burn) % thin == 0)
            samples[:, div(i-burn, thin)] = x[:]
        end
    end

    samples
end
