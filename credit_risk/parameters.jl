module CreditRisk

" Checks array's size is equal to expected, raise ArgumentError otherwise "
macro checksize(expected, array)
    # use caller's context
    expected = esc(expected)
    array = esc(array)
    return quote
        name = summary($array)
        actual = size($array)
        expected = $expected
        if expected != actual
            throw(ArgumentError("""Incorrect dimension for $name:
                expected $expected != actual $actual"""))
        end
    end
end

" Parameters related to the credit risk problem "
struct Parameter
    # Number of creditors
    N::Int32
    # Number of credit states
    C::Int32
    # Dimension for Systematic Risk Factor Z, i.e. Z ∼  N(0, I_S)
    S::Int32
    # Tail for computing P(L ⩾ l)
    l::Float32

    # Computed
    # Credit Migration Matrix,
    # CMM[n, c] has probability of n-th creditor moving to credit state c
    cmm::Array{Float32, 2}
    # Exposure at Default
    # Value lost if n-th creditor is in credit state c=1, i.e. defaults
    ead::Array{Float32, 1}
    # Percentage loss/gain
    # lgc[n, c] represent loss from creditor n move to credit state c
    lgc::Array{Float32, 2}
    # β[n, c] indicates n-th creditor's sensitivity to systematic risk factor Z
    β::Array{Float32, 2}
    # Initial credit state for each creditor
    # In case of binary credit states, initially at c=2, i.e. non-default state
    cn::Array{Int32, 1}
    # Threshold for creditor n migrating from current credit state `cn` to any `c`
    # H[n, c] represent threshold for migrating to c from current state
    H::Array{Float32, 2}

    function Parameter(N, C, S, l, cmm, ead, lgc, β, cn, H)
        N >= 0 || throw(ArgumentError("Invalid number of creditors: $N"))
        C >= 0 || throw(ArgumentError("Invalid number of credit states: $C"))
        S >= 0 || throw(ArgumentError("Invalid Dimension for systematic risk factor: $S"))
        0 <= l <= 1 || throw(ArgumentError("Invalid tail probability $l"))
        @checksize (N, C)   cmm
        @checksize (N,)     ead
        @checksize (N, C)   lgc
        @checksize (N, S)   β
        @checksize (N,)     cn
        @checksize (N, C)   H
        new(N, C, S, l, cmm, ead, lgc, β, cn, H)
    end
end

n = 20
c = 4
s = 10

x = ones(n, c)
y = ones(n, s)
z = ones(n)


p = Parameter(n,c,s,0.2,x,z,x,y,z,x)
print(summary(p))

end
