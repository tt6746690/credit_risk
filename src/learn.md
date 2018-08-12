

#### Resources

+ [julia](https://github.com/JuliaLang/julia)
+ [julia doc](https://docs.julialang.org/en/latest/)
+ [distribution](https://juliastats.github.io/Distributions.jl/)

#### Starter

```julia
module Learn

function simulate_π(n)
    in_quad_circle(x, y) = x^2 + y^2 <= 1
    τ = 0
    for i = 1:n
        x, y = rand(), rand()
        if in_quad_circle(x, y)
            τ += 1
        end
    end
    π = 4*τ / n
    return π
end

π = simulate_π(10000)
println(π)  # 3.142

end
```


##### Random Numbers

```julia
rand(Int, 2)
rand(MersenneTwister(0), Dict(1=>2, 3=>4)) # 1=>2
randn((2, 3))
```
- uses Mersenne twister library
    - uniform random numbers
- `rand([rng=GLOBAL_RNG], [S], [dims...])`
    - `rand!([rng], A, [S=eltype(A)])`
- `randn([rng], [T=Float64], [dims...])`
    - normally distributed random number of type `T` with mean 0 and stddev 1
    - `randn!([rng=GLOBAL_RNG], A::AbstractArray) -> A`
- `randstring([rng=GLOBAL_RNG], [chars], [len=8])`
- `seed!([rng=GLOBAL_RNG], seed) -> rng`
    - reseed random number generator
    - if seed provided, sequence of number is reproducible



##### Modules

```julia
module MyModule
using Lib

using BigLib: thing1, thing2

import Base.show

importall OtherLib

export MyType, foo

struct MyType
    x
end

bar(x) = 2x
foo(a::MyType) = bar(a.x) + 1

show(io::IO, a::MyType) = print(io, "MyType $(a.x)")
end
```
- module
    - introduce global scope with `module Name ... end`



##### Constructor
```julia
type Foo
    bar
    baz
end
foo = Foo(1,2)
```
- outer constructor methods
    - convenience in constructing objects
        - can only call inner constructor
    - create `Foo` with 0 or 1 argument
    - `Foo(x) = Foo(x, x)`
    - note this once is called `Foo()=  Foo(0)` instead of default arg
- inner constructor methods
    - goal
        - enforce invariants and
        - allow construction of self-referential objects
    - syntax
        - declared inside `type` block
        - access to locally existent function called `new`, which creates objects of the block's type
        - if defined, no default constructor methods
            - equivalent to `Foo(bar, baz) = new(bar, baz)`
    ```julia
    struct OrderedPair
        x::Real
        y::Real
        OrderedPair(x,y) = x > y ? error("out of order") : new(x, y)
    end
    OrderedPair(2,1) # gives error
    ```
- best practice
    - as few as possible inner constructor, taking all arguments explicitly and enforce error checking and transformation
    - outer function used to supplied default arguments or auxiliary transformations
- incomplete initialization
    ```julia
    mutable struct SelfReferential
       obj::SelfReferential
    end
    b = SelfReferential(a)  # where do obj come from
    ```
    - `new`
        - allowed to be called with fewer than number of fields the type has
        - return an object with unspecified fields uninitialized
        - inner constructor method use this incomplete object to finish initialization and assign it to `b`
    ```julia
    mutable struct SelfReferential
        obj::SelfReferential
        SelfReferential() = (x = new(); x.obj = x)
    end
    ```
    - plain data type
        - primitive types and immutable struct of other plain data types is uninitialized to random value
- parametric constructor
    ```julia
    struct Point{T<:Real}
        x::T
        y::T
    end

    Point(1,2)          # Point{Int64}(1, 2)
    Point(1.0,2.0)      # Point{Float64}(1.0, 2.5)
    Point{Int64}(1, 2)  # Point{Int64}(1, 2)
    ```
    - `Point{T}` is a distinct constructor function for each type `T`
        - `Point{T<:Real}` provides an inner constructor for each `T`
    ```julia
    # Automatically generate 1 inner constructor for each possible T<:Real
    # and a general outer constructor that takes in a pair of real arguments
    struct Point{T<:Real}
        x::T
        y::T
        Point{T}(x, y) where {T<:Real} = new(x, y)
    end
    Point(x::T, y::T) where {T<:Real} = Point{T}(x, y)
    ```
    - examples
        - `Point{Int64}(1,2)` invokes inner constructor
        - `Point(1,2)` invokes outer constructor,
            - works for args of same type
            - raise no method error if args of different type
    - A constructor that promotes integer to float
        - `Point(x::Int64, y::Float64) = Point(convert(Float64,x),y)`
    - Or just use `promote`
        - promote to a common type
        - `Point(x::Real, y::Real) = Point(promote(x,y)...)`
- case study
```julia
struct OurRational{T<:Integer} <:Real
    num::T
    den::T
    function ourRational{T}(num::T, den::T) where T<:Integer
        # Invariant: checks they are all non-zero
        if num == 0 && den == 0
            error("invalid rational 0/0")
        end
        # Transformation: Convert to lowest term ...
        g = gcd(den, num)
        num = div(num, g)
        den = div(den, g)
        new(num, den)
    end
end

# for OurRational(1,1)
OurRational(n::T, d::T) where {T<:Integer} = OurRational{T}(n,d)
# for promoting n,d to a common type
OurRational(n::Integer, d::Integer) = OurRational(promote(n,d)...)
# for default argument for `d=1`
OurRational(n::Integer) = OurRational(n,one(n))
# define an operator ⊘
⊘(n::Integer, d::Integer) = OurRational(n,d)
# custom for rational y
⊘(x::Integer, y::OurRational) = (x*y.den) ⊘ y.num
# for complex values
⊘(x::Complex, y::Real) = complex(real(x) ⊘ y, imag(x) ⊘ y)
⊘(x::Real, y::Complex) = (x*y') ⊘ real(y*y')
function ⊘(x::Complex, y::Complex)
   xy = x*y'
   yy = real(y*y')
   complex(real(xy) ⊘ yy, imag(xy) ⊘ yy)
end
```
- outer-only constructors
    - to suppress





##### Type System
- type system
    - concrete types
        - may not subtype each other
        - final, and only have abstract type as their super types
    - no division between object and non-object values
        - all values are true objects having a type belonging to a fully connected type tree
    - no compile-time type
        - only run-time type
    - only values, not variable has type
        - variable are names bound to values
    - abstract/concrete types can be paramterized by other types, or by symbols
- annotation with `::`,
    - is an instance of
    - type assertion
        - assert value of expression on the left is an instance of type on the right
            - if type on right is concrete, then value on left must have that type as impl
            - if type on left is abstract, then value on left must be implemented by a concrete type that is a subtype of the abstract type
- abstract types
    - describe a set of related concrete types
    - define new abstract type
        - default parent type is `Any`
    ```julia
    abstract type <<name>> end
    abstract type <<name>> <: <<supertype>> end
    ```
    - `Union{}`
        - bottom of type graph
        - no object is an instance of `Union{}`
    - `<:`
        - is a subtype of
- primitive types
    ```
    primitive type «name» «bits» end
    primitive type «name» <: «supertype» «bits» end
    ```
- composite types
    - a collection of anmed fields, an instance of which can be treated as a single value
    - all values are objects, but functions are not bundled with associated methods
    ```julia
    struct Foo
        bar
        baz::Int
        qux::Float64
    end
    foo = Foo("Hello, world.", 23, 1.5)
    fieldnames(foo)
    ```
    - constructor
        - automatically generated
            - accepts any, and call `convert` to convert them to types of the field
            - accepts arguments that match the field types exactly
    - immutable, cannot be modified after construction
        - can contain mutable fields, i.e. `array`
        - just the fields of object cannot be changed to point to different objects
        - use `mutable` if necessary, object allocated on heap always ...
    - immutability
        - immutable object pass around by copy
        - mutable type pass around by reference
- declared types
    - properties of abstract,primitive,composite types
        - are explicitly declared
        - have names
        - have supertypes
        - may have parameters
    - all an instance ofr `DataType`, which is the type of any of these types
        - `typeof(Real)`
        - can be abstract or concrete.
        - can have
            - specified size
            - storage layout
            - field names (optionally)
- Type Unions
    - a special abstract type, includes object of all instances of any of its argument types
    - `IntOrString = Union{Int, AbstractString}`
- parametric types
    - types can take parameters
    - All declared types (`DataType`) can be paramterized
    - parametric composite types
        ```
        struct Point{T}
            x::T
            y::T
        end
        Point{Float64} <: DataType
        ```
- operation on types
    - `isa(obj, type)`
    - `typeof()`
    - `supertype()`



##### Multi-dimensional Arrays
- array
    - a collection of objects stored in a multi-dimensional grid
    - function argument pass by reference, so mutation to input argument visible in parent function
- basics
    - eltype, length, ndims, size, indices, eachindex, stride
- construction
    - `A = Array{T, N}` construct array containing type `T` elements with `N` dimensions
    - `reshape(A, dims...)`, share same underlying data, `:` is inferred
    - `copy` and `deepcopy`
    - `similar(A, T, dims ...)`
    - `rand(T, dims...)`, an array with random i.i.d. uniformly distributed in [0,1]
    - `eye(T,n)` identity
    - `fill(x, dims...)`
- concatenation
    - `cat(k, A...)` concat along dimension k, dimension not in `dims` should have same size,
    - `vcat` i.e. [A; B; C; ...]   equiv `cat(1, A, ...Tokenize)`
    - `hcat` i.e. [A B C ...]  equiv `cat(2, A, ...)`
- typed array initializers
    - `T[A, B, C, ...]`
- comprehensions
    - `A = [F(x, y, ...) for x=rx, y=ry, ...]`
- generator expressions
    - comprehension written without brackets, yield a generator
    - `sum(1/n^2 for n=1:1000)`
- indexing
    - `X = A[I_1, I_2, ..., I_n]`
        - `:` select all indices in the dimension
        - `from:to:stride` select contiguous or strided sub selections
        - `end` used to represent last index of each dimension
            - `X[1:end-1]`, all but last item in 1d array
- assignment
    - `A[I_1, I_2, .., I_n] = X`
- indices types
    - scalar index, non-boolean
    - array of scalar index,
        - []
        - ranges
        - array <: AbstractArray
    - object representing array of scalar indices
- logical indexing
    - index by boolean array select elements where value of mask is `true`
    - `X[[false, true, true, false], :]`
- iteration
    - `for a in A ... end`
    - `for i in eachindex(A) ... end`
        - `i` here is `Int` of `A` is an array type with fast linear indicing. otherwise it is a `CartesianIndex`
- array and vectorized operation and functions
    - -, +, * , /, \ , ^, ==, !=, ≈,
        - element wise
    - dot syntax for vectorized operation for EVERY binary operator
        - `f.(args...)` , i.e. `sin.(x)`
- broadcasting
    - expands singleton dimensions in array arguments to match corresponding dimension without extra memory
    - applies functions element-wise
    ```julia
    a, A = rand(2,1), rand(2,3);
    repmat(a,1,3)+A  # too wasteful

    broadcast(+, a, A) # efficient
    ```
- implementation
    - `AbstractArray{T,N}`
        - `AbstractVector` and `AbstractMatrix` for 1-D and 2-D case
        - may have different undelrying structure,
        - but generally have `size` `getindex` `setindex!` ...
    - `DenseArray`
        - abstract subtype of `AbstractArray`
        - arrays that store elements contiguously in column major order
        - concrete types like `Vector` and `Matrix` are alises for 1-D and 2-D case
    - `SubArray`
        - performs indexing by sharing memory with original array, rather than copy
            - just store the index, and not data
        - created from `view`, same as calling `getindex`
            - but data left in place
            - `@views` will convert expression, code block, or any `array[...]` to a view instead
    - `BitArray`
        - space-efficient packed boolean arrays, 1 bit per boolean value
    - strided array
        - elements laid out in regular offsets
        - defines `stride(A)`
        - note
            - `DenseArray` is a strided array with stride=1
    - QR decomposition of small section of large array without temporaries and calling LAPACK function
    ```julia
    a = rand(10,10)
    b = view(a, 2:2:8, 2:2:4)
    (q, r) = qr(b)
    q
    r
    ```
