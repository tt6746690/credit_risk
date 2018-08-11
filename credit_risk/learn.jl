
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



"""
Constructor
```
type Foo
    bar
    baz
end
foo = Foo(1,2)
```
- outer constructor methods
    - convenience in constructing objects
    - create `Foo` with 0 or 1 argument
    - `Foo(x) = Foo(x, x)`
    - note this once is called `Foo()=  Foo(0)` instead of default arg
- inner constructor methods
    -



"""

"""
Type System
- type system
    - concrete types
        - may not subtype each other
        - final, and only have abstract type as their supertypes
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
    ```
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
    ```
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
"""



"""
Multi-dimensional Arrays
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
        - `from:to:stride` select contiguous or strided subselections
        - `end` used to represnet last index of each dimension
            - `X[1:end-1]`, all but last item in 1d array
- assignment
    - `A[I_1, I_2, .., I_n] = X`
- indices types
    - scalaar index, nonboolean
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
    - -, +, *, /, \, ^, ==, !=, ≈,
        - elementwise
    - dot syntax for vecotrized operation for EVERY binary operator
        - `f.(args...)` , i.e. `sin.(x)`
- broadcasting
    - expands singleton dimensions in array arguments to match corresponding dimension wihtout extra memory
    - applies functions elementwise
    ```
    a, A = rand(2,1), rand(2,3);
    repmat(a,1,3)+A  # too wasteful

    broadcast(+, a, A) # effifcient

    ```

"""
