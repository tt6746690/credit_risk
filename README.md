# credit_risk


#### Setup 

```
# 1. Install julia v0.7 @ https://julialang.org/downloads/platform.html

# 2. Open julia repl
$ julia

# 3. Enter package mode
julia> ]

# 4. Activate package 
(v0.7) pkg> activate .

# 5. Download dependencies
(CreditRisk) pkg> instantiate

# 6. Update package and precompile modules
(CreditRisk) pkg> update; precompile

# 7. Back to julia repl, and start hacking!
julia> using CreditRisk
```

#### Layout

```
.
├── data            # pre-saved `Parameter
├── docs            # documentation 
├── notebook        # jupyter notebook for easier coding
├── paper           # papers
├── scripts         # scripts that runs julia code
├── src             # source directory  
└── test            # test
```


#### notes

+ https://github.com/JuliaIO/JLD.jl/pull/227
