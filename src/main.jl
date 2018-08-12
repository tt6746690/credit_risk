include("CreditRisk.jl")

module tst
import Main.CreditRisk: Parameter


n = 20
c = 4
s = 10
l = 0.2

x = ones(n, c)
y = ones(n, s)
z = ones(n)

# p = Parameter(n,c,s,0.2,x,z,x,y,z,x)
p = Parameter(n,c,s,l)
print(summary(p))
print(p)

end
