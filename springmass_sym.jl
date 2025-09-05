using Symbolics

@variables m::Real,c::Real, k::Real
@variables x::Complex, t::Real
@variables λ1::Complex, λ2::Complex
@variables x0::Real

λ1 = -c + im*sqrt(4*m*k-c^2)
λ2 = -c - im*sqrt(4*m*k-c^2)
C1 = -λ2/(λ1-λ2)*x0
C2 = λ1/(λ1-λ2)*x0
x = C1*exp(λ1*t) + C2*exp(λ2*t)
println(simplify(x))
