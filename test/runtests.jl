using ApproxFun, RatFun, Test


x = Fun()
f = cos(x) // exp(x)

@test f(0.1) ≈ cos(0.1)/exp(0.1)
@test 2f isa RationalFun
@test (2f)(0.1) ≈ 2f(0.1)
@test inv(f)(0.1) ≈ exp(0.1)/cos(0.1)
@test (f*f)(0.1) ≈ f(0.1)^2
@test (f/f)(0.1) ≈ 1
