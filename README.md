# RatFun.jl
RatFun is a package for working with functions expressible as a Fun divided by another Fun. The numerator and denominator may be different types of Funs.

## Example usage

Define two functions on an interval using ApproxFun:

```using ApproxFun
using RatFun
p = Fun(x->exp(-10x^2))
q = Fun(x->-1+.5x*sin(100x^2)
```

You can define a RationalFun, representing p/q, as follows.

```r = p // q
```

We use Plots.jl to plot it.

```using Plots
plot(plot(p,title="p"),plot(q,title="q"),plot(r,title="r=p/q"),layout=grid(3,1))
```

<img src=images/ExampleRatFun1.png>

You can also do the same for Funs with more exotic types, such as those with Dirac deltas, as follows.

```p2 = .7*DiracDelta(-1.25) - .9*DiracDelta(.5) + p
q2 = .5*KroneckerDelta(-1.25) - .7*KroneckerDelta(.5) + q
r2 = p2 // q2
plot(plot(p2,title="p2"),plot(q2,title="q2"),plot(r2,title="r2=p2/q2"),layout=grid(3,1))
```

<img src=images/ExampleRatFun2.png>

RationalFuns with Dirac deltas in the numerator are used in the SpectralMeasures package.
