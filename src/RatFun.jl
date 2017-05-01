module RatFun
using Base,Compat,ApproxFun,RecipesBase

import ApproxFun:evaluate,dimension,domain,setdomain,PiecewiseSpace,DiracSpace,PointSpace,∞
import Base:+,-,*,/,.*,.-,./,.+,getindex

export RationalFun, inv

immutable RationalFun{S1,T1,S2,T2}
    p::Fun{S1,T1}
    q::Fun{S2,T2}
end

domain(r::RationalFun) = domain(r.p)
setdomain(r::RationalFun,d) = RationalFun(setdomain(r.p,d),setdomain(r.q,d))

function evaluate(r::RationalFun,x)
    (r.p)(x)/(r.q)(x)
end

@compat (r::RationalFun)(x) = evaluate(r,x)

for op = (:*,:.*)
    @eval $op(r1::RationalFun,r2::RationalFun)=RationalFun($op(r1.p,r2.p),$op(r1.q,r2.q))
    @eval $op(r::RationalFun,a::Union{Number,Fun}) = RationalFun($op(r.p,a),r.q)
    @eval $op(a::Union{Number,Fun},r::RationalFun) = RationalFun($op(a,r.p),r.q)
end

Base.inv(r::RationalFun) = RationalFun(r.q,r.p)

Base.vec(r::RationalFun) = RationalFun.(vec(r.p),vec(r.q))

(./)(r1::RationalFun,r2::RationalFun)=r1.*inv(r2)
(./)(a,r::RationalFun)=a.*inv(r)

(/)(r1::RationalFun,r2::RationalFun)=r1*inv(r2)
(/)(a,r::RationalFun)=a*inv(r)

(/)(r::RationalFun,a)=RationalFun(r.p,r.q*a)
(./)(r::RationalFun,a)=r./a

for op = (:+,:.+,:-,:.-)
  @eval $op(r1::RationalFun,r2::RationalFun) = RationalFun($op((r1.p*r2.q),(r2.p*r1.q)),r1.q*r2.q)
end

Base.convert(::Type{Fun},r::RationalFun) = r.p/r.q

# The padding in this function can be improved
# No support for functions with poles within the domain
function plotptsvals(r::RationalFun)
    p = r.p
    q = r.q
    plen = ncoefficients(p)
    qlen = ncoefficients(q)
    if dimension(space(p)) == ∞ && dimension(space(q)) == ∞
        p=pad(p,10plen+10*qlen+1000)
        q=pad(q,10plen+10*qlen+1000)
        r = RationalFun(p,q)
    elseif dimension(space(p)) == ∞ && dimension(space(q)) < ∞
        p=pad(p,10plen+10*dimension(space(q))+1000)
        q=pad(q,10plen+10*dimension(space(q))+1000)
        r = RationalFun(p,q)
    elseif dimension(space(q)) == ∞ && dimension(space(p)) < ∞
        q=pad(q,dimension(space(p))+10qlen+1000)
        p=pad(p,dimension(space(p))+10qlen+1000)
        r = RationalFun(p,q)
    else
        p=pad(p,dimension(space(p))+dimension(space(q)))
        r = RationalFun(p,q)
    end
    return points(r.p),values(r.p)./values(r.q)
end

@recipe function f{S,T<:Real}(g::RationalFun{S,T})
    plotptsvals(g)
end

@recipe function f{S,T<:Real}(x::AbstractVector{T},g::RationalFun{S,T})
    x,g(x)
end


@recipe function f{S1<:PiecewiseSpace,S2<:PiecewiseSpace,T1<:Real,T2<:Real}(r::RationalFun{S1,T1,S2,T2})
    vp = vec(r.p)
    vq = vec(r.q)

    if !isempty(vp)
        @series begin
            primary --> true
            RationalFun(vp[1],vq[1])
        end
    end

    for k=2:length(vp)
        @series begin
            primary := false
            RationalFun(vp[k],vq[k])
        end
    end
end

# For dirac space, we draw a dotted line extending to infinity
@recipe function f{S1<:DiracSpace,T1<:Real,S2<:PointSpace,T2<:Real}(r::RationalFun{S1,T1,S2,T2})
    p = r.p
    q = r.q
    pts=space(p).points
    n=length(pts)
    ws=pad(p.coefficients./q.coefficients,length(pts))
    @series begin
        primary --> true
        ones(2)*pts[1],[0,1]*ws[1]
    end

    if length(ws) > 1
        @series begin
            primary := false
            ones(2)*pts[2:end]',[0,1]*ws[2:end]'
        end
    end

    @series begin
        primary := false
        linestyle := :dot
        ones(2)*pts',[1,2]*ws'
    end
end

# for PointSpace, we draw just a line
@recipe function f{S1<:PointSpace,T1<:Real,S2<:PointSpace,T2<:Real}(r::RationalFun{S1,T1,S2,T2})
    p = r.p
    q = r.q
    pts=space(p).points
    n=length(pts)
    ws=pad(p.coefficients./q.coefficients,length(pts))

    @series begin
        primary --> true
        ones(2)*pts[1],[0,1]*ws[1]
    end

    @series begin
        primary := false
        ones(2)*pts[2:end]',[0,1]*ws[2:end]'
    end
end

end #module
