"""Exact checks for the two Region-B signed-kernel comparisons."""
import sympy as s

if not __debug__:
    raise RuntimeError("Run without -O: these exact checks require assertions.")

r,u,v,A,B,h=s.symbols('r u v A B h', positive=True)
ea,ed,em=s.symbols('ea ed em',positive=True)
def H(e):return ((e+1/e)-r*(e-1/e))/2
K=-H(ea)*H(ed)*H(em)
def polynomial(d,face):
    d=list(map(s.Rational,d))
    b=d[0]+(d[1]-d[3])/2
    rates=(d[2]+(d[1]+d[3])/2,(d[3]-d[1])/2,(d[1]+d[3])/2)
    def op(f):return sum(rate*e*s.diff(f,e) for rate,e in zip(rates,(ea,ed,em)))-2*r*b*f
    p=K-op(op(K))
    values=(u*v,v,v) if face=='F' else (u*v,v,u*v)
    return s.factor(p.subs(dict(zip((ea,ed,em),values))))
directions={3:('1/2',0,0,-1),7:('-1/2',0,0,-1),
            4:('1/2',0,-1,1),8:('-1/2',0,-1,1)}
eta=(r-1)*B-r-1
Q3=A*B*(B*(2*r**3-r**2-4*r+3)-2*r**3+r**2+2*r-1)+B*(-2*r**3-r**2+2*r+1)+2*r**3+r**2-4*r-3
Q4=A*B*(B*(2*r**3-3*r**2+1)-2*r**3+3*r**2+2*r-3)+B*(-2*r**3-3*r**2+2*r+3)+2*r**3+3*r**2-1
diffs={j:s.factor(polynomial(directions[j],'F')-polynomial(directions[j+4],'F')) for j in (3,4)}
for j,Q,sign in [(3,Q3,-1),(4,Q4,1)]:
    target=sign*r*eta*Q
    assert s.factor(4*u*v**3*diffs[j]-target.subs({A:u**2,B:v**2}))==0
astar=(r+1)/((r-1)*B)
map_b=1+2*h/(r-1)
assert s.factor(Q3.subs(A,1).subs(B,map_b)-4*(h-1)*((2*r+3)*h+r))==0
assert s.factor(Q3.subs(A,astar).subs(B,map_b)-4*(h-1)*(r+1))==0
primitive=s.integrate(r*eta*Q4,r)
J=s.factor(primitive-primitive.subs(r,1))
k=1-h
edge1=-2*(r-1)/3*(1+k+k**2+(2*r+1)**2*k**3)
edgemax=-2*(r-1)/15*(5+(4*r**2+12*r+14)*k+(16*r**2+28*r+11)*k**2)
assert s.factor(J.subs(A,1).subs(B,map_b)-edge1)==0
assert s.factor(J.subs(A,astar).subs(B,map_b)-edgemax)==0
eta_u=(r-1)*u**2-r-1
diag3=r*eta_u*((r-1)*(r+2)*u**2-(r-2)*(r+1))/(2*u**2)
diag4=-r**2*eta_u**2/(2*u**2)
for j,target in [(3,diag3),(4,diag4)]:
    p=(polynomial(directions[j],'G')-polynomial(directions[j+4],'G')).subs(v,1)
    assert s.factor(p-target)==0
print('PASS: direct directional kernels, four endpoint certificates, and both diagonal differences.')
