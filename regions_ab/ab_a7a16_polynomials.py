"""Check the face factorizations directly from the common kernel."""
import sympy as s

if not __debug__:
    raise RuntimeError("Run without -O: these exact checks require assertions.")

a,b,c,r,D=s.symbols('a b c r D', positive=True)
sh=lambda x:(x-1/x)/2
Y=sh(a); S=sh(a*b); h=sh(b); M=(a*b+1/(a*b))/(2*S)
Kf=-Y**2/S*(M-r)-2*Y*h*(M-r)**2-S*h**2*(M-r)**3
Kr=-S*Y*(M-r)**2-S**2*h*(M-r)**3
K=sh(a/c)/Y*Kf+sh(c)/Y*Kr
cancel=lambda v:s.cancel(s.together(v))
def Q(v,d):
    ky=d[1]+d[2]; kz=d[3]-d[1]; ke=d[1]
    return cancel(ky*a*s.diff(v,a)+s.Rational(kz,2)*b*s.diff(v,b)+ke*c*s.diff(v,c))
for j,d in [(11,(1,0,-1,0)),(14,(0,1,0,-1)),(16,(0,0,1,0))]:
    kl=s.Rational(d[0])-(s.Rational(d[3])-d[1])/2
    qK=Q(K,d)
    P=cancel(K-Q(qK,d)+4*r*kl*qK-4*r*r*kl*kl*K)
    EP=cancel(-qK.subs(r,M)*Q(M,d))
    beta=-Y/(S*h)
    base=(D-M)*(D-M+1)*(D-M-1)
    expectedF={11:4*S*h*h*base*(D-beta)**2,14:4*S*h*h*D*base*(D-beta),16:0}[j]
    expectedR={11:4*S*S*h*D*base*(D-beta),14:4*S*S*h*D*D*base,16:0}[j]
    for name,subs,expected in [('F',{c:1},expectedF),('R',{c:a},expectedR)]:
        actual=cancel(P.subs(subs).subs(r,M-D))
        assert cancel(actual-expected)==0,(j,name)
        eactual=cancel(EP.subs(subs))
        etarget=Y*Y/S**5 if j in (11,16) and name=='F' else 0
        assert cancel(eactual-etarget)==0,(j,name,'endpoint')
        print(j,name,'PASS',flush=True)
