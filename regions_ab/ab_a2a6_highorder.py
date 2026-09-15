"""Exact high-order certificate for A6, without numerical sampling."""
import sympy as s
from ab_a2a6_verify import check

if not __debug__:
    raise RuntimeError("Run without -O: these exact checks require assertions.")

t,z,X=s.symbols('t z X',positive=True)
rt=s.sqrt(2)
sh=lambda n:(z**n-z**(-n))/2
ch=lambda n:(z**n+z**(-n))/2
S=sh(1);C=ch(1);M=C/S
Dt=lambda f:s.diff(f,t)+z*s.diff(f,z)
G=(8-25*M**2+15*M**4-15*M*(M*M-1)**2*t)*S/(16*rt*M)
G1=s.factor(-S*S*Dt(G));G2=s.factor(-S*S*Dt(G1))
P=193-430*M*M+240*M**4+t*(-315*M+510*M**3-240*M**5)
Wder=(-45*t/S+Dt(P/S))/(384*rt)
check('A2 endpoint primitive derivative',Wder+1/(2*rt*C)+(4*M*M-1)*G/S**2)
E=s.factor(((25-32*X**2-40*ch(2)+15*ch(4)+24*X*sh(2))*G
    +4*(8*X-3*sh(2))*G1-8*G2)/(2*S**6))
a,b=s.symbols('a b')
mom=s.symbols('m0:4')
ee=s.exp(-2*X*M)
def Dfull(f):
    return (Dt(f)-s.diff(f,a)/(2*C)-s.diff(f,b)/(2*S)
        -sum(s.diff(f,mi)*M**i*ee*G/S**2 for i,mi in enumerate(mom)))
U=2*(-C*(mom[0]+2*mom[2])+S*(mom[1]+2*mom[3]))+2*ee*G/S**3+(2*a*C-4*b*S+S/C)/rt
BD=S*(-2*Dfull(Dfull(U))-10*M*Dfull(U)+12*U)
BE=-20*(mom[1]+2*mom[3])+ee*E-2*rt*(-10*b+5/C+1/C**3)
check('A6 differential operator identity',BD-BE)
Nt=-2*rt*(5/S-5*S/C**2-3*S/C**4)
J=rt*(Dt(E)+2*X*E/S**2+20*M*(1+2*M*M)*G/S**2+s.exp(2*X*M)*Nt)
p0=120*t*C-80*S-15*sh(3)+sh(5)
p1=120*t*(26+15*ch(2))-2415*sh(2)-12*sh(4)-7*sh(6)
p2=3240*t+4320*t*ch(2)+1080*t*ch(4)+256*C*(4+ch(2))*S**7+3*(-584*sh(2)-394*sh(4)-24*sh(6)+sh(8))
p3=-360*t+256*sh(2)-40*sh(4)+sh(8)
expected={0:-5*p0/(2*S**8),1:p1/(4*rt*S**9),2:-p2/(8*S**10*C**3),3:-rt*p3/S**11}
for k in range(4):
    check('A6 J derivative boundary '+str(k),rt**k*s.diff(J,X,k).subs(X,0)-expected[k])
check('A6 fourth derivative',rt**4*s.diff(J,X,4)+256*s.exp(2*X*M)*(4+ch(2))/S**5)
def deriv(f,n):
    for _ in range(n):f=s.expand(Dt(f))
    return f
check('p0 fifth',deriv(p0,5)-10*S*(12*t-104*sh(2)+625*sh(4)))
check('p1 fifth',deriv(p1,5)+192*S*(-600*t*C+695*sh(3)+567*sh(5)))
check('p2 eighth',deriv(p2,8)-1280*(864*t*(ch(2)+64*ch(4))+3086*sh(2)+54784*sh(4)-137781*sh(6)+65536*sh(8)+78125*sh(10)))
check('p3 first',Dt(p3)-512*(3+ch(2))*S**6)
for name,p,n in [('p0',p0,5),('p1',p1,5),('p2',p2,8),('p3',p3,1)]:
    initial=[s.factor(deriv(p,k).subs({t:0,z:1})) for k in range(n)]
    expected_init=[0]*n
    if name=='p2':expected_init[-1]=5898240
    assert initial==expected_init,(name,initial)
    print('PASS',name,'initial derivatives',initial,flush=True)
