"""Exact verification of the one-variable and A9/A10 certificates."""
import sympy as s

if not __debug__:
    raise RuntimeError("Run without -O: these exact checks require assertions.")

a,t,xi,r=s.symbols('a t xi r', positive=True)
A,B=s.symbols('A B')
sh=(a-a**-1)/2; ch=(a+a**-1)/2; M=ch/sh
shn=lambda n:(a**n-a**-n)/2
chn=lambda n:(a**n+a**-n)/2
J=s.symbols('J0:4')
gM=(8-25*M**2+15*M**4-15*M*(M**2-1)**2*t)/(16*s.sqrt(2)*M/sh)
em=s.symbols('em', positive=True)
def simp(v): return s.cancel(s.together(v))
def Dt(v):
    v=s.diff(v,t)+a*s.diff(v,a)-s.diff(v,A)/(2*ch)-s.diff(v,B)/(2*sh)+s.diff(v,em)*2*xi/sh**2*em
    return simp(v)
def D(v):
    ans=Dt(v)
    for k,j in enumerate(J): ans+=s.diff(v,j)*(-1/sh**2)*em*gM*M**k
    return simp(ans)
p=2*((-1+6*r*r)*ch+r*(9-14*r*r)*sh)
F=sum(s.expand(p).coeff(r,k)*J[k] for k in range(4))+2/sh**3*em*gM+(2*A*ch-4*B*sh+sh/ch)/s.sqrt(2)
LF=simp(-D(D(F))-5*M*D(F)+6*F)
E=simp(D(sh*LF)/em)
print('Constructed the exact one-variable differential expression.',flush=True)
def DX(v): return simp(s.diff(v,xi)-2*M*em*s.diff(v,em))
def DXn(v,n):
    for _ in range(n): v=DX(v)
    return v
assert all(s.diff(E,j)==0 for j in J)
B3=-360*t+256*shn(2)-40*shn(4)+shn(8)
B2=-600*t*ch+435*sh+54*shn(3)+2*shn(5)-shn(7)
B1=-120*t*(2+chn(2))+185*shn(2)-4*shn(4)+shn(6)
B12=-120*t*ch+120*sh-5*shn(3)+3*shn(5)
C2=300*t-165*shn(2)-3*shn(4)+7*shn(6)
checks={
 'E0': E.subs({xi:0,em:1}),
 'E1':DXn(E,1).subs({xi:0,em:1})+3*B1/(16*s.sqrt(2)*sh**9),
 'E2':DXn(E,2).subs({xi:0,em:1})-B2/(4*s.sqrt(2)*sh**10),
 'E3':DXn(E,3).subs({xi:0,em:1})+B3/(4*s.sqrt(2)*sh**11),
 'E4':DXn(E,4)+16*s.sqrt(2)/em*(4+chn(2))/sh**5,
 'B3':Dt(B3)-512*(3+chn(2))*sh**6,
 'B2':Dt(B2)+2*sh*C2,
 'C2':Dt(C2)-96*(13+7*chn(2))*sh**4,
 'B1':Dt(B1)-4*sh*B12,
 'B12':Dt(B12)-30*sh*(-4*t+shn(4)),
}
for k,v in checks.items():
    assert simp(v)==0,k
    print(k,'PASS',flush=True)

# A9 transport forcing H=-G_t/2+coth(t)G-csch(t)F.
p9=1+7*r*r-8*r**4+(-3+8*r*r)*((1+r*r)*chn(2)-2*r*shn(2))
# Its integral coefficients reduce to 12 r(r^2-1), and the rest is explicit.
assert simp(-Dt(p9)/2+M*p9-p/sh-12*r*(r*r-1))==0
Hrest=(4*B-2/ch)/s.sqrt(2)-4/sh**4*em*gM
Ht=simp(Dt(Hrest)-12*M*(M*M-1)/sh**2*em*gM)
N=-5*(19*xi+24*t)*ch-30*(xi+4*t)*chn(2)*ch+30*sh+65*shn(3)+3*shn(5)+xi/sh*(120*t*ch**2+ch*shn(5))-32*sh**7/em
assert simp(Ht-em*N/(16*s.sqrt(2)*sh**8*ch**2))==0
Q=simp(sh*DX(N).subs({xi:0,em:1}))
V=simp(120*t*ch-80*sh-15*shn(3)+shn(5))
assert simp(N.subs({xi:0,em:1})+2*ch**2*V)==0
assert simp(DXn(N,2)+128*ch**2*sh**5/em)==0
Q7=Q
V7=V
for _ in range(7): Q7=Dt(Q7); V7=Dt(V7)
assert simp(Q7-64*(339*chn(2)-3584*chn(4)+15309*chn(6)-16384*chn(8)+120*t*shn(2)))==0
assert simp(V7-5*(152*ch-6561*chn(3)+15625*chn(5)+24*t*sh))==0
for k in range(7):
    assert simp(Q.subs({a:1,t:0}))==0
    assert simp(V.subs({a:1,t:0}))==0
    Q=Dt(Q); V=Dt(V)
print('A9 transport certificates PASS',flush=True)

# A16 and the integrated A7/A8 endpoint certificate.
K=-60*t+45*shn(2)-9*shn(4)+shn(6)
assert simp(Dt(K)-192*sh**6)==0
reduced=gM+sh**4/(2*s.sqrt(2)*ch)
V16=120*t*ch-115*sh+6*shn(3)-6*shn(5)+shn(7)
C16=60*t-5*shn(2)-23*shn(4)+7*shn(6)
assert simp(Dt(reduced)-3*V16/(128*s.sqrt(2)*sh**4))==0
assert simp(Dt(V16)-2*sh*C16)==0
assert simp(Dt(C16)-32*(19+21*chn(2))*sh**4)==0
L=s.symbols('L') # Li_2(e^-2t)-4[t b(t)+Li_2(e^-t)], L'=2t/sinh t.
J78=691835/sh**7+50400*L+98304*24*B*M+(15*t*(-5053*ch+2681*chn(3)-805*chn(5)+105*chn(7))-377077*shn(3)+120695*shn(5)-16857*shn(7))/sh**8
DJ78=Dt(J78)+s.diff(J78,L)*2*t/sh
DDJ78=Dt(DJ78)
target=786432*a*a*(1+a*a)/(a*a-1)**10*((a*a-1)**7*24*B+36*a**7*K)
assert simp(DDJ78-target)==0
print('A16 and A7/A8 endpoint derivatives PASS',flush=True)
