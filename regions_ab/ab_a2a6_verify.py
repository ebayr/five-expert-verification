"""Exact symbolic identities for Region A, controls 2--6.

This checks the residual kernels, boundary terms, face comparisons, and
hyperbolic trace identities used in the proof.  No numerical sampling is used.
"""
import sympy as s

if not __debug__:
    raise RuntimeError("Run without -O: these exact checks require assertions.")

A,B,C,r=s.symbols('A B C r', positive=True)
rt=s.sqrt(2)
sh=lambda z:(z-1/z)/2
ch=lambda z:(z+1/z)/2
T=A*B**2*C
M=ch(T)/sh(T)
H=lambda z:ch(z)-r*sh(z)
K=-H(T)*H(C/A)*H(A*C)
V=(s.atan(1/T)*ch(C/A)*ch(T)*ch(A*C)
   +s.atanh(1/T)*sh(C/A)*sh(T)*sh(A*C)-sh(C**2)/2)/rt
directions={2:(-1,0,0,0),3:(1,-1,0,0),4:(0,1,-1,0),5:(0,0,1,-1),6:(0,0,0,1)}
def diff(expr,d):
    return sum(di*z*s.diff(expr,z) for di,z in zip(d[1:],(A,B,C)))
def parts(j):
    d=directions[j]; L=2*d[0]+d[1]-d[3]
    DK=diff(K,d)
    P=s.factor(K-diff(DK,d)/4+L*r*DK-L*L*r*r*K)
    E=s.factor(-diff(M,d)*DK.subs(r,M)/4)
    N=s.expand(V-diff(diff(V,d),d)/4)
    return P,E,N
def check(name,expr):
    result=s.simplify(s.factor(expr.rewrite(s.exp)))
    assert result==0,(name,result)
    print('PASS',name,flush=True)

if __name__=='__main__':
    p={j:parts(j) for j in directions}
    check('A3/A4 face integral equality',(p[3][0]-p[4][0]).subs(A,1))
    check('A3/A4 face endpoint equality',(p[3][1]-p[4][1]).subs(A,1))
    check('A4/A6 endpoint equality',p[4][1]-p[6][1])
    check('A5/A6 endpoint equality',p[5][1]-p[6][1])
    t=s.symbols('t',positive=True)
    face={A:1,B:s.exp(t/2),C:1}
    check('A6 first face polynomial',p[6][0].subs(face)-(1+2*r*r)*(-s.cosh(t)+r*s.sinh(t))/4)
    check('A6 first face endpoint',p[6][1].subs(face)-1/(4*s.sinh(t)**3))
    Nfirst=(2*s.atan(s.exp(-t))*s.cosh(t)-4*s.atanh(s.exp(-t))*s.sinh(t)+s.tanh(t))/(8*rt)
    check('A6 first face N4',p[6][2].subs(face)-Nfirst)
    diag={A:s.exp(t/2),B:1,C:s.exp(t/2)}
    expected=(3*(r*r-1)+(1-7*r*r)*s.cosh(2*t)+2*(r+2*r**3)*s.sinh(2*t))/8
    # Replace hyperbolic functions by exponentials for the exact comparison.
    check('A6 second face polynomial',s.expand((p[6][0].subs(diag)-expected).rewrite(s.exp)))
    Ndiag=(s.atan(s.exp(-t))*(3-s.cosh(2*t))+5*s.sinh(t)-4*s.atanh(s.exp(-t))*s.sinh(2*t))/(8*rt)
    check('A6 second face N4',p[6][2].subs(diag)-Ndiag)
    for j in (4,5):
        expected=-(r*r-1)*(1-s.cosh(2*t)+r*s.sinh(2*t))/2
        check('A%d/A6 diagonal polynomial'%j,s.expand(((p[j][0]-p[6][0]).subs(diag)-expected).rewrite(s.exp)))
        h=(s.exp(2*t)-1)*s.atan(s.exp(-t))+(s.exp(2*t)+1)*s.atanh(s.exp(-t))-2*s.exp(t)
        check('A%d/A6 diagonal N4'%j,(p[j][2]-p[6][2]).subs(diag)-(1-s.exp(-2*t))*h/(4*rt))
    a,b,c=s.symbols('a b c',positive=True)
    F=(8*a**4*c+(b+c)**2*(8*b*(1+b)*(2+b)+(8+b*(32+19*b))*c+2*(4+7*b)*c**2+3*c**3)
       +2*a*(b+c)*(4*b*(2+3*b*(2+b))+(16+b*(52+35*b))*c+2*(14+17*b)*c**2+11*c**3)
       +8*a**3*(2*c+(b+c)*(b+4*c))+a*a*(8*c+(b+c)*(24*b*(1+b)+(64+67*b)*c+43*c*c)))
    # A=1, C^2=1+a, B^4*C^2=1+a+b, and r=(2+a+b+c)/(a+b+c).
    trans={A:1,C:s.sqrt(1+a),B:((1+a+b)/(1+a))**s.Rational(1,4),r:(2+a+b+c)/(a+b+c)}
    scale=1/(32*s.sqrt(1+a+b)*(1+a))
    check('A3 quintic certificate',p[3][0].subs(trans)+scale*8*F/(a+b+c)**5)
    check('A3 endpoint certificate',p[3][1].subs(trans)-scale*64*b*b*(1+a+b)**2/(a+b)**5)
    integral=sum(s.expand(16*F).coeff(c,k)*s.factorial(k)*s.factorial(5-k)/(s.factorial(6)*(a+b)**(6-k)) for k in range(6))
    NN=(165*a**4+6*a**3*(28+115*b)+10*b*b*(4+b*(8+9*b))+4*a*b*(44+5*b*(26+27*b))+a*a*(56+b*(608+975*b)))
    check('A3 endpoint domination',integral-64*b*b*(1+a+b)**2/(a+b)**5-4*NN/(15*(a+b)**5))
    D=s.factor(((p[5][0]-p[6][0])/(r*r-1)).subs(A,1))
    check('A5 comparison endpoint 1',D.subs(r,1)+(C*C-1)/(B*B*C**3))
    check('A5 comparison endpoint M',D.subs(r,M.subs(A,1))+2*sh(B*B)*sh(C)/sh(B*B*C)**2)
    check('A5 comparison convexity',s.diff(D,r,2)-(C*C-1)**2*(1+B**4*C*C)/(2*B*B*C**3))
    HH=(T*T-1)*s.atan(1/T)+(T*T+1)*s.atanh(1/T)-2*T
    check('A5 N4 comparison',(p[5][2]-p[6][2]).subs(A,1)-(C**4-1)*HH.subs(A,1)/(4*rt*B*B*C**3))
    u,v=s.symbols('u v',positive=True)
    cot=ch(u*u)/sh(u*u); csc=1/sh(u*u)
    upper=ch(u*u*v)/sh(u*u*v)
    DY=lambda f:u*s.diff(f,u)/2
    DZ=lambda f:v*s.diff(f,v)/2
    for j in (2,3,6):
        PF,EF,NF=[s.factor(f.subs({A:1,B:u,C:v})) for f in p[j]]
        PG,EG,NG=[s.factor(f.subs({A:u,B:1,C:u*v},simultaneous=True)) for f in p[j]]
        check('A%d second face E zero'%j,EG)
        check('A%d second hyperbolic polynomial'%j,DY(PG)-2*cot*PG+2*csc*PF)
        check('A%d second hyperbolic endpoint'%j,DY(upper)*PG.subs(r,upper)+2*csc*EF)
        check('A%d second hyperbolic N4'%j,DY(NG)-2*cot*NG+2*csc*NF)
        if j in (2,6):
            check('A%d first hyperbolic polynomial'%j,DY(PF)-2*DZ(PF)-2*cot*PF+2*csc*PG)
            check('A%d first hyperbolic endpoint'%j,DY(EF)-2*DZ(EF)-2*cot*EF+2*csc*EG)
            check('A%d first hyperbolic N4'%j,DY(NF)-2*DZ(NF)-2*cot*NF+2*csc*NG)
