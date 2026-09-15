"""Exact correction-kernel definitions used by ab_bother_verify.py."""
import sympy as s
r,u,v=s.symbols('r u v',positive=True)
# u=e^y, v=e^(z/2)
def H(e): return ((e+1/e)-r*(e-1/e))/2
def J(e): return ((e-1/e)-r*(e+1/e))/2
def P(d,face='F'):
    d=list(map(s.Rational,d)); b=d[0]+(d[1]-d[3])/2
    tt=d[2]+(d[1]+d[3])/2; dd=(d[3]-d[1])/2; mm=(d[1]+d[3])/2
    exps=(u*v,v,v) if face=='F' else (u*v,v,u*v)
    h0,h1,h2=map(H,exps); j0,j1,j2=map(J,exps)
    k=h0*h1*h2
    dk=tt*j0*h1*h2+dd*h0*j1*h2+mm*h0*h1*j2
    return s.expand((tt*tt+dd*dd+mm*mm-1+4*r*r*b*b)*k+2*tt*dd*j0*j1*h2+2*tt*mm*j0*h1*j2+2*dd*mm*h0*j1*j2-4*r*b*dk)
dirs={'A7':(0,-1,0,0),'B7':('-1/2',0,0,-1),'A8':(-1,1,-1,0),'B8':('-1/2',0,-1,1),
'A6':(0,0,0,1),'B6':('-1/2',1,0,0),'A5':(0,0,1,-1),'B5':('1/2',-1,1,0),
'A10':(-1,0,0,1),'B10':('-3/2',1,0,0),'A9':(-1,0,1,-1),'B9':('-1/2',-1,1,0)}
