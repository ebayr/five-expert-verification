"""Exact four-expert residual definitions used by ab_bother_verify.py."""
import sympy as s
T,M,L,A,B=s.symbols('T M L A B',positive=True)
def co(x):return (x+1/x)/2
def si(x):return (x-1/x)/2
W=(A*co(L)*co(T)*co(M)+B*si(L)*si(T)*si(M)-si(M*L)/2)/s.sqrt(2)
V=W.subs(L,1/L)
def dt(f): return T*s.diff(f,T)-s.diff(f,A)/(2*co(T))-s.diff(f,B)/(2*si(T))
def dm(f): return M*s.diff(f,M)
def dl(f): return L*s.diff(f,L)
def der(f,e):
 a,b,c=map(s.Rational,e)
 return ((a+c)/2+b)*dt(f)+(a+c)*dm(f)/2+(c-a)*dl(f)/2
def res(e):return s.cancel(V-der(der(V,e),e))
es={'B4':(0,-1,1),'B5':(-1,1,0),'B6':(1,0,0),'B7':(0,0,-1),'B8':(0,-1,1),'B9':(-1,1,0),'B10':(1,0,0)}
