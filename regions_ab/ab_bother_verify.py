"""Exact Region B trace, paired-control, and reversal identities.

The companion modules are imported from this script's directory.  No notebook,
manuscript source, numerical sampling, or working-directory data is required.
"""
from ab_bother_kernel import s, r, u, v, P, dirs

if not __debug__:
    raise RuntimeError("Run without -O: these exact checks require assertions.")

def same(a,b=0):
    assert s.cancel(a-b)==0, s.factor(a-b)
pf={j:P(dirs[j]) for j in ('B5','B6','B7','B8','B9','B10')}
pg={j:P(dirs[j],'G') for j in pf}
cy=(u+1/u)/(u-1/u); sy=2/(u-1/u)
e1={j:u*s.diff(pf[j],u)-v*s.diff(pf[j],v)-2*cy*pf[j]+2*sy*pg[j] for j in pf}
e2={j:u*s.diff(pg[j],u)-2*cy*pg[j]+2*sy*pf[j] for j in pf}
same(pf['B5'],pf['B6']); same(pf['B9'],pf['B10']);same(pg['B8'],pg['B9'])
same(e1['B5']+e1['B6']);same(e1['B9']+e1['B10'])
same(e1['B7']);same(e1['B8'])
mt=(u*u*v*v+1)/(u*u*v*v-1);st=2*u*v/(u*u*v*v-1)
for j in ['B6','B7','B10']:
    same(e2[j])
    same(-st**2*pg[j].subs(r,mt)+s.Rational(1,2)*st**5*(u-1/u)/2)
print('PASS correction polynomials, endpoint terms, and trace identities')

from ab_bother_fourrational import T, M, L, A, B, co, si, W, V, dt, der, res, es
ff={j:res(e).subs(M,L) for j,e in es.items()}
gg={j:res(e).subs(M,T) for j,e in es.items()}
ee1={j:-L*s.diff(ff[j],L)-2*co(T/L)/si(T/L)*ff[j]+2*gg[j]/si(T/L) for j in ff}
ee2={j:dt(gg[j])-2*co(T/L)/si(T/L)*gg[j]+2*ff[j]/si(T/L) for j in ff}
same(ff['B5'],ff['B6']);same(ff['B9'],ff['B10']);same(gg['B8'],gg['B9'])
same(ee1['B5']+ee1['B6']);same(ee1['B9']+ee1['B10'])
same(ee1['B7'],-1/s.sqrt(2));same(ee1['B8'],-1/s.sqrt(2))
for j in ['B6','B7','B10']:same(ee2[j])
print('PASS four-expert errors and trace identities')

Delta=V-W
expected=si(L)*(co(M)-2*B*si(T)*si(M))/s.sqrt(2)
same(Delta,expected)
same(Delta-der(der(Delta,(0,1,0)),(0,1,0)),si(L)*si(T/M)/(s.sqrt(2)*si(T)))
for dlt in [(1,0,-1),(-1,1,-1)]:same(Delta-der(der(Delta,dlt),dlt))
print('PASS all direct reversal identities')
