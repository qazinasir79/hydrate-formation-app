import numpy as np, json, sys, math
sys.path.insert(0,'.')
import fitv13 as F
d = json.load(open('params_v8.json'))
p = np.array([d['dK0'][g] for g in F.G5] + [d['dK1'][g] for g in F.G5])
r = F.resid(p)
n = 0
print('== anchors ==')
for (g,T,P,w,st) in F.ANCHS:
    print(f'{g:5s} T={T:7.2f} st={st} w={w} r={r[n]:+.3f}'); n+=1
print('== mix ==')
for (T,gf,P,st) in F.MIX:
    print(f'mix T={T:.2f} r={r[n]:+.3f}'); n+=1
print('== sol ==')
for (g,T,P,x,w) in F.SOL:
    print(f'{g:5s} T={T:.2f} P={P/1e6:.2f} r={r[n]:+.3f}'); n+=1
print('cost', 0.5*float(np.sum(r**2)))
