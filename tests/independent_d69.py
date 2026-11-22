#!/usr/bin/env python3
"""Independent D69 rational-distance verifier, no import of primary Gram generator."""
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
from pathlib import Path
import json,hashlib
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'results'/'certificate_d69.json'
cert=json.loads(p.read_text(encoding='utf-8'))
parent_file=ROOT/'prior'/'d68'/'results'/'certificate_d68.json'
parent=json.loads(parent_file.read_text(encoding='utf-8'))
archive=ROOT/'prior'/'d68'/'frozen_scientific_tag_v0.1.0.zip'
assert hashlib.sha256(archive.read_bytes()).hexdigest()=='a9f58e6f88d057158b1ce5ba4f4e5e9fb3f2aaeddae2493447e396724e36e63a'
from zipfile import ZipFile
with ZipFile(archive) as archived_parent:
    assert archived_parent.testzip() is None
    assert archived_parent.read('results/certificate_d68.json')==parent_file.read_bytes()
assert hashlib.sha256(parent_file.read_bytes()).hexdigest()=='16bdab4b4166e40817dc38561f1cc14c4f7cba994313988bc85879d6c7584b13'
core=[
    (F(0),F(0),F(1)),(F(1,2),F(0),F(1,2)),
    (F(-1,4),F(3,4),F(1,2)),(F(-1,4),F(-3,4),F(1,2)),
    (F(0),F(0),F(-1)),(F(1,2),F(0),F(-1,2)),(F(-1,2),F(0),F(-1,2))
]
def d2(i,j):
    if i==j:return F(0)
    if i<7 and j<7:
        a,b,c=(core[i][k]-core[j][k] for k in range(3))
        return 3*a*a+b*b+c*c
    if 7<=i<12 and 7<=j<12:return F(5,2)
    return F(2)

N=74
pairs=list(combinations(range(N),2))
computed={(i+1,j+1):d2(i,j) for i,j in pairs}
supplied={(int(e['i']),int(e['j'])):F(e['d2']) for e in cert['labelled_squared_distances']}
assert cert['schema']=='POLYBRUTUS-D69-EXACT-CERTIFICATE-v1' and cert['status']=='PASS'
assert len(computed)==len(supplied)==2701 and computed==supplied
old={(int(e['i']),int(e['j'])):F(e['d2']) for e in parent['labelled_squared_distances']}
assert len(old)==2628 and all(computed[k]==v for k,v in old.items())
assert all(d2(i,73)==F(2) for i in range(73))
expected=Counter({F(4):6,F(7):2,F(8):2670,F(9):3,F(10):10,F(12):6,F(13):2,F(16):2})
assert Counter(4*d2(i,j) for i,j in pairs)==expected
assert sum(computed.values(),F(0))==F(21655,4)

D=s.Matrix(N,N,lambda i,j:s.Rational(d2(i,j).numerator,d2(i,j).denominator))
H=s.eye(N)-s.ones(N)/N
assert DomainMatrix.from_Matrix(-H*D*H/2).to_field().rank()==69
oldD=D[:73,:73];Ho=s.eye(73)-s.ones(73)/73
assert DomainMatrix.from_Matrix(-Ho*oldD*Ho/2).to_field().rank()==68

L=s.zeros(N)
for i,j in pairs:
    w=1/D[i,j]
    L[i,i]+=w;L[j,j]+=w;L[i,j]-=w;L[j,i]-=w
assert DomainMatrix.from_Matrix(L-s.Rational(37)*s.eye(N)).to_field().rank()==11
assert DomainMatrix.from_Matrix(L-s.Rational(73,2)*s.eye(N)).to_field().rank()==70
assert cert['affine_rank']==cert['centered_gram_rank']==69
assert cert['parent_distances_preserved']==2628 and cert['new_distances']==73
# Negative tests on labelled metric and histogram-insensitive relabelling.
wrong=computed.copy();wrong[(1,74)]=F(9,4)
assert wrong!=supplied
swap=lambda i:8 if i==1 else 1 if i==8 else i
permuted={tuple(sorted((swap(i),swap(j)))):value for (i,j),value in computed.items()}
assert permuted!=computed and Counter(permuted.values())==Counter(computed.values())

print('INDEPENDENT D69 PASS: 2701/2701 rational pairs, 2628/2628 parent pairs, Gram rank 69')
print('INDEPENDENT IDEAL GRAPH PASS: eigenvalue 37 x63, eigenvalue 73/2 x4')
print('CERTIFICATE_SHA256='+hashlib.sha256(p.read_bytes()).hexdigest())
