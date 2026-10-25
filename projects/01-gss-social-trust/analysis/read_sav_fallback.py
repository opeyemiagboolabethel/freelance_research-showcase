"""Minimal audited reader for bytecode-compressed SPSS $FL2 numeric fields.
Used only where the standard pyreadstat library is unavailable; never a replacement for
metadata verification or an independent read with IBM SPSS / R haven.
"""
from __future__ import annotations
import struct, io, math
from pathlib import Path
import pandas as pd

U32 = struct.Struct('<i'); F64=struct.Struct('<d')
def get_i(f):return U32.unpack(f.read(4))[0]
def read_sav_numeric(path,cols):
    f=open(path,'rb'); b=f.read(176)
    assert b[:4]==b'$FL2','SPSS SAV compression type unsupported'
    layout,nslots,compression,weight,ncases=struct.unpack_from('<5i',b,64)
    bias=struct.unpack_from('<d',b,84)[0]
    assert (layout,compression)==(2,1),(layout,compression)
    names=[];types=[];labels=[]; extended={}
    while True:
        offset=f.tell();r=get_i(f)
        if r==2:
            typ,haslabel,nmissing,printfmt,writefmt=struct.unpack('<5i',f.read(20))
            name=f.read(8).decode('latin-1').strip()
            desc=''
            if haslabel:
                labelsize=get_i(f);desc=f.read((labelsize+3)//4*4)[:labelsize].decode('utf-8','replace')
            if nmissing: f.seek(abs(nmissing)*8,1)
            names.append(name);types.append(typ);labels.append(desc)
        elif r==3:
            n=get_i(f)
            for j in range(n):
                f.seek(8,1);l=f.read(1)[0];f.seek(((l+1+7)//8)*8-1,1)
        elif r==4:
            n=get_i(f);f.seek(4*n,1)
        elif r==6:
            n=get_i(f);f.seek(80*n,1)
        elif r==7:
            subtype,size,n=struct.unpack('<iii',f.read(12))
            if size<0 or n<0 or size*n>200000000:raise RuntimeError(f'bad type7 {subtype} {size} {n} offset {offset}')
            payload=f.read(size*n)
            if subtype==13:
                try:
                    for fragment in payload.decode('utf-8','replace').split('\t'):
                        if '=' in fragment:
                            short,long=fragment.split('=',1);extended[short]=long
                except Exception:pass
        elif r==999:
            reserved=get_i(f);break
        else:raise RuntimeError(f'unknown dictionary record {r} at offset {offset}')
    print('dictionary',len(names),'slots',nslots,'case count',ncases,'data offset',f.tell(),'compression',compression)
    assert len(names)==nslots,(len(names),nslots)
    wanted=[name.upper() for name in cols]
    ix={x:[i for i,n in enumerate(names) if n.upper()==x or extended.get(n,'').upper()==x] for x in wanted}
    print('found',ix)
    absent=[x for x,y in ix.items() if not y]
    if absent:print('missing requested cols',absent)
    pos={i:x for x,y in ix.items() for i in y}
    results={x:[] for x in wanted}
    nonnumeric={x:[] for x in wanted if ix[x] and types[ix[x][0]]!=0}
    buff=[];idx=0;cases=0
    while cases<ncases:
        codebytes=f.read(8)
        if len(codebytes)!=8:raise RuntimeError(f'EOF at case {cases} element {idx}')
        for code in codebytes:
            if code==0:continue
            if code==252:break
            if code==253:
                raw=f.read(8)
                if len(raw)!=8: raise RuntimeError('Unexpected EOF reading literal')
            elif code==254:raw=b' '*8
            elif code==255:raw=None
            else:raw=None
            if idx in pos:
                nm=pos[idx]
                if types[idx]==0:
                    value=F64.unpack(raw)[0] if code==253 else (code-bias if 1<=code<=251 else math.nan)
                    if value < -1e250: value=math.nan
                else:
                    value=raw.decode('utf-8','replace').strip() if raw else None
                results[nm].append(value)
            idx+=1
            if idx==nslots:
                idx=0;cases+=1
                if cases==ncases:break
        if code==252 and cases<ncases:raise RuntimeError(f'Premature end at {cases} cases')
    print('parsed cases',cases,'bytes read',f.tell())
    out=pd.DataFrame({x:results[x] for x in wanted if ix[x]})
    return out,{'record_count':ncases,'slots':nslots,'compression':compression,'metadata':{x:{'position':ix[x][0],'label':labels[ix[x][0]],'type':types[ix[x][0]]} for x in ix if ix[x]},'reader':'Independent narrow-purpose Python SPSS bytecode parser'}

if __name__=='__main__':
    import sys,json
    fields='ID YEAR HAPPY TRUST CONFINAN CONEDUC AGE SEX DEGREE MARITAL WTSSNRPS WTSSPS WTSSNRPS_AS VPSU VSTRAT SAMPLE BALLOT MODE PARTIAL GSSMODE'.split()
    df,meta=read_sav_numeric(sys.argv[1],fields)
    print(df.head().T.to_string())
    print(df.describe().T[['count','mean','min','max']].to_string())
    print(meta)
