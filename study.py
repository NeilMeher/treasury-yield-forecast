import argparse
import csv
from datetime import date
from hashlib import sha256
import json
from pathlib import Path
from subprocess import run

import numpy as np
from figure import Plot

ROOT = Path(__file__).resolve().parent
URL = "https://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv"
MATS = np.array([1,2,3,5,7,10])
HORIZONS = (1,3,12)
NAMES = ("No change","Maturity AR(1)","DNS")
STRESS = {"GFC":("2008-09","2009-06"),"2020":("2020-03","2020-12"),"2022":("2022-01","2022-12")}


def load_months(path):
    monthly = {}
    with path.open(encoding="utf-8-sig",newline="") as file:
        for line in file:
            if line.startswith("Date,"):
                break
        else:
            raise ValueError("Missing GSW header")
        for row in csv.DictReader(file,fieldnames=line.strip().split(",")):
            day = date.fromisoformat(row["Date"])
            if not 1985<=day.year<=2025:
                continue
            raw = [row[f"SVENY{m:02d}"] for m in MATS]
            if any(val in ("","NA",".") for val in raw):
                continue
            month = (day.year,day.month)
            if month not in monthly or day > monthly[month][0]:
                monthly[month] = day,[float(val) for val in raw]
    expected = [(y,m) for y in range(1985,2026) for m in range(1,13)]
    if sorted(monthly) != expected:
        raise ValueError(f"Missing months: {sorted(set(expected)-monthly.keys())}")
    dates,values = zip(*(monthly[month] for month in expected))
    payload = "\n".join(",".join([str(day),*(f"{x:.6f}" for x in v)]) for day,v in zip(dates,values))
    return dates,np.array(values),sha256(payload.encode()).hexdigest()


def loadings(maturities):
    x = .0609*np.asarray(maturities)*12
    slope = np.divide(-np.expm1(-x), x, out=np.ones_like(x), where=x!=0)
    return np.column_stack((np.ones_like(x),slope,slope-np.exp(-x)))


def ar_coefficients(values):
    x,y = values[:-1],values[1:]
    xmean,ymean = x.mean(axis=0),y.mean(axis=0)
    centered = x-xmean
    scale = np.sum(centered**2,axis=0)
    beta = np.divide(np.sum(centered*(y-ymean),axis=0),scale,out=np.zeros_like(scale),where=scale>1e-18)
    return ymean-beta*xmean,beta


def forecast(values,horizon):
    a,b=ar_coefficients(values)
    state=values[-1]
    for _ in range(horizon):
        state=a+b*state
    return state


def walk(values,initial):
    matrix=loadings(MATS)
    factors=np.linalg.lstsq(matrix,values.T,rcond=None)[0].T
    forecasts = {}
    for h in HORIZONS:
        predictions,actual,origins=[],[],[]
        for i in range(initial-1,len(values)-h):
            window=values[:i+1]
            predictions.append([window[-1],forecast(window,h),matrix@forecast(factors[:i+1],h)])
            actual.append(values[i+h])
            origins.append(i)
        forecasts[h]=(np.array(predictions),np.array(actual),np.array(origins))
    return forecasts,factors


def block_ci(gaps,reps,seed_value):
    n=len(gaps)
    length=12
    start=np.random.default_rng(seed_value).integers(0,n,size=(reps,(n+length-1)//length,1))
    mask=((start+np.arange(length))%n).reshape(reps,-1)[:,:n]
    return np.quantile(gaps[mask].mean(axis=1),[.025,.975],axis=0)


def evaluate(dates,runs,quick):
    rows=[]
    for h,(forecasts,actual,origins) in runs.items():
        squared=((forecasts-actual[:,None,:])*100)**2
        score=np.sqrt(squared.mean(axis=0))
        intervals=block_ci(squared[:,1:,:]-squared[:,0:1,:],100 if quick else 1000,707+h)
        targets=[str(dates[x+h])[:7] for x in origins]
        masks={"All":np.ones(len(origins),dtype=bool),
               **{key:np.array([lo<=t<=hi for t in targets]) for key,(lo,hi) in STRESS.items()}}
        for key,mask in masks.items():
            if not any(mask):
                continue
            rmses=np.sqrt(squared[mask].mean(axis=0))
            for j,name in enumerate(NAMES):
                for i,mat in enumerate(MATS):
                    ratio = rmses[j,i]/rmses[0,i] if rmses[0,i] else None
                    bounds = intervals[:,j-1,i] if j and key=="All" else (None,None)
                    rows.append(dict(period=key,horizon=h,maturity=int(mat),model=name,n=int(mask.sum()),
                                     rmse_bp=round(rmses[j,i],6),relative_rmse=round(ratio,6) if ratio is not None else None,
                                     mse_gap_lo=round(bounds[0],6) if bounds[0] is not None else None,
                                     mse_gap_hi=round(bounds[1],6) if bounds[1] is not None else None))
        print(f"{h}m pooled RMSE (bp), {len(origins)} origins: {np.sqrt(squared.mean(axis=(0,2))).round(3)}")
    return rows


def heatmap(out,rows,label):
    cells=[r for r in rows if r["period"]=="All" and r["model"]!="No change"]
    parts=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 780 402" role="img">',
           '<title>Relative out-of-sample Treasury forecast error</title>',
           '<rect width="780" height="402" fill="#fbfaf6"/>',
           '<g font-family="Arial,sans-serif" fill="#163043">',
           '<text x="24" y="38" font-size="24">Treasury forecast tournament</text>',
           f'<text x="24" y="62" font-size="13">{label} · RMSE ratio vs no change · below 1.00 improves</text>']
    for i,mat in enumerate(MATS):
        parts.append(f'<text x="{209+i*95}" y="94" text-anchor="middle">{mat}Y</text>')
    for row,(model,h) in enumerate((m,h) for m in NAMES[1:] for h in HORIZONS):
        y=106+row*44
        parts.append(f'<text x="25" y="{y+25}" font-size="14">{model} · {h}m</text>')
        for col,mat in enumerate(MATS):
            value=next(c["relative_rmse"] for c in cells if c["model"]==model and c["horizon"]==h and c["maturity"]==mat)
            color="#e4e4d9" if value is None or abs(value-1)<.02 else "#b8ded7" if value<1 else "#efd0bb"
            x=163+col*95
            value="—" if value is None else f"{value:.2f}"
            parts.append(f'<rect x="{x}" y="{y}" width="92" height="39" rx="3" fill="{color}"/>')
            parts.append(f'<text x="{x+46}" y="{y+25}" text-anchor="middle" font-size="15">{value}</text>')
    parts.append('</g></svg>')
    (out/"relative-error.svg").write_text("\n".join(parts)+"\n",encoding="utf-8")


def charts(out,dates,runs,factors,rows,label):
    heatmap(out,rows,label)
    years=np.array([d.year+(d.month-1)/12 for d in dates])
    colors=("#176876","#93406b","#9c6137")
    p=Plot("Estimated Nelson–Siegel factors",label,(min(years),max(years)),
           (float(factors.min())-.3,float(factors.max())+.3),"Observation year","Factor (percentage points)")
    for i,name in enumerate(("Level","Slope","Curvature")):
        p.series(zip(years,factors[:,i]),colors[i],name,dots=False)
    p.save(out/"factors.svg")

    predictions,actual,origins=runs[3]
    target=years[origins+3]
    actual=actual[:,1]
    no_change=predictions[:,0,1]
    dns=predictions[:,2,1]
    all_y=np.r_[actual,no_change,dns]
    p=Plot("Three-month forecast of the 2-year yield",label,
           (min(target),max(target)),(float(all_y.min())-.2,float(all_y.max())+.2),
           "Realization year","Yield (%)")
    for line,name,color in zip((actual,no_change,dns),("Realized","No change","DNS"),colors):
        p.series(zip(target,line),color,name,dots=False)
    p.save(out/"forecast.svg")


def synthetic():
    rng=np.random.default_rng(707)
    days=[date(1990+i//12,i%12+1,25) for i in range(120)]
    factors=np.cumsum(rng.normal(0,.08,(len(days),3)),axis=0)+[4,-1,.5]
    return days,factors@loadings(MATS).T


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument("--quick",action="store_true")
    cli.add_argument("--download",action="store_true")
    cli.add_argument("--input",type=Path,default=ROOT/"data"/"feds200628.csv")
    cli.add_argument("--output",type=Path,default=ROOT/"report")
    args=cli.parse_args()
    if args.quick:
        dates,values=synthetic()
    else:
        if args.download:
            args.input.parent.mkdir(exist_ok=True,parents=True)
            run(["curl","--fail","--silent","--show-error","--location","--max-time","60",
                 URL,"--output",str(args.input)],check=True)
        dates,values,fingerprint=load_months(args.input)
        if fingerprint!=(ROOT/"source-sha256.txt").read_text().strip():
            raise ValueError(f"Source monthly snapshot differs: {fingerprint}")
    initial=60 if args.quick else 180
    runs,factors=walk(values,initial)
    rows=evaluate(dates,runs,args.quick)
    args.output.mkdir(exist_ok=True,parents=True)
    with (args.output/"accuracy.csv").open("w",newline="",encoding="utf-8") as file:
        writer=csv.DictWriter(file,fieldnames=rows[0],lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    charts(args.output,dates,runs,factors,rows,"Synthetic fixture" if args.quick else "Fed fitted zero-coupon yields")
    a,b=ar_coefficients(values[:initial])
    audit=dict(origin=str(dates[initial-1]),target=str(dates[initial]),maturity=2,
               value=values[initial-1,1],intercept=a[1],slope=b[1],
               maturity_ar_1m=runs[1][0][0,1,1],realization=values[initial,1])
    (args.output/"one-origin.json").write_text(json.dumps({k:round(v,10) if isinstance(v,float) else v for k,v in audit.items()},indent=2)+"\n")
    print(audit)


if __name__=="__main__":
    main()
