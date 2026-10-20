"""Step 8 - representation of NUTS 2 regions and NUTS 3 counties (share of samples vs share of area).

Input : data/input/nuts_romania_count.xlsx (sheets "NUTS 2" and "NUTS 3": sample counts per unit and
        Eurostat NUTS typology codes). Area of each county (km2) is fixed in the dictionary below.
Output: outputs/tables/nuts2_table.csv, nuts3_table.csv; summary statistics are printed to the console.
"""
import _bootstrap  # noqa: F401
from pathlib import Path
import pandas as pd, numpy as np
ROOT = Path(__file__).resolve().parents[1]
X = ROOT / "data" / "input" / "nuts_romania_count.xlsx"
n2=pd.read_excel(X,sheet_name="NUTS 2"); n3=pd.read_excel(X,sheet_name="NUTS 3")
area={"RO121":6250,"RO421":7746,"RO311":6822,"RO211":6622,"RO111":7539,"RO112":5358,"RO212":4987,"RO122":5361,"RO221":4766,"RO321":240,"RO222":6101,"RO422":8532,"RO312":5087,"RO113":6672,"RO223":7104,"RO123":3707,"RO313":4056,"RO411":7425,"RO224":4465,"RO314":3544,"RO412":5572,"RO124":6637,"RO423":7072,"RO315":4455,"RO213":5477,"RO322":1564,"RO114":6303,"RO413":4942,"RO125":6705,"RO214":5897,"RO414":5503,"RO316":4715,"RO115":4420,"RO116":3867,"RO126":5432,"RO215":8553,"RO317":5788,"RO424":8692,"RO225":8484,"RO216":5317,"RO415":5764,"RO226":4854}
n3["area"]=n3.NUTS_ID.map(area)
print(n3.area.isna().sum(), len(n3), n3.numar_prob.sum(), n3.area.sum())
n3["reg"]=n3.NUTS_ID.str[:4]
n2a=n3.groupby("reg").agg(area=("area","sum"),n3sum=("numar_prob","sum")).reset_index()
n2=n2.merge(n2a,left_on="NUTS_ID",right_on="reg")
assert (n2.n3sum==n2.numar_prob).all()
def stat(df,key):
    df=df.copy(); N=df.numar_prob.sum(); A=df.area.sum()
    df["pct_samples"]=100*df.numar_prob/N; df["pct_area"]=100*df.area/A
    df["ratio"]=df.pct_samples/df.pct_area; df["dens100"]=100*df.numar_prob/df.area
    return df
n2=stat(n2,"n2"); n3=stat(n3,"n3")
print(n2[["NUTS_NAME","numar_prob","pct_samples","area","pct_area","ratio","dens100"]].round(2).to_string())
print(n3.sort_values("dens100",ascending=False)[["NUTS_NAME","numar_prob","area","ratio","dens100"]].round(2).to_string())
# Gini of county densities (unweighted) and sample-based concentration
def gini(x):
    x=np.sort(np.asarray(x,float)); k=len(x); c=np.cumsum(x); return (k+1-2*c.sum()/c[-1])/k
print("Gini of counts over counties", gini(n3.numar_prob), "Gini densities", gini(n3.dens100))
# Lorenz area-weighted: share of samples against share of area
o=n3.sort_values("dens100"); cs=np.cumsum(o.numar_prob)/o.numar_prob.sum(); ca=np.cumsum(o.area)/o.area.sum()
# area-weighted Gini
xs=np.r_[0,ca]; ys=np.r_[0,cs]; G=1-2*np.trapezoid(ys,xs); print("area-weighted Gini", G)
top=n3.sort_values("numar_prob",ascending=False)
print("top1..5 share", (top.numar_prob.cumsum()/1091).head(6).round(3).tolist(), top.NUTS_NAME.head(6).tolist())
print("counties with ratio<0.5:", (n3.ratio<0.5).sum(), " ratio>2:", (n3.ratio>2).sum(), "<=10 samples", (n3.numar_prob<=10).sum())
print("median dens", n3.dens100.median(), "range", n3.dens100.min(), n3.dens100.max())
# typology
for col,lab in [("MOUNT_TYPE",{1:"mountain (>50% area)",4:"non-mountain"}),("URBN_TYPE",{1:"predominantly urban",2:"intermediate",3:"predominantly rural"}),("COAST_TYPE",{1:"coastal",3:"non-coastal"})]:
    g=n3.groupby(col).agg(k=("NUTS_ID","size"),n=("numar_prob","sum"),area=("area","sum")).reset_index()
    g["pct_s"]=100*g.n/g.n.sum(); g["pct_a"]=100*g.area/g.area.sum(); g["ratio"]=g.pct_s/g.pct_a; g["d"]=100*g.n/g.area; g[col]=g[col].map(lab)
    print(g.round(2).to_string())
out = ROOT / "outputs" / "tables"; out.mkdir(parents=True, exist_ok=True)
n2.to_csv(out / "nuts2_table.csv", index=False); n3.to_csv(out / "nuts3_table.csv", index=False)
