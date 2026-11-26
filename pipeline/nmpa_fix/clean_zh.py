import pandas as pd, re
from pathlib import Path
HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "data"
d=pd.read_csv(DATA / "nmpa_candidates.csv")
pc=pd.read_csv(DATA / "pubchem_results.csv")
d=d.merge(pc[["query","status"]],left_on="key2",right_on="query",how="left")
FORMS=["注射用","注射液","缓释片","控释片","肠溶片","肠溶胶囊","分散片","口腔崩解片","咀嚼片","泡腾片","含片","薄膜衣片","软胶囊","胶囊","颗粒","干混悬剂","混悬剂","口服溶液","口服液","溶液","滴眼液","眼膏","滴耳液","滴鼻液","喷雾剂","气雾剂","吸入粉雾剂","吸入","乳膏","软膏","凝胶贴膏","凝胶","贴剂","贴膏","巴布膏","栓","酊","搽剂","洗剂","糖浆","散","丸","乳剂","干乳剂","脂质体","微球","植入剂","粉针","冻干粉","片","注射","原料药","缓释","控释","肠溶","(Ⅱ)","(Ⅰ)","（Ⅱ）","（Ⅰ）","Ⅱ","Ⅰ"]
SALT_PRE=["盐酸","氢溴酸","硫酸","甲磺酸","苯磺酸","对甲苯磺酸","马来酸","富马酸","酒石酸","枸橼酸","磷酸","醋酸","乳酸","琥珀酸","丁二磺酸","乙磺酸","门冬氨酸","赖氨酸","精氨酸","葡甲胺","硝酸","草酸","苹果酸","萘磺酸","双羟萘酸","依地酸","水合","无水"]
SALT_SUF=["钠","钾","钙","镁","锌","铝","二钠","氯化物","碘化物","溴化物"]
def clean(z):
    z=str(z); z=re.sub(r"[\(（][^)）]*[\)）]","",z)
    for f in sorted(FORMS,key=len,reverse=True): z=z.replace(f,"")
    return z.strip()
def strip_salt(z):
    for s in sorted(SALT_PRE,key=len,reverse=True):
        if z.startswith(s) and len(z)>len(s)+1: z=z[len(s):]; break
    for s in sorted(SALT_SUF,key=len,reverse=True):
        if z.endswith(s) and len(z)>len(s)+1: z=z[:-len(s)]; break
    return z
d["zh_clean"]=d.zh.map(clean); d["zh_base"]=d.zh_clean.map(strip_salt)
d.to_csv(HERE / "nmpa_candidates_clean.csv", index=False)
open(HERE / "zh_names.txt", "w").write("\n".join(sorted(x for x in set(d.zh_clean.dropna()) | set(d.zh_base.dropna()) if 1 < len(x) < 40)))
nf=d[d.status!="ok"]
print(len(d),len(nf), nf.zh_clean.nunique(), nf.zh_base.nunique())
print(nf[["zh","zh_clean","zh_base","year"]].sample(25,random_state=2).to_string())
print(d[d.key2.str.contains("carrimycin")][["zh","zh_clean","zh_base"]])
