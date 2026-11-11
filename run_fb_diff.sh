#!/usr/bin/env bash
# run_fb_diff.sh -- Nov 11 meeting deliverables 1 and 2.
#
#   bash run_fb_diff.sh          submit everything + run the fast part now
#   bash run_fb_diff.sh analyze  (when queue empty) fold changes + full-model fronts
#
# (1) FEEDBACK OFF: IFN-g secretion feedback removed. Deck arg 28
#     (max_ifng_production, 2.0) -> 0.5, so secretion = 0.5 x k (the basal level,
#     same as an unstimulated cell); nothing else changed. Male + female, NSEED seeds.
# (2) EFFECTIVE DIFFUSION, male vs female, two ways:
#   (a) FAST (minutes, login node): thesis simplified model + AR+ T cells that
#       internalize bound IFN-g at the deck rate (0.0016667/s). Sweep AR+ density:
#       0 (female-like) up to dense. -> figures/Deff_vs_ARp.png + table.
#   (b) FULL MODEL: one seed per sex, standard params, run to day 17.5 (504) with
#       free IFN-g dumps every 6 h; IFN-g penetration depth into the tumour vs time.
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
SRC=rsc_out_2026-07-03; DECK="$PWD/$SRC/in.full"; BIN="$PWD/spk_sexbias"
NSEED="${NSEED:-20}"; TIME="${TIME:-30:00:00}"
MS=$SRC/sex_bias_initialization_files/d14_M_4_SPPARKS_init
FS=$SRC/sex_bias_initialization_files/d14_F_4_SPPARKS_init
FEMAWK='{ if ($0 ~ /recruitment/ && $0 !~ /app_style/){ if($0 ~ /AR\+/)$6=0; else if($0 ~ /Tp recruitment/||$0 ~ /Tm recruitment/)$6=$6/0.63 }
          if ($0 ~ /Tp differentiation/ && $0 !~ /app_style/){ if($0 ~ /AR\+ Tm/)$6=0; else $6=$6*2 } print }'
mkdir -p figures

if [ "${1:-}" = "analyze" ]; then
python3 - <<'PY'
import glob, math, numpy as np
def fold(f):
    hdr=None; first=last=t=None
    for line in open(f):
        s=line.split()
        if "sum(i1)" in s:
            if hdr: break
            hdr=s; k=s.index("sum(i1)"); continue
        if hdr and len(s)==len(hdr):
            try: tt=float(s[0]); v=float(s[k])
            except ValueError: continue
            if first is None: first=v
            last=v; t=tt
    return last/first if first and t is not None and t>=1000 else None
print("=== (1) feedback off: day-21 fold change ===")
res={}
for sex in ("male","female"):
    v=[x for x in (fold(f) for f in glob.glob(f"runs_big/fboff/{sex}/s*/log.spparks")) if x]
    if v:
        m=np.mean(v); se=np.std(v,ddof=1)/math.sqrt(len(v)) if len(v)>1 else 0; res[sex]=m
        print(f"  {sex:6s} n={len(v):3d} fold={m:.3f} +/- {se:.3f}")
if len(res)==2: print(f"  ratio male/female = {res['male']/res['female']:.2f}   (standard 2.38/1.55 = 1.54)")

print("\n=== (2b) full model: IFN-g penetration into tumour ===")
L=600; N=200
def read(p):
    cols=None; data=[]; inv=False
    for line in open(p):
        s=line.split()
        if not s: continue
        if s[0]=="Values": inv=True; continue
        if cols is None and s[0]=="id": cols=[c for c in s if c!="values"]; continue
        if inv and cols and len(s)>=len(cols):
            try: data.append([float(x) for x in s[:len(cols)]])
            except ValueError: pass
    a=np.array(data); return {c:a[:,i] for i,c in enumerate(cols)}
def depth_map(inside):
    # distance (in IFN-g pixels, 30 um) from the tumour edge, inside the tumour
    d=np.where(inside,np.inf,0.0); changed=True; k=0
    while changed and k<N:
        k+=1; P=np.pad(d,1,mode="constant",constant_values=0.0)
        nb=np.minimum.reduce([P[:-2,1:-1],P[2:,1:-1],P[1:-1,:-2],P[1:-1,2:]])+1
        new=np.where(inside,np.minimum(d,nb),0.0); changed=not np.array_equal(new,d); d=new
    return d
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(6.4,5)); out=[]
for sex,col in (("male","#1f4db8"),("female","#b81f5b")):
    files=sorted(glob.glob(f"runs_front/{sex}/dump.front.*"),key=lambda p:int(p.split(".")[-1]))
    if len(files)<3: print(f"  {sex}: not enough dumps"); continue
    d0=read(files[0]); ix=np.clip(d0["x"].astype(int),0,L-1)//3; iy=np.clip(d0["y"].astype(int),0,L-1)//3
    can=np.zeros((N,N)); np.add.at(can,(iy,ix),d0["i1"])
    inside=can>0; dep=depth_map(inside); maxd=int(dep[inside].max())
    T=[];X2=[]
    for p in files[1:]:
        d=read(p); g=np.zeros((N,N)); np.add.at(g,(iy,ix),d["d1"])
        prof=np.array([g[(dep==r)].mean() if (dep==r).any() else 0 for r in range(1,maxd+1)])
        if prof.max()<=0: continue
        thr=0.1*prof[0] if prof[0]>0 else 0.1*prof.max()
        below=np.where(prof<thr)[0]; r=(below[0] if len(below) else len(prof))*30.0
        t=int(p.split(".")[-1])/6.0  # hours
        T.append(t); X2.append(r*r)
    T=np.array(T); X2=np.array(X2)
    if len(T)>=3:
        sl=np.polyfit(T,X2,1)[0]; out.append((sex,sl))
        ax.plot(T,X2,"o-",color=col,label=f"{sex} (slope {sl:.0f} um^2/h)")
        print(f"  {sex:6s} penetration^2 slope = {sl:.0f} um^2/h  (~D_eff {sl/3600/4:.3f} um^2/s, 1D estimate)")
if len(out)==2: print(f"  male/female slope ratio = {out[0][1]/out[1][1]:.2f}")
ax.set_xlabel("Time (hours)"); ax.set_ylabel("IFN-γ penetration depth² (µm²)")
ax.set_title("IFN-γ spread into the tumour, full model"); ax.legend(frameon=False); fig.tight_layout()
fig.savefig("figures/Front_full_model_m_f.png",dpi=200); print("  wrote figures/Front_full_model_m_f.png")
PY
exit 0
fi

[ -x "$BIN" ] || { echo "binary $BIN missing"; exit 1; }

# ---------------- (1) feedback off: submit ----------------
FB='/^app_style/ s/ 4 4 2\.0 1000 1000 1000 / 4 4 0.5 1000 1000 1000 /'
sed "$FB" "$DECK" | grep -q " 4 4 0.5 1000 1000 1000 " || { echo "!! feedback edit did not match app_style line"; exit 1; }
JOB=$PWD/runs_big/_job_generic.sh
[ -f "$JOB" ] || { echo "missing $JOB (created by do_tasks.sh)"; exit 1; }
for sex in male female; do
  base=runs_big/fboff/$sex; [ -f $base/cfg.sh ] && { echo "  $base already submitted"; continue; }
  mkdir -p $base; printf '%s\n' "$FB" > $base/extra.sed
  [ $sex = male ] && { sl=$MS; fem=0; } || { sl=$FS; fem=1; }
  cat > $base/cfg.sh <<EOF
BASE=$PWD/$base
SEEDBASE=9000
SLIDEPATH=$PWD/$sl
SLIDENAME=$(basename $sl)
FEM=$fem
DECK=$DECK
BIN=$BIN
EXTRASED=$PWD/$base/extra.sed
EOF
  sbatch --job-name=fb_${sex:0:1} --array=1-$NSEED --time=$TIME --ntasks=1 \
    --output=$PWD/$base/slurm_%A_%a.out --export=ALL,CFG=$PWD/$base/cfg.sh "$JOB" >/dev/null \
    && echo "  submitted feedback-off $sex x$NSEED"
done

# ---------------- (2b) full-model front tracking: submit ----------------
for sex in male female; do
  d=runs_front/$sex; [ -f $d/job.sh ] && { echo "  $d already submitted"; continue; }
  mkdir -p $d/sex_bias_initialization_files
  [ $sex = male ] && sl=$MS || sl=$FS
  cp -f $sl $d/sex_bias_initialization_files/
  sed -e "/^seed/ s/[0-9][0-9]*/8377/" -e "s|read_sites.*|read_sites sex_bias_initialization_files/$(basename $sl)|" \
      -e 's|^#dump .*|dump 1 sites 36.0 dump.front.* id x y i1 d1|' -e 's|^run .*|run 504|' "$DECK" > $d/in.simulation
  grep -q "^dump 1 sites 36.0" $d/in.simulation || { echo "!! no #dump line in deck; front runs not submitted"; break; }
  [ $sex = female ] && { awk "$FEMAWK" $d/in.simulation > $d/t && mv $d/t $d/in.simulation; }
  cat > $d/job.sh <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=front_${sex:0:1}
#SBATCH --time=$TIME
#SBATCH --ntasks=1
#SBATCH --output=slurm_%j.out
module load foss/2023a
cd "\$SLURM_SUBMIT_DIR"
mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
  ( cd $d && sbatch job.sh >/dev/null ) && echo "  submitted full-model front run $sex (to day 17.5)"
done

# ---------------- (2a) simplified model with AR+ cells: run now ----------------
CPP=~/CD16_NK92_project/ifng_eff_diff_arp.cpp
module load foss/2023a 2>/dev/null || true
g++ -O3 -std=c++17 $CPP -o ifng_eff_diff_arp || { echo "compile failed"; exit 1; }
mkdir -p diff_arp
for a in 0 0.03 0.1 0.3 1; do ./ifng_eff_diff_arp diff_arp/arp_$a 5 0.01 $a 2>/dev/null & done; wait
python3 - <<'PY'
import csv, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def fit(t,x):
    mx=x.max(); m=(x>0.05*mx)&(x<0.55*mx)
    if m.sum()<3: m=np.ones_like(x,bool)
    return np.polyfit(t[m],x[m],1)[0]
rows=[]
for a in ("0","0.03","0.1","0.3","1"):
    R={}
    for r in csv.DictReader(open(f"diff_arp/arp_{a}_front.csv")):
        R.setdefault(r["run"],[]).append((float(r["time_h"]),float(r["front_um2"])))
    b=np.array(R["bound"]); z=np.array(R["bare"])
    d=20.0*fit(b[:,0],b[:,1])/fit(z[:,0],z[:,1]); rows.append((float(a),d))
print("\n=== (2a) simplified model: D_eff vs AR+ T-cell density ===")
print("  AR+ cells per 30um pixel   D_eff (um^2/s)")
for a,d in rows: print(f"  {a:>8g}                  {d:.3f}")
print("  (0 = female-like; male slide average is ~0.02-0.05 per pixel)")
fig,ax=plt.subplots(figsize=(6,4.6))
ax.plot([a for a,_ in rows][1:],[d for _,d in rows][1:],"o-",color="#1f4db8",label="with AR+ internalization")
ax.axhline(rows[0][1],ls="--",color="#b81f5b",label="no AR+ cells (female-like)")
ax.set_xscale("log"); ax.set_xlabel("AR+ T cells per IFN-γ pixel"); ax.set_ylabel("D_eff (µm²/s)")
ax.set_title("Effective IFN-γ diffusion vs AR+ density"); ax.legend(frameon=False); fig.tight_layout()
fig.savefig("figures/Deff_vs_ARp.png",dpi=200); print("  wrote figures/Deff_vs_ARp.png")
PY
echo
echo "Queue: squeue -u \$USER     When fb_* and front_* are done:  bash run_fb_diff.sh analyze"
