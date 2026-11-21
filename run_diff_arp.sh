#!/usr/bin/env bash
# run_diff_arp.sh -- Jay's follow-up: same simplified model (thesis 4.3), original D = 20 um^2/s,
# with 1x, 2x and 4x the AR+ density (0.03, 0.06, 0.12 AR+ cells per 30 um pixel = 33, 67, 133 per mm^2),
# plus male D = 40 and female (no AR+) for reference. Prints where x^2 levels off and the effective
# diffusion constant in the linear (early) part, D_eff = slope/2 from x^2 = 2 D_eff t (the thesis convention,
# which gives 0.75 um^2/s for the female case vs 0.76 in thesis Fig 4.8).
# Saves figures/Diffusion_male_ARp_density.png. ~5 min.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash run_diff_arp.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
mkdir -p diff_arp figures && cd diff_arp
g++ -O3 -std=c++17 ../../ifng_eff_diff_arp.cpp -o ifd || { echo "compile failed"; exit 1; }
DAYS=${DAYS:-15}
#            name AR+/pixel D multiplier
for c in "m1 0.03 1" "a2 0.06 1" "a4 0.12 1" "m2 0.03 2" "f1 0 1"; do
  set -- $c; ./ifd $1 $DAYS 0.01 $2 $3 2> $1.log &
done
echo "running 5 simulations for $DAYS days of model time (a few minutes)..."; wait
SRC=../diff_arp python3 - <<'PY'
import csv, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def load(p):
    r=[x for x in csv.DictReader(open(p+"_front.csv")) if x["run"]=="bound"]
    return np.array([float(x["time_h"]) for x in r]), np.array([float(x["front_um2"]) for x in r])
def deff(t,x,h=24):
    m=t<=h; return np.polyfit(t[m]*3600,x[m],1)[0]/2
runs=[("m1","male, D = 20, AR+ ×1 (33/mm²)","#1f4db8","-"),("a2","male, D = 20, AR+ ×2 (67/mm²)","#1f4db8","--"),
      ("a4","male, D = 20, AR+ ×4 (133/mm²)","#1f4db8",":"),("m2","male, D = 40, AR+ ×1","#6a8fe0","-."),
      ("f1","female, D = 20, no AR+","#b81f5b","-")]
print(f"{'run':38s} {'x^2 level-off (um^2)':>21s} {'front (um)':>10s} {'D_eff 0-24 h (um^2/s)':>22s}")
fig,ax=plt.subplots(figsize=(6.8,5))
for p,lab,col,ls in runs:
    t,x=load(p); end=x[-1] if p!="f1" else float("nan")
    print(f"{lab:38s} {x[-1]:21.3g} {x[-1]**.5:10.0f} {deff(t,x):22.2f}" + ("   (still rising, no level-off)" if p=="f1" else ""))
    ax.plot(t,x,ls,color=col,lw=2,label=lab)
ax.set_xlabel("Time (hours)"); ax.set_ylabel(r"Front position$^2$, $x^2$ ($\mu m^2$)"); ax.set_ylim(0,8e5)
ax.set_title("IFN-γ spread: AR+ density and D"); ax.legend(frameon=False,fontsize=8.5); ax.grid(True,color="#e6e6e6",lw=.6)
fig.tight_layout(); fig.savefig("../figures/Diffusion_male_ARp_density.png",dpi=200); print("wrote figures/Diffusion_male_ARp_density.png")
PY
