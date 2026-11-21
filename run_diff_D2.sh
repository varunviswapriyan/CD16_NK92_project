#!/usr/bin/env bash
# run_diff_D2.sh -- Jay's question: male IFN-g front with D and 2D (simplified model, thesis 4.3).
# Compiles ifng_eff_diff_arp.cpp, runs male (0.03 AR+ cells per pixel) and female (0) at D = 20 and 40
# for 15 days in parallel (~3 min), prints where x^2 levels off and saves figures/Diffusion_male_D_vs_2D.png.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash run_diff_D2.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
mkdir -p diff_D2 figures && cd diff_D2
g++ -O3 -std=c++17 ../../ifng_eff_diff_arp.cpp -o ifd || { echo "compile failed"; exit 1; }
DAYS=${DAYS:-15}
#            name  AR+/pixel  D multiplier
for c in "m1 0.03 1" "m2 0.03 2" "f1 0 1" "f2 0 2"; do
  set -- $c; ./ifd $1 $DAYS 0.01 $2 $3 2> $1.log &
done
echo "running 4 simulations for $DAYS days of model time (a few minutes)..."; wait
python3 - <<'PY'
import csv, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def load(p):
    r=[x for x in csv.DictReader(open(p+"_front.csv")) if x["run"]=="bound"]
    return [float(x["time_h"]) for x in r],[float(x["front_um2"]) for x in r]
for p,lab in (("m1","male   D=20"),("m2","male   D=40"),("f1","female D=20"),("f2","female D=40")):
    t,x=load(p); print(f"{lab}: x^2 at end {x[-1]:.3g} um^2 (front {x[-1]**.5:.0f} um); change over last 5 days {100*(x[-1]/x[-241]-1):.1f}%")
a=load("m1")[1][-1]; b=load("m2")[1][-1]; print(f"male level-off ratio 2D/D = {b/a:.2f}")
fig,ax=plt.subplots(figsize=(6.4,5))
for p,lab,ls in (("m1","male, D = 20 µm²/s","-"),("m2","male, D = 40 µm²/s","--")):
    t,x=load(p); ax.plot(t,x,ls,color="#1f4db8",lw=2,label=lab); ax.axhline(x[-1],color="grey",lw=.8,ls=":")
    ax.text(5,x[-1]*1.02,f"{x[-1]/1e5:.2f}×10⁵ µm²",fontsize=9,color="grey")
ax.set_xlabel("Time (hours)"); ax.set_ylabel(r"Front position$^2$, $x^2$ ($\mu m^2$)"); ax.set_title("IFN-γ spread in males: D vs 2D")
ax.legend(frameon=False); ax.grid(True,color="#e6e6e6",lw=.6); fig.tight_layout(); fig.savefig("../figures/Diffusion_male_D_vs_2D.png",dpi=200)
print("wrote figures/Diffusion_male_D_vs_2D.png")
PY
