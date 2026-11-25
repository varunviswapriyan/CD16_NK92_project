#!/usr/bin/env bash
# run_diff_nounbind.sh -- Das's request: the simplified IFN-g model (thesis 4.3) with and without unbinding.
# K = 100 is the original model (IFN-g binds and unbinds, equilibrium each step); K = 0 means IFN-g never unbinds,
# so free IFN-g binds until the free receptors in a pixel are used up. Male = 0.03 AR+ cells per 30 um pixel, female = 0.
# D = 20 um^2/s, 15 days, ~2-3 min. Saves figures/Diffusion_unbinding_vs_none.png.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash run_diff_nounbind.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
mkdir -p diff_nounbind figures && cd diff_nounbind
g++ -O3 -std=c++17 ../../ifng_eff_diff_arp.cpp -o ifd || { echo "compile failed"; exit 1; }
#            name AR+/pixel Dmult K
for c in "m1 0.03 1 100" "f1 0 1 100" "m0 0.03 1 0" "f0 0 1 0"; do
  set -- $c; ./ifd $1 15 0.01 $2 $3 $4 2> $1.log &
done
echo "running 4 simulations for 15 days of model time..."; wait
python3 - <<'PY'
import csv, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def load(p):
    r=[x for x in csv.DictReader(open(p+"_front.csv")) if x["run"]=="bound"]
    return np.array([float(x["time_h"]) for x in r]), np.array([float(x["front_um2"]) for x in r])
fig,axs=plt.subplots(1,2,figsize=(11,4.4))
print(f"{'run':42s} {'x^2 at 15 d (um^2)':>19s} {'front (um)':>10s}")
for ax,(m,f,title,ymax,fit) in zip(axs,(("m1","f1","IFN-γ unbinds (K = 100, as before)",2e6,False),
                                       ("m0","f0","IFN-γ never unbinds (K = 0)",2e5,True))):
    for p,lab,col in ((f,"female, no AR+","#c0143c"),(m,"male, AR+ 0.03 per pixel","#1f4db8")):
        t,x=load(p); ax.plot(t,x,color=col,lw=1.6 if fit else 2,label=lab)
        extra="   levels off"
        if p.startswith("f"):
            k=np.polyfit(t,x,1); deff=k[0]/3600/2
            extra=f"   still rising, D_eff over 15 d = {deff:.2f} um^2/s"
            if fit: ax.plot(t,np.polyval(k,t),"--",color="k",lw=1,label=f"straight-line fit (D_eff = {deff:.2f} µm²/s)")
        print(f"{title.split(' (')[0]+', '+lab:42s} {x[-1]:19.3g} {x[-1]**.5:10.0f}{extra}")
    ax.set_ylim(0,ymax); ax.set_title(title); ax.set_xlabel("Time (hours)"); ax.grid(True,color="#e6e6e6",lw=.6)
    ax.set_ylabel(r"Front position$^2$, $x^2$ ($\mu m^2$)"); ax.legend(frameon=False,loc="upper left",fontsize=9)
axs[1].text(0.98,0.03,"steps: the front moves one 30 µm pixel at a time",transform=axs[1].transAxes,ha="right",fontsize=8,color="#5a6475")
fig.tight_layout(); fig.savefig("../figures/Diffusion_unbinding_vs_none.png",dpi=200); print("wrote figures/Diffusion_unbinding_vs_none.png")
PY
