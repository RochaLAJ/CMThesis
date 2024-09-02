import MDAnalysis as mda
import warnings
import pandas as pd
import MDAnalysis.analysis.align
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

from google.colab import drive
drive.mount("/content/drive", force_remount=True)
warnings.filterwarnings("ignore")
%matplotlib inline


def reorderLegend(ax=None,order=None,unique=False):
    if ax is None: ax=plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    if order is not None:
        keys=dict(zip(order,range(len(order))))
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t,keys=keys: keys.get(t[0],np.inf)))
    if unique:  labels, handles= zip(*unique_everseen(zip(labels,handles), key = labels))
    ax.legend(handles, labels)
    return(handles, labels)
  
i = mda.Universe('/content/drive/MyDrive/wt/leap/RETwt_leap.pdb','/content/drive/MyDrive/wt/MD/prod_imgWT.nc')
t = mda.Universe('/content/drive/MyDrive/M918T/leap/RET_M918T_leap.pdb','/content/drive/MyDrive/M918T/MD/prod_imgM918T.nc')
v = mda.Universe('/content/drive/MyDrive/M918V/leap/RET_M918V_leap.pdb','/content/drive/MyDrive/M918V/MD/prod_imgM918V.nc')

R_I = rms.RMSD(i,
             i,
             select='backbone',
             ref_frame=1).run()

R_T = rms.RMSD(t,
             t,
             select='backbone',
             ref_frame=1).run()

R_V = rms.RMSD(v,
               v,
               select='backbone',
               ref_frame=1).run()



df = pd.DataFrame(R_T.rmsd,
                  columns=['Frame_M918T', 'Time (ns)_M918T',
                           'p.Met918Thr'])
df2 = pd.DataFrame(R_V.rmsd,
                   columns=['Frame_M918V', 'Time (ns)_M918V',
                            'p.Met918Val'])

df3 = pd.DataFrame(R_I.rmsd, columns=['Frame_WT', 'Time (ns)_WT', 'Wildtype'])

final_df = pd.concat([df, df2, df3], axis=1)



final_df['Frame'] = final_df['Frame_M918T'] * 0.1
ax = final_df.plot(x='Frame', y=['p.Met918Val', 'p.Met918Thr', 'Wildtype'], color=["orange", "red", "darkgreen"])
ax.set_title("RMSD Comparison Between \n Wildtype, p.Met918Val and p.Met918Thr",fontweight="bold")
ax.set_ylabel('RMSD (Å)')
ax.set_xlabel('Nanoseconds')
reorderLegend(ax,['Wildtype', 'p.Met918Val', 'p.Met918Thr'])
plt.savefig('/content/drive/MyDrive/FinalThesisMTC/RET/v2/RMSDRenamed.png', dpi=600)

plt.hist(final_df['Wildtype'], bins=100, alpha=0.5, label="Wildtype", color='darkgreen')
plt.hist(final_df['M918V'], bins=100, alpha=0.5, label="p.Met918Val", color='orange')
plt.hist(final_df['M918T'], bins=100, alpha=0.5, label="p.Met918Thr", color='red')
plt.xlabel("RMSD (Å)")
plt.ylabel("Nanoseconds")
plt.title("Distribution of RMSD Values", fontweight="bold")
plt.legend(loc='upper right')
plt.savefig('/content/drive/MyDrive/FinalThesisMTC/RET/v2/RMSD_HistogramRenamed.png', dpi=600)
