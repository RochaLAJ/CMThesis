import MDAnalysis as mda
import warnings
import pandas as pd
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import align
from pmda.rms import RMSF
from google.colab import drive
#drive.mount("/content/drive", force_remount=True)
warnings.filterwarnings("ignore")
%matplotlib inline

data = pd.read_csv("/content/drive/MyDrive/FinalThesisMTC/RET/v2/Tables/rmsf.csv")
data['Residue'] = data['Residue'] + 698
data.head()

plt.plot(data["Residue"], data["Wildtype"], color="darkgreen", label="Wildtype")
plt.plot(data["Residue"], data["M918V"], color="orange", label="M918V")
plt.plot(data["Residue"], data["M918T"], color="red", label="M918T")
plt.xlabel('Residue number')
plt.ylabel('RMSF ($\AA$)')
plt.title("RMSF Comparison Between \n Wildtype, p.Met918Val and p.Met918Thr", fontweight="bold")
#plt.legend(loc='upper center', ncol=3)

plt.axvspan(700, 729, zorder=0, alpha=0.2, color='red')
plt.axvspan(730, 750, zorder=0, alpha=0.2, color='orange')
plt.axvspan(751, 800, zorder=0, alpha=0.2, color='red')


plt.axvspan(801, 870, zorder=0, alpha=0.2, color='red')
plt.axvspan(870, 894, zorder=0, alpha=0.2, color='red')


plt.axvspan(895, 905, zorder=0, alpha=0.2, color='darkgreen')
plt.axvspan(906, 934, zorder=0, alpha=0.2, color='red')
plt.axvspan(934, 950, zorder=0, alpha=0.2, color='red')

plt.axvspan(951, 954, zorder=0, alpha=0.2, color='darkgreen')

plt.axvspan(955, 1010, zorder=0, alpha=0.2, color='red')
plt.axvspan(1010, 1015, zorder=0, alpha=0.2, color='darkgreen')
plt.savefig('/content/drive/MyDrive/FinalThesisMTC/RET/v2/Figures/Fig3_RMSFRenamed', dpi=600)
