import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyemma
import re 
import argparse 
from pyemma.util.contexts import settings
from misc.cv_err import MiscFunctions
from pathlib import Path


def vamp_bars(pdb, files):


    file = MiscFunctions.get_filename(pdb)
    outputname = file + "_vampBars.jpg"
    
    torsions_feat = pyemma.coordinates.featurizer(pdb)
    torsions_feat.add_backbone_torsions(cossin=True, periodic=False)
    torsions_data = pyemma.coordinates.load(files, features=torsions_feat)
    labels = ['backbone\ntorsions']

    positions_feat = pyemma.coordinates.featurizer(pdb)
    positions_feat.add_selection(positions_feat.select_Backbone())
    positions_data = pyemma.coordinates.load(files, features=positions_feat)
    labels += ['backbone atom\npositions']

    dim = 1

    fig, axes = plt.subplots(1, 3, figsize=(12, 3), sharey=True)

    for ax, lag in zip(axes.flat, [5, 10, 20]):
      torsions_scores = MiscFunctions.score_cv(torsions_data, lag=lag, dim=dim)
      scores = [torsions_scores.mean()]
      errors = [torsions_scores.std()]
      positions_scores = MiscFunctions.score_cv(positions_data, lag=lag, dim=dim)
      scores += [positions_scores.mean()]
      errors += [positions_scores.std()]
      ax.bar(labels, scores, yerr=errors, color=['C0', 'C1'])
      ax.set_title(r'lag time $\tau$={:.1f}ns'.format(lag * 0.1))
      if lag == 5:
        vamp_bars_plot = dict(
            labels=labels, scores=scores, errors=errors, dim=dim, lag=lag)
    axes[0].set_ylabel('VAMP2 score')
    fig.tight_layout()
    fig.savefig(outputname, dpi=600)

    return positions_data 


def cluster_kmeans(pdata, output):
    
    tica = pyemma.coordinates.tica(pdata, lag=100)
    tica_output = tica.get_output()
    tica_concatenated = np.concatenate(tica_output)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    pyemma.plots.plot_density(*tica_concatenated[:, :2].T, ax=axes[1], logscale=True)
    axes[1].set_xlabel('IC 1')
    axes[1].set_ylabel('IC 2')
    fig.tight_layout()
    fig.savefig(output, dpi=600)

    return tica_output, tica_concatenated


def profile(ticaData, output):
    
    fig, axes = plt.subplots(4, 1, figsize=(12, 5), sharex=True)
    x = 0.1 * np.arange(ticaData[0].shape[0])
    for i, (ax, tic) in enumerate(zip(axes.flat, ticaData[0].T)):
        ax.plot(x, tic)
        ax.set_ylabel('IC {}'.format(i + 1))
    axes[-1].set_xlabel('time / ns')
    fig.tight_layout()
    fig.savefig(output, dpi=600)


def free_energy(tica, label1, label2, output):

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    cluster = pyemma.coordinates.cluster_kmeans(
    tica, k=50, max_iter=50, stride=10, fixed_seed=1)
    msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=100, dt_traj='0.1 ns')
    pyemma.plots.plot_contour(
        *tica[:, :2].T,
        msm.pi[np.concatenate(msm.dtrajs_active)],
        ax=axes[0],
        mask=True,
        cbar_label='stationary distribution')
    pyemma.plots.plot_free_energy(
        *tica[:, :2].T,
        weights=np.concatenate(msm.trajectory_weights()),
        ax=axes[1],
        legacy=False)
    for ax in axes.flat:
        ax.set_xlabel('IC 1')
    axes[0].set_ylabel('IC 2')
    axes[0].set_title(label1, fontweight='bold')
    axes[1].set_title(label2, fontweight='bold')
    fig.tight_layout()
    fig.savefig(output, dpi=600)


if __name__ == "__main__":
     parser = argparse.ArgumentParser()
     parser.add_argument("-pdb")
     parser.add_argument("-trajectory")

     args = parser.parse_args()


     file = MiscFunctions.get_filename(args.pdb)
     outputname = file + "_clusterKmeans.jpg"
     
     pdb_path = args.pdb
     trajectory_path = args.trajectory
     pdata = vamp_bars(pdb_path, trajectory_path)
     tdata, tconc = cluster_kmeans(pdata, outputname)

     outputname = file + "_cProfile.jpg"
     profile(tdata, outputname)

     outputname = file + "_fEnergy.jpg"
     free_energy(tconc, "M918T - Stationary Distribution", "M918T - Reweighted free energy surface", outputname)