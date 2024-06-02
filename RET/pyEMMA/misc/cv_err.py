import pyemma
import numpy as np 
from pathlib import Path

class MiscFunctions:

	def score_cv(data, dim, lag, number_of_splits=10, validation_fraction=0.5):

		with pyemma.util.contexts.settings(show_progress_bars=False):
			nval = int(len(data) * validation_fraction)
			scores = np.zeros(number_of_splits)
			for n in range(number_of_splits):
				ival = np.random.choice(len(data), size=nval, replace=False)
				vamp = pyemma.coordinates.vamp(
                       [d for i, d in enumerate(data) if i not in ival], lag=lag, dim=dim)
				scores[n] = vamp.score([d for i, d in enumerate(data) if i in ival])
		return scores


	def its_separation_err(ts, ts_err):
		return ts[:-1] / ts[1:] * np.sqrt(
		        (ts_err[:-1] / ts[:-1])**2 + (ts_err[1:] / ts[1:])**2)


	def get_filename(pdb):
		return Path(pdb).stem




