"""Thin wrapper around a released ChromBPNet (bias-corrected) Keras model.

This is the ONE module that needs the GPU/TF environment (chrombpnet's tensorflow) and a
downloaded model, so it cannot be unit-tested in a CPU-only env. It is written against the
documented ChromBPNet model interface; the two things to CONFIRM against an actual extracted
ENCODE tar are marked CONFIRM below:

  1. the filename of the bias-corrected model inside the tar (glob '*nobias*.h5' by default), and
  2. the output ordering [profile_logits, logcounts] (ChromBPNet's standard head order).

The bias-corrected ("nobias") model is the right one for variant effects: it predicts accessibility
with the enzyme/Tn5 bias already removed, so ref-vs-alt differences reflect regulatory sequence.
Input: (N, 2114, 4) one-hot. Outputs: profile logits (N, 1000) and log-counts (N, 1).
"""
from __future__ import annotations
from pathlib import Path
import glob
import numpy as np


class ChromBPNetModel:
    INPUT_LEN = 2114
    PROFILE_LEN = 1000

    def __init__(self, model_dir, model_glob="*nobias*.h5"):
        self.model_dir = Path(model_dir)
        matches = sorted(glob.glob(str(self.model_dir / "**" / model_glob), recursive=True))
        if not matches:
            raise FileNotFoundError(
                f"No bias-corrected model matching {model_glob!r} under {self.model_dir}. "
                f"CONFIRM the filename inside the extracted ENCODE tar and pass model_glob.")
        self.model_path = matches[0]
        self._model = None

    def _load(self):
        if self._model is None:
            # compile=False avoids needing ChromBPNet's custom loss objects just for inference.
            from tensorflow import keras
            self._model = keras.models.load_model(self.model_path, compile=False)
        return self._model

    def predict(self, onehot: np.ndarray, batch_size: int = 256):
        """onehot: (N, 2114, 4) -> (profile_logits (N,1000), logcounts (N,)).

        CONFIRM head order: ChromBPNet returns [profile, counts]. If a given release swaps them,
        flip the two lines below.
        """
        onehot = np.asarray(onehot, dtype=np.float32)
        if onehot.shape[1] != self.INPUT_LEN:
            raise ValueError(f"expected input length {self.INPUT_LEN}, got {onehot.shape[1]}")
        out = self._load().predict(onehot, batch_size=batch_size, verbose=0)
        profile_logits, logcounts = out[0], out[1]
        return np.asarray(profile_logits), np.asarray(logcounts).reshape(-1)
