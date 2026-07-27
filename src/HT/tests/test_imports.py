"""Tests for optional plotting dependencies."""

import os
from pathlib import Path
import subprocess
import sys
import textwrap


def test_numerical_core_imports_without_matplotlib():
    """The model and Hamiltonian are usable in a NumPy-only environment."""

    repository_root = Path(__file__).resolve().parents[3]
    subprocess.run(
        [
            sys.executable,
            "-c",
            textwrap.dedent(
                """
                import sys

                sys.modules["matplotlib"] = None

                import numpy as np
                from HT.models.hofstadter import Hofstadter

                model = Hofstadter(1, 3, lat="square")
                num_bands, *_ = model.unit_cell()
                hamiltonian = model.hamiltonian(np.array([0.0, 0.0]))

                assert num_bands == 3
                assert hamiltonian.shape == (3, 3)
                assert np.isfinite(hamiltonian).all()
                assert np.allclose(hamiltonian, hamiltonian.conj().T)
                """
            ),
        ],
        check=True,
        cwd=repository_root,
        env={
            **os.environ,
            "PYTHONPATH": str(repository_root / "src"),
        },
    )
