import tarfile
import os
import collections
import numpy as np
import pandas as pd

class DunbrackLibrary:
    """
    Parses and stores the Dunbrack rotamer library.
    """
    def __init__(self):
        self.lib = self._load_rotamer_library()

    def _load_rotamer_library(self):
        """
        Loads the rotamer library from the tarball.
        """
        amino_acids = [
            "arg", "asn", "asp", "cys", "gln", "glu", "his", "ile", "leu",
            "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val"
        ]
        db = {}

        columns = collections.OrderedDict()
        columns["T"] = np.str_
        columns["Phi"] = np.int64
        columns["Psi"] = np.int64
        columns["Count"] = np.int64
        columns["r1"] = np.int64
        columns["r2"] = np.int64
        columns["r3"] = np.int64
        columns["r4"] = np.int64
        columns["Probabil"] = np.float64
        columns["chi1Val"] = np.float64
        columns["chi2Val"] = np.float64
        columns["chi3Val"] = np.float64
        columns["chi4Val"] = np.float64
        columns["chi1Sig"] = np.float64
        columns["chi2Sig"] = np.float64
        columns["chi3Sig"] = np.float64
        columns["chi4Sig"] = np.float64

        data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
        tar_path = os.path.join(data_dir, "dunbrack_rotamer.tar.gz")

        with tarfile.open(tar_path, "r:gz") as tar:
            lib_path = "dunbrack-rotamer/original/SimpleOpt1-0/ALL.bbdep.rotamers.lib"
            try:
                f = tar.extractfile(lib_path)
                if f is not None:
                    df = pd.read_csv(
                        f,
                        names=list(columns.keys()),
                        dtype=columns,
                        comment="#",
                        sep='\s+',
                        engine="c",
                    )
                    for amino_acid in amino_acids:
                        db[amino_acid.upper()] = df[df['T'] == amino_acid.upper()]
            except KeyError:
                print(f"Warning: Could not find {lib_path} in the rotamer library tarball.")
        return db

    def get_rotamers(self, three_letter, phi, psi, dphi=10, dpsi=10):
        """
        Return list of rotamers (rotamer dicts) whose φ,ψ bins are within dphi, dpsi of given angles.
        """
        df = self.lib.get(three_letter.upper())
        if df is None:
            return []

        # Find the closest phi/psi bin
        phi_dist = np.abs(((phi - df['Phi'] + 180) % 360) - 180)
        psi_dist = np.abs(((psi - df['Psi'] + 180) % 360) - 180)

        filtered_df = df[(phi_dist <= dphi) & (psi_dist <= dpsi)]

        if filtered_df.empty:
            # If no rotamers are found in the bin, find the single closest one
            dist_sq = phi_dist**2 + psi_dist**2
            closest_idx = dist_sq.idxmin()
            filtered_df = df.loc[[closest_idx]]

        rotamers = []
        for _, row in filtered_df.iterrows():
            chi_vals = [row['chi1Val'], row['chi2Val'], row['chi3Val'], row['chi4Val']]
            chi_sigs = [row['chi1Sig'], row['chi2Sig'], row['chi3Sig'], row['chi4Sig']]

            num_chi = 0
            for c in chi_vals:
                if not np.isnan(c):
                    num_chi += 1

            rotamers.append({
                'phi': row['Phi'],
                'psi': row['Psi'],
                'prob': row['Probabil'],
                'chi': tuple(chi_vals[:num_chi]),
                'sigma': tuple(chi_sigs[:num_chi])
            })
        return rotamers