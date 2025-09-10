"""A module for loading the Dunbrack rotamer library."""
import tarfile
import os
import collections
import numpy as np
import pandas as pd


class DunbrackLibrary:
    """
    Parses and stores the Dunbrack rotamer library.
    """

    _COLUMNS = collections.OrderedDict(
        [
            ("T", np.str_),
            ("Phi", np.int64),
            ("Psi", np.int64),
            ("Count", np.int64),
            ("r1", np.int64),
            ("r2", np.int64),
            ("r3", np.int64),
            ("r4", np.int64),
            ("Probabil", np.float64),
            ("chi1Val", np.float64),
            ("chi2Val", np.float64),
            ("chi3Val", np.float64),
            ("chi4Val", np.float64),
            ("chi1Sig", np.float64),
            ("chi2Sig", np.float64),
            ("chi3Sig", np.float64),
            ("chi4Sig", np.float64),
        ]
    )

    def __init__(self):
        self.lib = self._load_rotamer_library()

    def _load_rotamer_library(self):
        """
        Loads the rotamer library from the tarball.
        """
        amino_acids = [
            "arg", "asn", "asp", "cys", "gln", "glu", "his", "ile",
            "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp",
            "tyr", "val",
        ]
        db = {}

        data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
        tar_path = os.path.join(data_dir, "dunbrack_rotamer.tar.gz")

        try:
            with tarfile.open(tar_path, "r:gz") as tar:
                lib_path = "dunbrack-rotamer/original/SimpleOpt1-0/ALL.bbdep.rotamers.lib"
                f = tar.extractfile(lib_path)
                if f is not None:
                    df = pd.read_csv(
                        f,
                        names=list(self._COLUMNS.keys()),
                        dtype=self._COLUMNS,
                        comment="#",
                        sep=r"\s+",
                        engine="c",
                    )
                    for amino_acid in amino_acids:
                        db[amino_acid.upper()] = df[df["T"] == amino_acid.upper()]
        except (KeyError, FileNotFoundError, tarfile.ReadError) as e:
            print(f"Warning: Could not load rotamer library from {tar_path}. Error: {e}")
        return db

    def get_rotamers(  # pylint: disable=too-many-arguments
        self, three_letter, phi, psi, dphi=10, dpsi=10
    ):
        """
        Return list of rotamers whose phi,psi bins are within dphi, dpsi of given angles.
        """
        df = self.lib.get(three_letter.upper())
        if df is None:
            return []

        # Find the closest phi/psi bin
        phi_dist = np.abs(((phi - df["Phi"] + 180) % 360) - 180)
        psi_dist = np.abs(((psi - df["Psi"] + 180) % 360) - 180)

        filtered_df = df[(phi_dist <= dphi) & (psi_dist <= dpsi)]

        if filtered_df.empty:
            # If no rotamers are found in the bin, find the single closest one
            dist_sq = phi_dist**2 + psi_dist**2
            closest_idx = dist_sq.idxmin()
            filtered_df = df.loc[[closest_idx]]

        return self._dataframe_to_rotamer_dicts(filtered_df)

    def _dataframe_to_rotamer_dicts(self, df):
        """Convert a dataframe of rotamers to a list of dictionaries."""
        rotamers = []
        for _, row in df.iterrows():
            chi_vals = [row["chi1Val"], row["chi2Val"], row["chi3Val"], row["chi4Val"]]
            chi_sigs = [row["chi1Sig"], row["chi2Sig"], row["chi3Sig"], row["chi4Sig"]]

            num_chi = 0
            for c in chi_vals:
                if not np.isnan(c):
                    num_chi += 1

            rotamers.append(
                {
                    "phi": row["Phi"],
                    "psi": row["Psi"],
                    "prob": row["Probabil"],
                    "chi": tuple(chi_vals[:num_chi]),
                    "sigma": tuple(chi_sigs[:num_chi]),
                }
            )
        return rotamers
