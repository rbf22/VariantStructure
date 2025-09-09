import math
from bisect import bisect_left
from collections import defaultdict

class DunbrackLibrary:
    """
    Parses backbone-dependent Dunbrack rotamer library file.
    Expects lines: RES φ_bin ψ_bin rotamer_name prob chi1 chi2... chiN sigma1 sigma2... sigmaN
    """
    def __init__(self, filepath):
        self.lib = defaultdict(list)
        with open(filepath) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                res, phi, psi, rot, prob = parts[0], float(parts[1]), float(parts[2]), parts[3], float(parts[4])
                # assume fixed number of chi and sigma following
                # Example: chi1 chi2 chi3 sigma1 sigma2 sigma3 ...
                chi = tuple(float(x) for x in parts[5:5+4])  # adapt length
                sig = tuple(float(x) for x in parts[5+4:5+8])
                self.lib[res].append({'phi': phi, 'psi': psi, 'rot': rot, 'prob': prob, 'chi': chi, 'sigma': sig})
        # Build index for fast lookup if needed

    def get_rotamers(self, three_letter, phi, psi, dphi=5, dpsi=5):
        """
        Return list of rotamers (rotamer dicts) whose φ,ψ bins are within dphi, dpsi of given angles.
        """
        candidates = self.lib.get(three_letter, [])
        res = []
        for entry in candidates:
            if abs(((phi - entry['phi'] + 180) % 360) - 180) <= dphi and \
               abs(((psi - entry['psi'] + 180) % 360) - 180) <= dpsi:
                res.append(entry)
        return res