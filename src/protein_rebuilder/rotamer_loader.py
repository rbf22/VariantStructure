from collections import defaultdict

class DunbrackLibrary:
    """
    A dummy Dunbrack rotamer library.
    """
    def __init__(self):
        self.lib = defaultdict(list)
        # Create a dummy library for a few amino acids
        self.lib['VAL'] = [
            {'phi': -60, 'psi': 140, 'rot': '1', 'prob': 0.5, 'chi': (60,), 'sigma': (5,)},
            {'phi': -140, 'psi': 140, 'rot': '2', 'prob': 0.3, 'chi': (180,), 'sigma': (5,)},
            {'phi': -70, 'psi': -40, 'rot': '3', 'prob': 0.2, 'chi': (-60,), 'sigma': (5,)},
        ]
        self.lib['LEU'] = [
            {'phi': -60, 'psi': 140, 'rot': '1', 'prob': 0.5, 'chi': (60, 60), 'sigma': (5, 5)},
            {'phi': -140, 'psi': 140, 'rot': '2', 'prob': 0.3, 'chi': (180, 60), 'sigma': (5, 5)},
            {'phi': -70, 'psi': -40, 'rot': '3', 'prob': 0.2, 'chi': (-60, 180), 'sigma': (5, 5)},
        ]

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