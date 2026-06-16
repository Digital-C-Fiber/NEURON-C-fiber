# I used this script to iteratively find candidates with closer alignment after pertubating only one conductance in analysis-COVID.ipynb.

import neuron
import main_CMi
from multiprocessing import Pool

neuron.load_mechanisms('./MOD_Tigerholm')


params_orig = {
    'gPump': -0.00485891709589456,
    'gNav17Parent': 0.22831806252579703,
    'gNav17': 0.33198932295899997,
    'gNav18Parent': 0.2411879110860178,
    'gNav18': 0.04520936848183172,
    'gNav19': 0.0032208034421239012,
    'gKs': 0.0030962055552702606,
    'gKf': 0.025994684685553986,
    'gH': 0.009167424690232623,
    'gKdr': 0.013875365047132446,
    'gKna': 0.0014444412024341238,
    'vRest': -55.0
}

# gNav17/gNav18 scaling also to parents
experiments = [
    {"gNav17": 0.3, "gNav18": 0.30},
    {"gNav17": 0.2, "gNav18": 0.20, "gH": 0.20},
    {"gNav17": 0.2, "gNav18": 0.30, "gH": 0.30},
    {"gNav17": 0.3, "gNav18": 0.30, "gH": 0.40},
    {"gNav17": 0.2, "gNav18": 0.30, "gKdr": 0.20},
    {"gNav17": 0.3, "gNav18": 0.30, "gKs": -0.20},
    {"gNav17": 0.3, "gNav18": 0.30, "gKs": -0.20, "gKdr": 0.20},
    {"gNav17": 0.3, "gNav18": 0.30, "gKs": -0.20, "gKdr": 0.20, "gH": 0.40},
    {"gNav17": 0.3, "gNav18": 0.30, "gKs": -0.20, "gKdr": 0.20, "gH": 0.30},
    {"gNav17": 0.3, "gNav18": 0.30, "gKdr": 0.20, "gH": 0.30},
    {"gNav17": 0.3, "gNav18": 0.30, "gKs": -0.20, "gH": 0.30},
    {"gNav17": 0.3, "gNav18": 0.30, "gKs": -0.20, "gKdr": 0.20, "gH": 0.40},
    {"gNav17": 0.3, "gNav18": 0.30, "gKs": -0.10, "gKdr": 0.10, "gH": 0.20}
]

PARENT_COUPLED = {"gNav17": "gNav17Parent", "gNav18": "gNav18Parent"}

SIM_KWARGS = {
    "prot":          42,
    "sine":          False,
    "scalingFactor": 0.1,
}


def apply_scales(base, scales):
    """Scales ionic conductances while accounting for the fact that some conductances differ in parent branch"""
    p = base.copy()
    for name, dg in scales.items():
        p[name] = base[name] * (1 + dg)
        if name in PARENT_COUPLED:
            parent = PARENT_COUPLED[name]
            p[parent] = base[parent] * (1 + dg)
    return p



params = [apply_scales(params_orig, exp) for exp in experiments]
print(f"Running {len(params)} simulations.")


if __name__ == "__main__":
    with Pool(processes=5) as pool:
        async_results = [
            pool.apply_async(main_CMi.run, kwds={**SIM_KWARGS, **p})
            for p in params
        ]
        output = [r.get() for r in async_results]
