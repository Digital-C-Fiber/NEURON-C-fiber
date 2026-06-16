import neuron
import main_CM
from multiprocessing import Pool

neuron.load_mechanisms('./MOD_Tigerholm')


params_orig = {
    "gPump":        -0.0025,
    "gNav17":        0.24686453257354574,
    "gNav17Parent":  0.13115152763095123,
    "gNav18":        0.37673567973121774,
    "gNav18Parent":  0.2343920300579889,
    "gNav19":        0.00017254238997420438,
    "gKs":           0.008865226128662577,
    "gKf":           0.02709394494148292,
    "gH":            0.014140202887083592,
    "gKdr":          0.008469950837206652,
    "gKna":          0.001398204170298818,
    "vRest":        -55,
}


# gNav17/gNav18 scaling is automatically applied to their Parent counterparts
# for doing a few simulations to test out
experiments = [
    {"gNav17": -0.40, "gNav18": -0.30, "gH": -0.40},
    {"gNav17": -0.40, "gNav18": -0.40, "gH": -0.20},
    {"gNav17": -0.30, "gNav18": -0.30, "gH": -0.40},
    {"gNav17": -0.30, "gNav18": -0.40, "gH": -0.20},
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
            pool.apply_async(main_CM.run, kwds={**SIM_KWARGS, **p})
            for p in params
        ]
        output = [r.get() for r in async_results]
