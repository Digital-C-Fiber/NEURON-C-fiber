import neuron as neuron
import main_CM
from multiprocessing import Pool

neuron.load_mechanisms('./MOD_Tigerholm')

params = []


params_orig = {
    "gPump": -0.0025,
    "gNav17": 0.24686453257354574,
    "gNav17Parent": 0.13115152763095123,
    "gNav18": 0.37673567973121774,
    "gNav18Parent": 0.2343920300579889,
    "gNav19": 0.00017254238997420438,
    "gKs": 0.008865226128662577,
    "gKf": 0.02709394494148292,
    "gH": 0.014140202887083592,
    "gKdr": 0.008469950837206652,
    "gKna": 0.001398204170298818,
    "vRest": -55
}

changes = [-0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.6, 0.8]

changes_v_rest  = [-0.3, -0.2, -0.1, 0.1, 0.2, 0.4]

changes_gNav18  = [-0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.3]

params_dict = {}

# assign parameter sets
for name, val in params_orig.items():
    params_dict[name] = {}

    if name in ['gNav17Parent', 'gNav18Parent', 'gKf']:
        continue
    elif name in ['gNav17', 'gNav18']:
        parent_name = name + 'Parent'
        if name == 'gNav18':
            for dg in changes_gNav18:
                new_p = params_orig.copy()
                new_p[name] = params_orig[name] * (1 + dg)
                new_p[parent_name] = params_orig[parent_name] * (1 + dg)
                params_dict[name][dg] = new_p
                params.append(list(new_p.values()))
        else:
            for dg in changes:
                new_p = params_orig.copy()
                new_p[name] = params_orig[name] * (1 + dg)
                new_p[parent_name] = params_orig[parent_name] * (1 + dg)
                params_dict[name][dg] = new_p
                params.append(list(new_p.values()))
    elif name == 'vRest':
        for dg in changes_v_rest:
            new_p = params_orig.copy()
            new_p[name] = params_orig[name] * (1 + dg)
            params_dict[name][dg] = new_p
            params.append(list(new_p.values()))
    else:
        for dg in changes:
            new_p = params_orig.copy()
            new_p[name] = params_orig[name] * (1 + dg)
            params_dict[name][dg] = new_p
            params.append(list(new_p.values()))

print(f"Generated {len(params)} unique parameter combinations.")

if __name__ == '__main__':
    # Parallelize
    with Pool(processes=5) as pool: # 6 core cpu -> 5 processes so 1 is empty (hopefully?)
        results = []
        for param in params:
            # Create kwargs dictionary with varied parameters
            kwargs = {
                'prot': 42,
                'sine': False,
                'scalingFactor': 0.1,
                'gPump': param[0],
                'gNav17Parent': param[1],
                'gNav17': param[2],
                'gNav18Parent': param[3],
                'gNav18': param[4],
                'gNav19': param[5],
                'gKs': param[6],
                'gKf': param[7],
                'gH': param[8],
                'gKdr': param[9],
                'gKna': param[10],
                'vRest': param[11]
            }

            # Submit job to pool with keyword arguments
            result = pool.apply_async(main_CM.run, kwds=kwargs)
            results.append(result)

        # Wait for all results to complete
        output = [r.get() for r in results]