def get_delta_for_para(key):
    if key.endswith('r0'):
        delta = 0.05
    elif key.endswith('e0'):
        delta = 0.001
    elif key.endswith('bi'):
        delta = 0.005

    # temperature dependence
    elif key.endswith('dr'):
        delta = 0.01
    elif key.endswith('de'):
        delta = 0.02
    elif key.endswith('dl'):
        delta = 0.01
    elif key.endswith('d2'):
        delta = 0.01

    else:
        raise Exception('Unknown parameter: ' + key)

    return delta

def get_bound_for_para(key):
    if key.endswith('r0'):
        bound = (3, 5)
        if key.startswith('h_'):
            bound = (2, 3)
    elif key.endswith('e0'):
        bound = (0.01, 0.5)
        if key.startswith('h_'):
            bound = (0.01, 0.05)
    elif key.endswith('bi'):
        bound = (-0.6, 0.6)

    # temperature dependence
    elif key.endswith('dr'):
        bound = (-0.1, 0.01)
    elif key.endswith('de'):
        bound = (-0.01, 0.2)
    elif key.endswith('dl'):
        bound = (-0.01, 0.1)
    elif key.endswith('d2'):
        bound = (-0.01, 0.1)

    else:
        raise Exception('Unknown parameter: ' + key)

    return bound


