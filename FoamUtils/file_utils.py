import os


def get_sorted_times(dir):
    """
    Get the sorted time directories in a given directory,
    e.g. postProcessing/wallHeatFlux/
    """
    # return [float(string) for string in os.listdir(dir)]
    times = os.listdir(dir)
    return sorted(times, key=lambda x: float(x))

def read_parameters(path="parameters"):
    """
    Read the parameters file, in OpenFOAM format.
    """
    p = {}
    with open(path, "r") as f:
        for line in f:
            if line.startswith("//") or line.startswith("/*"):
                continue
            if line.startswith("}"):
                break
            key, value = line.split(" ", 1)
            try:
                value = float(value.strip()[:-1])
            except ValueError:
                value = value.strip()[:-1]
            p[key] = value
    return p

