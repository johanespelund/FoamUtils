import os
import numpy as np


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

def load_wallHeatFlux(filename, start=0, end=-1):
    wall_heat_flows = {}
    times = set()

    # Find unique patch names
    lines = open(filename, "r").readlines()
    for line in lines:
        if line[0] != "#":
            line = line.rstrip().split("\t")
            patch_name = line[1]
            if patch_name not in wall_heat_flows:
                wall_heat_flows[patch_name] = []
            else:
                break

    # Read wall heat flows
    lines = open(filename, "r").readlines()
    for line in lines[
        : -len(wall_heat_flows)
    ]:  # Last time is for some reason written twice, so don't read last entries
        if line[0] != "#":
            line = line.rstrip().split("\t")
            time = float(line[0])
            patch_name = line[1]
            times.add(time)
            wall_heat_flows[patch_name].append(float(line[-1]))

    times = np.array(sorted(list(times)))[start:end]
    for key in wall_heat_flows:
        wall_heat_flows[key] = np.array(wall_heat_flows[key])[start:end]

    return times, wall_heat_flows

