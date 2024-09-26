import os
import numpy as np
import pandas as pd


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


def load(filename):
    data = np.loadtxt(filename, comments="#")
    times = data[start:end, 0]
    values = [data[start:end, i] for i in range(1, data.shape[1])]
    return (times, *values)


def load_all_times(times_dir, use_latest_run=True):
    """
    Load data from all time directories in a given directory.
    - times_dir: Directory containing time directories, e.g. "postProcessing/forces"
    - use_latest_run: If True, only load the latest run in each time directory, e.g. "wallHeatFlux_0.dat"
    Returns a merged pandas DataFrame with all data.
    """
    times = get_sorted_times(times_dir)
    data = pd.DataFrame()

    for time in times:
        if use_latest_run:
            files = os.listdir(f"{times_dir}/{time}")
        else:
            files = os.listdir(f"{times_dir}/{time}")
            files = [f for f in files if not f.endswith(f"_{time}")]
        for file in files:
            filename = f"{times_dir}/{time}/{file}"

            # The last row with comments are the column names
            # Need to read this separately to get the column names

            with open(filename, "r") as f:
                for line in f:
                    if line[0] != "#":
                        break
            column_names = line[1:].split()
            df = pd.read_csv(filename, comment="#", names=column_names, delim_whitespace=True)

            df["time"] = float(time)
            data = pd.concat([data, df], axis=0)
    return data


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


def load_sampling_set(sampling_dir, time, set_name):
    """
    Load data from OpenFOAM sample sets.
    Parameters:
    - sampling_dir: Name of the sampling directory, e.g. "postProcessing/linesampling"
    - time: Time directory, e.g. "0.1"
    - set_name: Name of the set, e.g. "centerline"

    Returns:
    - pandas DataFrame with the data
    """
    files = os.listdir(f"{sampling_dir}/{time}/")
    samples = [f.split("_")[0] for f in files]
    file_type = ''


    for file in files:
        if f"{set_name}_" in file:
            filename = f"{sampling_dir}/{time}/{file}"
            filetype = file.split(".")[1]
            
            if filetype == "csv":
                df = pd.read_csv(filename)
                df.columns = [str(col) for col in df.columns]
                df.columns = [col.replace("_0", "x").replace("_1", "y").replace("_2", "z") for col in df.columns]


            elif filetype == "xy":
                fields_string = file.replace(f"{set_name}_", "").split(".")[0]
                fields = ["z"] + fields_string.split("_")
                raw_data = np.loadtxt(filename)
                data = {}

                # Split vector fields into components
                offset = 0
                for i in range(len(fields)):
                    j = i + offset
                    field_name = fields[i]
                    if "U" in field_name:
                        data[f"{field_name}x"] = raw_data[:, j]
                        offset += 1
                        j = i + offset
                        data[f"{field_name}y"] = raw_data[:, j]
                        offset += 1
                        j = i + offset
                        data[f"{field_name}z"] = raw_data[:, j]
                    else:
                        data[field_name] = raw_data[:, j]

                df = pd.DataFrame(data)

            # Remove rows that have duplicated values in the first column
            df.drop_duplicates(subset=df.columns[0], inplace=True)

            return df

    return data

