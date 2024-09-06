import os


def get_sorted_times(dir):
    """
    Get the sorted time directories in a given directory,
    e.g. postProcessing/wallHeatFlux/
    """
    # return [float(string) for string in os.listdir(dir)]
    times = os.listdir(dir)
    return sorted(times, key=lambda x: float(x))
