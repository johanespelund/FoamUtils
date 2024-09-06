from setuptools import setup


def readfile(filename):
    with open(filename, "r") as f:
        return f.read()

def parse_requirements(filename):
    """Load requirements from a pip requirements file."""
    with open(filename, 'r') as f:
        return f.read().splitlines()

setup(
    name="FoamUtils",
    version="0.1",
    description="Utility functions for OpenFOAM",
    long_description=readfile('README.md'),
    author="Johan R. Espelund",
    author_email="johan.espelund@sintef.no",
    url="https://github.com/johanespelund/FoamUtils",
    py_modules=[
        "plot_residuals",
        "ThermophysicalProperties",
        "file_utils",
        "wall_profiles",
    ],
    install_requires=parse_requirements('requirements.txt'),
    entry_points={
        "console_scripts": [
            "plot_res = FoamUtils.plot_residuals:main",
            "plot_thermo = FoamUtils.ThermophysicalProperties:main",
        ]
    },
)
