from setuptools import setup


def readfile(filename):
    with open(filename, "r+") as f:
        return f.read()


setup(
    name="FoamUtils",
    version="0.1",
    description="Utility functions for OpenFOAM",
    # long_description=readfile('README.md'),
    author="Johan R. Espelund",
    author_email="johan.espelund@sintef.no",
    url="",
    py_modules=[
        "plot_residuals",
        "ThermophysicalProperties",
        "file_utils",
        "wall_profiles",
    ],
    # license=readfile('LICENSE'),
    entry_points={
        "console_scripts": [
            "plot_res = FoamUtils.plot_residuals:main",
            "plot_thermo = FoamUtils.ThermophysicalProperties:main",
        ]
    },
)
