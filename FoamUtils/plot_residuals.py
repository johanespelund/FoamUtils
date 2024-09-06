import click
import matplotlib.pyplot as plt
import pandas as pd
import scienceplots

plt.style.use(["science", "nature"])

### Definition of the plot function. Options are passed as arguments to the main function ###


def plot(dir, res_type, save, dpi, ncol, legend_box, time_unit, time_folder):
    # Define the file path
    file_path = f"{dir}/solverInfo/{time_folder}/solverInfo.dat"

    # Read the first line to get the column names
    with open(file_path, "r") as f:
        f.readline()
        columns = f.readline().strip().split()
        # Remove the '#' from the column names
        columns.pop(0)

    # Read the data file
    data = pd.read_csv(file_path, skiprows=2, delim_whitespace=True, names=columns)

    # Convert time to desired unit
    if time_unit == "min":
        data["Time"] /= 60
    elif time_unit == "h":
        data["Time"] /= 3600

    # Select the desired columns
    selected_columns = [
        col for col in data.columns if col == "Time" or col.endswith(res_type)
    ]
    selected_data = data[selected_columns]

    # Plot the initial residuals
    ax = selected_data.plot(x="Time", logy=True, lw=0.5)
    plt.legend(ncol=ncol)
    if legend_box:
        plt.legend(frameon=True, fancybox=True, ncol=ncol)

    plt.ylabel("Scaled residuals [-]")
    plt.xlabel(f"Time [{time_unit}]")
    plt.tight_layout()
    plt.grid()

    if save:
        plt.savefig(f"{dir}/residuals_{res_type}.{save}", dpi=dpi)
        print(f"Saved plot to {dir}/residuals_{res_type}.{save}")

    plt.show()


### Definition of the command line options ###


@click.command()
@click.option("--dir", default="postProcessing", help="Which directory to plot from.")
@click.option(
    "--res-type",
    default="initial",
    help='Which residual type to plot: "initial" or "final".',
)
@click.option("--save", default="", help="Save the plot to a specified filetype.")
@click.option("--dpi", default=600, help="Resolution of image (DPI).")
@click.option("--ncol", default=1, help="Number of columns in the plot.")
@click.option("--legend-box", is_flag=True, help="Draw box around legend.")
@click.option(
    "--time-unit", type=click.Choice(["s", "min", "h"]), default="s", help="Time unit."
)
@click.option("--time-folder", default=0, help="From which time to plot.")
def main(dir, res_type, save, dpi, ncol, legend_box, time_unit, time_folder):
    plot(**locals())


if __name__ == "__main__":
    main()
