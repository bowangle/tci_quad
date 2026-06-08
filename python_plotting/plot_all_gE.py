from extern.tci_quad.python_plotting.plotting_error import load_dat, plot_tt_comparison

from pathlib import Path
import matplotlib.pyplot as plt

for file in Path(".").rglob("*_error.dat"):
    a = load_dat(str(file))

    filename = str(file).removesuffix("_error.dat")

    plot_tt_comparison(a, filename=filename)

    plt.close()