import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

def activateNicePlots():
    nice_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    }
    mpl.rcParams.update(nice_fonts)
    plt.rcParams['image.cmap'] = 'viridis'
    sns.set_palette(sns.color_palette())
