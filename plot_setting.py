import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

# Use the pgf backend (must be done before import pyplot interface)
inches_per_pt = 1 / 72.27
textwidth = 448.1309 * inches_per_pt
golden_ratio = (5**.5 - 1) / 2
# print(f'Golden ratio = {golden_ratio:.2f}')
slidewidth = 13.33  # inches 
slideheight = 7.5  # inches

settings = {
    'axes.grid': True,
    'axes.grid.which': 'both',
    'xtick.color': '#DDDDDD',
    'ytick.color': '#DDDDDD',
    "axes.labelsize": "medium",
    "axes.titlesize": "medium",
    "figure.labelsize": "medium",
    "figure.titlesize": "medium",
    # Make the legend/label fonts a little smaller
    "legend.fontsize": "small",
    "legend.title_fontsize": "small",
    "xtick.labelsize": "small",
    "ytick.labelsize": "small",
    "lines.linewidth": 0.8,
}

report_fast = settings | {
    "font.size": 9,
}

report_tex = settings | {
    "text.usetex": True,
    "font.family": "Helvetica",
    "font.size": 12,
}

# Following doesnt allow previewing
report_pgf_tex = settings | {
    "font.family": "sans-serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
}

pressentation_tex = report_tex | {
    "font.size": 14,
}


# \textwidth of the document can be determine by \showthe\textwidth and checking the logs
# for our document it las was 448.1309 pts
def set_size(width=textwidth, height=textwidth*golden_ratio, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
        Document width in inches, or string of predined document type
    height: float, optional
        Height of the figure in inches. Default is ``width * golden_ratio``
    subplots: array-like, optional
        The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
        Dimensions of figure in inches

    Source
    ------
    https://jwalton.info/Embed-Publication-Matplotlib-Latex/
    """


    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2
    # golden_ratio = (1 + golden_ratio) * 0.5

    # Figure height in inches
    if height is None:
        fig_height = width * golden_ratio * (subplots[0] / subplots[1])
    else:
        fig_height = height    

    return (width, fig_height)
