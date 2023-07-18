import matplotlib.pyplot as plt

fig, axs = plt.subplots(4, 2)
gs = axs[3, 0].get_gridspec()
# remove the underlying axes
for ax in axs[3, 0:]:
    ax.remove()
axbig = fig.add_subplot(gs[3, 0:])
#axbig.annotate('Big Axes \nGridSpec[1:, -1]', (0.1, 0.5),
#               xycoords='axes fraction', va='center')

fig.tight_layout()

plt.show()
