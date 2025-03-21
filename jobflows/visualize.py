import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

def _figure_settings():
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=26)
    plt.rc('ytick', labelsize=26)

    fig = plt.figure(figsize=(22, 12))
    ax = fig.add_subplot(1, 1, 1)

    border_width = 1.5
    ax.spines['top'].set_linewidth(border_width)
    ax.spines['right'].set_linewidth(border_width)
    ax.spines['bottom'].set_linewidth(border_width)
    ax.spines['left'].set_linewidth(border_width)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction='in', length=8, width=border_width, which='major',
                   top=True, right=True, pad=10)
    ax.tick_params(direction='in', length=4, width=border_width, which='minor',
                   top=True, right=True, pad=10)
    return fig, ax

def plot_msd(lagtimes, msd_values, annotate_times=None):
    
    fig, ax = _figure_settings()
    
    ax.plot(lagtimes, msd_values, color='black', linewidth=3.0, label=r'MSD ($\mathrm{\AA}^2$)')

    ax.set_xlabel(r"$\rm \tau$ (ps)", fontsize=26, labelpad=15)
    ax.set_ylabel(r"MSD ($\mathrm{\AA}^2$)", fontsize=26, color='black', labelpad=15)
    ax.tick_params(axis='y', colors='black')
    ax.yaxis.label.set_color('black')

    # Ensure lower y-limit is 0 (MSD should be nonnegative)
    ax.set_ylim(bottom=0)

    # Create a second y-axis for sqrt(MSD) using secondary_yaxis
    def msd_to_r(x):
        # Avoid negative inputs by ensuring x is nonnegative
        return np.sqrt(np.maximum(x, 0))

    def r_to_msd(x):
        return x**2

    ax_sqrt = ax.secondary_yaxis('right', functions=(msd_to_r, r_to_msd))
    ax_sqrt.set_ylabel(r'$\sqrt{\mathrm{MSD}}$ ($\mathrm{\AA}$)', fontsize=26, 
                       color='red', labelpad=15)
    ax_sqrt.tick_params(axis='y', colors='red')
    ax_sqrt.spines['right'].set_color('red')

    # Set secondary axis limits based on primary axis limits explicitly
    primary_ylim = ax.get_ylim()
    ax_sqrt.set_ylim(msd_to_r(primary_ylim[0]), msd_to_r(primary_ylim[1]))

    x_max = ax.get_xlim()[1]
    if annotate_times:
        for t_interest in annotate_times:
            idx = np.argmin(np.abs(lagtimes - t_interest))
            t_val = lagtimes[idx]
            msd_val = msd_values[idx]
            r_val   = np.sqrt(msd_val)

            # Mark the point on the MSD curve with enhanced aesthetics:
            ax.plot(
                t_val, msd_val, 
                marker='o', markersize=12, 
                markerfacecolor='red', markeredgecolor='white', 
                markeredgewidth=2, linestyle='None'
            )
            
            # Draw dashed lines to axes

            # Draw vertical line: from (t_val, msd_val) down to (t_val, 0)
            ax.plot([t_val, t_val], [msd_val, 0], '--', color='black')

            # Draw horizontal line: from (t_val, msd_val) to (x_max, msd_val)
            ax.plot([t_val, x_max], [msd_val, msd_val], '--', color='black')

            # Annotate the point with its time and sqrt(MSD)
            ax.annotate(
                rf"(τ={t_val:.0f} ps, $\langle \Delta r\rangle$={r_val:.2f} Å)",
                xy=(t_val, msd_val),
                xytext=(25, 10),
                textcoords="offset points",
                fontsize=22,
                color='black',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.8)
            )

    plt.tight_layout()
    plt.savefig("msd.png")

