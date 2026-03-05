import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation, FFMpegWriter

def main():
    # --- load data ---
    out = np.genfromtxt("trajectory.dat", dtype=float, delimiter=" ", skip_header=1)
    masses = np.genfromtxt("trajectory.dat", dtype=float, delimiter=" ", max_rows=1)
    n_steps = out.shape[0]
    dt = 0.01
    time = np.arange(n_steps) * dt
    
    # positions
    r1 = out[:, 0:3]
    r2 = out[:, 3:6]
    
    # velocities
    v1 = out[:, 6:9]
    v2 = out[:, 9:12]
    
    # reduced particle
    reduced = r1 - r2
    redvel = v1 - v2
    
    # center of mass (constant for binary ICs)
    com = (masses[0] * r1[0] + masses[1] * r2[0]) / np.sum(masses)
    
    # energy of reduced particle
    h = 0.5 * np.linalg.norm(redvel, axis=1)**2 \
        - np.sum(masses) / np.linalg.norm(reduced, axis=1)
    energy_err = np.abs((h - h[0]) / h[0])
    
    # --- figure layout ---
    fig = plt.figure(figsize=(14, 7))
    gs = GridSpec(2, 2, width_ratios=[2.2, 1], height_ratios=[1, 1])
    
    ax_bin = fig.add_subplot(gs[:, 0])      # big left panel
    ax_red = fig.add_subplot(gs[0, 1])      # top-right
    ax_eng = fig.add_subplot(gs[1, 1])      # bottom-right
    
    # --- binary plot ---
    ax_bin.set_title("Binary system")
    ax_bin.set_aspect("equal")
    ax_bin.grid(ls="--", alpha=0.5)
    
    # Fix: Use min/max of both r1 and r2 to set limits
    all_x = np.concatenate([r1[:, 0], r2[:, 0]])
    all_y = np.concatenate([r1[:, 1], r2[:, 1]])
    ax_bin.set_xlim(1.2 * all_x.min(), 1.2 * all_x.max())
    ax_bin.set_ylim(1.2 * all_y.min(), 1.2 * all_y.max())
    ax_bin.set_xlabel(r"$x$")
    ax_bin.set_ylabel(r"$y$")
    
    line1, = ax_bin.plot([], [], "b-", lw=1, label="Body 1")
    line2, = ax_bin.plot([], [], "r--", lw=1, label="Body 2")
    p1, = ax_bin.plot([], [], "bo", ms=8)
    p2, = ax_bin.plot([], [], "ro", ms=8)
    com_pt, = ax_bin.plot(com[0], com[1], "kx", ms=10, mew=2, label="CoM")
    ax_bin.legend(loc="best")
    
    # --- reduced particle plot ---
    ax_red.set_title("Reduced particle")
    ax_red.set_aspect("equal")
    ax_red.grid(ls="--", alpha=0.5)
    ax_red.set_xlabel(r"$\Delta x$")
    ax_red.set_ylabel(r"$\Delta y$")
    ax_red.plot(reduced[:, 0], reduced[:, 1], "k-", lw=1, alpha=0.3)
    red_pt, = ax_red.plot([], [], "bo", ms=8)
    
    # --- energy plot ---
    ax_eng.set_title("Energy conservation")
    ax_eng.set_yscale("log")
    ax_eng.grid(ls="--", alpha=0.5)
    ax_eng.plot(time, energy_err, "b-", lw=1)
    t_marker = ax_eng.axvline(0, color="r", ls="--", lw=2)
    ax_eng.set_xlabel(r"$t$")
    ax_eng.set_ylabel(r"$|\Delta h / h_0|$")
    
    plt.tight_layout()
    
    # --- animation update ---
    def update(i):
        # binary trajectories
        line1.set_data(r1[:i, 0], r1[:i, 1])
        line2.set_data(r2[:i, 0], r2[:i, 1])
        p1.set_data([r1[i, 0]], [r1[i, 1]])
        p2.set_data([r2[i, 0]], [r2[i, 1]])
        
        # reduced particle
        red_pt.set_data([reduced[i, 0]], [reduced[i, 1]])
        
        # energy marker
        t_marker.set_xdata([time[i], time[i]])
        
        return line1, line2, p1, p2, red_pt, t_marker
    
    ani = FuncAnimation(
        fig,
        update,
        frames=n_steps,
        interval=20,
        blit=True  # Changed to True for better performance
    )
    
    print("Saving animation... this may take a while")
    writer = FFMpegWriter(fps=30, bitrate=1800)
    ani.save("binary.mp4", writer=writer, dpi=300)
    print("Animation saved as binary.mp4")
    
    plt.show()

if __name__ == "__main__":
    main()
