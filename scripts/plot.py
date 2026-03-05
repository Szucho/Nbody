import numpy as np
import matplotlib.pyplot as plt


def main():

    #loading the output file of the cpp code
    out = np.genfromtxt("trajectory.dat", dtype=float, delimiter=" ", skip_header=1)
    masses = np.genfromtxt("trajectory.dat", dtype=float, delimiter=" ", max_rows=1)

    #check if number of bodies & dimension (output file) is correct
    n_bodies = int(len(out[0])//6) #3D -> number of bodies is number of times (x,y,z)
    print(f"Number of bodies based on dimension and output vectors: {n_bodies}")
    print(f"Actual number of bodies: {len(masses)}")
    #we expect binaries, therefore we can check the reduced particle
    reduced = out[:,0:3]-out[:,3:6]

    #centre of mass of the binary
    com = (masses[0]*out[0,0:3]+masses[1]*out[0,3:6])/np.sum(masses)

    #energy of reduced particle
    redvel = out[:,6:9]-out[:,9:]
    h = 0.5*np.linalg.norm(redvel, axis=1)**2 - np.sum(masses)/np.linalg.norm(reduced, axis=1)

    #plot of the reduced particle
    plt.figure(figsize=(8,8))
    plt.plot(reduced[:,0], reduced[:,1], "b-", label="reduced body")
    plt.grid(ls="--", alpha=0.5)
    plt.xlabel(r"$x$[AU]")
    plt.ylabel(r"$y$[AU]")
    plt.legend(loc="best")
    plt.axis("equal")
    plt.show()

    #plot showing the 2 bodies 
    plt.figure(figsize=(8,8))
    plt.plot(out[:,0], out[:,1], "b-", label="body 1")
    plt.plot(out[:,3], out[:,4], "r--", label="body 2")
    plt.scatter(com[0], com[1], marker="x", color="k", label="CoM")
    plt.grid(ls="--", alpha=0.5)
    plt.xlabel(r"$x$[AU]")
    plt.ylabel(r"$y$[AU]")
    plt.axis("equal")
    plt.legend(loc="best")
    plt.show()
    
    #plot showing energy of the reduced particle
    plt.figure(figsize=(8,8))
    plt.plot(np.arange(0,50.01,0.01), abs((h-h[0])/h[0]), "b-", label="integrated energy")
    plt.axhline(h[0], ls="--", color="black", alpha=0.5, label="initial energy")
    plt.grid(ls="--", alpha=0.5)
    plt.legend(loc="best")
    plt.xlabel(r"$t$")
    plt.yscale("log")
    plt.ylabel(r"$|\frac{h-h_0}{h_0}|$")
    plt.show()




if __name__=="__main__":
    main()
