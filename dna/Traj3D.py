import math

# For drawing
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from dna.RotTable import RotTable


class Traj3D:
    """Represents a 3D trajectory"""

    # Vertical translation (elevation) between two di-nucleotides
    __MATRIX_T = np.array(
        [[1, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 1, -3.38/2],
         [0, 0, 0, 1]]
    )

    def __init__(self):
        self.__Traj3D = {}

    def getTraj(self) -> dict:
        return self.__Traj3D

    def compute(self, dna_seq: str, rot_table: RotTable):

        # Matrice cumulant l'ensemble des transformations géométriques engendrées par la séquence d'ADN
        total_matrix = np.eye(4)  # Identity matrix

        # On enregistre la position du premier nucléotide
        self.__Traj3D = [np.array([0.0, 0.0, 0.0, 1.0])]

        matrices_Rz = {}
        matrices_Q = {}
        # On parcourt la sequence, nucléotide par nucléotide
        for i in range(1, len(dna_seq)):
            # On lit le dinucleotide courant
            dinucleotide = dna_seq[i-1]+dna_seq[i]
            # On remplit au fur et à mesure les matrices de rotation
            if dinucleotide not in matrices_Rz:
                matrices_Rz[dinucleotide], matrices_Q[dinucleotide] = \
                    self.__compute_matrices(rot_table, dinucleotide)

            # On calcule les transformations géométriques
            # selon le dinucleotide courant,
            # et on les ajoute à la matrice totale
            total_matrix = \
                total_matrix \
                @ self.__MATRIX_T \
                @ matrices_Rz[dinucleotide] \
                @ matrices_Q[dinucleotide] \
                @ matrices_Rz[dinucleotide] \
                @ self.__MATRIX_T

            # On calcule la position du nucléotide courant
            # en appliquant toutes les transformations géométriques
            # à la position du premier nucléotide
            self.__Traj3D.append(total_matrix @ self.__Traj3D[0])

    def __compute_matrices(self, rot_table: RotTable, dinucleotide: str):

        Omega = math.radians(rot_table.getTwist(dinucleotide))
        # Create rotation matrix of theta on Z axis
        cO = math.cos(Omega/2)
        sO = math.sin(Omega/2)
        matrices_Rz = \
            np.array([[cO, sO, 0, 0],
                      [-sO, cO, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

        sigma = rot_table.getWedge(dinucleotide)
        delta = rot_table.getDirection(dinucleotide)
        alpha = math.radians(sigma)
        beta = math.radians(delta - 90)
        # Rotate of -beta on Z axis
        # Rotate of -alpha on X axis
        # Rotate of beta on Z axis
        cb = math.cos(beta)
        sb = math.sin(beta)
        ca = math.cos(alpha)
        sa = math.sin(alpha)
        matrices_Q = \
            np.array([[cb, -sb, 0, 0],
                      [sb, cb, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]]) \
            @ np.array([[1, 0, 0, 0],
                        [0, ca, -sa, 0],
                        [0, sa, ca, 0],
                        [0, 0, 0, 1]]) \
            @ np.array([[cb, sb, 0, 0],
                        [-sb, cb, 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])

        return matrices_Rz, matrices_Q

    def draw(self):
        self.fig = plt.figure()
        self.ax = plt.axes(projection='3d')
        xyz = np.array(self.__Traj3D)
        x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
        self.ax.plot(x[1:-1], y[1:-1], z[1:-1])
        self.ax.scatter(x[0], y[0], z[0], c='red')
        self.ax.scatter(x[-1], y[-1], z[-1], c='green')
        plt.show()

    def getDistance(self):
        return math.sqrt(self.energy())

    def write(self, filename: str):
        self.fig.savefig(filename)

    def energy(self):
        xyz = np.array(self.__Traj3D)
        x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
        return (x[0]-x[-1])**2 + (y[0]-y[-1])**2 + (z[0]-z[-1])**2
