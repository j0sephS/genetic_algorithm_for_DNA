from dna.RotTable import RotTable
from dna.Traj3D import Traj3D


def traditionnal_main(seq, filename, JSON_filename):
    """Fonction appelée par __main__.py (pour éviter le bouclage des imports)"""
    rot_table = RotTable(JSON_filename)
    traj = Traj3D()
    traj.compute(seq, rot_table)

    # print(traj.getTraj())
    print("Distance:", traj.getDistance())
    traj.draw()
    traj.write(filename+".png")
