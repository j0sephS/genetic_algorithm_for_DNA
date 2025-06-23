from dna.Recuit import Recuit
from dna.RotTable import RotTable
from dna.Traj3D import Traj3D


def test_initial_state():
    seqs = ["AGCT", "CGTA"]
    initial_state = RotTable()
    recuit = Recuit(seqs, initial_state, 100, 10)
    assert recuit.state == initial_state
    assert recuit.e == recuit.energy(initial_state)
    assert recuit.k == 0


def test_generate_new_state():
    seqs = ["AGCT", "CGTA"]
    initial_state = RotTable()
    recuit = Recuit(seqs, initial_state, 100, 10)
    new_state = recuit.generateNewState()
    assert new_state != recuit.state


def test_energy_calculation():
    seqs = ["AGCT", "CGTA"]
    initial_state = RotTable()
    # Pour faire une ligne droite -> plus long qu'avec la table
    initial_state2 = RotTable('test_table.json')
    recuit = Recuit(seqs, initial_state, 100, 10)
    recuit2 = Recuit(seqs, initial_state2, 100, 10)
    assert isinstance(recuit.e, float)
    assert isinstance(recuit2.e, float)
    assert recuit2.e > recuit.e


def test_probability():
    seqs = ["AGCT", "CGTA"]
    initial_state = RotTable()
    recuit = Recuit(seqs, initial_state, 100, 10)
    T1, T2 = 100, 1000
    dE1, dE2 = 100, 1000
    P11 = recuit.probability(dE1, T1)
    P21 = recuit.probability(dE2, T1)
    P12 = recuit.probability(dE1, T2)
    P22 = recuit.probability(dE2, T2)
    assert 0 <= P11 <= 1 and 0 <= P21 <= 1 and 0 <= P12 <= 1 and 0 <= P22 <= 1
    assert P21 < P11 and P21 < P12  # Plus grande différence d'énergie -> plus petite probabilité
    assert P11 < P12 and P21 < P22  # Plus grande température -> plus grande probabilité


def test_calculate_temp():
    seqs = ["AGCT", "CGTA"]
    initial_state = RotTable()
    recuit = Recuit(seqs, initial_state, 100, 10)
    old_temp = recuit.temp
    temp = recuit.calculateTemp(1)
    assert recuit.temp < old_temp


def test_iterate():
    seqs = ["AGCT", "CGTA"]
    initial_state = RotTable()
    recuit = Recuit(seqs, initial_state, 100, 10)
    recuit.iterate()
    assert recuit.k == 1


def test_run():
    seqs = ["AGCT", "CGTA"]
    initial_state = RotTable()
    recuit = Recuit(seqs, initial_state, 10, 10)
    recuit.run()
    assert (recuit.k == recuit.k_max and recuit.e > recuit.e_max) or (
        recuit.k < recuit.k_max and recuit.e <= recuit.e_max)
