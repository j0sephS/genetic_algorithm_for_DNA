import math
import random
import numpy as np
from dna.Traj3D import Traj3D
from dna.RotTable import RotTable
import os
import json
import copy


class Recuit:
    """Recuit simulé

    Args:
        seqs (list): Liste des séquences à comparer
        initial_state (RotTable): Modèle de conformation initial
        k_max (int): Nombre maximal d'itérations
        e_max (int): Energie seuil pour arrêter l'algorithme
    """

    def __init__(self, seqs, initial_state, k_max, e_max):
        self.seqs = seqs
        self.initial_state = initial_state
        self.state = initial_state
        self.e = self.energy(initial_state)
        self.k = 0
        self.k_max = k_max
        self.e_max = e_max
        self.temp = 150000
        """self.initial_delta_temp = 1000
        self.delta_temp = 1000
        self.stuck = 0"""

    def generateNewState(self):
        """Génère un nouvel état voisin en modifiant légèrement l'état actuel
        return: RotTable -- Nouvel état
        """

        new_state = copy.deepcopy(self.state)

        # Modifier légèrement l'état
        for key in new_state.rot_table:
            ranges = new_state.getRanges(key)
            delta_Twist, delta_Wedge, delta_Direction = random.uniform(-min(ranges[0])/3, min(
                ranges[0])/3), random.uniform(-min(ranges[1])/3, min(ranges[1])/3), 0
            new_state.updateRangesAndValues(key, [delta_Twist, delta_Wedge, delta_Direction])

        return new_state

    def energy(self, state):
        """Calcule l'énergie d'un état

        Args:
            state (RotTable): Etat à évaluer

        Returns:
            float -- Energie de l'état
        """

        diff = 0

        # Pour chaque séquence, on calcule l'énergie de la trajectoire (distance entre le départ et l'arrivée)^2 (sqrt prend beaucoup de temps et est inutile dans une fonction d'évaluation)
        traj = Traj3D()

        # dists = [] #Pour le terme de normalisation
        for seq in self.seqs:
            traj.compute(seq, state)
            # dists.append(traj.getDistance())
            diff += traj.energy()

        # On ajoute un terme qui vise à trouver des valeurs valables pour toute séquence -> on normalise par la longueur de la séquence:
        # diff += abs(dists[0]/len(self.seqs[0]) - dists[1] / len(self.seqs[1])) * 10000000

        return diff

    def probability(self, energy_diff, temperature):
        """Calcule la probabilité d'accepter un nouvel état, même si son énergie est plus élevée

        Args:
            energy_diff (float): Différence d'énergie entre le nouvel état et l'actuel
            temperature (float): Température actuelle

        Returns:
            float -- Probabilité d'accepter le nouvel état
        """

        # On veut une fonction qui :
        # * se rapproche de 0 quand energy_diff augmente (on rechigne à accepter un état vraiment 'moins bon' (énergie beaucoup plus élevée))
        # * se rapproche de 1 quand energy_diff diminue (on accepte plus facilement des états 'moins bons' s'ils sont proches)
        # * la température joue un rôle important: plus elle est élevée, plus on accepte facilement des états 'moins bons'

        return math.exp(-energy_diff/temperature)

    def calculateTemp(self, k):
        """Calcule la température actuelle

        Args:
            k (int): Itération actuelle

        Returns:
            float -- Température actuelle
        """

        # La température doit être très élevée au début puis décroître (on veut converger vers un minimum global)
        self.temp = 0.992*self.temp
        return self.temp

    def iterate(self):
        """Itère une fois dans l'algorithme de recuit simulé"""
        new_state = self.generateNewState()
        new_energy = self.energy(new_state)
        if new_energy < self.e or random.random() < self.probability(
                new_energy - self.e, self.calculateTemp(self.k)):
            self.state = new_state
            self.e = new_energy
        self.k += 1

    def run(self):
        """Lance l'algorithme de recuit simulé"""
        while self.k < self.k_max and self.e > self.e_max:
            self.iterate()
            # Affichage toutes les 10 itérations
            if not self.k % 10:
                print(f"iteration:{self.k}       energy:{self.e:.2f}       temp:{self.temp:.2f}")
        return self.state

    def write(self, filename="results/recuit_result"):
        """Enregistre l'état final dans un fichier JSON"""
        i = 1
        while os.path.exists(f"{filename}{i}.json"):
            i += 1
        with open(f"{filename}{i}.json", "w") as file:
            json.dump(self.state.rot_table, file, indent=4)
        print("Result saved in", f"{filename}{i}.json")


def recuit_main(seqs, JSON_filename, max_iters=100):
    """Fonction appelée par __main__.py (pour éviter le bouclage des imports)"""
    recuit = Recuit(seqs, RotTable(JSON_filename), max_iters, 10)
    print("---- Lancement de l'algorithme du recuit simulé ----")
    recuit.run()
    traj = Traj3D()
    dist = []
    for seq in seqs:
        traj.compute(seq, recuit.state)
        dist.append(traj.getDistance())
        print("Distance:", traj.getDistance())
        # traj.draw()
    recuit.write()
