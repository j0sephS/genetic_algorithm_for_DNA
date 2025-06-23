# =============================================================================
# Ce module implémente un algorithme génétique pour optimiser des paramètres
# de rotation (twist, wedge, direction) appliqués à une trajectoire 3D.
# =============================================================================

from dna.RotTable import RotTable
from dna.Traj3D import *
from math import *
import random
from copy import deepcopy
import time
import matplotlib.pyplot as plt

# =============================================================================
# Classe Individu
# Représente un individu au sein d'une population dans le cadre de l'algorithme
# génétique. Chaque individu contient :
#   - un objet RotTable (stockant les paramètres de rotation pour chaque
#     dinucléotide)
#   - une Traj3D (pour calculer un score via une fonction fitnesse (distance))
#   - un score (distance, dans l'exemple)
#   - un dictionnaire de plages de "bruit" autorisées, utilisées lors de la mutation à savoir ici les bornes décrites dans le format json
# =============================================================================


class Individu:

    def __init__(self):
        self.data = RotTable()  # Paramètres de rotation (Twist, Wedge, Direction)
        self.traj = Traj3D()
        self.score = None       # Score de l'individu (calculé via calcul_dist)
        self.bruit = {}         # Dictionnaire pour stocker les seuils min et max de bruit/ne pas sortir des bornes pour
        # les paramètres

        # Récupère la table par défaut puis initialise les plages de bruit
        table = self.data.getTable()
        for key in table:
            # Pour chaque dinucléotide, on stocke trois tuples (min, max):
            # C'est à dire les bornes max et min possibles pour faire varier
            # le twist, wedge et la direction dans les bornes convenues
            self.bruit[key] = [
                (-table[key][3], table[key][3]),
                (-table[key][4], table[key][4]),
                (-table[key][5], table[key][5])
            ]

    def __str__(self):
        """
        Output : String 

        Méthode d'overriding pour un affichage plus propre de chaque individu:
        - Sa dernière coordonée
        - Le score obtenu via la fonction fitness
        """

        s = "\nLast (x,y,z) : " + str(self.getLastPoint()) \
            + "\nScore : " + str(self.score)
        return s

    def copy(self):
        """
        Output : Individu

        Méthode d'overriding pour un clonage en profondeur d'un individu
        Le but est d'éviter d'associer des adresses mémoire directement aux attributs d'un individu,
        afin d'empêcher qu'un individu copié, en modifiant ses valeurs, ne répercute ces changements sur son clone
        """

        new = Individu()
        # Copie des paramètres de rotation
        for di in self.data.getTable():
            new.setDinucleotide(
                di,
                self.data.getTwist(di),
                self.data.getWedge(di),
                self.data.getDirection(di)
            )
        # Copie de la trajectoire, du score et du bruit
        new.traj = deepcopy(self.traj)
        new.score = self.score
        new.bruit = deepcopy(self.bruit)
        return new

    # -------------------------------------------------------------------------
    # Getters et setters
    # -------------------------------------------------------------------------
    def getTraj(self):
        return self.traj

    def getScore(self):
        return self.score

    def getData(self):
        return self.data

    def getBruit(self):
        return self.bruit

    # Récupère le dernier point de la trajectoire 3D (x, y, z)
    def getLastPoint(self):
        return self.traj.getTraj()[-1][:3]

    def setData(self, rot):
        self.data = rot

    def setScore(self, n):
        self.score = n

    # Met à jour les trois paramètres (Twist, Wedge, Direction) d'un dinucléotide
    def setDinucleotide(self, di, T, W, D):
        self.data.setTwist(di, T)
        self.data.setWedge(di, W)
        self.data.setDirection(di, D)

    # -------------------------------------------------------------------------
    # Méthode principale
    # -------------------------------------------------------------------------

    def add_bruit(self, dinucleotide: str):
        """
        Input : 
        - self : Individu
        - dinucleotide : String représentant le di_nucléatide dont les valeurs vont être changées

        Output : None

        Applique un bruit aléatoire à un dinucléotide donné, c'est-à-dire qu'on va
        piocher une valeur aléatoire dans les plages disponibles (cf attribut bruit), l'ajouter au
        paramètre (Twist / Wedge / Direction), puis mettre à jour ces plages (cf attribut bruit) pour
        éviter que le bruit ne s'accumule de manière illimitée (en dehors des bornes données)
        """

        # On récupère les bornes max et min des différents paramètres stockés dans l'attribut bruit
        twist_min, twist_max = self.bruit[dinucleotide][0][0], self.bruit[dinucleotide][0][1]
        wedge_min, wedge_max = self.bruit[dinucleotide][1][0], self.bruit[dinucleotide][1][1]
        direction_min, direction_max = self.bruit[dinucleotide][2][0], self.bruit[dinucleotide][
            2][1]

        # On tire trois valeurs uniformes indépendantes dans nos bornes possibles
        unif_t = random.uniform(twist_min, twist_max)
        unif_w = random.uniform(wedge_min, wedge_max)
        unif_d = random.uniform(direction_min, direction_max)

        # Met à jour le dinucléotide avec ces nouveaux paramètres
        self.setDinucleotide(
            dinucleotide,
            self.data.getTwist(dinucleotide) + unif_t,
            self.data.getWedge(dinucleotide) + unif_w,
            self.data.getDirection(dinucleotide) + unif_d
        )

        # Adapte les plages de bruit (décalage en fonction de la valeur prise)
        twist_min, twist_max = twist_min - unif_t, twist_max - unif_t
        wedge_min, wedge_max = wedge_min - unif_w, wedge_max - unif_w
        direction_min, direction_max = direction_min - unif_d, direction_max - unif_d

        # Met à jour le dictionnaire de bruit
        self.bruit[dinucleotide] = [
            (twist_min, twist_max),
            (wedge_min, wedge_max),
            (direction_min, direction_max)
        ]


# =============================================================================
# Classe Genetique
# Gère la population d'individus dans l'algorithme génétique. Elle inclut :
#   - la création initiale de la population
#   - les méthodes de sélection (élitisme, roulette, tournoi)
#   - les méthodes de croisement (1 point ou n points)
#   - la mutation
#   - la mise à jour des scores
# =============================================================================
class Genetique:

    def __init__(self, len_pop):
        # Crée une liste d'individus de taille len_pop
        self.population = [Individu() for _ in range(len_pop)]

        # Pour chaque individu, on ajoute un bruit initial sur tous les dinucléotides (1ere genération aléatoire)
        for ind in self.population:
            for dinucleotide in ind.data.getTable():
                ind.add_bruit(dinucleotide)

        self.len_pop = len_pop
        self.best_individu = None  # On stock le meilleur individu (c'est a dire distance minimale)

    def __str__(self):
        """
        Output : String 

        Méthode d'overriding pour un affichage plus propre de la population :
        Elle s'appuie sur la méthode __str__ recréée dans la classe Individu
        """

        s = "-----------------------------------\nFst (x,y,z) : [0. 0. 0.]\n"
        for i in range(self.len_pop):
            s += "ind" + str(i) + " : " + self.population[i].__str__()
            if (i != self.len_pop - 1):
                s += "\n\n\n"
            else:
                s += "\n-----------------------------------"
        return s

    # -------------------------------------------------------------------------
    # Getter
    # -------------------------------------------------------------------------

    def getBest_individu(self):
        return self.best_individu

    # =============================================================================
    # Méthodes principales
    # =============================================================================

    def selection(self, method='elitisme', rate=0.5):
        """
        Input : 
        - Genetique
        - method : str, prédéfini à 'elitisme'
        - rate : float, permettant de définir la proportion de "bons individus" à garder,
                prédéfini à 50%

        Sélectionne les individus selon la méthode choisie :
            - 'elitisme' : on garde les meilleurs individus
            - 'roulette' : on tire au sort selon une probabilité proportionnelle à
                        l'inverse du score
            - 'tournoi'  : on choisit deux individus au hasard et on prend le meilleur
        """

        if method == 'elitisme':
            self.selection_elitisme(rate)
        elif method == 'roulette':
            self.selection_roulette(rate)
        elif method == 'tournoi':
            self.selection_tournoi(rate)
        else:
            raise ValueError(f"Unknown selection method: {method}")

    # -------------------------------------------------------------------------
    # Types de selection
    # -------------------------------------------------------------------------

    def selection_elitisme(self, rate=0.5):
        """
        Input : 
        - Genetique
        - rate : float, permettant de définir la proportion de "bons individus" à garder,
            prédéfini à 50%

        Output : None -> Met a jour les attributs population et len_pop

        Selection visant à garder rate*100% des meilleurs individus (ceux avec un score minimal)

        """

        # On coupe la taille de la population à rate * len_pop
        self.len_pop = int(rate * self.len_pop)

        # On trie la population par ordre croissant de distance (calcul_dist),
        # puis on ne garde que les meilleurs (c'est à dire les rate*100% premiers)
        self.population = sorted(
            self.population,
            key=lambda x: calcul_dist(x),
            reverse=False
        )[: self.len_pop]

        # Si le meilleur individu global n'est pas défini ou si le meilleur
        # actuel est mieux noté, on le met à jour
        if self.best_individu == None or self.population[0].score < self.best_individu.score:
            # appel de copy de individu (cf ligne 59)
            self.best_individu = self.population[0].copy()

    def selection_roulette(self, rate=0.5):
        """
        Input : 
        - Genetique
        - rate : float, permettant de définir la proportion de "bons individus" à garder,
            prédéfini à 50%

        Output : None -> Met a jour les attributs population et len_pop

        Selection visant à garder rate*100% des individus tirés au sort selon une probabilité proportionnelle 
        à leur score (plus petit=plus probable)
        """

        proba = [1 / x.score for x in self.population]
        # Normalisation des probabilités. On inverse les probabilités proportionnelles
        # au score afin qu'un score plus petit soit plus probable qu'un score plus gros
        proba = [p / sum(proba) for p in proba]

        selected = []
        pop = int(rate * self.len_pop)
        for _ in range(pop):
            # On choisit un individu selon les probabilités calculées
            ind: Individu = np.random.choice(self.population, p=proba)
            selected.append(ind)
            # Mise à jour du meilleur individu, si nécessaire
            if self.best_individu == None or ind.score < self.best_individu.score:
                self.best_individu = ind.copy()  # appel de copy de individu (cf ligne 59)

        # On met à jour nos attributs
        self.population = selected
        self.len_pop = pop

    def selection_tournoi(self, rate=0.5):
        """
        Input : 
        - Genetique
        - rate : float, permettant de définir la proportion de "bons individus" à garder,
            prédéfini à 50%

        Output : None -> Met a jour les attributs population et len_pop

        Selection où on choisit deux individus au hasard et on prend le meilleur.
        Toujours garder rate*100% des individus
        """

        selected = []
        pop = int(rate * self.len_pop)
        for _ in range(pop):
            # On tire 2 individus au hasard et on prend le meilleur (score minimal)
            tournament = random.sample(self.population, 2)
            winner = min(tournament, key=lambda x: x.score)
            selected.append(winner)
            # Mise à jour du meilleur si nécessaire
            if self.best_individu is None or winner.score < self.best_individu.score:
                self.best_individu = winner.copy()  # appel de copy de individu (cf ligne 59)

        # On met à jour nos attributs
        self.population = selected
        self.len_pop = pop

    # -------------------------------------------------------------------------
    # Méthode pour la fonction fitness
    # -------------------------------------------------------------------------
    def refresh_score(self, seq):
        """
        Input : 
        - Genetique
        - seq : str, séquence d'adn traitées

        Output : None -> Recalcul le score et met a jour son attribut (score) de
        chaque individu dans la population

        """

        for individu in self.population:
            individu.traj.compute(seq, individu.data)  # Calcule la trajectoire
            score = calcul_dist(individu)              # Calcule la distance finale
            individu.setScore(score)                   # Met a jour l'attribut

    # -------------------------------------------------------------------------
    # Méthode pour le croisement
    # -------------------------------------------------------------------------

    def croisement_n_point(self, n=2):
        """
        Input : 
        - Genetique
        - n : int, nombre de croisement 

        Output : None -> met a jour les attributs

        Méthode du croisement en n point (de 1 jusqu'à len_pop)
        """

        # On reajuste les valeurs de n pour ne pas sortir des bornes possibles
        if n > self.len_pop:
            n = self.len_pop
        elif n <= 0:
            n = 1

        # Chaque pair de parent va créer 2 fils (il est donc néccéssaire de parcourir
        # len_pop//2 pour garder une population constante)

        for _ in range(0, self.len_pop, 2):

            # Selection de 2 parents au hasards pour le croisement
            parent1, parent2 = random.sample(self.population, 2)
            child1, child2 = Individu(), Individu()
            table = parent1.data.getTable().keys()

            # On récupère n dinucléotides (points de croisement)
            crossover_points = random.sample(list(table), n)

            # booléen qui permettra d'échanger les parties des chromosomes
            between_cross_point = False

            # On parcourt l'ensemble des dinucléotides
            for dinucleotide in table:

                # On récupére les paramètres des 2 parents
                p1_T, p1_W, p1_D = (
                    parent1.data.getTwist(dinucleotide),
                    parent1.data.getWedge(dinucleotide),
                    parent1.data.getDirection(dinucleotide)
                )
                p2_T, p2_W, p2_D = (
                    parent2.data.getTwist(dinucleotide),
                    parent2.data.getWedge(dinucleotide),
                    parent2.data.getDirection(dinucleotide)
                )

                # Chaque fois qu'on rencontre un point de croisement,
                # on inverse la logique de copie
                if dinucleotide in crossover_points:
                    between_cross_point = not between_cross_point

                # Par défaut, child1 <- parent1 et child2 <- parent2
                # mais si on est entre deux points, on échange
                c_fst, c_snd = child1, child2
                if between_cross_point:
                    c_fst, c_snd = child2, child1

                # On ajoute les chromosomes
                c_fst.setDinucleotide(dinucleotide, p1_T, p1_W, p1_D)
                c_fst.bruit[dinucleotide] = parent1.bruit[dinucleotide].copy()
                c_snd.setDinucleotide(dinucleotide, p2_T, p2_W, p2_D)
                c_snd.bruit[dinucleotide] = parent2.bruit[dinucleotide].copy()

            # ajout des 2 fils croisés
            self.population.extend([child1, child2])

        # Il se peut qu'en fonction de la parité de la taille de départ de notre population
        # prendre la motié de cette dernière rajoute un individu en trop à notre population
        # lors du croisement.
        # Exemple : 50//2 = 25
        # 25//2 (cf ligne 381) = 12
        # Le croisement donne donc 51 éléments
        # cette condition permet d'enlever donc le fils en trop
        if (len(self.population) % 2 != 0):
            self.population = self.population[:len(self.population) - 1]

        self.len_pop = len(self.population)

    # -------------------------------------------------------------------------
    # Méthode pour la mutation
    # -------------------------------------------------------------------------

    def mutation(self, seuil, rate=0.5):
        """
        Input : 
        - seuil : borne max d'un intervalle où une proba pourra tombée dedans pour vérifier une condition
        - rate : % de la population qui va être mutée

        Output : None -> Met a jour attibuts

        Avec une probabilité seuil, on prend un certain nombre (rate*taille_pop)
        d'individus au hasard et on leur applique une perturbation (add_bruit)
        sur un dinucléotide aléatoire.
        """

        len_m = int(rate * self.len_pop)
        P_m = random.uniform(0, 1)
        if P_m < seuil:
            for _ in range(len_m):
                individu = random.choice(self.population)  # choix individu au hasard
                mutation_point = random.choice(
                    list(individu.data.getTable().keys()))  # dinucléotide au hasard
                individu.add_bruit(mutation_point)  # bruit sur ce dinucléotide


# =============================================================================
# Fin des classes / Début des fonctions pour l'algorithme
# =============================================================================

# -------------------------------------------------------------------------
# Fonction fitness
# -------------------------------------------------------------------------
def calcul_dist(individu: Individu):
    """
        Input : 
        - seuil : borne max d'un intervalle où une proba pourra tombée dedans pour vérifier une condition
        - rate : % de la population qui va être mutée

        Output : None -> Met a jour attibuts

        Avec une probabilité seuil, on prend un certain nombre (rate*taille_pop)
        d'individus au hasard et on leur applique une perturbation (add_bruit)
        sur un dinucléotide aléatoire.
    """
    x, y, z = individu.getLastPoint()[0], individu.getLastPoint()[1], individu.getLastPoint()[2]
    # Distance de l'origine : sqrt(x^2 + y^2 + z^2)
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))


#  -----------------------------------------------------------------------------
# Algorithme
# -----------------------------------------------------------------------------
def algo_genetique(seq, taille,istest=False, n=2, algorithme_selection='elitisme', rate=0.5) -> Individu:
    """
    Input : 
        - seq : chaine d'adn à calculer
        - taille : int pour le nombre d'individu -> attribut len_pop de Genetique
        - n : int pour le croisement_n_point (init=2)
        - rate : float % de la population qui va être mutée

    Output :
        - Indidivu qui minimise la distance pour notre problème

    # Fonction principale mettant en œuvre l'algorithme génétique :

    #   - Crée la population de taille 'taille'

    #   - Boucle jusqu'à un certain critère d'arrêt (acc != 40) Tant que la génération ne donne pas 40 fois le meme meilleur score
    #   - À chaque itération :
    #       1. Sélection
    #       2. Croisement (n points)
    #       3. Mutation
    #       4. Mise à jour des scores
    #       5. On conserve le meilleur individu 
    """

    pop = Genetique(taille)
    pop.refresh_score(seq)
    best = pop.getBest_individu()
    acc = 0
    seuil = 0.5

    print("---- Lancement de l'algorithme génétique ----")
    while (acc != 40):  # Tant que la génération ne donne pas 40 fois le meme meilleur score
        # Sélection
        pop.selection(algorithme_selection, rate)
        # Croisement à n points
        pop.croisement_n_point(n)
        # Mutation
        pop.mutation(seuil)
        # Mise à jour des scores
        pop.refresh_score(seq)

        tmp = pop.getBest_individu()
        print(str(acc) + " :")
        # Si on trouve un meilleur individu, on reset acc
        if (best == None or tmp.score < best.score):
            acc = 0
            best = tmp
        else:
            acc += 1
        # Tous les 5 itérations sans nouveau changement, on augmente les probas de mutation
        # pour augmenter les chances de trouver une nouvelle solution
        if (acc % 5 == 0):
            seuil += 0.3
        print(pop.getBest_individu())

    # On teste le meilleur individu final pour vérifier qu'il respecte
    # les contraintes de la RotTable d'origine (pas de dépassement)
    table = pop.getBest_individu().getData().getTable()
    if isInBounds(table):
        print("\033[92mVrai\033[0m")
    else:
        print("\033[91mFaux\033[0m")
    print(pop.getBest_individu().getData().getTable())
    if not istest: pop.getBest_individu().traj.draw()
    return pop.getBest_individu()


#  -----------------------------------------------------------------------------
# Fonctions secondaires
# -----------------------------------------------------------------------------

def isInBounds(table):
    """
    Input : 
        - table : dict des valeurs de notre solution

    Output : Boolean

    Permet de savoir si la solution respecte
    les contraintes de la RotTable d'origine (pas de dépassement des bornes min et max)
    """

    tdep = RotTable().getTable()

    # on parcourt les nucléotides notre solution
    for di in tdep:
        # On regarde si Twist est bien dans les bornes min et max grace aux valeurs originales avec tdep
        if table[di][0] > tdep[di][0] + tdep[di][3] or table[di][0] < tdep[di][0] - tdep[di][3]:
            print(di)
            return False
        # On regarde si Wedge est bien dans les bornes min et max grace aux valeurs originales avec tdep
        if table[di][1] > tdep[di][1] + tdep[di][4] or table[di][1] < tdep[di][1] - tdep[di][4]:
            print(di)
            return False
        # On regarde si Direction est bien dans les bornes min et max grace aux valeurs originales avec tdep
        if table[di][2] > tdep[di][2] + tdep[di][5] or table[di][2] < tdep[di][2] - tdep[di][5]:
            print(di)
            return False
    # si tous les nucléotides sont validés
    return True


def stats(seq):
    """
    Permet de comparer l'impact de la taille de la population (populations)
    et de la méthode de sélection (selection_methods) sur le score final.
    Affiche un graphe score vs population.
    """

    populations = [20, 50, 100, 200]
    selection_methods = ['elitisme', 'roulette', 'tournoi']
    plt.figure(figsize=(10, 6))
    for method in selection_methods:
        scores = []
        for pop in populations:
            scores.append(algo_genetique(seq, pop, method).getScore())
        plt.plot(populations, scores, label=f'Score & Population ({method})')

    plt.xlabel('Population')
    plt.ylabel('Score')
    plt.legend()
    plt.show()
    plt.savefig("score_vs_population.png")
