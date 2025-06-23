import pytest
from dna.Genetic import *
from dna.RotTable import *
from dna.Traj3D import *
from unittest.mock import patch

# ====== TESTS POUR LA CLASSE Individu ======

def test_individu_init():
    t = RotTable().getTable()
    individu = Individu()
    assert type(individu.data) is type(RotTable())
    assert type(individu.traj) is type(Traj3D())
    assert individu.score is None
    
    for key in t:
        assert individu.bruit[key]==[(-t[key][3],t[key][3]),(-t[key][4],t[key][4]),(-t[key][5],t[key][5])]

def test_individu_copy():
    individu1 = Individu()
    individu2 = individu1.copy()
    individu1.getTraj().compute("AA",individu1.getData())
    individu2.getTraj().compute("AA",individu2.getData())
    assert isinstance(individu2, Individu)
    assert individu1 is not individu2
    assert individu1.traj is not individu2.traj
    assert individu1.traj.getDistance() == individu2.traj.getDistance()
    assert individu1.data is not individu2.data
    assert individu1.data.getTable() == individu2.data.getTable()
    assert individu1.getBruit() == individu2.getBruit()
    assert individu1.bruit is not individu2.bruit
    individu1.add_bruit('AA')
    assert individu1.data.getTable() != individu2.data.getTable()
    assert individu1.getBruit() != individu2.getBruit()
    assert individu1.score==individu2.score
    

def test_individu_getters():
    individu = Individu()
    assert individu.getTraj() == individu.traj
    assert individu.getScore() == individu.score
    assert individu.getData() == individu.data
    assert individu.getBruit() == individu.bruit

def test_individu_setters():
    individu = Individu()
    individu.setScore(100)
    assert individu.getScore() == 100
    r = RotTable("dna/table.json")
    r.getRanges('AA')
    individu.setData(r)
    assert individu.getData() == r
    individu.setDinucleotide('AA',1,1,1)
    assert individu.getData().getTable()['AA'] == [1,1,1] + r.getRanges('AA')

def test_individu_add_bruit():
    individu = Individu()
    individu.add_bruit()
    assert isinstance(individu.bruit, dict)

r = RotTable()
m = r.getRanges('AA')[0][0]
@patch("random.uniform", return_value=m) # forcer la valeur de random unif
def test_individu_add_bruit(n):
    individu = Individu()
    print(n)
    original_data = deepcopy(individu.data)
    individu.add_bruit('AA')
    
    modified_data = individu.data
    assert original_data.getTable() != modified_data.getTable()
    assert modified_data.getTwist('AA') == original_data.getTwist('AA') + m
    assert modified_data.getWedge('AA') == original_data.getWedge('AA') + m
    assert modified_data.getDirection('AA') == original_data.getDirection('AA') + m

    assert individu.getBruit()['AA'] == [(-r.getRanges('AA')[0][0]-m,r.getRanges('AA')[0][1]-m),
                                         (-r.getRanges('AA')[1][0]-m,r.getRanges('AA')[1][1]-m),
                                         (-r.getRanges('AA')[2][0]-m,r.getRanges('AA')[2][1]-m)
                                        ]



# ====== TESTS POUR LA CLASSE Genetique ======

def test_genetique_init():
    genetique = Genetique(10)
    assert isinstance(genetique, Genetique)
    for ind in genetique.population: assert isinstance(ind,Individu)
    assert genetique.len_pop == 10

    

def test_genetique_getter():
    genetique = Genetique(10)
    assert genetique.getBest_individu()==None
    genetique.refresh_score("AA")
    genetique.selection_elitisme()
    assert genetique.getBest_individu().score==genetique.population[0].score




def test_selection_elitisme():
    genetique = Genetique(10)
    genetique.refresh_score("AA")
    l = [x.score for x in genetique.population]
    l.sort()
    r = 0.5
    len = genetique.len_pop
    genetique.selection_elitisme(r)
    assert genetique.len_pop == int(len*r)
    l = l[:genetique.len_pop]
    l2 = [x.score for x in genetique.population]
    assert l==l2

def test_selection_elitisme():
    genetique = Genetique(10)
    genetique.refresh_score("AA")
    l = [x.score for x in genetique.population]
    l.sort()
    r = 0.5
    len = genetique.len_pop
    genetique.selection_elitisme(r)
    assert genetique.len_pop == int(len*r)
    l = l[:genetique.len_pop]
    l2 = [x.score for x in genetique.population]
    assert l==l2

def test_selection_roulette_bias_towards_low_score():
    # Vérifie que les individus à faible score sont sélectionnés plus souvent.
    for _ in range (200):
        genetique_instance = Genetique(200)
        genetique_instance.refresh_score("AA")
        genetique_instance.selection_roulette(rate=1.0) 
        selected_scores = [ind.score for ind in genetique_instance.population]

        low_score_count = sum(1 for score in selected_scores if score <= 5)
        high_score_count = sum(1 for score in selected_scores if score > 5)

    assert low_score_count >= high_score_count 

    genetique_instance.selection_roulette(rate=0.5)
    for individu in genetique_instance.population:
        assert individu in genetique_instance.population

    best_score = min(ind.score for ind in genetique_instance.population)

    assert genetique_instance.best_individu is not None 
    assert round(genetique_instance.best_individu.score,2) == round(best_score,2)


def test_selection_tournois_bias_towards_low_score():
    # Vérifie que les individus à faible score sont sélectionnés plus souvent.
    for _ in range (200):
        genetique_instance = Genetique(200)
        genetique_instance.refresh_score("AA")
        genetique_instance.selection_tournoi(rate=1.0)
        selected_scores = [ind.score for ind in genetique_instance.population]

        low_score_count = sum(1 for score in selected_scores if score <= 5)
        high_score_count = sum(1 for score in selected_scores if score > 5)

    assert low_score_count >= high_score_count 

    genetique_instance.selection_tournoi(rate=1.0)
    for individu in genetique_instance.population:
        assert individu in genetique_instance.population
    
    best_score = min(ind.score for ind in genetique_instance.population)

    assert genetique_instance.best_individu is not None 
    assert round(genetique_instance.best_individu.score,2) == round(best_score,2)


def test_refresh_score():
    genetique_instance = Genetique(200)
    genetique_instance.refresh_score("ATCG")
    genetique_instance.selection_elitisme(1.0)
    for ind in genetique_instance.population: assert ind.score!=None


def test_croisement_n_point_population():
    genetique_instance = Genetique(200)
    initial_size = genetique_instance.len_pop
    genetique_instance.croisement_n_point()
    
    expected_size = 2 * initial_size
    if expected_size % 2 != 0:  
        expected_size -= 1

    assert genetique_instance.len_pop == expected_size

    table_reference = Individu().getData().getTable()

    for individu in genetique_instance.population:
        assert all(dinucleotide in table_reference for dinucleotide in individu.data.getTable())

def test_algo_genetique():
    ind = algo_genetique("AA",10,True)
    assert isInBounds(ind.getData().getTable())
