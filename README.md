# EI Algorithme génétique

## Installation

Un environnement Python est nécessaire. Ce repo a été testé sous python 3.12, mais il devrait fonctionner sous la majorité des versions 3.

- `pip install -r requirements.txt` permet d'installer **tous les paquets** requis.

## Utilisation

Tout se lance à partir d'une commande dans le terminal : `python -m dna`
Cette commande va construire et afficher le brin d'adn extrait de _data/plasmid_8k.fasta_ en utilisant le modèle de conformation de 1993 : _dna/table.json_.

Plusieurs paramètres permettent de changer ce comportement. Ils sont explicités en tapant `python -m dna --help`, et détaillés ici :

- `-m [traditional | recuit | genetic] (par défaut : traditional)` permet de choisir le mode de fonctionnement. Se référer au paragraphe "Modes".
- `-d [path_to_file] (par défaut : data/plasmid_8k.fasta)` permet de choisir une séquence d'ADN. Inutile pour le mode `recuit` qui s'entraîne sur les deux séquences automatiquement.
- `-j [path_to_file] (par défaut : dna/table.json)` permet de choisir un modèle de conformation. Dans le mode `traditional`, c'est celui qui sera utilisé pour le tracé. Dans les autres, ce sera le modèle à partir duquel évoluer.
- `-i [positive integer] (par défaut : 100)`permet de choisir le nombre d'itérations maximum. Utile pour le mode `recuit`
- `-p [positive integer] (par défaut : 10)`permet de choisir la taille de la population. Utile pour le mode `genetic`
- `-s`permet de tracer des courbes comparatives des différents modes de séléction. Nécessite le mode `genetic` pour fonctionner

## Tests

Une couverture du code est réalisée grâce à `pytest`.
La version du rapport la plus récente est directement accessible par le raccourci _Coverage.lnk_.
Pour regénérer le rapport, il faut éxécuter la commande suivante : `pytest --cov=dna --cov-report html test_*.py`.

## Modes

- **traditional :** Ce mode calcule la trajectoire dans l'espace d'une séquence ADN et l'affiche.
- **recuit :** Ce mode essaie d'affiner le modèle de conformation d'entrée afin de faire boucler les séquences ADN. Un algorithme de recuit est utilisé, avec une énergie dépendant principalement de la distance entre le début de la séquence et sa fin.
- **genetic :** Ce mode essaie d'affiner le modèle de conformation d'entrée afin de faire boucler les séquences ADN. Il utilise un algorithme génétique.
