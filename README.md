<h1 align="center"> 💻 Guides pour Cage Moléculaire ⚛️</h1>
<p>
</p>

<!-- badges: start -->
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Génération de guides de construction de cage moléculaire spécifique à un substrat. La génération se fait par ajout de motifs liants intéragissants avec le substrat puis par création de chemins moléculaires reliant les motifs liants.

En cours : Création des chemins entre les motifs liants et amélioration des motifs liants.

>Projet de recherche de l'équipe ALMOST du [laboratoire DAVID](https://www.david.uvsq.fr/accueil/).

## Compilation

Pour compiler le projet utilisez :
```sh
make
```

## Utilisation

### Démo
Pour lancer une démo sur le substrat YILLAG :
```sh
make demo
```

### Exécution
Pour exécuter le programme entrez :
```sh
./bin/cageMol.exe -i [fichier_substrat.xyz]
```
Puis les paramètres alpha et sizemax peuvent être aussi modifiés.
Alpha est utilisé pour la génération d'une enveloppe concave et sizemax correspond au nombre d'atomes maximum que l'on veut dans un chemin qu'on génère.
```sh
alpha (défaut 3) : -a [double]

sizemax (défaut 5) : -s [entier]
```
Pour avoir de l'aide : 
 ```sh
-h
```

### Nettoyage des fichiers

Pour supprimer l'éxécutable et les fichiers objets :
```sh
make clean
```
Pour supprimer en plus les résultats : 
```sh
make mrproper
```

### Visualisation

Pour visualiser les résultats vous pouvez utiliser Pymol ou tout autres logiciels de visualisation moléculaire.

## Documentation 

Ouvrir avec un navigateur web le fichier index.html, se trouvant dans doc/html. [Documentation](./doc/html/index.html)

## À faire

Pathfinding : implémenter l'algorithme JPS (article : The Jump Point Search Pathfinding System in 3D) pour essayer de remplacer A*

Pour convenir aux chimistes : avoir des chemins qui se ressemble au maximum "symétrie" (améliore la création de la cage), avoir des motifs avec 3 chemins pour tridimensionner la cage

Optimisation de l'ensemble

## publications

### Thèse

**Marie Bricage**. Modélisation et Algorithmique de graphes pour la construction de structures moléculaires.. Géométrie algorithmique [cs.CG]. Université Paris Saclay (COmUE), 2018. Français. [⟨NNT : 2018SACLV031⟩](https://www.theses.fr/2018SACLV031). [⟨tel-01955838⟩](https://theses.hal.science/tel-01955838)

## Auteurs

<div align="center">

👤 **Anne FERNET** : Github: [@uvsq21915170](https://github.com/uvsq21915170)

👤 **Noé DEMANGE** : Github: [@NoeDemange](https://github.com/NoeDemange)

👤 **Priscille DAOULAS** :  Github: [@priscdls](https://github.com/priscdls)

👤 **Marie BRICAGE** : Github: [@Isima](https://github.com/Isima)

</div>

## Encadrants
 <div align="center">
  <b>Sandrine VIAL &emsp; Yann STROZECKI &emsp; Dominique BARTH</b>
</div>

## Contributeurs
<div align="center">
  <b>Anne DUVEAU &emsp; Rafael PREAULT &emsp; Alexis GUIGAL</b>
</div>

### Chimiste

<div align="center">
  <b>Olivier DAVID</b>
</div>

***
<div>
    <table align="center">
        <tr>
            <td>
                <div>
                    <a href="https://www.david.uvsq.fr">
                    <img width=150px src="https://www.david.uvsq.fr/wp-content/themes/david/src/img/logo_david.svg" alt="Logo Laboratoire david"</a>
                </div>
            </td>
            <td>
                <a href="http://www.uvsq.fr/">
                <img src="https://www.david.uvsq.fr/wp-content/themes/david/src/img/partners/logo-universite-versailles.svg" alt="Université de Versailles"></a>
            </td>
            <td>
                <a href="https://www.universite-paris-saclay.fr/fr">
                <img src="https://www.david.uvsq.fr/wp-content/themes/david/src/img/partners/logo-universite-paris-saclay.svg" alt="Université de Paris Saclay"></a>
            </td>
        </tr>
    </table>
</div>