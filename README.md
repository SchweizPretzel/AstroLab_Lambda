# AstroLab : Simulateur de projectiles dans l'espace

*Par Demba Cissé, Liam Feller, Marc-Antoine Gauthier et Levy Loïk Enok Koutoua dans le cadre du cours de Développement d'Applications dans un Environnement Graphique*

Une application graphique permettant de simuler des mouvements de projectile dans l’espace selon diverses conditions. L’application est destinée aux étudiants et professeurs du domaine de la physique de l’espace et de l’astronomie. Par la simulation, elle aide à comprendre et visualiser les mécanismes des mouvements de projectiles météoriques selon la gravité d’astres environnants, tout en offrant une flexibilité dans le positionnement d’astres.

## Installation et Lancement

Prérequis : Python 3.10 ou plus récent.

```bash
git clone https://github.com/Enok01/ProjetGroupe.git
cd ProjetGroupe
pip install -r requirements.txt
python -m src.main
```

> Le programme se lance depuis la racine du dépôt : tous les imports passent par `src.` (ex. `from src.astre import AstreFixe`).

## Notions élémentaires de physique

- **Astres fixes** : chaque astre tourne en cercle autour du centre (0, 0) qui est le soleil, ou autour de son `parent` pour la lune (Demba implémente la fonctionnalité en extra). Sa position est donnée directement par une formule, sans intégration :
  θ(t) = 2π·t / T + φ, puis r(t) = d · (cos θ, sin θ).
  Les astres ne s'attirent pas entre eux.
- **Projectiles** : seuls les projectiles subissent la gravité de tous les astres **actifs** (loi de Newton) :
  a = Σᵢ G·Mᵢ·(rᵢ − r) / ‖rᵢ − r‖³.
  Leur trajectoire est calculée pas à pas (méthode de Verlet a confirmer).
- **Fin d'un projectile** : collision quand sa distance à un astre est inférieure au rayon de cet astre, ou expiration quand il dépasse sa durée de vie.

Unités utilisées dans tout le programme :

| Grandeur           | Unité                                           |
|--------------------|-------------------------------------------------|
| Distance, position | UA (unité astronomique)                         |
| Temps              | jour                                            |
| Masse              | masse terrestre (Terre = 1)                     |
| Vitesse            | UA / jour                                       |
| G                  | ≈ 8,888 × 10⁻¹⁰ UA³ / (masse terrestre · jour²) |

## Format du fichier CSV  

Un seul fichier décrit les astres et les projectiles ;
la colonne `type` indique lequel, et les colonnes qui ne s'appliquent pas restent vides.


## Pseudo-Preuve de la méthode de détection des collisions

Posons le projectile $p$ et l'astre le plus proche de celui-ci $a$,
ainsi que leurs vecteurs positions respectifs $\vec{r_p}$ et $\vec{r_a}$
et leurs vitesses respectives $\vec{v_p}$ et $\vec{v_a}$. Soit :

$$\vec{r_p} = \begin{bmatrix} x_p & y_p & z_p \end{bmatrix}\hspace{1cm}
\vec{v_p} = \begin{bmatrix} v_{xp} & v_{yp} & v_{zp} \end{bmatrix}$$

$$\vec{r_a} = \begin{bmatrix} x_a & y_a & z_a \end{bmatrix}\hspace{1cm}
\vec{v_a} = \begin{bmatrix} v_{xa} & v_{ya} & v_{za} \end{bmatrix}$$

Posons également deux vecteurs, $\vec{s_p}$ et $\vec{s_a}$, dépendants de $t$,
un paramètre d'interpolation dans l'intervale $[0,1]$.

$$\vec{s_p} = \vec{r_p} + t\vec{v_p}$$

$$\vec{s_p} = \vec{r_a} + t\vec{v_a}$$

Finalement, posons $R_p$, le rayon du projectile, et $R_a$, le rayon de l'astre.

$\exists t \in \left[0,1\right] \hspace{.1cm} | \hspace{.1cm} \|\vec{s_p} - \vec{s_a}\| ≤ R_p + R_a \implies$ il y a une collision

Donc, si l'équation qui suit a une solution, il y a une collision :
$$\|\vec{s_p} - \vec{s_a}\|^2 - (R_p + R_a)^2 = 0$$

Développons $\|\vec{s_p} - \vec{s_a}\|^2$.

$\|\vec{s_p} - \vec{s_a}\|^2$
$= (s_{xp} - s_{xa})^2 + (s_{yp} - s_{ya})^2 + (s_{zp} - s_{za})^2$

$= (x_p + v_{xp}t - x_a - v_{xa}t)^2 + (y_p + v_{yp}t - y_a - v_{ya}t)^2 + (z_p + v_{zp}t - z_a - v_{za}t)^2$

Posons $\Delta \vec{r} = \vec{r}_p - \vec{r}_a$ et $\Delta \vec{v} = \vec{v}_p - \vec{v}_a$
ainsi que leurs composantes en utilisant la même convention que précédemment.

$= (\Delta x + \Delta v_{xp}t)^2 + (\Delta y + \Delta v_{yp}t)^2 + (\Delta z + \Delta v_{zp}t)^2$

$= \Delta x^2 + 2\Delta x \Delta v_{xp}t + \Delta v_{xp}^2t^2 + \Delta y^2 + 2\Delta y \Delta v_{yp}t + \Delta v_{yp}^2t^2 + \Delta z^2 + 2\Delta z \Delta v_{zp}t + \Delta v_{xp}^2t^2$

$= (\Delta v_{xp}^2 + \Delta v_{yp}^2 + \Delta v_{zp}^2)t^2 + 2(\Delta x \Delta v_{xp} + \Delta y \Delta v_{yp} + \Delta z \Delta v_{zp})t + (\Delta x^2 + \Delta y^2 + \Delta z^2)$

$= (\Delta \vec{v} \cdot \Delta \vec{v}) t^2 + 2(\Delta \vec{r} \cdot \Delta \vec{v}) t + (\Delta \vec{r} \cdot \Delta \vec{r})$

Remplaçant $\|\vec{s_p} - \vec{s_a}\|^2$ dans l'équation d'origine, on obtient :

$(\Delta \vec{v} \cdot \Delta \vec{v}) t^2 + 2(\Delta \vec{r} \cdot \Delta \vec{v}) t + (\Delta \vec{r} \cdot \Delta \vec{r}) - (R_p + R_a)^2$ = 0

... ce qui est une équation quadratique, donc nous pouvons utiliser le discriminant afin de qualifier ses solutions :

$D = b^2-4ac$

$D = 4(\Delta \vec{r} \cdot \Delta \vec{v})^2 - 4(\Delta \vec{v} \cdot \Delta \vec{v}) \left((\Delta \vec{r} \cdot \Delta \vec{r}) - R_p - R_a\right)$

Si $D<0$, il n'y a pas de collision. Sinon, continuons.

$t = \frac{-b \pm \sqrt{D}}{2a}$

$t_1 = \frac{-2(\Delta \vec{r} \cdot \Delta \vec{v}) - \sqrt{D}}{2\Delta \vec{v} \cdot \Delta \vec{v}}$ et 
$t_2 = \frac{-2(\Delta \vec{r} \cdot \Delta \vec{v}) + \sqrt{D}}{2\Delta \vec{v} \cdot \Delta \vec{v}}$

Si $\left(t_1 \in [0,1]\right) \vee \left(t_2 \in [0,1]\right)$, alors il y a une collision.

$\blacksquare$

À noter que la preuve est valide en 2D également, il faut simplement omettre la composante en z pour chaque vecteur.
