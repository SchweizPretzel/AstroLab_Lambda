import math

from astre import Astre
from src.astre import AstreFixe, Projectile, Vecteur2, G
from src.exceptions import HorsSystemeError, TempsMaximumAtteinte

import numpy as np

DISTANCE_MAX_AJOUTER = 2  # Distance en UA ajoutées pour les bornes maximales du projectile pour le considérer comme 'hors du système'


def calculer_trajectoire(projectile: Projectile, liste_astre: list[AstreFixe], dt: float = 0.1, temps_maximum: int = 3653) -> AstreFixe:
    """_Crée la trajectoire du projectile en calculant sa position et sa vitesse sous les gravités des astres du système.
        Les états calculés sont enregistrés dans l'historique du projectile.
        Retourne l'astre fixe toucher ou Souleve un Exception selon ce qui ce passe à la place_

    Args:
        projectile (Projectile): _Projectille utiliser pour la trajectoire_
        liste_astre (list[AstreFixe]): _Système des astres fixes dans une liste_
        dt (float, optional): _pas de temps entre les positions (en jour)_. Initialiser à 0.1 (recommandé pour la rapidité) si plus précis et un peux moins rapide -> 0.01.
        temps_maximum (int, optional): _nombre de position calculée maximum_. Initialisé à 3653 jours (10 ans).

    Raises:
        HorsSystemeError: _si le projectile est hors du système._
        TempsMaximumAtteinte: _si le nombre de temps maximum a été atteint._

    Returns:
        AstreFixe: _retourne l'astre fixe seulement s'il a eu contacte avec un astre fixe_
    """

    # Vérifi si toutes les valeurs données sont correct
    if not isinstance(projectile, Projectile):
        raise ValueError("Physique : p doit etre un Projectile.")
    if not all(isinstance(astre, AstreFixe) for astre in liste_astre):
        raise ValueError(
            "Physique : Il faut que toutes les astres dans liste_astre soit des AstreFixe."
        )
    if not (isinstance(dt, (int, float)) and isinstance(temps_maximum, (int, float))):
        raise ValueError("Physique : dt et temps_maximum doivent être des nombres.")
    if dt < 0 or temps_maximum < 0:
        raise ValueError("Physique : dt et temps_maximum doivent être positives.")

    Nb_iterations_enregistrement = int(1 / dt)

    # détermine la longueur maximale pour arrêter la création de la liste
    astre_loin = max(liste_astre, key=lambda astre: astre.distance)
    distance_max = astre_loin.distance + DISTANCE_MAX_AJOUTER

    # créer une liste de position que le projectile a faite
    pos = projectile.pos0
    vit = projectile.vit0
    t = projectile.temps_apparition  # en jour : les astres partent de leur position au lancement
    nombre_iteration_max = int(temps_maximum / dt)

    for i in range(nombre_iteration_max):

        new_pos, new_vit = prochaine_pos_vit(pos, liste_astre, dt, vit, t)
        contacte, new_pos = est_en_contacte(pos, new_pos, projectile, liste_astre, t, dt)

        if contacte is not None: # s'il y a un contacte
            projectile.ajouter_etat(new_pos, new_vit)  # ajoute la dernière position du projectile (sur l’astre)
            return contacte

        if new_pos.norme2() > distance_max**2:
            raise HorsSystemeError("Le projectile est hors du système.")

        # Enregistre seulement les positions une fois par jour
        if i % Nb_iterations_enregistrement == 0 and i != 0:
            projectile.ajouter_etat(new_pos, new_vit)

        pos, vit = new_pos, new_vit
        t += dt

    raise TempsMaximumAtteinte("Le nombre de temps maximum a été atteint.")


def acceleration(pos: Vecteur2, astres: list[AstreFixe], t: float) -> Vecteur2:
    """Calcule l'accélération (en UA/j²) que les astres actifs exercent sur le projectile a la position `pos`

    Args:
        t (float): le temps auquel on calcule l'acceleration
        pos (Vecteur2): la position du projectile (en UA)
        astres (list[AstreFixe]): liste qui contient tous les astres fixes

    Returns:
        Vecteur2: accélération du projectile (en UA/j²)
    """
    # Un astre inactif n'attire pas
    actifs = [a for a in astres if a.actif]
    if not actifs:
        return Vecteur2.nul()

    # Tous les astres d'un coup avec numpy (une boucle Python par astre était le plus lent du programme)
    masses = np.array([a.masse for a in actifs])
    vecteurs_astre_projectile = np.array([(a.position(t) - pos).unpack() for a in actifs])
    distances2 = (vecteurs_astre_projectile ** 2).sum(axis=1)
    masses[distances2 < 1e-12] = 0  # Pour éviter les divisions par 0 : un astre au même point n'attire pas
    distances2 = np.maximum(distances2, 1e-12)

    # Fg = Gm1m2/d**2 et a = F/m1, donc a = G * m2/d**2 dans la direction vecteur / d, soit G * m2 * vecteur / d**3.
    # @ (produit matriciel) additionne les contributions de tous les astres : https://numpy.org/doc/stable/reference/generated/numpy.matmul.html
    somme_fg = (G * masses / distances2 ** 1.5) @ vecteurs_astre_projectile
    return Vecteur2(somme_fg[0], somme_fg[1])


def prochaine_pos_vit(pos: Vecteur2, astres: list[AstreFixe], dt: float, v: Vecteur2, t: float) -> tuple[Vecteur2, Vecteur2]:
    """Calcule la prochaine vitesse et la prochaine position du projectile en utilisant
    l’algorithme de Verlet pour la conservation d’énergie à long terme.
    (Pour plus d’information, consulter la vidéo suivante : https://www.youtube.com/watch?v=ii-N61PeeQU&t=747s)

    Args:
        pos (Vecteur2): la position x et y du projectile à ce moment (en UA)
        astres (list[AstreFixe]): liste qui contient toutes les positions des astres et les astres fixes du systeme
        dt (float): pas de temps (en j)
        v (Vecteur2): vitesse du projectile (en UA/j)

    Returns:
        tuple (new_pos, new_vit) :
            new_pos -> nouvelle position x et y du projectile (en UA)

            new_vit -> nouvelle vitesse x et y du projectile (en UA/j)
    """
    # 1- Trouvons a(t) -> a
    a = acceleration(pos, astres, t)

    # 2- Trouvons la position au temps t + dt -> new pos
    new_pos = pos + v * dt + 0.5 * a * dt**2

    # 3- Trouvons a(t + dt) -> new_a
    new_a = acceleration(new_pos, astres, t + dt)

    # 4- Trouvons v(t+dt) -> new_v
    new_v = v + 1 / 2 * (new_a + a) * dt

    return new_pos, new_v


def est_en_contacte(pos: Vecteur2, new_pos: Vecteur2, p: Projectile, astres: list[AstreFixe],
                    t: float, dt: float) -> tuple[AstreFixe | None, Vecteur2]:
    """_Vérifie s'il est en contact ou non en calculant si la trajectoire des deux points a des points avec l’astre le plus proche_

    Args:
        dt (float): différentielle du temps
        t (float): temps
        pos (Vecteur2): _Position x et y du projectile au moment précédant (en UA)_
        new_pos (Vecteur2): _Nouvelle position x et y du projectile maintenant (en UA)_
        p (Projectile): _Projectile lancé_
        astres (list[AstreFixe]): _Liste qui contient touts les astres fixes du système_

    Returns:
        AstreFixe: _l'astre touché et la position où elle a toucher OU None s'il ne touche rien_
    """

    # Trouver l'astre le plus proche du projectile
    astre_proche = None
    astre_proche_position = None
    distance_min = float("inf")
    for a in astres:
        pos_astre = a.position(t)
        d2 = (pos_astre - new_pos).norme2()

        if d2 < distance_min:
            distance_min = d2
            astre_proche = a
            astre_proche_position = pos_astre

    if astre_proche is None:
        return None, new_pos
    astre_proche_position: Vecteur2 # Redéclare pour que PyCharm ferme sa gueule

    rayon_max = astre_proche.rayon + p.rayon

    if abs(new_pos.x - pos.x) < 1e-12:  # éviter une division par zéro
        return None, new_pos

    # Détection des collisions par méthode vectorielle
    v_proj = 1/dt * (new_pos - pos)

    delta_r = pos - astre_proche_position
    delta_v = v_proj - astre_proche.vitesse(t)

    A = delta_v.norme2()
    B = delta_r.dotprod(delta_v)
    C = delta_r.norme2() - rayon_max * rayon_max


    discriminant = B * B - 4 * A * C

    if discriminant < 0:
        return None, new_pos

    t1 = (-B - math.sqrt(discriminant)) / (2 * A)
    t2 = (-B + math.sqrt(discriminant)) / (2 * A)

    if 0 < t1 < 1:
        return astre_proche, astre_proche_position

    elif 0 < t2 < 1:
        return astre_proche, astre_proche_position
    return None, new_pos
