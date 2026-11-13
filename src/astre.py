"""
Classes des objets de la simulation (Levy) : Astre (classe abstraite) et AstreFixe (orbite circulaire).

Documentation officielle :
    numpy : https://numpy.org/doc/stable/reference/
    classes abstraites (abc) : https://docs.python.org/3/library/abc.html
    pandas Series.get (lecture d'une ligne du CSV) : https://pandas.pydata.org/docs/reference/api/pandas.Series.get.html
"""
from abc import ABC, abstractmethod

import numpy as np

from src.exceptions import DebutHistorique, FinHistorique

#Constante gravitationnelle convertie dans les unités du programme : UA³ / (masse terrestre · jour²) ≈ 8,888e-10
G_SI = 6.67430e-11          # m³ / (kg · s²)
MT_EN_KG = 5.9722e24        # masse terrestre en kg
JOUR_EN_S = 86400.0         # jour en secondes
UA_EN_M = 149597870700.0    # unité astronomique en mètres
G = G_SI * MT_EN_KG * JOUR_EN_S ** 2 / UA_EN_M ** 3


def en_booleen(valeur) -> bool:
    """_summary_

    Args:
        valeur (_type_): _description_

    Raises:
        ValueError: _description_

    Returns:
        bool: _description_
    """
    #Convertit la valeur « actif » du CSV ou d'un formulaire en booléen.
    #Attention : bool("false") vaut True en Python, d'où cette fonction.
    #np.bool_ : le type booléen de numpy, que pandas donne quand il lit true/false dans le CSV
    #https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.bool
    if isinstance(valeur, (bool, np.bool_)):
        return bool(valeur)
    texte = str(valeur).strip().lower()
    if texte in ("true", "vrai", "oui", "1", "1.0"):
        return True
    if texte in ("false", "faux", "non", "0", "0.0"):
        return False
    raise ValueError(f"« {valeur} » n'est pas une valeur valide pour actif (attendu true ou false)")


class Astre(ABC):
    
    #Classe abstraite commune aux astres fixes et aux projectiles.
    """
    Attributs :
        nom : nom unique de l'objet.
        masse : masse en masses terrestres (Terre = 1), comme dans systeme_solaire.csv.
        couleur : couleur d'affichage dans le graphique matplotlib (ex. "orange" ou "#FDB813").
        taille : taille d'affichage (rayon relatif à celui de la Terre dans le CSV).
        rayon : rayon de collision en UA.
        actif : un objet inactif n'attire pas et n'est pas simulé.
    """

    #Constructeur de la classe Astre
    def __init__(self, nom: str, masse: float, couleur: str, taille: float, rayon: float, actif: bool = True):
        self.nom = nom
        self.couleur = couleur
        self.actif = en_booleen(actif)
        #Les valeurs numériques peuvent venir du CSV ou d'un formulaire : on vérifie leur conversion
        try:
            self.masse = float(masse)
            self.taille = float(taille)
            self.rayon = float(rayon)
        except (TypeError, ValueError):
            raise ValueError(f"Astre « {nom} » : masse, taille et rayon doivent être des nombres") from None
        if self.masse < 0 or self.rayon < 0 or self.taille <= 0:
            raise ValueError(f"Astre « {nom} » : masse et rayon doivent être positifs, taille strictement positive")

    @abstractmethod
    def position(self, t: float) -> np.ndarray:
        #Renvoie la position (x, y) en UA au temps t (jours), tableau numpy de forme (2,).
        pass


class AstreFixe(Astre):
    #Astre en orbite circulaire autour du centre (0, 0) ou autour de son parent (lune).
    """
    Attributs :
        distance : rayon de l'orbite en UA (0 = immobile au centre, comme le Soleil).
        periode : période de révolution en jours (0 = immobile).
        phase : angle de l'astre à t = 0, en degrés (comme dans systeme_solaire.csv).
        parent : astre autour duquel il tourne, ou None pour le centre.
    """

    #Constructeur de la classe AstreFixe
    def __init__(self, nom: str, masse: float, couleur: str, taille: float, rayon: float, actif: bool = True,
                 distance: float = 0.0, periode: float = 0.0, phase: float = 0.0, parent=None):
        """_summary_

        Args:
            nom (str): _description_
            masse (float): _description_
            couleur (str): _description_
            taille (float): _description_
            rayon (float): _description_
            actif (bool, optional): _description_. Defaults to True.
            distance (float, optional): _description_. Defaults to 0.0.
            periode (float, optional): _description_. Defaults to 0.0.
            phase (float, optional): _description_. Defaults to 0.0.
            parent (_type_, optional): _description_. Defaults to None.

        Raises:
            ValueError: _description_
            ValueError: _description_
            ValueError: _description_
        """
        super().__init__(nom, masse, couleur, taille, rayon, actif)
        try:
            self.distance = float(distance)
            self.periode = float(periode)
            self.phase = float(phase)
        except (TypeError, ValueError):
            raise ValueError(f"AstreFixe « {nom} » : distance, période et phase doivent être des nombres") from None
        if self.distance < 0 or self.periode < 0:
            raise ValueError(f"AstreFixe « {nom} » : distance et période doivent être positives")
        if parent is self:
            raise ValueError(f"AstreFixe « {nom} » ne peut pas être son propre parent")
        self.parent = parent

    def _temps(self, t) -> float:
        #Convertit le temps en nombre de jours, avec un message clair si t n'est pas un nombre.
        try:
            return float(t)
        except (TypeError, ValueError):
            raise ValueError(f"AstreFixe « {self.nom} » : le temps doit être un nombre de jours, reçu {t!r}") from None

    def _angle(self, t: float) -> float:
        #Angle en radians au temps t : theta = 2*pi*t/T + phase (l'angle reste égal à la phase si T = 0).
        #La phase est donnée en degrés, on la convertit en radians pour cos et sin.
        t = self._temps(t)
        #np.radians convertit des degrés en radians : https://numpy.org/doc/stable/reference/generated/numpy.radians.html
        #np.pi vaut 3.14159... : https://numpy.org/doc/stable/reference/constants.html
        phase = np.radians(self.phase)
        if self.periode > 0:
            return 2 * np.pi * t / self.periode + phase
        return phase

    def position(self, t: float) -> np.ndarray:
        #Position r(t) = centre + d * (cos(theta), sin(theta)), le centre étant (0, 0) ou le parent.
        theta = self._angle(t)
        centre = np.zeros(2) if self.parent is None else self.parent.position(t)
        #np.array crée le vecteur [x, y] ; np.cos et np.sin prennent un angle en radians
        #https://numpy.org/doc/stable/reference/generated/numpy.array.html
        #https://numpy.org/doc/stable/reference/generated/numpy.cos.html
        return centre + self.distance * np.array([np.cos(theta), np.sin(theta)])

    @classmethod
    def depuis_ligne(cls, ligne):
        #Crée un AstreFixe à partir d'une ligne du CSV (dictionnaire ou ligne pandas).
        #Les cases vides valent la valeur par défaut ; le parent (lune) est relié ensuite par le chargeur CSV.
        try:
            nom = str(ligne["nom"]).strip()
        except (KeyError, TypeError):
            raise ValueError("Ligne CSV invalide : la colonne « nom » est manquante") from None

        def nombre(cle, defaut=0.0):
            #ligne.get(cle, defaut) : la valeur de la colonne, ou defaut si la colonne n'existe pas
            #(marche pour un dictionnaire et pour une ligne pandas) : https://pandas.pydata.org/docs/reference/api/pandas.Series.get.html
            valeur = ligne.get(cle, defaut)
            if valeur is None or str(valeur).strip() == "":
                return defaut
            try:
                valeur = float(valeur)
            except (TypeError, ValueError):
                raise ValueError(f"Astre « {nom} », colonne « {cle} » : « {valeur} » n'est pas un nombre") from None
            #np.isnan : vrai si la valeur est NaN, ce que pandas met dans une case vide du CSV
            #https://numpy.org/doc/stable/reference/generated/numpy.isnan.html
            return defaut if np.isnan(valeur) else valeur

        return cls(
            nom=nom,
            masse=nombre("masse"),
            couleur=str(ligne.get("couleur") or "#FFFFFF").strip(),
            taille=nombre("taille", 5.0),
            rayon=nombre("rayon"),
            actif=en_booleen(ligne.get("actif", True)),
            distance=nombre("distance"),
            periode=nombre("periode"),
            phase=nombre("phase"),
        )


class Projectile(Astre):
    def __init__(self, nom, masse, couleur, taille, rayon, pos0,
                 vit0, temps_apparition, duree_vie, actif=True):
        """_summary_

        Args:
            nom (_type_): _description_
            masse (_type_): _description_
            couleur (_type_): _description_
            taille (_type_): _description_
            rayon (_type_): _description_
            pos0 (_type_): _description_
            vit0 (_type_): _description_
            temps_apparition (_type_): _description_
            duree_vie (_type_): _description_
            actif (bool, optional): _description_. Defaults to True.
        """

        super().__init__(nom, masse, couleur, taille, rayon, actif)
        self.pos0 = np.array(pos0, dtype=float)  # Position initiale du projectile, liste de coordonnées [x, y]
        self.vit0 = np.array(vit0, dtype=float)  # Vitesse initiale du projectile, liste de coordonnées [vx, vy]
        self.temps_apparition = temps_apparition  # Temps d'apparition du projectile : question a poser...
        self.duree_vie = duree_vie  # Durée de vie du projectile: question a poser...

        self.historique = np.array([[*self.pos0,
                                     *self.vit0]])  # Liste pour stocker l'historique des positions du projectile : unpacking de pos0 et vit0 pour créer un tableau 2D avec les coordonnées et vitesses initiales par rappre a un N
        self.etat = "actif"  # État du projectile : peut être actif ou inactif, mais obligatoirement actifs sa création, question a poser...: Que faire si on réinitialise mais que il avait été désactivé par le gestionnaire d'objets ?
        self.cible = None  # Cible du projectile?: question a poser...

    def reinitialiser(self):
        """
         hypothèses: lors de la réinitialisation on recommence la simulation depuis
         le début, donc je vais supprimer l'historique, remettre son état à l'état de base
         et la cible (je sais toujours pas le ça à quoi) mais on va la réinitialiser
         parce qu'elle peut changer, les autres paramètres ne peuvent pas changer,
         donc on n'a pas besoin de le réinitialiser
        """

        self.historique = np.array([[*self.pos0, *self.vit0]])  # on oublie les pas calculés
        self.etat = "actif"  # Réinitialiser l'état à actif
        self.cible = None  # Réinitialiser la cible à None

    def ajouter_etat(self, pos, vit):
        """
        Ajoute un nouvel état du projectile à l'historique.
            Args:
                position (list): Position actuelle du projectile [x, y].
                vitesse (list): Vitesse actuelle du projectile [vx, vy].
            """

        nouvelle_ligne = np.array([[*pos, *vit]])  # Créer une nouvelle ligne avec la position et la vitesse
        self.historique = np.vstack((self.historique, nouvelle_ligne))  # Ajouter la ligne à l'historique

        """ 
            est ce que je dois mettre à jour l'état du projectile en fonction de la durée de vie?
            ou la simulation va gérer ça? je pense que c'est la simulation qui va gérer ça, 
            donc je ne vais pas le faire ici:C’est la Simulation qui compare le temps actuel à t_apparition  + duree_vie, 
            qui met etat = "expire" et qui arrête ensuite d’appeler ajouter_etat

            """

    def position(self, t):
        """
        Calcule la position du projectile à un instant donné t.

        Args:
            t (float): Temps de la simulation, en jours (1 pas = 1 journée)

        Returns:
            np.array: Position du projectile à l'instant t.
        """
        try:
            return self.etat_a(t)[:2]  # Retourne seulement les coordonnées x et y de l'état complet à l'instant t
        except DebutHistorique as e:
            raise DebutHistorique(str(e))
        except FinHistorique as e:
            raise FinHistorique(str(e))

    def etat_a(self, pas):
        """
        Retourne l'état complet du projectile à un pas donné de la simulation.

        Args:
        pas (int): Pas de la simulation (1 pas = 1 journée).

        Returns:
        np.ndarray: [x, y, vx, vy]. Retourne l'état au pas demandé
        Raises:
            DebutHistorique: Si le pas est avant l'apparition du projectile.
            FinHistorique: Si le pas dépasse les pas déjà calculés.
        """
        if pas < self.temps_apparition:
            raise DebutHistorique("Le projectile n'a pas encore été créé")  # Retourne l'état initial avant l'apparition

        ligne = int(pas - self.temps_apparition)  # Nombre de pas depuis l'apparition

        if ligne < len(self.historique):
            return self.historique[ligne]  # Retourne ligne conplete
        raise FinHistorique(
            "Le projectile ne s'est pas rendu là")  # Retourne le dernier état connu si le pas dépasse les pas calculés
