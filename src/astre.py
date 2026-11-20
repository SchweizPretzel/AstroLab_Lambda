"""
Classes des objets de la simulation (Levy) : Astre (classe abstraite) et AstreFixe (orbite circulaire).

Documentation officielle :
    numpy : https://numpy.org/doc/stable/reference/
    classes abstraites (abc) : https://docs.python.org/3/library/abc.html
    pandas Series.get (lecture d'une ligne du CSV) : https://pandas.pydata.org/docs/reference/api/pandas.Series.get.html
"""
from abc import ABC, abstractmethod

import math
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
    """
    Class abstraite pour AstreFixe et Projectile
    """
    def __init__(self, nom: str, masse: float, couleur: str, taille: float, rayon: float, actif: bool = True):
        """

        Args:
            nom (str): Le nom de l'astre
            masse (float): La masse de l'astre, en MT.
            couleur (str): La couleur d'affichage de l'astre. Doit être lisible par matplotlib.
            taille (float): La taille d'affichage de l'astre.
            rayon (float): Le rayon de l'astre, en UA.
            actif (bool): Si l'astre devrait être simulé ou affiché.
        """
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
    def position(self, t: float) -> Vecteur2:
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

    def position(self, t: float) -> Vecteur2:
        # Position r(t) = centre + d * (cos(theta), sin(theta)), le centre étant (0, 0) ou le parent.
        theta = self._angle(t)
        if self.parent is None:
            centre = Vecteur2.nul()
        else:
            centre = self.parent.position(t)

        return centre + Vecteur2.v_polaire(self.distance, theta)

    def vitesse(self, t: float) -> Vecteur2:
        theta = self._angle(t)
        if self.parent is None:
            v_parent = Vecteur2.nul()
        else:
            v_parent = self.parent.vitesse()

        if self.periode > 0:
            return v_parent + Vecteur2.v_polaire(2 * np.pi * self.distance, theta)
        else:
            return Vecteur2.nul()


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
    def __init__(self, nom: str, masse: float, couleur: str, taille: float, rayon: float, pos0: Vecteur2,
                 vit0: Vecteur2, temps_apparition: int, duree_vie: int, actif=True):
        """Créé un nouveau projectile qui peut être simulé.

        Args:
            nom (str): Le nom du projectile.
            masse (float): La masse du projectile, en MT.
            couleur (str): La couleur d'affichage du projectile. Doit être lisible par matplotlib.
            taille (float): La taille d'affichage du projectile.
            rayon (float): Le rayon du projectile, en UA.
            pos0 (Vecteur2): Le vecteur de position initiale du projectile, en UA.
            vit0 (Vecteur2): Le vecteur de vitesse initiale du projectile, en UA/jour
            temps_apparition (int): La valeur de temps où le projectile apparaît.
            duree_vie (int): La durée de vie du projectile.
            actif (bool): Si l'astre devrait être simulé ou affiché.
        """

        super().__init__(nom, masse, couleur, taille, rayon, actif)
        self.pos0 = pos0  # Position initiale du projectile, liste de coordonnées [x, y]
        self.vit0 = vit0  # Vitesse initiale du projectile, liste de coordonnées [vx, vy]
        self.temps_apparition = temps_apparition  # Temps d'apparition du projectile : question a poser...
        self.duree_vie = duree_vie  # Durée de vie du projectile: question à poser...

        self.historique = np.array([[*(self.pos0.unpack()), *(self.vit0.unpack())]])
        self.etat = "actif"
        self.cible = None

    def reinitialiser(self):
        self.historique = np.array([[*(self.pos0.unpack()), *(self.vit0.unpack())]])  # on oublie les pas calculés
        self.etat = "actif"  # Réinitialiser l'état à actif
        self.cible = None  # Réinitialiser la cible à None

    def ajouter_etat(self, pos, vit):
        """
        Ajoute un nouvel état du projectile à l'historique.
            Args:
                pos (Vecteur2): Position actuelle du projectile.
                vit (Vecteur2): Vitesse actuelle du projectile.
            """

        nouvelle_ligne = np.array([[pos.x, pos.y, vit.x, vit.y]])  # Créer une nouvelle ligne avec la position et la vitesse
        self.historique = np.vstack((self.historique, nouvelle_ligne))  # Ajouter la ligne à l'historique

    def position(self, t: float):
        """
        Calcule la position du projectile à un instant donné t.

        Args:
            t (float): Temps de la simulation, en jours (1 pas = 1 journée)

        Returns:
            np.array: Position du projectile à l'instant t.
        """
        try:
            etat = self.etat_a(int(t))
            return etat[0]  # Retourne seulement les coordonnées x et y de l'état complet à l'instant t
        except DebutHistorique as e:
            raise DebutHistorique(str(e))
        except FinHistorique as e:
            raise FinHistorique(str(e))

    def etat_a(self, pas: int) -> tuple[Vecteur2, Vecteur2]:
        """
        Args:
            pas (int): L'étape de la simulation pour laquelle on cherche l'état du projectile, en jours.

        Returns:
            Le vecteur de position et le vecteur de vitesse au pas spécifié.
        """
        if pas < self.temps_apparition:
            raise DebutHistorique("Le projectile n'a pas encore été créé à ce temps.")

        i_ligne = int(pas - self.temps_apparition)  # Nombre de pas depuis l'apparition

        if i_ligne >= len(self.historique):
            raise FinHistorique("Le projectile n'existe plus à ce temps.")

        ligne = self.historique[i_ligne]

        return Vecteur2(ligne[0], ligne[1]), Vecteur2(ligne[2], ligne[3])


class Vecteur2:
    """
    Un vecteur 2-dimensionnel qui supporte l'addition (+), la soustraction (-), la multiplication par un scalaire (*)
    et le produit scalaire (dotprod).
    """
    def __init__(self, x, y):
        self.__x = x
        self.__y = y

    @property
    def x(self):
        return self.__x

    @property
    def y(self):
        return self.__y

    def norme2(self) -> float:
        return float(self.__x * self.__x + self.__y * self.__y)

    def norme(self) -> float:
        return math.sqrt(self.norme2())

    def __add__(self, other: Vecteur2) -> Vecteur2:
        return Vecteur2(self.__x + other.__x, self.__y + other.__y)

    def __neg__(self) -> Vecteur2:
        return Vecteur2(-self.__x, -self.__y)

    def __sub__(self, other: Vecteur2) -> Vecteur2:
        return self + (-other)

    def __mul__(self, other: float) -> Vecteur2:
        return Vecteur2(other * self.__x, other * self.__y)

    def __rmul__(self, other: float) -> Vecteur2:
        return Vecteur2(other * self.__x, other * self.__y)

    def unpack(self) -> tuple[float, float]:
        return self.__x, self.__y

    def dotprod(self, other: Vecteur2) -> float:
        return self.__x * other.__x + self.__y * other.__y

    @staticmethod
    def nul() -> Vecteur2:
        return Vecteur2(0, 0)

    @staticmethod
    def pack_t(arg: tuple[float, float]) -> Vecteur2:
        return Vecteur2(arg[0], arg[1])

    @staticmethod
    def v_polaire(r: float, theta: float):
        return Vecteur2(r * np.cos(theta), r * np.sin(theta))
