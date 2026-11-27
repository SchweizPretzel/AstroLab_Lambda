from math import sqrt, cos, sin, radians, pi


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
        return sqrt(self.norme2())

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
        """
        Returns:
            Le vecteur nul
        """
        return Vecteur2(0, 0)

    @staticmethod
    def pack_t(arg: tuple[float, float]) -> Vecteur2:
        """
        Génère un vecteur à partir d'un tuple contenant ses composantes
        Args:
            arg: tuple contenant les composantes du vecteur

        Returns:
            Le vecteur résultant
        """
        return Vecteur2(arg[0], arg[1])

    @staticmethod
    def v_polaire(r: float, theta: float):
        """
        Génère un vecteur 2d avec rayon r et orientation theta
        Args:
            r: rayon
            theta: orientation; angle entre le vecteur et l'axe des x

        Returns:
            Le vecteur résultant
        """
        return Vecteur2(r * cos(theta), r * sin(theta))


class Vecteur3:
    """
    Un vecteur 3-dimensionnel qui supporte l'addition (+), la soustraction (-), la multiplication par un scalaire (*),
    le produit scalaire (dotprod) et le produit vectoriel (crossprod).
    """
    def __init__(self, x, y, z):
        self.__x = x
        self.__y = y
        self.__z = z

    @property
    def x(self):
        return self.__x

    @property
    def y(self):
        return self.__y

    @property
    def z(self):
        return self.__z

    def norme2(self) -> float:
        return float(self.__x * self.__x + self.__y * self.__y + self.__z * self.__z)

    def norme(self) -> float:
        return sqrt(self.norme2())

    def __add__(self, other: Vecteur3) -> Vecteur3:
        return Vecteur3(self.__x + other.__x, self.__y + other.__y, self.__z + other.__z)

    def __neg__(self) -> Vecteur3:
        return Vecteur3(-self.__x, -self.__y, -self.__z)

    def __sub__(self, other: Vecteur3) -> Vecteur3:
        return self + (-other)

    def __mul__(self, other: float) -> Vecteur3:
        return Vecteur3(other * self.__x, other * self.__y, other * self.__z)

    def __rmul__(self, other: float) -> Vecteur3:
        return Vecteur3(other * self.__x, other * self.__y, other * self.__z)

    def unpack(self) -> tuple[float, float, float]:
        return self.__x, self.__y, self.__z

    def dotprod(self, other: Vecteur3) -> float:
        return self.__x * other.__x + self.__y * other.__y + self.__z + other.__z

    @staticmethod
    def nul() -> Vecteur3:
        """
        Returns:
            Le vecteur nul
        """
        return Vecteur3(0, 0, 0)

    @staticmethod
    def pack_t(arg: tuple[float, float, float]) -> Vecteur3:
        """
        Génère un vecteur à partir d'un tuple contenant ses composantes
        Args:
            arg: tuple contenant les composantes du vecteur

        Returns:
            Le vecteur résultant
        """
        return Vecteur3(arg[0], arg[1], arg[2])

    @staticmethod
    def v_cyl(r: float, theta: float, h: float):
        """
        Génère un vecteur avec longueur horizontale r, orientation theta, et hauteur h
        Args:
            r: rayon de la projection du vecteur su le plan horizontal
            theta: orientation; angle entre la projection du vecteur sur le plan horizontal et l'axe des
            h: hauteur; composante verticale du vecteur

        Returns:
            Le vecteur résultant
        """
        return Vecteur3(r * cos(theta), h, r * sin(theta))

    @staticmethod
    def v_polaire(r: float, theta: float, phi: float):
        """
        Génère un vecteur 3d avec rayon r, orientation theta, et inclinaison phi
        Args:
            r: rayon
            theta: orientation; angle entre la projection du vecteur sur le plan horizontal et l'axe des x
            phi: inclinaison; angle entre le vecteur et le plan horizontal

        Returns:
            Le vecteur résultant
        """
        return Vecteur3(r * cos(theta) * cos(phi), r * sin(theta) * cos(phi), r * sin(phi))


class Ellipse3:
    def __init__(self, orientation: float, inclinaison: float, a: float, b: float, angles="rad"):
        if a == 0 or b == 0:
            raise ValueError("Une ellipse ne peut avoir une distance nulle")

        if angles == "deg":
            self.__orientation = radians(orientation)
            self.__inclinaison = radians(inclinaison)
        else:
            self.__orientation = orientation
            self.__inclinaison = inclinaison

        self.__a = Vecteur3.v_polaire(a, orientation, inclinaison)
        self.__b = Vecteur3.v_polaire(b, orientation + pi/2, 0)

    def __call__(self, theta: float):
        return
