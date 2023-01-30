from typing import Tuple
from algebra import Vector

from math import sin, cos

class Pendulum:
    def __init__(self, L: float = 1.0, m: float = 1.0) -> None:
        """
        :param L: длина стержня
        :param m: масса шара
        """
        self.L = L
        self.m = m

    def calculate_path(self, theta: float, dtheta: float, x0: float = 0, y0: float = 0) -> None:
        """Вычисляет x,y координаты маятника"""
        self.theta = theta
        self.dtheta = dtheta
        self.x = self.L * Vector(map(sin, theta)) + x0
        self.y = self.L * Vector(map(cos, theta)) + y0

    def get_max_x(self) -> float:
        """Максимальное значение x координаты, которой достигает этот маятник"""
        return max(self.x)

    def get_max_y(self) -> float:
        """Максимальное значение y координаты, которой достигает этот маятник"""
        return max(self.y)

    def get_max_coordinates(self) -> Tuple[float, float]:
        """Максимальные координаты, которых достигает эта система"""
        return (self.get_max_x(), self.get_max_y())
