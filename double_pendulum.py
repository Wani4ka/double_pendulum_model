"""Model of a double pendulum"""
# Author: Chris Greening
# Date: 2019-07-15
from math import radians
from copy import copy
from typing import Tuple, List

from scipy import constants

import odeint
from pendulum import Pendulum
from equations import derivative
from algebra import Matrix, Vector

class DoublePendulum:
    tmax = 30.0
    dt = 0.05

    def __init__(self, L1: float = 1, L2: float = 1, m1: float = 1, m2: float = 1,
                 y0: List[float] = None, g: float = -1 * constants.g, odeint_method = odeint.rk4) -> None:
        """Создает двойной маятник с заданными параметрами и начальными условиями

        Параметры
        ----------
        L1 : float = 1
            Длина стержня первого маятника
        L2 : float = 1
            Длина стержня второго маятника
        m1 : float = 1
            Масса шара первого маятника
        m2 : float = 1
            Масса шара второго маятника
        y0 : List[float] = [90, 0, 90, 0]
            Начальный угол и угловая скорость маятников (в градусах): [theta1, w1, theta2, w2]
        g : float = scipy.constants.g
            Ускорение свободного падения
        odeint_method : function = odeint.rk4
            Метод интегрирования
        """
        if y0 is None:
            y0 = [90, 0, 90, 0]
        for i in range(len(y0)):
            y0[i] = radians(y0[i])
        self.pendulum1 = Pendulum(L1, m1)
        self.pendulum2 = Pendulum(L2, m2)
        self.y0 = Vector(y0)
        self.g = g
        self.odeint_method = odeint_method

        # Интегрируем
        self._calculate_system()

        self.max_length = self.pendulum1.L + self.pendulum2.L

    def get_frame_x(self, i: int) -> Tuple[int]:
        return (0, self.pendulum1.x[i], self.pendulum2.x[i])

    def get_frame_y(self, i: int) -> Tuple[int]:
        return (0, self.pendulum1.y[i], self.pendulum2.y[i])

    def get_frame_coordinates(self, i: int) -> Tuple[Tuple[int]]:
        return (self.get_frame_x(i), self.get_frame_y(i))

    def get_max_x(self) -> float:
        return self.pendulum2.get_max_x()

    def get_max_y(self) -> float:
        return self.pendulum2.get_max_y()

    def get_max_coordinates(self) -> float:
        return self.pendulum2.get_max_coordinates()

    def _calculate_system(self) -> None:
        """Решение ОДУ и вычисление траекторий для обоих маятников"""
        args = (self.pendulum1.L, self.pendulum2.L, self.pendulum1.m, self.pendulum2.m, self.g)
        self.t, self.y = self.odeint_method(derivative, self.y0, 0, self.tmax, self.dt, args=args)
        y_matrix = Matrix(self.y)

        # Вычисление траектории первого маятника
        self.pendulum1.calculate_path(
            theta=y_matrix.column(0),
            dtheta=y_matrix.column(1)
        )
        # Вычисление траектории второго маятника
        self.pendulum2.calculate_path(
            theta=y_matrix.column(2),
            dtheta=y_matrix.column(3),
            x0=self.pendulum1.x,
            y0=self.pendulum1.y
        )

        self.w = y_matrix.column(1)

    @classmethod
    def create_multiple_double_pendula(
            cls, num_pendula: int = 1, L1: float = 1.0,
            L2: float = 1.0, m1: float = 1.0, m2: float = 1.0,
            y0: List[float] = None, dtheta: float = .05, odeint_method = odeint.rk4) -> List["DoublePendulum"]:
        """Создает несколько двойных маятников, каждый из которых на dtheta градусов ниже предыдущего"""
        pendula = []
        if y0 is None:
            y0 = [90, 0, 90, 0]
        for _ in range(num_pendula):
            double_pendulum = cls(
                L1=L1,
                L2=L2,
                m1=m1,
                m2=m2,
                y0=copy(y0),
                odeint_method=odeint_method
            )
            pendula.append(double_pendulum)
            y0[0] += dtheta
        return pendula
