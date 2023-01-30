from math import sin, cos
from algebra import Vector

# первая производная от y = theta1, z1, theta2, z2
def derivative(t, y, L1, L2, m1, m2, g):
    theta1, z1, theta2, z2 = y
    theta1dot, theta2dot = z1, z2

    c, s = cos(theta1 - theta2), sin(theta1 - theta2)

    z1dot = (
                    m2 * g * sin(theta2) * c -
                    m2 * s * (L1 * z1 ** 2 * c + L2 * z2 ** 2) -
                    (m1 + m2) * g * sin(theta1)
            ) / (L1 * (m1 + m2 * s ** 2))
    z2dot = (
                    (m1 + m2) *
                    (L1 * z1 ** 2 * s - g * sin(theta2) + g * sin(theta1) * c) +
                    m2 * L2 * z2 ** 2 * s * c
            ) / (L2 * (m1 + m2 * s ** 2))

    return Vector([theta1dot, z1dot, theta2dot, z2dot])
