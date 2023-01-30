from algebra import arrange, Vector


# Метод эйлера
def euler(f, y0, t_start, t_max, h, args=()):
    t = arrange(t_start, t_max, h)
    y = [y0]
    for i in range(1, len(t)):
        y.append(y[-1] + h * f(t[i], y[-1], *args))
    return t, y


# Классический РК4
def rk4(f, y0, t_start, t_max, h, args=()):
    t = arrange(t_start, t_max, h)
    y = [y0]
    for i in range(1, len(t)):
        prev = y[-1]
        k1 = f(t[i], prev, *args)
        k2 = f(t[i] + h / 2, prev + (h / 2) * k1, *args)
        k3 = f(t[i] + h / 2, prev + (h / 2) * k2, *args)
        k4 = f(t[i] + h, prev + h * k3, *args)
        y.append(prev + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4))
    return t, y


# Адамса-Башфорта 3-шаговый
def adams_bashforth3(f, y0, t_start, t_max, h, args=()):
    t = arrange(t_start, t_max, h)
    # Разгон
    y = rk4(f, y0, t_start, t_start + h * 3, h, args=args)[1][:]

    for i in range(3, len(t)):
        k3 = f(t[i - 2], y[-3], *args)
        k2 = f(t[i - 1], y[-2], *args)
        k1 = f(t[i], y[-1], *args)
        y.append(y[-1] + h * (23 * k1 - 16 * k2 + 5 * k3) / 12)
    return t, y


# Рунге-Кутты-Фельберга с подбором шага
def rkf(f, y0, t_start, t_max, h, args=()):
    # https://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
    def rkf_step(y, t, step):
        tolerance = 1e-5
        k1 = f(t, y, *args)
        k2 = f(t, y + step * k1 / 4, *args)
        k3 = f(t, y + k1 * step * (3 / 32) + k2 * step * 9 / 32, *args)
        k4 = f(t, y + k1 * step * (1932 / 2197) + k2 * (-7200) * step / 2197 + k3 * 7296 * step / 2197, *args)
        k5 = f(t, y + k1 * 439 * step / 216 + k2 * (-8) * step / 1 + k3 * 3680 * step / 513 + k4 * (-845) * step / 4104,
               *args)
        k6 = f(
            t,
            y + k1 * step * (-8) / 27 + k2 * step * 2 + k3 * step * (-3544) / 2565 + k4 * step * 1859 / 4104 +
            k5 * step * (-11) / 40,
            *args
        )

        y1 = y + k1 * step * 16 / 135 + k3 * step * 6656 / 12825 + k4 * step * 28561 / 56430 + k5 * step * (
            -9) / 50 + k6 * step * 2 / 55
        y2 = y + k1 * step * 25 / 216 + k3 * step * 1408 / 2565 + k4 * step * 2197 / 4104 + k5 * step * (-1) / 5
        epsilon = (y1 - y2).abs().norm()
        if epsilon > tolerance:
            step *= (tolerance / (2 * epsilon)) ** 0.25
            return rkf_step(y, t, step)
        return y2, t + step

    t = t_start
    tlist = [t]
    answer = [y0]
    while t < t_max:
        step = rkf_step(answer[-1], t, h)
        answer.append(step[0])
        t = step[1]
        tlist.append(t)
    return tlist, answer


# Дормана-Принса 5(4) с подбором шага
def rkdp(f, y0, t_start, t_max, h, args=()):
    a21 = 1 / 5
    a31, a32 = 3 / 40, 9 / 40
    a41, a42, a43 = 44 / 55, -56 / 15, 32 / 9
    a51, a52, a53, a54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
    a61, a62, a63, a64, a65 = 9017 / 3186, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656

    c2, c3, c4, c5, c6 = 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1

    b1 = (35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84)
    b2 = (5197 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40)

    def rkdp_step(y, t, step):
        tolerance = 1e-4
        k1 = f(tlist[-1], y, *args)
        k2 = f(tlist[-1] + c2 * step, y + step * a21 * k1, *args)
        k3 = f(tlist[-1] + c3 * step, y + step * (a31 * k1 + a32 * k2), *args)
        k4 = f(tlist[-1] + c4 * step, y + step * (a41 * k1 + a42 * k2 + a43 * k3), *args)
        k5 = f(tlist[-1] + c5 * step, y + step * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), *args)
        k6 = f(tlist[-1] + c6 * step, y + step * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5), *args)
        k = [k1, k2, k3, k4, k5, k6]
        fsal = Vector([0 for _ in range(len(k1))])
        for i in range(6):
            fsal += k[i] * b1[i] * step
        k7 = f(tlist[-1] + step, answer[-1] + fsal, *args)
        k.append(k7)
        y1 = y + fsal
        y2 = y
        for i in range(7):
            y2 += k[i] * b2[i] * step
        epsilon = (y1 - y2).abs().norm()
        if epsilon > tolerance:
            step *= (tolerance / (2 * epsilon)) ** 0.25
            return rkdp_step(y, t, step)
        return y2, t + step


    t = t_start
    tlist = [t]
    answer = [y0]
    while t < t_max:
        step = rkdp_step(answer[-1], t, h)
        answer.append(step[0])
        t = step[1]
        tlist.append(t)
    return tlist, answer


summary = [
    (euler, 'Метод Эйлера'),
    (rk4, 'Классический метод Рунге-Кутты 4 порядка'),
    (adams_bashforth3, '3-шаговый метод Адамса-Башфорта'),
    (rkf, 'Метод Рунге-Кутты Фельберга с подбором шага'),
    (rkdp, 'Метод Дормана-Принса 5(4)7 с подбором шага')
]
