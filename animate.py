import random
from typing import List

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from colorsys import hls_to_rgb
from random import random

from double_pendulum import DoublePendulum

from odeint import summary

colors = {}
dt = 0.05


def rand_color() -> str:
    h,s,l = random(), 0.5 + random()/2.0, 0.4 + random()/5.0
    rgb = [int(256*i) for i in hls_to_rgb(h,l,s)]
    return '#' + ''.join("%02X" % i for i in rgb)

def animate(i, pendula_axes):
    """Annimation frame"""
    time_template = 'time = %.1fs'
    return_arr = []
    for double_pendulum, ax_data in pendula_axes:
        _, line, time_text = ax_data
        idx = i
        frame_x, frame_y = double_pendulum.get_frame_coordinates(idx)
        line.set_data(frame_x, frame_y)
        if double_pendulum.t[idx] >= 30 - dt:
            plt.close()
            return []
        time_text.set_text(time_template % (double_pendulum.t[idx]))
        return_arr.extend([
            line,
            time_text,
        ])
    return return_arr


def create_axes(
        fig: "matplotlib.figure.Figure",
        pendula: List["DoublePendulum"]
):
    axes = []
    longest_double_pendulum = max(pendula, key=lambda x: x.max_length)
    for i in range(len(pendula)):
        color = rand_color()
        ax = _create_individual_axis(
            longest_double_pendulum=longest_double_pendulum,
            fig=fig
        )
        line, = ax.plot([], [], 'o-', lw=2, color=color)
        colors[color] = pendula[i]
        time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
        axes.append((ax, line, time_text))
    return axes


def _create_individual_axis(
        longest_double_pendulum: "DoublePendulum",
        fig: matplotlib.figure.Figure,
):
    ax = fig.add_subplot(
        111,
        autoscale_on=False,
        xlim=(
            -longest_double_pendulum.max_length,
            longest_double_pendulum.max_length
        ),
        ylim=(
            -longest_double_pendulum.max_length,
            longest_double_pendulum.max_length
        ),
    )
    ax.set_aspect('equal')
    ax.grid()
    return ax


def make_figure(odeint_method, method_name):
    """Создает маятники"""
    fig = plt.figure()
    fig.suptitle(method_name)
    pendula = DoublePendulum.create_multiple_double_pendula(num_pendula=5, y0=[90, 0, 90, 0], L2=2, m1=3, odeint_method=odeint_method)
    # pendula = DoublePendulum.create_multiple_double_pendula(num_pendula=10, y0=[100, 0, 90, 0], L1=2, m1=3)
    # pendula = DoublePendulum.create_multiple_double_pendula(num_pendula=10, y0=[0, 0, 90, 0], L1=2, m1=3)

    axes = create_axes(fig=fig, pendula=pendula)
    pendula_axes = list(zip(pendula, axes))
    animation.FuncAnimation(
        fig,
        lambda i: animate(i, pendula_axes),
        list(range(1, len(pendula[0].t))),
        interval=25,
        blit=True,
    )
    # mngr = plt.get_current_fig_manager()
    # mngr.window.wm_geometry("+100+100")
    plt.show()


if __name__ == "__main__":
	for method in summary:
		make_figure(*method)
