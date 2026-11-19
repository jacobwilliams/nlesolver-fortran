#!/usr/bin/env python3
"""Data for the 'minimum-norm step vs. closest solution' figure of the technical note.

One equation in two unknowns, F(x, y) = y - x^2, starting from x0 = (1, 0):
  * plain:    minimum-norm Newton steps,  p = -J^+ F
  * anchored: steps that minimize ||x_k + p - x0|| subject to J p = -F
  * closest:  the point of the parabola closest to x0

Usage (from the repository root): python3 docs/scripts/parabola.py
"""
import math
import os

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')


def F(x, y):
    return y - x * x


def J(x, y):
    return (-2.0 * x, 1.0)


X0 = (1.0, 0.0)


def iterate(anchored, iters=12):
    x, y = X0
    pts = [(x, y)]
    for _ in range(iters):
        a, b = J(x, y)
        jj = a * a + b * b
        if anchored:
            dx, dy = x - X0[0], y - X0[1]
            r = -F(x, y) + (a * dx + b * dy)       # J q = -F + J d,  q = p + d
            q = (r * a / jj, r * b / jj)           # minimum-norm q
            p = (q[0] - dx, q[1] - dy)
        else:
            r = -F(x, y)
            p = (r * a / jj, r * b / jj)           # minimum-norm p
        x, y = x + p[0], y + p[1]
        pts.append((x, y))
    return pts


def main():
    os.makedirs(DATA, exist_ok=True)
    for name, anchored in [('plain', False), ('anchored', True)]:
        with open(os.path.join(DATA, f'parabola_{name}.csv'), 'w') as f:
            f.write('x,y\n')
            for p in iterate(anchored):
                f.write('%.10f,%.10f\n' % p)
    # closest point: minimize (x-1)^2 + x^4  ->  2(x-1) + 4x^3 = 0
    x = 0.6
    for _ in range(50):
        x -= (2 * (x - 1) + 4 * x ** 3) / (2 + 12 * x * x)
    with open(os.path.join(DATA, 'parabola_closest.csv'), 'w') as f:
        f.write('x,y\n%.10f,%.10f\n' % (x, x * x))
    p = iterate(False)[-1]
    print('minimum-norm Newton: (%.6f, %.6f), distance %.6f' % (p[0], p[1], math.hypot(p[0] - 1, p[1])))
    print('closest point:       (%.6f, %.6f), distance %.6f' % (x, x * x, math.hypot(x - 1, x * x)))


if __name__ == '__main__':
    main()
