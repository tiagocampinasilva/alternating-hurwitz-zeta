"""Algebraic checks for equations (4), (5), (7)--(9), (12) and (14)."""
import sympy as sp


def residuals():
    i, pi = sp.I, sp.pi
    x, y, h, r, u = sp.symbols("x y h r u", real=True)
    z, f = sp.symbols("z f")
    a, w, v = sp.symbols("a w v", nonzero=True)
    b = sp.symbols("b", positive=True)
    half = sp.Rational(1, 2)

    yield "Abel-Plana, logistic weight", (1 + sp.tanh(b/2)).rewrite(sp.exp)/2 - 1/(1+sp.exp(-b))
    kernel = (a*w - v)/(1-a)
    yield "(8)-(9), splitting, a=exp(-2*pi*y)", (w-v)/(1-a) - w - kernel
    yield "(8), exponential form, exp(D)=v/(a*w)", kernel - w*(1-v/(a*w))/(1/a-1)

    theta = sp.atan(y*h**2)
    height = pi/h**2 + x/h
    phase_difference = 2*pi*y - 2*height*theta + i*theta
    yield "(8), continuous value at zero", -sp.diff(phase_difference, y).subs(y, 0)/(2*pi) - x*h/pi + i*h**2/(2*pi)
    power_difference = (1+i*y*h**2)**(-half-i*height) - (1-i*y*h**2)**(-half-i*height)
    yield "(7), quotient at the origin", sp.diff(power_difference, y).subs(y, 0)/(2*pi) + i*(half+i*height)*h**2/pi

    psi = -pi*z/h - (half+i*height)*sp.log(1+i*h*z)
    phi = x*z - i*pi*z**2/2
    p = -i*z/2 - i*x*z**2/2 - pi*z**3/3
    taylor = sp.series(psi, h, 0, 2).removeO().expand()
    yield "Lemma 5.1, constant phase coefficient", taylor.coeff(h, 0) - phi
    yield "Lemma 5.1, linear phase coefficient", taylor.coeff(h, 1) - p

    angle = sp.atan(r/(sp.sqrt(2)+r))
    delta = r/sp.sqrt(2) - angle
    denominator = 1 + sp.sqrt(2)*r + r**2
    yield "Lemma 5.2, derivative of the rotated angle", sp.diff(angle, r) - 1/(sp.sqrt(2)*denominator)
    yield "Lemma 5.2, delta prime divided by r", sp.diff(delta, r)/r - (1+r/sp.sqrt(2))/denominator
    yield "Lemma 5.2, derivative of delta prime divided by r", sp.diff(sp.diff(delta, r)/r, r) + (1/sp.sqrt(2)+2*r+r**2/sp.sqrt(2))/denominator**2

    profile = (1+i)/2 * sp.exp(-i*x**2/(2*pi)) * sp.erfc(-(1-i)*x/(2*sp.sqrt(pi)))
    yield "(4), F(0)", profile.subs(x, 0) - (1+i)/2
    yield "(4), differential identity", pi*sp.diff(profile, x) - 1 + i*x*profile
    f1 = (1-i*x*f)/pi

    def differentiate(expr):
        return sp.diff(expr, x) + f1*sp.diff(expr, f)

    f2 = differentiate(f1)
    f3 = differentiate(f2)
    yield "Second derivative of F", pi*f2 + i*f + i*x*f1
    yield "Third derivative of F", pi*f3 + 2*i*f1 + i*x*f2
    j1 = -i*f1/2 - i*x*f2/2 - pi*f3/3
    yield "(5), integration by parts", j1 - i*(f1-x*f2)/6
    yield "(5), G minus its Gaussian-moment expression", half+i/(6*pi)-x**2*(1-i*x*f)/(6*pi**2) - (half+j1)

    upper_pi = sp.Rational(22, 7)
    majorant = (u/2 + u**2/2 + upper_pi*u**3/3)**2 + u**2/2 + 2*u**3/3 + upper_pi*u**4/2
    printed = 3*u**2/4 + 7*u**3/6 + 241*u**4/84 + 22*u**5/21 + 484*u**6/441
    yield "(14), majorant polynomial", majorant - printed

    for j, full, odd in ((1, sp.Rational(1,24), sp.Rational(1,8)),
                         (3, sp.Rational(1,240), sp.Rational(1,16))):
        standard = sp.factorial(j)*sp.zeta(j+1)/(2*pi)**(j+1)
        odd_moment = sp.factorial(j)*(1-sp.Integer(2)**(-j-1))*sp.zeta(j+1)/pi**(j+1)
        yield f"Lemma 4.1, moment {j}", standard - full
        yield f"Lemma 4.2, odd-mode moment {j}", odd_moment - odd


def main():
    total = failed = 0
    for name, expression in residuals():
        total += 1
        try:
            difference = sp.simplify(expression)
        except (ValueError, TypeError, NotImplementedError) as error:
            print(f"UNRESOLVED  {name}: {error}")
            failed += 1
            continue
        if difference == 0:
            print(f"{name}: 0")
        else:
            print(f"UNCONFIRMED  {name}: {difference}")
            failed += 1
    print(f"{total-failed}/{total} symbolic residuals reduced to zero.")
    return int(failed != 0)


if __name__ == "__main__":
    raise SystemExit(main())
