"""Exact rational checks for Lemmas 4.1--5.2, equations (14)--(17) and Appendix 1."""
from fractions import Fraction as Q
from math import factorial


def equal(label, observed, expected):
    if observed != expected:
        raise ArithmeticError(f"{label}: obtained {observed}, expected {expected}")
    print(f"{label} = {observed}")


def positive(label, margin):
    if margin <= 0:
        raise ArithmeticError(f"{label}: nonpositive margin {margin}")
    print(f"{label}: margin {margin}")


def exp_bracket(z, degree):
    lower = sum(z**k/factorial(k) for k in range(degree+1))
    upper = lower + z**(degree+1)/factorial(degree+1)/(1-z/(degree+2))
    return lower, upper


def transition_margins(label, d):
    r0, root_lower = Q(3,10), Q(140,99)
    positive(label + " Gaussian", Q(157,50)*d/r0**2 - Q(9,8) - Q(5,22))
    positive(label + " exponential", 10*Q(157,50)*d/r0 - 1/root_lower - Q(67,20))


def transition_checks():
    r0, lower, upper = Q(3,10), Q(140,99), Q(99,70)
    positive("sqrt(2), lower square", 2-lower**2)
    positive("sqrt(2), upper square", upper**2-2)
    v = r0/(lower+r0)
    transition_margins("Arctangent", r0/upper - v + v**3/3 - v**5/5)
    lower, upper = Q(1414213,10**6), Q(1414214,10**6)
    positive("Appendix lower square", 2-lower**2)
    positive("Appendix upper square", upper**2-2)
    step = r0/4096
    total = 0
    for j in range(4096):
        r = j*step
        b = lower*r+r*r
        total += (10**12*b) // (upper*(1+b))
    d = step*total/10**12
    equal("Appendix lower sum", d, Q(1592131749274149,40960000000000000))
    transition_margins("Lower sum", d)


def main():
    pi_upper, a, kappa = Q(22,7), Q(9,8), Q(67,20)
    h = Q(1,10)
    transition_checks()
    equal("Completed-square constant", Q(1,2)/(4*Q(5,22)), Q(11,20))
    positive("exp(11/20) < 7/4", Q(7,4)-exp_bracket(Q(11,20), 6)[1])
    cutoff = Q(461,100)
    positive("exp(4.61) > 100", exp_bracket(cutoff, 16)[0]-100)
    positive("exp(0.461) < 5/3", Q(5,3)-exp_bracket(cutoff*h, 6)[1])
    positive("Local multiplier < 3", 3-Q(5,3)*(exp_bracket(Q(1), 6)[1]-1))
    phase_bound = (2*h+h*h)*cutoff + 2*pi_upper*cutoff**3*h**4/3
    equal("Short-interval phase", phase_bound, Q(51902943991,52500000000))
    positive("Phase < 0.989", Q(989,1000)-phase_bound)
    pi_lower = Q(157,50)
    positive("Tail prefactor 3/pi < 1", 1-3/pi_lower)
    positive("Tail exponent below -5/2", pi_lower-3)
    equal("Tail comparison at q=100, squared", Q(100)**(-5), Q(1,10**5)**2)
    positive("Tail denominator exceeds 2/3", Q(1,3)-Q(100)**(-6))
    positive("2*pi-h > pi", pi_lower-h)

    # Appendix 1 counterexample: A=100*pi*(1-pi/4), |D| <= 2*A+pi/4.
    positive("Counterexample A > 64", 100*pi_lower*(1-pi_upper/4)-64)
    positive("Counterexample 5|D| < 3150", 3150-5*(200*pi_upper+pi_upper/4))
    positive("Counterexample numerator > 3150", (Q(64)**3/6-1)/2-3150)

    # Square the first derivative bound and add the second derivative bound.
    first = (Q(0), Q(1,2), Q(1,2), pi_upper/3)
    p = [Q(0), Q(0), Q(1,2), Q(2,3), pi_upper/2, Q(0), Q(0)]
    for i, left in enumerate(first):
        for j, right in enumerate(first):
            p[i+j] += left*right
    equal("P_b coefficients (degrees 2--6)", tuple(p[2:]),
          (Q(3,4), Q(7,6), Q(241,84), Q(22,21), Q(484,441)))

    positive("Squared g0 bound", Q(67,80)**2-pi_upper/(4*a))
    gaussian = [Q(67,80), 1/(2*a)]  # Even entries propagate the bound on g0.
    for j in range(2, 7):
        gaussian.append((j-1)*gaussian[j-2]/(2*a))
    exponential = [Q(factorial(j))/kappa**(j+1) for j in range(7)]
    for j in range(2,7):
        print(f"j={j}: Gaussian bound {gaussian[j]}, exponential moment {exponential[j]}")
    jbar = sum(p[j]*(Q(7,4)*gaussian[j]+exponential[j]) for j in range(2,7))
    equal("Jbar", jbar, Q(533948559110614545113,69278297361197722560))

    # The two full-kernel moments are evaluated symbolically in the companion script.
    m1, m3 = Q(1,24), Q(1,240)
    odd1, odd3 = 3*m1, 15*m3
    equal("Odd-mode moments", (odd1,odd3), (Q(1,8),Q(1,16)))
    local = (2+h)*m1 + 2*pi_upper*h**3*m3/3
    tail = Q(1,10**5)
    split, b3 = 5*local+tail, 3*local+tail
    whole = (2+h)*odd1 + 2*pi_upper*h**3*odd3/3
    equal("B_s", split, Q(689147,1575000))
    equal("B_g", whole, Q(22061,84000))
    equal("B_3", b3, Q(3,12)+Q(3,240)+Q(22,7*120000)+tail)
    main_sum, whole_sum, appendix_sum = (jbar/2+b for b in (split,whole,b3))
    equal("Equation (16)", main_sum, Q(371609064358084350147257,86597871701497153200000))
    equal("Equation (17)", whole_sum, Q(14258445238530472531037,3463914868059886128000))
    equal("Appendix 1", appendix_sum, Q(356452924783991052355217,86597871701497153200000))
    positive("43/10 minus (16)", Q(43,10)-main_sum)
    positive("4117/1000 minus (17)", Q(4117,1000)-whole_sum)
    positive("4117/1000 minus Appendix 1", Q(4117,1000)-appendix_sum)


if __name__ == "__main__":
    main()
