global s_values, feigenbaum_ratio
s_values = [-2.944271909999, 2, 3.23606797749979]  # known values of s_0, s_1, s_2
feigenbaum_ratio = []


# To find the iterates and their derivatives
def iteration_function(mu, k):
    x = 0.5
    xd = 0
    for i in range(k):
        xd = x * (1 - x) + mu * (1 - 2 * x) * xd
        x = mu * x * (1 - x)
    return x, xd


def newton_method(k, j):
    mu = initial_guess(j)
    feigenbaum_ratio.append(delta(j - 1))
    for i in range(4):
        h = (iteration_function(mu, k)[0] - 0.5) / iteration_function(mu, k)[1]
        mu = mu - h
    return mu


# To find the initial guess based on s_n
def initial_guess(n):
    return s_values[n - 1] + (s_values[n - 1] - s_values[n - 2]) / delta(n - 1)


# To find the Feigenbaum ratios
def delta(n):
    if n == 2:
        return 4
    else:
        return (s_values[n - 1] - s_values[n - 2]) / (s_values[n] - s_values[n - 1])


# To update the list of s_n
for i in range(3, 20):
    s_values.append(newton_method(2 ** (i - 1), i))

print(" n        ", "s_n", "             Feigenbaum ratio")
for i in range(len(feigenbaum_ratio)):
    print('{:2d} {:.19f} {:.19f}'.format(i + 2, s_values[i + 3], feigenbaum_ratio[i]))
