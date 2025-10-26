import math

x = input("Degree of the taylor polynom: ")
degree = int(x)
a = math.pi / 4
coeffs = []

# sin and cos derivatives cycle every 4 terms:
# 0: sin
# 1: cos
# 2: -sin
# 3: -cos
# 4: sin ...
for n in range(degree):
    if n % 4 == 0:
        deriv = math.sin(a)
    elif n % 4 == 1:
        deriv = math.cos(a)
    elif n % 4 == 2:
        deriv = -math.sin(a)
    else:  # n % 4 == 3
        deriv = -math.cos(a)
    
    coeff = deriv / math.factorial(n)
    coeffs.append(coeff)

for c in coeffs:
    print(f"{c:.17g}")
