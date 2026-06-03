import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import scipy.special as special
from scipy.integrate import solve_ivp
import pandas as pd

x = np.linspace(0.01, 10, 100)
j0 = special.jv(0, x)
j1 = special.jv(1, x)
j2 = special.jv(2, x)
j3 = special.jv(3, x)
j4 = special.jv(4, x)
j5 = special.jv(5, x)
j6 = special.jv(6, x)

pd.set_option("display.max_rows", 100)

def bessel_function(x, y, u):
    y, v = y
    dydx = v, -v/x - (x**2-u**2)*y/x**2
    return dydx

def initial_ju(u, x):
    k = np.arange(0, 100)
    return np.sum(((-1)**k / (special.factorial(k) * special.factorial(k + u))) * ((x/2)**(u+2*k))) 

def initial_dju_dx(u):
    k = np.arange(0, 100)
    X = sp.symbols('X')
    ju = np.sum(((-1)**k / (special.factorial(k) * special.factorial(k + u))) * ((X/2)**(u+2*k)))
    df = sp.diff(ju, X)
    return df.subs(X, 0.01)


u = 0
x_values = np.linspace(0.01, 10, 100)
sol = solve_ivp(bessel_function, (0.01, 10), [initial_ju(u, 0.01), initial_dju_dx(u)], args=(u,), t_eval=x_values)
pd.compare=pd.DataFrame({
    "x": x_values,
    "Built-in j0": j0,
    "Computed j0": sol.y[0]
})
print(pd.compare)
plt.plot(x, j0, label="Built-in j0", linestyle="dashed", color="skyblue")
plt.plot(x_values, sol.y[0], label="Computed j0",color="blue")
plt.title("Bessel Function of the First Kind (j0)")
plt.xlabel("x")
plt.ylabel("j0(x)")
plt.legend()
plt.grid()
plt.show()

print("---"*30)

u = 1
sol = solve_ivp(bessel_function, (0.01, 10), [initial_ju(u, 0.01), initial_dju_dx(u)], args=(u,), t_eval=x_values)
pd.compare1=pd.DataFrame({
    "x": x_values,
    "Built-in j1": j1,
    "Computed j1": sol.y[0]
})
print(pd.compare1)
plt.plot(x, j1, label="Built-in j1", linestyle="dashed", color="lightcoral")
plt.plot(x_values, sol.y[0], label="Computed j1",color="red")
plt.title("Bessel Function of the First Kind (j1)")
plt.xlabel("x")
plt.ylabel("j1(x)")
plt.legend()
plt.grid()
plt.show()

print("---"*30)

u = 2
sol = solve_ivp(bessel_function, (0.01, 10), [initial_ju(u, 0.01), initial_dju_dx(u)], args=(u,), t_eval=x_values)
pd.compare2=pd.DataFrame({
    "x": x_values,
    "Built-in j2": j2,
    "Computed j2": sol.y[0]
})
print(pd.compare2)
plt.plot(x, j2, label="Built-in j2", linestyle="dashed", color="lightgreen")
plt.plot(x_values, sol.y[0], label="Computed j2", color="green")
plt.title("Bessel Function of the First Kind (j2)")
plt.xlabel("x")
plt.ylabel("j2(x)")
plt.legend()
plt.grid()
plt.show()

print("---"*30)

u = 3
sol = solve_ivp(bessel_function, (0.01, 10), [initial_ju(u, 0.01), initial_dju_dx(u)], args=(u,), t_eval=x_values)
pd.compare3=pd.DataFrame({
    "x": x_values,
    "Built-in j3": j3,
    "Computed j3": sol.y[0]
})
print(pd.compare3)
plt.plot(x, j3, label="Built-in j3", linestyle="dashed", color="lightyellow")
plt.plot(x_values, sol.y[0], label="Computed j3", color="orange")
plt.title("Bessel Function of the First Kind (j3)")
plt.xlabel("x")
plt.ylabel("j3(x)")
plt.legend()
plt.grid()
plt.show()

print("---"*30)

u = 4
sol = solve_ivp(bessel_function, (0.01, 10), [initial_ju(u, 0.01), initial_dju_dx(u)], args=(u,), t_eval=x_values)
pd.compare4=pd.DataFrame({
    "x": x_values,
    "Built-in j4": j4,
    "Computed j4": sol.y[0]
})
print(pd.compare4)
plt.plot(x, j4, label="Built-in j4", linestyle="dashed", color="lightgray")
plt.plot(x_values, sol.y[0], label="Computed j4", color="gray")
plt.title("Bessel Function of the First Kind (j4)")
plt.xlabel("x")
plt.ylabel("j4(x)")
plt.legend()
plt.grid()
plt.show()

print("---"*30)

u = 5
sol = solve_ivp(bessel_function, (0.01, 10), [initial_ju(u, 0.01), initial_dju_dx(u)], args=(u,), t_eval=x_values)
pd.compare5=pd.DataFrame({
    "x": x_values,
    "Built-in j5": j5,
    "Computed j5": sol.y[0]
})
print(pd.compare5)
plt.plot(x, j5, label="Built-in j5", linestyle="dashed", color="lightpink")
plt.plot(x_values, sol.y[0], label="Computed j5", color="magenta")
plt.title("Bessel Function of the First Kind (j5)")
plt.xlabel("x")
plt.ylabel("j5(x)")
plt.legend()
plt.grid()
plt.show()

print("---"*30)

u = 6
sol = solve_ivp(bessel_function, (0.01, 10), [initial_ju(u, 0.01), initial_dju_dx(u)], args=(u,), t_eval=x_values)
pd.compare6=pd.DataFrame({
    "x": x_values,
    "Built-in j6": j6,
    "Computed j6": sol.y[0]
})
print(pd.compare6)
plt.plot(x, j6, label="Built-in j6", linestyle="dashed", color="lightseagreen")
plt.plot(x_values, sol.y[0], label="Computed j6", color="teal")
plt.title("Bessel Function of the First Kind (j6)")
plt.xlabel("x")
plt.ylabel("j6(x)")
plt.legend()
plt.grid()
plt.show()








