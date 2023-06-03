import numpy as np
import sympy as sp
import re
import math

def main():
    global coefficients
    choice = eval(input("\n(1) Root Find \n(2) Solving Ordinary Differential Equation(O.D.E) \nChoose A Solution For The Appriopriate Problem to be solve :"))
    while choice not in [1, 2]:
        print("Guess should be btw 1 and 2")
        choice = eval(input("\n(1) Root Find \n(2) Solving Ordinary Differential Equation(O.D.E) \nChoose A Solution For The Appriopriate Problem to be solve :"))
    else:
        if choice == 1:
            coefficients = []
            while not coefficients:
                coefficients = input("Enter your function f(x) in terms of coefficient of powers X's seperated by a space in form(a1, a2, a3, ....an): ").strip().split()
            else:
                initial_intervals =initialApproximation()
                method_answer = eval(input("\n\nRoot Finding\n(1)Bisecton \n(2)Newton-Raphson \n(3)Secant or LineraInterpolaton \nSelect Prefferable method's to solving problem numerical: "))
                while not  1<=method_answer<=3:
                    method_answer = input("\nRoot Finding\n(1)Bisecton \n(2)Newton-Raphson \n(3)Secant or LineraInterpolaton \nSelect Prefferable method's to solving problem numerical: ")
                else:
                    print("\nThe Root lies in the interval(s): ", initial_intervals, end=" ")
                    if method_answer == 1:
                        print("The root(s) is/are: ", end="")
                        [print(myBisectionMethod((x0, x1)), end="  ") for x0, x1 in initial_intervals]
                    elif method_answer == 2:
                        print("The root(s) is/are: ", end="")
                        [print(myNewtonRaphson((x0, x1)), end="  ") for x0, x1 in initial_intervals]
                    else:
                        print("The root(s) is/are: ", end="")
                        [print(myLinearInterpolation(x0, x1), end="  ") for x0, x1 in initial_intervals]
        else:
            method_of_solution = int(input("\nSolution To O.D.E \n(1)Taylor series (2)Picard's theorem (3)Euler's theorem \n Select Prefferable method's to solving problem numerical:"))
            while method_of_solution not in range(1, 4):
                method_of_solution = input("\n(1)Taylor series (2)Picard's theorem (3)Euler's theorem \n Select Prefferable method's to solving problem numerical:")
            else:
                if method_of_solution == 1:
                    taylorSeries_method()
                elif method_of_solution == 2:
                    picards_method()
                else:
                    eulers_method()

def f(X):
    n= len(coefficients)-1
    total =0
    for x in coefficients:
        total += eval(x)*X**n
        n -= 1
    return total

def initialApproximation():
    array_of_X_coords = np.arange(-4, 5, 1)
    list_of_interval = []
    for i in range(len(array_of_X_coords)-1):
        f_x0 = f(array_of_X_coords[i])
        f_x1 = f(array_of_X_coords[i+1])

        if f_x0 * f_x1 < 0:
            list_of_interval.append((array_of_X_coords[i], array_of_X_coords[i+1]))

    return list_of_interval


def myBisectionMethod(interval, error= 0.5e-4):
    x0, x1 = interval
    if f(x0) * f(x1) >= 0:
        raise ValueError("The function must have opposite signs at the endpoints of the interval for the root to lie in the interval")
    
    while abs(x1-x0)> error:
        c = (x1 + x0)/2

        if f(x0)* f(c)< 0:
            x1 = c

        else:
            x0 = c
        
    return c

def myLinearInterpolation(interval, error= 0.5e-4):
    x0, x1 = interval
    if f(x0) * f(x1) >= 0:
        raise ValueError("The function must have opposite signs at the endpoints of the interval for the root to lie in the interval")
    
    while abs(x1-x0)>error:
        x2 = (x0*f(x1) - x1*f(x0))/(f(x1) -f(x0))

        if f(x0)* f(x2)< 0:
            x1 = x2

        else:
            x0 = x2
        
    return x2

def numerical_first_derivative(x, h=0.1):
    return (f(x+h) -f(x-h))/(2*h)

def numerical_intergration():
    pass


def myNewtonRaphson(interval, error=0.5e-3):
    x0,  x1 = interval
    xk = np.random.choice(np.linspace(x0, x1, 50))
    if f(x0) * f(x1) >= 0:
        raise ValueError("The function must have opposite signs at the endpoints of the interval for the root to lie in the interval")
 
    while abs(f(xk))>error:
        xk = xk - (f(xk)/numerical_first_derivative(xk))      
    return xk

def picards_method():
    print("-------------------Picard's-------------------------")
    n = eval(input("Enter nth order of iteration: ").lower())
    x0 = eval(input("Enter initial value of x0: "))
    y0 = eval(input("Enter initial value of y0 at x0=0: "))
    equation = input("Enter f'(x, y)-> dy_dx in terms x and y: ")
    while not re.match("[a-z]", equation ):
        equation = input("Enter f'(x, y)-> dy_dx in terms x and y: ").lower()
    x, y = sp.symbols('x, y')
    equation = sp.simplify(equation)
    expansion = y0
    for i in range(n):
        expansion = y0 + sp.integrate(sp.Subs(equation, y, expansion), x)
        print(f"(y{i+1}): {sp.simplify(expansion)}")

    x_value = eval(
        input("Enter an x value to substitute into the last iteration: "))
    result = sp.simplify(sp.Subs(expansion, x, x_value))
    print(f"y={result} at x = {x_value}")

def eulers_method():
    print("-------------------Euler's-------------------------")
    n = eval(input("Enter nth order of iteration: ").lower())
    x0 = eval(input("Enter initial value of x0: "))
    y0 = eval(input("Enter initial value of y0 at x0=0: "))
    h = eval(input("Enter h->(Xk-Xk_1) the small steps size: ").lower())
    equation = input("Enter f'(x, y)-> dy_dx in terms x and y: ")
    while not re.match("[a-z0-9]", equation):
        equation = input("Enter f'(x, y)-> dy_dx in terms x and y: ").lower()
    yk = y0
    xk = x0 
    x, y = sp.symbols('x, y')
    equation = sp.simplify(equation)
    for i in range(n): 
        yk = yk + h*sp.simplify(sp.Subs(equation, (x, y), (xk, yk)))
        xk = xk+h
        print(f"x{i+1}, y{i+1}: {xk}, {yk}")

def taylorSeries_method():
    print("-------------------Taylor's-------------------------")
    n = eval(input("Enter nth order of iteration: ").lower())
    x0 = eval(input("Enter initial value of x0: "))
    y0 = eval(input("Enter initial value of y0 at x0: "))
    xk = eval(input("Enter value of b to get y(b): "))
    yk = y0
    x, y = sp.symbols('x, y')
    equation = input("Enter f'(x, y)-> dy_dx in terms x and y: ")
    while not re.match("[a-z0-9]", equation):
        equation = input("Enter f'(x, y)-> dy_dx in terms x and y: ").lower()
    equation = sp.simplify(equation)
    for i in range(n):
        yk += ((x**(i+1))/math.factorial(i+1))*sp.Subs(sp.idiff(equation, y, x), (x, y), (x0, y0))
    print(f"y(x) = {sp.simplify(yk)} \n", f"x ={xk}, y={sp.simplify(sp.Subs(yk, x, xk))}")

main()