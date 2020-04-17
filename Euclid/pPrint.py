#this code was taken from Ferbruary 2020 Math240 project designed by Doug Baldwin, edited by Joseph Leva April 2020

#Its use for this project is to do polynomial multiplcation, subtraction and printing


# A class that represents polynomials for use in a simple "polynomial pocket calculator" application.


from operator import add, sub

from math import isclose
from sympy.parsing.sympy_parser import parse_expr
from sympy import *
from sympy import poly,Poly,QQ,ZZ,GF,sympify
from sympy.abc import x




class Polynomial :



    # Internally, polynomials are represented by two pieces of information:
    #   o coefficients, a list of the polynomial's coefficients
    #   o degree, the polynomial's degree.
    # These attributes obey several class invariants, including...
    #   o degree = len( coefficients ) - 1
    #   o If the polynomial is a_n x^n + a_{n-1} x^{n-1} + ... + a_1 x + a_0, then coefficients[0] = a_0,
    #     coefficients[1] = a_1, and so forth through coefficients[n] = a_n
    #   o coefficients[n] != 0, unless degree = 0.




    # Initialize a polynomial from a list of its coefficients. Preconditions for this constructor are...
    #   o If the intended polynomial is a_n x^n + a_{n-1} x^{n-1} + ... + a_1 x + a_0, then the coefficients list should
    #     contain a_0 through a_n, in that order.
    #   o The leading coefficient, a_n aka coefficients[n], must not be 0, unless that is the only term in the polynomial.

    def __init__( self, coefficients ) :
        self.set( coefficients )


    #Produce a string representation of this polynomial.

    def __str__(self) :

        poly = str(self.coefficients[0])

        for e in range(1, len(self.coefficients)) :
            if not isclose(self.coefficients[e], 0.0) :
                poly = poly + " + " + str((self.coefficients[e])) + "x^" + str(e)

        return poly




    # Change this polynomial, i.e., replace its value with a new value, specified by a new list of coefficients. This
    # has the same preconditions as the constructor.

    def set( self, coefficients ) :
        self.coefficients = coefficients
        self.degree = len( coefficients ) - 1




    # Get the list of coefficients for this polynomial.

    def getCoefficients( self ) :
        return self.coefficients




    # Get this polynomial's degree.

    def getDegree( self ) :
        return self.degree

    def normalizeCoefficients(coefficients):

        while len(coefficients) > 0 and coefficients[-1] == 0.0:
            coefficients = coefficients[: -1]

        if len(coefficients) == 0:
            coefficients = [0.0]

        return coefficients

    def poly2Sympy(self):
        poly = str(self.coefficients[0])

        for e in range(1, len(self.coefficients)):
            if not isclose(self.coefficients[e], 0.0):
                poly = poly + "+" + str((self.coefficients[e])) + "*x**" + str(e)

        return str(poly)

    def poly2SympyInt(self):
        array= []
        for i in range(0,len( self.coefficients)):
            array.append(int(self.coefficients[i]))
        poly = str(array[0])

        for e in range(1, len(array)):
            if not isclose(array[e], 0.0):
                poly = poly + "+" + str((array[e])) + "*x**" + str(e)

        return str(poly)







