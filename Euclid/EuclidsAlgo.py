from sympy.parsing.sympy_parser import parse_expr
from sympy import *
from sympy import poly,QQ,ZZ,GF,sympify
from sympy.abc import x
from pPrint import Polynomial


#The purpse of this program is to show the Extended Euclid's Algorithm
#Created by Joseph Leva for SUNY Geneseo Great Day April 2020
#files necessary to run: pPrint.py


#According to documentation, due to a bug in sympy galois fields represented GF() actually represented as modulus
#This helps us as we want modulus and not GF but if that bug were to be fixed this program would only work if m is prime.



#runDomain lets us call EuclidsAlgorithm  with whatever domain the user would like to input
def runDomain(a,b,m):
    if m == 'rat':
        EuclidsAlgorithm(a, b, m, 'QQ')
    else:
        EuclidsAlgorithm(a, b, m, GF(m))

#Euclids algorithm preforms gcd calculations takes polynomials, a and b, mod m(if applicable), and the domain of the
#polynomials as arguments
def EuclidsAlgorithm(a,b,m,dom):

    #the following statements set the polynomials in their domain and prepares them for further use
    a=poly(a)
    b= poly(b)
    a= a.as_poly(domain = dom)
    b= b.as_poly(domain = dom)
    ab = largesmall(a,b)
    a= ab[0]
    b= ab[1]
    #divide to obtain quotient q and remainder r
    q,r = div(a,b, domain= dom )

    #create empty arrays to accumulate quotients and remainders
    remainders=[ab[0],ab[1]]
    quotients= [[0],[0]]
    remainders.append(r)
    quotients.append(q)

    #loop collects remainders and quotients untill divide by 0 error, count used to determinne gcd

    x= True
    count = 2
    while x:
        try:
            count += 1
            q,r= div(remainders[count-2],remainders[count-1], domain = dom)
            remainders.append(r)
            quotients.append(q)


        except:
            x= False


    #the last value in remainders will always be 0, thus the gcd is the penultimate element in the list
    temp= remainders[count-2]
    gcd = temp

    #obtain leading coefficient->lc for use in obtaining monic equation Aa+Bb=gcd, made into its own sympy polynomial for div
    # to make the gcd a monic polynomial we divide it by its leading coefficient, luckily sympy has a function for that
    # we then take gcd as a polynomial argument to use for printing later

    lc = gcd.all_coeffs()[0]
    lc = sympify(float(lc))

    gcd = monic(gcd)
    gcd = Polynomial(gcd.all_coeffs()[::-1])

    #printData prints Recursive division, GCD, and rearanngement of equations so the user can see how the inverse euclid
    #algorithm s obtained
    printData(remainders,quotients,count,gcd)

    #use lc is used to circumvent sympy bug which doesnt allow for polynomials of degree 0 unless in function arguments
    useLC(remainders,quotients,count,gcd,lc,dom,m)

def EEA(remainders,quotients,count,gcd,dom, m,lcPoly,oneP):
    # Algorithm for obtaining A and B in linear combination: goes through the euclidalgorithm backwards
    A = oneP
    B = quotients[count - 2].neg()

    #pAb keeps track of string forms of A and B, ptrack tracks the amount of cycles A and A go through for output
    #pExplain explains to the user the calcualtion that is made
    pAB= [str(Polynomial(A.all_coeffs()[::-1])), str(Polynomial(B.all_coeffs()[::-1]))]
    pExplain=['']
    ptrack = 1


    for j in reversed(range(2, count - 2)):
        ptrack+= 1
        tmp = B
        # pAB.append statement is string form of B statement above
        pExplain.append( '( '+str(Polynomial(A.all_coeffs()[::-1]))+') - ( '+ str(Polynomial(quotients[j].all_coeffs()[::-1])) +' ) * ( '+ str(Polynomial(B.all_coeffs()[::-1]))+')')
        B = A - quotients[j] * B
        pAB.append(str(Polynomial(B.all_coeffs()[::-1])))
        A = tmp


    printReverseData(pAB,pExplain,ptrack,lcPoly.all_coeffs()[0])







#statements for printing linear combination, we print rationals also in float form for user readability
    if m == 'rat':
        # must divide by lc to maintain monic gcd
        A = div(A, lcPoly, domain= dom)[0].all_coeffs()[::-1]
        B = div(B, lcPoly, domain= dom)[0].all_coeffs()[::-1]
        Aq = Polynomial(A)
        Bq = Polynomial(B)
        a = Polynomial(remainders[0].all_coeffs()[::-1])
        b = Polynomial(remainders[1].all_coeffs()[::-1])
        Af = Polynomial(makeFloat(A))
        Bf = Polynomial(makeFloat(B))

        print('The linear Combination Yields:')
        print('\t( ' + str(Aq) + ' )( ' + str(a) + ' ) + ( ' + str(Bq) + ' )(' + str(b) + ' ) = ' + str(gcd) + '\n\n')
        print('In float form:')
        print('\t( ' + str(Af) + ' )( ' + str(a) + ' ) + ( ' + str(Bf) + ' )(' + str(b) + ' ) = ' + str(gcd) + '\n\n')
        print('___________________________________________________________________________________________________')

    else:
        #we use Posm to ensure the final output looks correct ie. -4 = 3 in mod 7, there is no need to print float values

        # must divide by lc to maintain monic gcd
        A = div(A, lcPoly, domain= dom)[0].all_coeffs()[::-1]
        B = div(B, lcPoly, domain= dom)[0].all_coeffs()[::-1]
        A = Polynomial(makePosm(makeFloat(A), m))
        B = Polynomial(makePosm(makeFloat(B), m))
        a = Polynomial(remainders[0].all_coeffs()[::-1])
        b = Polynomial(remainders[1].all_coeffs()[::-1])

        print('The linear Combination Yields:')
        print('\t( '+str(A)+' )( '+str(a)+' ) + ( '+ str(B)+' )('+str(b)+' ) = '+ str(gcd)+'\n')
        print('___________________________________________________________________________________________________')





#useLc are neccesary as sympy would not allow me to create polynomial pythons of degree 0

def useLC(remainders, quotients, count,gcd,lc,dom,m):

    return EEA(remainders, quotients, count, gcd,dom,m, lcPoly= poly(lc, x, domain= dom),oneP=poly(1,x,domain= dom))






#makeFloat used to make list out of sympy objects floats
def makeFloat(list):
    array=[]
    for i in range(0,len(list)):
        array.append(float(list[i]))

    return array

#makePosm used in EEAm to make all values positive for linear combination print
def makePosm(list,m):
    array = []
    for i in range(0, len(list)):
        array.append(list[i] % m)

    return array






#makePoly turns sympy object polynomials into 'Polynomial" polynomials objects from pPrint for printing
def makePoly(x):
    newarray=[]
    for i in range(0,len(x)):
        try:
            L = x[i].all_coeffs()
            newL= L[::-1]
            newarray.append(Polynomial(newL))
        except:
            newarray.append(Polynomial((x[i])))
    return newarray




#print data is meant to print out the recursive division, gcd and rearanging of equations for display to user
def printData(x,y,count,gcd):
    #We use makePoly to make remainders and quotients polynomial objects
    premainders= makePoly(x)
    pquotients= makePoly(y)
    print('___________________________________________________________________________________________________')
    print("Recursive Division:")
    print('\t( '+ str(premainders[0])+' )' + ' = '  + '( ' + str(premainders[1])+' )*( '+ str(pquotients[2]) + ' ) + ( '+ str(premainders[2]) +' )\n')
    #This loop prints equations from the while loop, starts at 3 so equations start with b
    for i in range(3,count):
        print('\t( '+str(premainders[i-2])+' )'+' = ('+ str(premainders[i-1])+' )*( '+ str(pquotients[i])+ ' ) + ( '+ str(premainders[i])+' )\n')

    print('___________________________________________________________________________________________________\n')
    print('Greatest Common Divisor: ' + '\n\tgcd((' + str(premainders[0]) + ' ),( ' + str(
        premainders[1]) + ' ))' + ' = ' + str(premainders[count-2]) + '\n')
    print('Monic Greatest Common Divisor: ' + '\n\tgcd((' + str(premainders[0]) + ' ),( ' + str(premainders[1]) + ' ))' + ' = ' + str(gcd) + '\n')
    print('___________________________________________________________________________________________________')
    print("Re-Arranging For Remainders:")
    print('\t( '+ str(premainders[2])+' )' + ' = '  + '( ' + str(premainders[0])+' ) - ( '+ str(premainders[1]) + ' )*( '+ str(pquotients[2]) +' )\n')
    #This loop prints equations from the while loop, starts at 3 so equations start with b
    for i in range(3,count):
        print('\t( '+str(premainders[i])+' )'+' = ('+ str(premainders[i-2])+' ) - ( '+ str(premainders[i-1])+ ' )*( '+ str(pquotients[i])+' )\n')
    print('___________________________________________________________________________________________________')


#printReverseData takes pAB list from EEA to print steps for the Inverse Euclidean Algorithm
def printReverseData(l,e,step,lc):
    print("Steps for the Inverse Euclidean Algorithm: ")
    print('Step 1:')
    print('\tA = ( ' + l[0] + ' )\n\tB = ( ' + l[1] + ' )')
    for i in range(1,step):
        print('Step '+ str(i+1)+ ':')
        print('\tA = ( '+ l[i]+ ' )\n\tB = ( '+ e[i]+' ) = ( '+ l[i+1]+' )')

    print('\nTo obtain monic linear combination, A and B must be divided by the leading coefficient of the gcd.')
    print('In this case the leading coefficient is ',lc)
    print('___________________________________________________________________________________________________\n')



#checks which polynomial should be dividend and divisor and sets them accordingly, if both are zero raises an error
def largesmall(a,b):
    if a== b and b == 0:
        raise ValueError("Can not divide zero by zero")
    dega= int(a.total_degree())
    degb= int(b.total_degree())
    if dega > degb:
        large = a
        small= b
    elif dega < degb:
        large = b
        small = a
    else:
        large = a
        small = b
        coeffa = a.all_coeffs()
        coeffb= b.all_coeffs()
        for i in range(0,len(coeffa)):
            if coeffa[i]< coeffb[i]:
                large = b
                small = a
                break
    return large,small









# The rest of code creates input environment

d=''

#inputPolynomial from Doug Baldwin 2020 math 240 Class
#the method takes string of number characters and outputs as a list of floats

def inputPolynomial( prompt ) :

    coefficients = [ float(c) for c in input( prompt ).split() ]
    return Polynomial( coefficients )


#check 1 and check2  allow the user to input polynomials as string of coefficients seperated by ' '
#it then turns these inputs into Polynomial  objects from pPrint to print to user
check1 = False
check2 = False

while check1 == False :

    print("\nInput polynomials as coefficients seperated by ' ' from degree= 0 to largest degree.")
    print(" For example: |input: 1 0 2 3 |output: 1 + 0x + 2x^2 + 3x^3 |\n")
    p = inputPolynomial("\tEnter first polynomial:")
    print("\nYou entered:", p)
    input1 = input("\tIf this is correct enter: 'yes' ")
    if input1 == 'yes' or input1== 'Yes':
        check1 = True



while check2 == False:

    q = inputPolynomial("\n\tEnter second polynomial: ")
    print("\nYou entered:", q)
    input2 = input("\tIf this is correct enter: 'yes' ")
    if input2.strip() == "yes" or input2== 'Yes' :
        check2 = True

#the following loop is used to check if the user wants to work in the rationals or integers ring m
#it sets polynomials p and q to be sympy expressions for the correct domain


loop = True
while loop == True :
    v= input("\nTo work with the Rationals enter 'Rationals' to work with integer rings enter 'Integers': ")
    if v == 'Rationals' or v == 'rationals':
        p = p.poly2Sympy()
        q = q.poly2Sympy()
        p = parse_expr(p)
        q = parse_expr(q)
        d='rat'
        loop = False

    elif v == 'Integers' or v == 'integers':
        p = p.poly2SympyInt()
        q = q.poly2SympyInt()
        p = parse_expr(p)
        q = parse_expr(q)
        d= int(input("Enter a natural number: "))
        loop = False




runDomain(p, q, d)

# Easy to read example: x^6+3x^4+2x^2+1 , x^3+2x+1
#  1 0 2 0 3 0 1  , 1 2 0 1
