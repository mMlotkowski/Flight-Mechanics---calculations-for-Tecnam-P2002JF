# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 10:25:02 2023

@author: mMlotkowski
"""

import math
import matplotlib.pyplot as mtp

"OPORY LOTU"
rpod = 1.35
rs = 1.1
rk = 1
rpod = 1.35
rust = 1.24

"PARAMETRY LOTU"
Vkmh = 56 #V MINIMALNA

Ns = 69000 #MOC

"WYMIARY SAMOLOTU"
m = 580  #MASA
S = 11.5  #POW SKRZYDŁĄ
lt = 6.61 #ROZPIETOSC
hk = 2.43 #WYS KADLUBA
bk = 1.1 #SREDNCA KADLUBA
ksi = 0
c = 18
b = 100

"PARAMETRY AERODYNAMICZNE"
MIp = 1.45*10**-5 #LEPKOSC
lamb = 6.43 #WYDLUZENIE
e = 1.78*(1-(0.045*(lamb**0.68)))-0.64  #OSWALD
Q = m*9.81 #CIEZAR
ro0 = 1.23  #GESTOSC

lista_Cx = []
lista_V = []
lista_Cx0 = []
lista_P = []
lista_Cz = []
lista_A = []


while Vkmh < 248.4:
    V = Vkmh/3.6
    lista_V.append(V)
    Re = (V*lt)/MIp
    rre = 47*(Re**-(0.2))
    CxSp = 0.0054*rs*(1+3*(c/b)*(math.cos(ksi))**2)*S
    CxSk = 0.031*rk*lt*(bk+hk)
    CxSsil = 0.15*bk*hk
    Cx0 = (rre*rpod*(rust*(CxSp+CxSk)+CxSsil))/S
    lista_Cx0.append(Cx0)
    Cz = (2*Q)/(ro0*(V**2)*S)
    lista_Cz.append(Cz)
    Cxi = (Cz**2)/(3.14*lamb*e)
    Cx = Cx0 + Cxi
    lista_Cx.append(Cx)
    P = Cx*((ro0*(V**2))/(2))*S
    lista_P.append(P)
    A = Cz/Cx
    lista_A.append(A)
    print(Cx)
    Vkmh = Vkmh+ 0.1

"P=f(v)"
mtp.figure(1)
mtp.plot(lista_V, lista_P, linewidth=1, color='black')
mtp.xlabel('Prędkosć')
mtp.ylabel('Moc')


"Cx = f(v)"
mtp.figure(2)
mtp.plot(lista_V, lista_Cx0, linewidth=1, color='black')
mtp.xlabel('Prędkosć')
mtp.ylabel('Cx0')


"A = f(v)"
mtp.figure(3)
mtp.plot(lista_V, lista_A, linewidth=1, color='black')
mtp.xlabel('Prędkosć')
mtp.ylabel('Doskonałosć aerodynamiczna')


mtp.show()