#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - Curso 2021/2022
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O'+str(i)+')0\t= '+str([round(j,3) for j in O[i]]))
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))

def muestra_robot(O,obj):
  # Muestra el robot graficamente
  plt.figure(1)
  plt.xlim(-L,L)
  plt.ylim(-L,L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
  plt.plot(obj[0], obj[1], '*')
  plt.show()
  input()
  plt.clf()

def matriz_T(d,th,a,al):
  # Calcula la matriz T (ángulos de entrada en radianes)
  
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# valores articulares arbitrarios para la cinemática directa inicial
th=[0.,0.,0.]
a =[5.,5.,5.]

prismatica = [False, True, False]

# # radianes
limMax = [90 * (pi/180), 10, 90 * (pi/180)]
limMin = [-90 * (pi/180), 0, -90 * (pi/180)]

# limMax = [90 * pi/180, 90 * pi/180, 90 * pi/180]
# limMin = [-90 * pi/180, -90 * pi/180, -90 * pi/180]

L = sum(a) # variable para representación gráfica
EPSILON = .01

plt.ion() # modo interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 3:
  sys.exit("python " + sys.argv[0] + " x y")
objetivo=[float(i) for i in sys.argv[1:]]

O=list(range(len(th)+1)) # Reservamos estructura en memoria
O[0]=cin_dir(th,a) # Calculamos la posicion inicial
print ("- Posicion inicial:")
muestra_origenes(O[0])

dist = float("inf")
prev = 0.
iteracion = 1
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist
  pos = len(th) - 1
  # Para cada combinación de articulaciones:
  for i in range(len(th)):
    # cálculo de la cinemática inversa:
    if not prismatica[pos - i]:
      # codigo para angulos
      xR = O[i][pos - i][0] # Punto del que se calcula th
      yR = O[i][pos - i][1] # Punto del que se calcula th
      xE = O[i][len(th)][0] # Punto final
      yE = O[i][len(th)][1] # Punto final
      
      x = objetivo[0] - xR
      y = objetivo[1] - yR
      
      alpha = atan2(y, x) # Incremento del angulo entre el objetivo y el punto del que estamos calculando th
      
      x = xE - xR
      y = yE - yR
      
      alpha2 = atan2(y, x) # Incremento del angulo entre el punto final y el punto del que estamos calculando th
      
      thita = alpha - alpha2
      th[pos - i] = thita + th[pos - i]
      
      # Comprobacion de limites

      while (th[pos - i] > pi) or (th[pos - i] < -pi):
        if th[pos - i] > pi:
          th[pos - i] -= 2 * pi
        if th[pos - i] < -pi:
          th[pos - i] += 2 * pi

      if limMax[pos - i] < th[pos - i]:
        th[pos - i] = limMax[pos - i]
      if limMin[pos - i] > th[pos - i]:
        th[pos - i] = limMin[pos - i]
      
      O[i+1] = cin_dir(th,a)
      
    else:
      # codigo para prismaticas
      w = 0
      for o in range(pos - i):
        w += th[o] # Angulos acumulados
      
      vector = [cos(w), sin(w)]
      rON = [objetivo[0] - O[i][len(th)][0], objetivo[1] - O[i][len(th)][1]]
      d = np.dot(vector, rON).tolist() 
      a[pos - i] += d
      # Limites
      if limMax[pos - i] < a[pos - i]:
        a[pos - i] = limMax[pos - i]
      if limMin[pos - i] > a[pos - i]:
        a[pos - i] = limMin[pos - i]

      O[i+1] = cin_dir(th,a)

  dist = np.linalg.norm(np.subtract(objetivo,O[-1][-1]))
  print ("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O[-1])
  muestra_robot(O,objetivo)
  print ("Distancia al objetivo = " + str(round(dist,5)))
  iteracion+=1
  O[0]=O[-1]

if dist <= EPSILON:
  print ("\n" + str(iteracion) + " iteraciones para converger.")
else:
  print ("\nNo hay convergencia tras " + str(iteracion) + " iteraciones.")
print ("- Umbral de convergencia epsilon: " + str(EPSILON))
print ("- Distancia al objetivo:          " + str(round(dist,5)))
print ("- Valores finales de las articulaciones:")
for i in range(len(th)):
  print ("  theta" + str(i+1) + " = " + str(round(th[i],3)))
for i in range(len(th)):
  print ("  L" + str(i+1) + "     = " + str(round(a[i],3)))
