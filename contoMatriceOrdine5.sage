## qui voglio vedere cosa succede di una cubica che viene
## stabilizzata da una matrice di ordine 5
## (dovrebbe avere almeno 10 punti di Eck)

## Se non ho sabgliato i conti, se una matrice di ordine 5 stabiliazza
## la generica cubica, allora deve madare L* in uno di questi 9 L set:
tuttiL2 = [("E2", "G4", "E5", "G3", "E3"), ("E2", "G4", "E5", "G3", "F14"),\
           ("E2", "G4", "E5", "G3", "F46"),("E2", "G4", "E6", "G3", "E3"), \
	   ("E2", "G4", "E6", "G3", "F14"), ("E2", "G4", "E6", "G3", "F45"),\
           ("E2", "G4", "F34", "G3", "F14"), ("E2", "G4", "F34", "G3", "F45"),\
	   ("E2", "G4", "F34", "G3", "F46")]

## procedo con il codice:
load('Imports/Utility.sage', 'Imports/Point.sage', 'Imports/Line.sage', 'Imports/TritangentPlane.sage', 'Imports/Group.sage', 'Imports/Cubic.sage')
import multiprocessing as mp
import pickle
Q = QQ
P.<x,y,z,t,b,c,d,e,f,l,m> = PolynomialRing(Q)
cubic_new = e*f*(2*x^2*y-2*x*y^2+x*z^2-x*z*t-y*t^2+y*z*t)+b*c*(x-t)*(x*z+y*t)+c*c*(z+t)*(y*t-x*z)+d*c*(y-z)*(x*z+y*t)+(e+f)*c*(x-y)*(y*t-x*z)
sing_cubics = (-1) * (-c + f) * (-c + e) * c * (c + f) * (c + e) * (-e + f)^2 * (-c*d + c*f + e*f) * (-c*d + c*e + e*f) * (-c^2 - c*d + e*f) * (b*c - c*f + e*f) * (b*c - c*e + e*f) * (b*c - c*d + 2*e*f) * (b*c - c^2 + e*f) * (b*c^2 + c^2*d + b*c*f - 2*c^2*f - c*d*f + 2*e*f^2) * (b*c^2 + c^2*d + b*c*e - 2*c^2*e - c*d*e + 2*e^2*f) * (-b*c^3 - 2*b*c^2*d + c^3*d + b*c^2*e + c^2*d*e + b*c^2*f + c^2*d*f + 3*b*c*e*f - 4*c^2*e*f - 3*c*d*e*f + 4*e^2*f^2)
line = Line([y, z])
general_cubic = Cubic(cubic_new, line, sing_cubics)
S = general_cubic

Lbase = ("E1", "G4", "E2", "G3", "E3")

## Guardo le condizioni che vengono imposte alla generica cubica S
## quando mando L* = Lbase in uno degli L-set di tuttiL2

mon = ((x+y+z+t)^3).monomials()

tutteCondizioni = []  ## qui accumulo l'L-set e le condizioni che impone
                      ## sui coefficienti
for L2 in tuttiL2:
  tt1 = cputime() 
  LbasePl = S.get_L_set_in_plucker(Lbase) 
  L2Pl = S.get_L_set_in_plucker(L2)
  M = find_projectivity(LbasePl, L2Pl)
  S2 = S.subs(change_coord(M))
  Mcf = matrix([[S.coefficient(mn) for mn in mon],\
                [S2.coefficient(mn) for mn in mon]])
  cond = list(Set(Mcf.minors(2)))
  print("Calcolo la dec primaria delle condizioni")
  print("Per accelerare i conti conviene prima saturare con c")
  J = P.ideal([remove_sing_factors(cc, S.sing_locus) for cc in cond if cc != 0])
  pd = J.saturation(ideal(c))[0].radical().primary_decomposition()
  tutteCondizioni.append([L2, pd])
  print("Tempo di calcolo:")
  print(cputime()-tt1)
  tt1 = cputime()
  print("condizioni:")
  print(pd)
  print("\n\n\n")
  sleep(1)
  