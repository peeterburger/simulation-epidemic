# Abschlussbeispiel Computational Statistics
# Gruppe [Namen bitte einfügen]

# Dependencies
import numpy
import random

# Parameter
k = [ 5, 10, 20 ] # Kontakte pro Person
m = [ 1, 5, 10 ] # Infizierte an Tag 1
p = [ 0.1, 0.25, 0.5 ]
n = 5000 # Anzahl der Personen im Modell
fpd = 5 # Erster möglicher Tag, an dem man sterben kann
isol = [ 7, 20 ] # Isolation ab Tag isol - 20 steht hierbei für keine Isolation (nachdem alle Kranken ab Tag 16 gesundet oder tot sind)
d = 0.001 # An jedem weiteren Tag beträgt die Sterbewahrscheinlichkeit 0.1%
# Das ergibt eine Sterbewahrscheinlichkeit von 1 - 0.999^5 bis 1 - 0.999^10 (im Durchschnitt ca. 0.75 %)

# Adjazenzmatrix initialisieren
def createAMatrix(k, n = 5000):
  A = numpy.full((n, n), False, dtype = bool)

  # iterate through rows
  for i in range(n):

    # calcualte index filter
    indexFilter = A.sum(axis = 0) < k
    for j in range(i + 1):
      indexFilter[j] = False

    rowsAllowed = sum(indexFilter)
    rowsRemaining = k - sum(A[i])

    if rowsRemaining > 0 and rowsAllowed > 0:
      numberNewContacts = min(rowsAllowed, rowsRemaining)

      newContacts = numpy.random.choice(numpy.concatenate([
          numpy.full(numberNewContacts, True),
          numpy.full(rowsAllowed - numberNewContacts, False)
      ]), rowsAllowed, replace=False)

      A[i, indexFilter] = newContacts
      A[indexFilter, i] = newContacts

  return A

# Berechne Matrizen und speichere
# A_k5 = createAMatrix(5)
# A_k10 = createAMatrix(10)
# A_k20 = createAMatrix(20)
# save(A_k5, A_k10, A_k20, file = "Adjazenzmatrizen.Rdata")

# Simulations-Funktion (Parameter k, m, p)
def simul_epidemic(n, k, m, p, isol = 20):

  # Adjazenzmatrix laden (je nach k)

  A = createAMatrix(k, n)

  # Krankheitsdauer (gleichverteilt zwischen 10 und 15 Tagen)
  # Tag 1 ist der Tag, an dem Person angesteckt wird - Genesung an Tag x+1 (Ansteckend bis Tag x)
  duration = numpy.random.randint(10, 16, n)

  # Zeithorizont der Simulation (100 reicht anscheinend, mehr ist aber problemlos möglich)
  T = 100
  # Matrix mit den Simulationsergebnissen
  results = numpy.full((T, n), "H")
  # Bisherige Dauer der Ansteckung
  daysSick = numpy.full(n, 0, dtype = int)
  # Neuinfektionen pro Tag
  newInfections = numpy.full(T, 0, dtype = int)

  # Initial infection at day 1
  results[0] = numpy.concatenate([numpy.full(m, "D"), numpy.full(n - m, "H")])

  for t in range(1, T):
    print(t)

    sickPeople = numpy.where(results[t - 1] == "D")[0]

    # stop the simulation, if further infection is impossible
    if sickPeople.size <= 0:
      break

    daysSick[sickPeople] += 1

    # print(sickPeople)

    for sickPerson in sickPeople:
      if daysSick[sickPerson] > duration[sickPerson]:
        # simualtion of recovery
        print(sickPerson, " cured")
        daysSick[sickPerson] = 0
        results[t:T, sickPerson] = "R"
      elif daysSick[sickPerson] >= fpd and random.random() < d:
        # simualtion of fatality
        print(sickPerson, " died")
        daysSick[sickPerson] = 0
        results[t:T, sickPerson] = "T"
      else:
        # current person remains sick
        results[t, sickPerson] = "D"

        # simulation of infection
        contacts = numpy.where(A[sickPerson])[0]
        for contact in contacts:
          if results[t, contact] == "H" and results[t - 1, contact] == "H" and random.random() < p:
            print("%s infected %s" % (sickPerson, contact))
            results[t, contact] = "D"
            newInfections[t] += 1

    # print(results[t])
    print("%d new infections" % newInfections[t])

  print(results[-1])

simul_epidemic(5000, k[0], m[0], p[0], 20)