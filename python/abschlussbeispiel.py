# Abschlussbeispiel Computational Statistics
# Gruppe [Namen bitte einfügen]

# Dependencies
import numpy as np

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
def createAMatrix(k, n):
  A = np.full((n, n), False)

  # iterate through rows
  for i in range(n):

    # calcualte the index filter
    indexFilter = A.sum(axis=0) < k
    for j in range(i + 1):
      indexFilter[j] = False

    nContactsPossible = sum(indexFilter)
    nContactsRemaining = k - sum(A[i])

    if nContactsRemaining > 0 and nContactsPossible > 0:
      nNewContacts = min(nContactsPossible, nContactsRemaining)

      # hypergeometric distribution of contacts
      newContacts = np.random.choice(np.concatenate([
          np.full(nNewContacts, True),
          np.full(nContactsPossible - nNewContacts, False)
      ]), nContactsPossible, replace=False)

      A[i, indexFilter] = newContacts
      A[indexFilter, i] = newContacts

  return A

# Simulations-Funktion (Parameter n, k, m, p)
def simul_epidemic(n, k, m, p, isol=20):

  A = createAMatrix(k, n)

  # Krankheitsdauer (gleichverteilt zwischen 10 und 15 Tagen)
  # Tag 1 ist der Tag, an dem Person angesteckt wird - Genesung an Tag x+1 (Ansteckend bis Tag x)
  duration = np.random.randint(10, 16, n)

  # Zeithorizont der Simulation (100 reicht anscheinend, mehr ist aber problemlos möglich)
  T = 100
  # Matrix mit den Simulationsergebnissen
  results = np.full((T, n), "H")
  # Bisherige Dauer der Ansteckung
  daysSick = np.full(n, 0)
  # Neuinfektionen pro Tag
  newInfections = np.full(T, 0)

  # Initial infection at day 1
  results[0] = np.concatenate([
      np.full(m, "D"),
      np.full(n - m, "H")
  ])

  for t in range(1, T):
    currentlySickPeople = np.where(results[t - 1] == "D")[0]
    daysSick[currentlySickPeople] += 1

    # stop the simulation, if further infection is impossible
    # if currentlySickPeople.size <= 0:
    #   break

    for currentlySickPerson in currentlySickPeople:
      if daysSick[currentlySickPerson] > duration[currentlySickPerson]:
        # simualtion of recovery
        daysSick[currentlySickPerson] = 0
        results[t:T, currentlySickPerson] = "R"
      elif daysSick[currentlySickPerson] >= fpd and np.random.random() < d:
        # simualtion of fatality
        daysSick[currentlySickPerson] = 0
        results[t:T, currentlySickPerson] = "T"
      else:
        # current person remains sick
        results[t, currentlySickPerson] = "D"

        # simulation of infection
        contacts = np.where(A[currentlySickPerson])[0]
        for contact in contacts:
          if results[t, contact] == "H" and results[t - 1, contact] == "H" and np.random.random() < p:
            results[t, contact] = "D"
            newInfections[t] += 1

  return results

print(simul_epidemic(100, k[0], m[0], p[0]))