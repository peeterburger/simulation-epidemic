# Abschlussbeispiel Computational Statistics
# Gruppe [Namen bitte einfügen]

# Dependencies
import numpy as np
import matplotlib.pyplot as plt

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

    nContactsRemaining = k - A[i].sum()

    # if no contacts have to be added, continue with the next person
    if nContactsRemaining <= 0:
      continue

    # all j < i are automatically excluded, since their column sum is
    # always greater or equal k
    indexFilter = A.sum(axis=0) < k

    # a person is not in contact with itself
    indexFilter[i] = False

    nContactsPossible = indexFilter.sum()
    
    # if it is not possible to add further contacts, continue with
    # the next person
    if nContactsPossible <= 0:
      continue

    nNewContacts = min(nContactsPossible, nContactsRemaining)

    possibleContacts = np.where(indexFilter)[0]
    chosenContacts = np.random.choice(possibleContacts, nNewContacts, replace=False)

    for contact in chosenContacts:
      A[i, contact] = A[contact, i] = True

  return A

# Simulations-Funktion (Parameter n, m, p)
def simul_epidemic(n, m, p, contactMatrix, isol=20):

  # Krankheitsdauer (gleichverteilt zwischen 10 und 15 Tagen)
  # Tag 1 ist der Tag, an dem Person angesteckt wird - Genesung an Tag x+1 (Ansteckend bis Tag x)
  duration = np.random.randint(10, 16, n)

  # Zeithorizont der Simulation (100 reicht anscheinend, mehr ist aber problemlos möglich)
  T = 100
  # Matrix mit den Simulationsergebnissen
  results = np.full((T, n), "H")
  # Bisherige Dauer der Ansteckung
  daysSick = np.full(n, 0)

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
        contacts = np.where(contactMatrix[currentlySickPerson])[0]
        for contact in contacts:
          if results[t, contact] == "H" and results[t - 1, contact] == "H" and np.random.random() < p:
            results[t, contact] = "D"

  return results

def evaluateResult(result):
  T = result.shape[0]

  currentlyInfected = np.empty(T, dtype=int)
  currentlyHealthy = np.empty(T, dtype=int)
  currentlyResistant = np.empty(T, dtype=int)
  currentlyDead = np.empty(T, dtype=int)
  newInfections = np.empty(T, dtype=int)

  # basic statistics
  for t, dailyStatistics in enumerate(result):
    currentlyInfected[t] = (dailyStatistics == 'D').sum()
    currentlyHealthy[t] = (dailyStatistics == 'H').sum()
    currentlyResistant[t] = (dailyStatistics == 'R').sum()
    currentlyDead[t] = (dailyStatistics == 'T').sum()

  # new infections
  for t, dailyStatistics in enumerate(result):
    if t == 0:
      newInfections[t] = (result[0] == 'D').sum()
    else:
      newInfections[t] = np.logical_and(
          result[t - 1] == 'H',
          result[t] == 'D'
      ).sum()

  return (currentlyInfected, currentlyHealthy,
      currentlyResistant, currentlyDead, newInfections)

def plotStatistics(statistics):
  plt.plot(np.arange(100), statistics[0])
  plt.plot(np.arange(100), statistics[1])
  plt.plot(np.arange(100), statistics[2])
  plt.plot(np.arange(100), statistics[3])
  plt.plot(np.arange(100), statistics[4])
  plt.xlabel("days")
  plt.ylabel("people")
  plt.savefig("plot3.pdf")

n = 2000

contactMatrix = createAMatrix(k[0], n)

result = simul_epidemic(n, m[0], p[0], contactMatrix)
statistics = evaluateResult(result)

plotStatistics(statistics)