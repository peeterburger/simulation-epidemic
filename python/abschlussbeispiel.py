import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# example parameters for a simulation
n = 5000                # total number of people
k = [ 5, 10, 20 ]       # number of contacts per person
m = [ 1, 5, 10 ]        # number of people infected at day one
p = [ 0.1, 0.25, 0.5 ]  # infection rate
i = [ 7, 20 ]           # number of days that pass before an infected person isolates itself
d = 0.001               # mortality rate per day
fpd = 5                 # minimum number of days before the desease might become fatal
minDaysInfected = 10    # minimum number of days a person might be infected
maxDaysInfected = 15    # maximum number of days a person might be infected

# creates a contact matrix based on the total number of people and the number
# of contacts per person
def createContactMatrix(n, k, usingFiles=False):

  filePath = "./data/%s_%s.npy" % (n, k)

  if usingFiles:
    try:
      return np.load(filePath)
    except IOError:
      print("File %s not found. Generating new matrix..." % filePath)

  contactMatrix = np.full((n, n), False)

  # iterate through rows
  for i in range(n):

    nContactsRemaining = k - contactMatrix[i].sum()

    # if no contacts have to be added, continue with the next person
    if nContactsRemaining <= 0:
      continue

    # all j < i are automatically excluded, since their column sum is
    # always greater or equal k
    indexFilter = contactMatrix.sum(axis=0) < k

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
      contactMatrix[i, contact] = contactMatrix[contact, i] = True

  if usingFiles:
    np.save(filePath, contactMatrix)

  return contactMatrix

# simulation of desease progression
def simulateEpidemic(n, m, p, i, contactMatrix):

  # time horizon for the simulation
  T = 100

  # duration of the desease (per-person)
  duration = np.random.randint(minDaysInfected, maxDaysInfected + 1, n)

  # matrix that holds the status of each person for each unit of time
  results = np.full((T, n), "H")

  # number of days each person is infected
  daysSick = np.full(n, 0)

  # initial 'm' infections at day 1
  results[0] = np.concatenate([
      np.full(m, "D"),
      np.full(n - m, "H")
  ])

  # simulate all following days
  for t in range(1, T):
    currentlySickPeople = np.where(results[t - 1] == "D")[0]
    daysSick[currentlySickPeople] += 1

    # stop the simulation, if further infection is impossible
    if currentlySickPeople.size <= 0:
      break

    for currentlySickPerson in currentlySickPeople:
      if daysSick[currentlySickPerson] > duration[currentlySickPerson]:
        # simualtion of recovery
        daysSick[currentlySickPerson] = 0
        results[t:T, currentlySickPerson] = "R"
      elif daysSick[currentlySickPerson] >= fpd and np.random.random() < d:
        # simualtion of fatality
        daysSick[currentlySickPerson] = 0
        results[t:T, currentlySickPerson] = "T"
      elif daysSick[currentlySickPerson] <= i:
        # person remains sick
        results[t, currentlySickPerson] = "D"

        # simulation of infection
        contacts = np.where(contactMatrix[currentlySickPerson])[0]
        for contact in contacts:
          if results[t, contact] == "H" and results[t - 1, contact] == "H" and np.random.random() < p:
            results[t, contact] = "D"

  return results

# returns a dataframe with usefull statistics extracted from the result of a simulation
def evaluateResult(result):
  T = result.shape[0]

  infected = np.empty(T, dtype=int)
  healthy = np.empty(T, dtype=int)
  resistant = np.empty(T, dtype=int)
  deaths = np.empty(T, dtype=int)
  newInfections = np.empty(T, dtype=int)

  # basic statistics
  for t, dailyStats in enumerate(result):
    infected[t] = (dailyStats == 'D').sum()
    healthy[t] = (dailyStats == 'H').sum()
    resistant[t] = (dailyStats == 'R').sum()
    deaths[t] = (dailyStats == 'T').sum()

  # new infections statistics
  for t in range(T):
    if t == 0:
      newInfections[t] = (result[0] == 'D').sum()
    else:
      newInfections[t] = np.logical_and(
        result[t - 1] == 'H',
        result[t] == 'D'
      ).sum()

  return pd.DataFrame.from_dict({
    "infected": infected,
    "healthy": healthy,
    "resistant": resistant,
    "number of deaths": deaths,
    "daily new cases": newInfections
  })

contactMatrix = createContactMatrix(n, k[0], True)

result = simulateEpidemic(n, m[0], p[0], 10, contactMatrix)

evaluateResult(result).plot()
plt.savefig('plot.pdf')