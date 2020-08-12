# Abschlussbeispiel Computational Statistics
# Gruppe [Namen bitte einfügen]

# Dependencies (Betrifft nur die Funktion rdunif für diskrete gleichverteilte Zufallszahlen)
library(purrr)

# Parameter
k = c(5, 10, 20) # Kontakte pro Person
m = c(1, 5, 10) # Infizierte an Tag 1
p = c(0.1, 0.25, 0.5)
n = 5000 # Anzahl der Personen im Modell
fpd = 5 # Erster möglicher Tag, an dem man sterben kann 
isol = c(7, 20) # Isolation ab Tag isol - 20 steht hierbei für keine Isolation (nachdem alle Kranken ab Tag 16 gesundet oder tot sind)
d = 0.001 # An jedem weiteren Tag beträgt die Sterbewahrscheinlichkeit 0.1%
# Das ergibt eine Sterbewahrscheinlichkeit von 1 - 0.999^5 bis 1 - 0.999^10 (im Durchschnitt ca. 0.75 %)

# Adjazenzmatrix initialisieren
createAMatrix = function(k, n = 5000) {
  nfails = 0
  while (T) {
    A = matrix(data = rep(F, n^2), ncol = n, nrow = n)
    
    # Indexing: Zuerst Zeile, dann Spalte
    
    for (i in 1:(n-1)) {
      # Zähle, wie viele Kontakte für Individuum i noch möglich sind
      remaininginrow = k - sum(A[1:(i-1),i])
      
      # Indexfilter: Zuerst werden alle ausgeschlossen für j <= i (obere Dreiecksmatrix bleibt übrig)
      indexfilter = c(rep(F, i),rep(T, n - i))
      for (j in (i+1):n) {
        # Dann: Ausgeschlossen, wenn Spalte schon hinreichend besetzt ist
        if (sum(A[1:(n-1), j]) == k) {
          indexfilter[j] = F
        }
      }
      
      for (j in (i+1):n) {
        
        if(remaininginrow > 0 & indexfilter[j] == T) {
          # Berechne Wahrscheinlichkeit für Kontakt (Anzahl der fehlenden Kontakte in Reihe / Anzahl der Spalten mit Index >= dieser Spalte, in die noch Kontakte gesetzt werden können)
          sampleprob = min(remaininginrow / sum(indexfilter[j:n]), 1)
          
          # Setze den Kontakt auf T, falls Zufall das ergibt
          A[i,j] = A[j, i] = sample(c(T, F), size = 1, prob = c(sampleprob, (1 - sampleprob)))
          
          if (A[i, j] == T) {
            # Wenn Wert T gesetzt wurde, fehlt in der Zeile einer weniger
            remaininginrow = remaininginrow - 1
          }
        }
      }
    }
    
    # Test: Sind alle Zeilen und Spalten mit k Werten besetzt? Falls nein: Warnung
    checksums = sum(c(rowSums(A), colSums(A)) != k)
    if(checksums != 0) {
      # warning("Adjazenzmatrix schlecht initialisiert. ", checksums, " Zeilen und/oder Spalten haben mehr oder weniger als ", k, " Kontakte.") # Warnung nicht notwendig, da sowieso wiederholt wird
      nfails = nfails + 1
      print(paste("Neuer Versuch (Nr. ", (nfails + 1), ")", sep = ""))
    } else {
      return(A)
      break
    }
  }
}

# Berechne Matrizen und speichere
# A_k5 = createAMatrix(5)
# A_k10 = createAMatrix(10)
# A_k20 = createAMatrix(20)
# save(A_k5, A_k10, A_k20, file = "Adjazenzmatrizen.Rdata")
# load("Adjazenzmatrizen.Rdata")

# Simulations-Funktion (Parameter k, m, p)
simul_epidemic = function(k, m, p, isol = 20){

  # Adjazenzmatrix laden (je nach k)
  
  if (k %in% c(5, 10, 20)) {
    A = eval(parse(text = paste("A_k", k, sep = ""))) # Wenn k gleich 5, 10 oder 20, wird eine der vorausberechneten Matrizen verwendet
  } else if (is.integer(k) & k < n & k > 0) {
    print("Adjazenzmatrix wird generiert. Das kann eine Weile dauern...")
    A = createAMatrix(k) 
  } else {
    simpleError("The value for n is not allowed.")
  }
  
  # Krankheitsdauer (gleichverteilt zwischen 10 und 15 Tagen)
  # Tag 1 ist der Tag, an dem Person angesteckt wird - Genesung an Tag x+1 (Ansteckend bis Tag x)
  duration = rdunif(n, 15, 10)
  
  # Kontakt zu infizierter Person (Funktion)
  riskyContact = function(t, i, p) {
    if (runif(1) < p) {
      # print(paste(t, i))
      results[t, i] <<- "D"
      daysSick[i] <<- 1
      newInfections[t] <<- newInfections[t] + 1
    }
  }
  
  # gesunden
  cure = function(t, i) {
    daysSick[i] <<- NA
    results[t:T, i] <<- "R"
  }
  
  # sterben
  die = function(t, i) {
    daysSick[i] <<- NA
    results[t:T, i] <<- "T"
  }
  
  # Zeithorizont der Simulation (100 reicht anscheinend, mehr ist aber problemlos möglich)
  T = 100
  # Matrix mit den Simulationsergebnissen
  results = matrix("H", nrow = T, ncol = n)
  # Bisherige Dauer der Ansteckung
  daysSick = rep(NA, n)
  # Neuinfektionen pro Tag
  newInfections = rep(0, T)
  
  # Anfängliche Ansteckung 
  initialInfection = function(m, n = 5000) sample(c(rep("D", m), rep("H", n - m)), size = n, replace = F)
  initial = initialInfection(m = m)
  sick = which(initial == "D")
  
  # Tag 1
  t = 1
  
  # Anfangs angesteckte sind noch krank
  results[1, sick] = "D"
  daysSick[sick] = 1
  
  # Anstecken
  for (i in sick) {
    contacts = which(A[i, ] == TRUE) # Kontakte von Person i 
    for (j in contacts) {
      riskyContact(t, j, p = p) # Leute, die noch nicht krank waren und mit Kranken einen Kontakt haben, haben einen riskanten Kontakt
    }
  }
  
  # Tag 2 und alle darauf folgenden Tage
  for (t in 2:T) {
    # Zuerst: Erfassen, wer am Tag davor noch krank war
    sick = which(results[t-1,] == "D")
    
    # Dauer-Counter steigt an
    daysSick = daysSick + 1
    
    for (i in sick) {
      if (daysSick[i] > duration[i]) {
        cure(t, i)   # Leute, die lang krank waren, werden gesund
      } else if (daysSick[i] >= fpd & runif(1) < d) {
        die(t, i) # Leute, die schon fpd Tage krank sind, können auch sterben
      } else {
        results[t, i] <- "D" # Leute, die nicht gesund werden oder sterben, sind weriterhin krank
      }
    }
    sick = which(results[t,] == "D") # Ansteckend ist nur mehr, wer nicht gesund wird, stirbt oder isoliert wird
    for (i in sick) {
      if (daysSick[i] <= isol) {
      contacts = which(A[i, ] == TRUE)
      for (j in contacts) {
        if(results[t-1, j] == "H") {
          riskyContact(t, j, p = p) # Leute, die noch nicht krank waren und mit Kranken einen Kontakt haben, haben einen riskanten Kontakt
        }
      }
      }
    }
    
  }
  
  nsick = function(v) sum(v == "D")
  nrecovered = function(v) sum(v == "R")
  ndead = function(v) sum(v == "T")
  
  # Tageweise Statistiken berechnen (krank, wieder gesund, tot, noch gesund)
  sick_stats <- apply(results, 1, nsick)
  recovered_stats <- apply(results, 1, nrecovered)
  dead_stats <- apply(results, 1, ndead)
  healthy_stats = n - sick_stats - recovered_stats - dead_stats
  
  summResult = data.frame(newInfections, sick_stats, recovered_stats, dead_stats, healthy_stats)
  colnames(summResult) = c("New Infections", "Sick", "Recovered", "Dead", "Not yet contracted Coronavirus")
  
  plot(summResult$Recovered, col = "forestgreen")
  points(summResult$Sick, col = "orange")
  points(summResult$Dead, col = "black")
  
  return(summResult)
}

# Kombinationen von Parametern (laut Vorschlag)
params = expand.grid(k, m, p)
colnames(params) <- c("k", "m", "p")
listnames = vector(mode = "character", length = nrow(params))

# für jede Kombination simulieren
simul_results = list()
for (r in 1:nrow(params)) {
  simul_results[[r]] = simul_epidemic(k = params[r, 1], m = params[r, 2], p = params[r, 3])
  listnames[r] = paste("k = ", params[r, 1], ", m = ", params[r, 2], ", p = ", params[r, 3], sep = "")
}
names(simul_results) = listnames

# Ergebnis: In allen Parameterkombinationen breitet sich das Virus sehr schnell aus
# Replikationsfaktor - Durchschnittliche Anzahl, die ein infiazaierter ansteckt ergibt sich aus k, p und der Krankheitslänge (im Durchschnitt 12.5 Tage)
# Berechnung für den Anfang einer Epidemie (da noch keine Immunität oder Tod), nachher niedriger
# Anzahl der Kontakte * (1 - Wahrscheinlichkeit, dass ein Kontakt über die gesamte Krankheit nicht angesteckt wird)
params$R0 = params$k * (1 - (1 - params$p)^max(isol, 12.5))

# Mit Isolation allerdings 
isol = 7
params$R0_isol = params$k * (1 - (1 - params$p)^min(isol, 12.5))
# R0 von 2.6 für k=5 und p = 0.1 - schon sehr nahe an den anfangs in den Medien kolportierten 2.5

simul_results_isol = list()
for (r in 1:nrow(params)) {
  simul_results_isol[[r]] = simul_epidemic(k = params[r, 1], m = params[r, 2], p = params[r, 3], isol = isol)
}


