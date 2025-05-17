from copy import copy
from multiprocessing import Manager, Pool
import time
from bacteria import bacteria
import numpy
import copy
from fastaReader import fastaReader
import random

if __name__ == "__main__":
    numeroDeBacterias = 4  
    iteraciones = 3      
    tumbo = 2  
    nado = 3
    secuencias = list()
    
    secuencias = fastaReader().seqs
    names = fastaReader().names

    numSec = len(secuencias)
    print("numSec: ", numSec)
    
    manager = Manager()
    poblacion = manager.list(range(numeroDeBacterias))
    names = manager.list(names)
    NFE = manager.list(range(numeroDeBacterias))

    # =========== INICIALIZACIÓN MEJORADA ===========
    # Más diversidad, pero pocos gaps al inicio (más chance de alineamientos buenos)
    def poblacionInicial():
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                sec = list(secuencias[j][:])
                # Solo 1 gap aleatorio por secuencia, fuerza diversidad pero NO saturación
                pos = random.randint(0, len(sec))
                sec = sec[:pos] + ['-'] + sec[pos:]
                bacterium.append(sec)
            poblacion[i] = list(bacterium)
    # -------- VERSIÓN ORIGINAL --------
    # def poblacionInicial():
    #     for i in range(numeroDeBacterias):
    #         bacterium = []
    #         for j in range(numSec):
    #             bacterium.append(secuencias[j])
    #         poblacion[i] = list(bacterium)
    # -----------------------------------

    def printPoblacion():
        for i in range(numeroDeBacterias):
            print(poblacion[i])

    # =========== GAP REMOVAL: Elimina columnas donde todas las secuencias son gaps ===========
    def eliminaGapsColumnasCompletas(poblacion):
        """
        Para cada bacteria, elimina columnas que son solo gaps ('-').
        Esto 'repara' alineamientos
        """
        for i in range(len(poblacion)):
            bacterium = list(poblacion[i])
            if not bacterium or not bacterium[0]:
                continue
            num_cols = len(bacterium[0])
            cols_to_remove = []
            for col in range(num_cols):
                if all(len(seq) > col and seq[col] == '-' for seq in bacterium):
                    cols_to_remove.append(col)
            # Elimina de atrás hacia adelante para no romper los índices
            for col in reversed(cols_to_remove):
                for seq in bacterium:
                    del seq[col]
            poblacion[i] = tuple(bacterium)
    # ========================================================================================

    operadorBacterial = bacteria(numeroDeBacterias)    
    veryBest = [None, None, None] #indice, fitness, secuencias
    
    # Registra el tiempo de inicio
    start_time = time.time()
    
    print("poblacion inicial ...")
    poblacionInicial() 

    for it in range(iteraciones):
        print(f"Iteración {it+1} de {iteraciones}")
        print("poblacion inicial creada - Tumbo ...")
        operadorBacterial.tumbo(numSec, poblacion, tumbo)    # Usar tu tumbo inteligente

        # ============ GAP REMOVAL EN CADA GENERACIÓN ============
        eliminaGapsColumnasCompletas(poblacion)
        # ========================================================

        print("Tumbo realizado - Cuadrando ...")
        operadorBacterial.cuadra(numSec, poblacion)
        print("poblacion inicial cuadrada - Creando granLista de Pares...")
        operadorBacterial.creaGranListaPares(poblacion)
        print("granList: creada - Evaluando Blosum Parallel")
        operadorBacterial.evaluaBlosum()  #paralelo
        print("blosum evaluado - creando Tablas Atract Parallel...")
        operadorBacterial.creaTablasAtractRepel(poblacion, 0.1, 0.002, 0.1, 0.001)
        operadorBacterial.creaTablaInteraction()
        print("tabla Interaction creada - creando tabla Fitness")
        operadorBacterial.creaTablaFitness()
        print("tabla Fitness creada ")
        bestIdx, bestFitness = operadorBacterial.obtieneBest(operadorBacterial.getNFE())
        if (veryBest[0] == None) or (bestFitness > veryBest[1]): #Remplaza el mejor 
            veryBest[0] = bestIdx
            veryBest[1] = bestFitness
            veryBest[2] = copy.deepcopy(poblacion[bestIdx])
        operadorBacterial.replaceWorst(poblacion, veryBest[0])
        operadorBacterial.resetListas(numeroDeBacterias)

    print("Very Best: ", veryBest)
    # Imprime el tiempo de ejecucion
    print("--- %s seconds ---" % (time.time() - start_time))
