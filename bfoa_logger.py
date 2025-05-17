import os
import csv
import time
import copy
import numpy as np

class BFOALogger:
    def __init__(self, num_runs=30, log_filename='bfoa_results.csv'):
        # Parámetros por defecto
        self.num_runs = num_runs
        self.log_filename = log_filename
        self.log_filepath = self._create_log_file()
        
        # Configuración de parámetros
        self.numeroDeBacterias = 4
        self.iteraciones = 3
        self.tumbo = 2
        self.dAttr = 0.1
        self.wAttr = 0.002
        self.hRep = 0.1
        self.wRep = 0.001
    
    def _create_log_file(self):
        """Crear archivo CSV para registro"""
        headers = [
            'Run', 
            'Total_Execution_Time', 
            'Best_Fitness', 
            'Best_BLOSUM_Score', 
            'Best_Interaction_Score', 
            'Global_NFE',
            'Parameters'
        ]
        
        os.makedirs('logs', exist_ok=True)
        filepath = os.path.join('logs', self.log_filename)
        
        with open(filepath, 'w', newline='') as csvfile:
            csv.writer(csvfile).writerow(headers)
        
        return filepath

    def _log_results(self, run, total_time, best_bacteria, veryBest, globalNFE):
        """Registrar resultados de cada corrida"""
        params = {
            'numeroDeBacterias': self.numeroDeBacterias,
            'iteraciones': self.iteraciones,
            'tumbo': self.tumbo,
            'dAttr': self.dAttr,
            'wAttr': self.wAttr,
            'hRep': self.hRep,
            'wRep': self.wRep
        }
        
        bestIdx, bestFitness = veryBest[0], veryBest[1]
        
        with open(self.log_filepath, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([
                run + 1,
                total_time,
                bestFitness,
                best_bacteria.blosumScore[bestIdx],
                best_bacteria.tablaInteraction[bestIdx],
                globalNFE,
                str(params)
            ])

    def _cargar_secuencias(self):
        """Cargar secuencias desde archivo FASTA"""
        with open('./multiFasta.fasta', 'r') as f:
            lineas = f.readlines()
        
        secuencias = []
        secuencia_actual = []
        
        for linea in lineas:
            if linea.startswith('>'):
                if secuencia_actual:
                    secuencias.append(''.join(secuencia_actual))
                    secuencia_actual = []
            else:
                secuencia_actual.append(linea.strip())
        
        if secuencia_actual:
            secuencias.append(''.join(secuencia_actual))
        
        return [list(seq) for seq in secuencias]

    def run_multiple_bfoa(self):
        """Ejecutar múltiples corridas de BFOA"""
        from bacteria import bacteria

        # Cargar secuencias
        secuencias = self._cargar_secuencias()
        print(f"Número de secuencias: {len(secuencias)}")
        print(f"Realizando {self.num_runs} corridas...")
        
        # Bucle de corridas
        for run in range(self.num_runs):
            print(f"\n--- Iniciando corrida {run + 1} de {self.num_runs} ---")
            
            # Inicialización de población
            poblacion = [list(bacterium) for bacterium in [secuencias] * self.numeroDeBacterias]
            
            # Inicialización del operador bacterial
            operadorBacterial = bacteria(self.numeroDeBacterias)
            veryBest = [None, None, None]
            globalNFE = 0
            
            # Tiempo de inicio
            start_time = time.time()
            
            # Bucle principal de optimización
            for it in range(self.iteraciones):
                operadorBacterial.tumbo(len(secuencias), poblacion, self.tumbo)
                operadorBacterial.cuadra(len(secuencias), poblacion)
                operadorBacterial.creaGranListaPares(poblacion)
                operadorBacterial.evaluaBlosum()
                operadorBacterial.creaTablasAtractRepel(poblacion, self.dAttr, self.wAttr, self.hRep, self.wRep)
                operadorBacterial.creaTablaInteraction()
                operadorBacterial.creaTablaFitness()
                
                globalNFE += operadorBacterial.getNFE()
                bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
                
                # Actualizar la mejor solución
                if (veryBest[0] is None) or (bestFitness > veryBest[1]):
                    veryBest[0] = bestIdx
                    veryBest[1] = bestFitness
                    veryBest[2] = copy.deepcopy(poblacion[bestIdx])
                
                # Reemplazar la peor solución con la mejor
                operadorBacterial.replaceWorst(poblacion, veryBest[0])
                operadorBacterial.resetListas(self.numeroDeBacterias)
            
            # Calcular tiempo total de ejecución
            total_time = time.time() - start_time
            
            # Registrar resultados
            self._log_results(run, total_time, operadorBacterial, veryBest, globalNFE)
            
            print(f"Corrida {run + 1} completada. Mejor Fitness: {veryBest[1]}")
            print(f"Tiempo de ejecución: {total_time} segundos")
        
        print(f"\nSe completaron {self.num_runs} corridas. Resultados guardados en {self.log_filepath}")

# Punto de entrada
if __name__ == "__main__":
    logger = BFOALogger(num_runs=30)
    logger.run_multiple_bfoa()