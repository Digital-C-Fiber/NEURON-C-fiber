import random
import numpy as np
from multiprocessing import Pool
import evaluate
import time
import csv
import math
import os
import copy
import time

gPump=0.0047891 
gNav17=0.10664 
gNav18=0.24271 
gNav19=9.4779e-05 
gKs=0.0069733 
gKf=0.012756 
gH=0.0025377 
gKdr=0.018002 
gKna=0.00042

#channel_scales=[0.005, 0.1, 0.2, 0.00009, 0.007, 0.01, 0.003, 0.02, 0.0004]
ions=[gPump, gNav17, gNav18, gNav19, gKs, gKf, gH, gKdr, gKna]

#particle class, initializes the particles
#dimension: number of ion channels to be optimized
class Particle:
    def __init__(self, dimension):
        
        #each channel is initialized with a random number in the boundaries set for the specific channel
        self.position = list(np.array([random.uniform(lower(ions[i]), upper(ions[i])) for i in range(dimension)]))
        
        #the velocity of each channel is set depending on the boundaries of the channel
        #self.velocity=np.array([random.uniform(-(upper(ions[i])-lower(ions[i])), upper(ions[i])-lower(ions[i])) for i in range(dimension)])
        self.velocity=np.array([0 for i in range(dimension)])
        
        self.best_position = self.position.copy()
        self.best_score = float('inf')

#get lower bound for channel
def lower(ion):
    return 0#because negative values don't make sense
    '''
    if ion==gKdr or ion==gNav18:
        return ion-2*ion
    else:
        return ion-4*ion
    '''
    #return ion-(0.5*ion)
    #return -4.5#beale
    #return -2*math.pi#mishra bird
    #return -10#holder
    #return -512#eggholder
    #return -100#schaffer, easom


#get upper bound for channel
def upper(ion):
    
    if ion==gKdr or ion==gNav18:
        return 2*ion
    else:
        return 4*ion
    
    #return ion+(0.5*ion)
    #return 4.5#beale
    #return 2*math.pi#mishra bird
    #return 10#holder
    #return 512#eggholder
    #return 100#schaffer, easom

#check if lower or upper bound is reached
def checkBound(pos, dimension):
    #print("###")
    #print("Particle position:"+str(pos))
    for i in range(dimension):
        if pos[i] < lower(ions[i]):
            #print("particle reached lower bound: "+str(pos[i]))
            pos[i]=lower(ions[i])
        elif pos[i] > upper(ions[i]):
            #print("particle reached upper bound: "+str(pos[i]))
            pos[i]=upper(ions[i])
        #pos[i]=round(pos[i],7)
    #print(pos)
    return list(pos)

def checkVelocity(velocity, dimension):
    for i in range(dimension):
        if abs(velocity[i]) > (upper(ions[i])-lower(ions[i]))*0.1:
            s=np.sign(velocity[i])
            velocity[i]=s*(upper(ions[i])-lower(ions[i]))*0.1
    return velocity
            

#update particle position
#particle: particle to be updated
#global_best_position: best position of all particles
#inertia_weight: weight for velocity
#cognitive_weight: weight for best position of particle
#social_weight: weight for best global position of all particles
#dimension: number of ion channels to be optimized
def update_particle(particle, global_best_position, randomParticle, inertia_weight, cognitive_weight, social_weight, c2, randPart_weight, c3, rand_velocity_weight, dimension, num_iterations, i, bestScore1, bestScore2):
    r1, r2, r3 = random.random(), random.random(), random.random()
    
    cognitive_component = np.array([(particle.best_position[i] - particle.position[i]) for i in range(dimension)])
    social_component = np.array([(global_best_position[i] - particle.position[i]) for i in range(dimension)])
    
    randPart_component = np.array([(randomParticle.best_position[i] - particle.position[i]) for i in range(dimension)])
    rand_velocity = np.array([random.uniform(-upper(ions[i])/2, upper(ions[i])/2) for i in range(dimension)])
    
    #inertia_weight is updated
    w_min=0.2
    w_max=1
    if i==0:
        inertia_weight=w_min
    else:
        inertia_weight = max(w_min, inertia_weight - ((w_max-w_min)/num_iterations)*i*(bestScore1-bestScore2))
        if inertia_weight==float('inf'):
            inertia_weight=w_max*i
        
        #print(bestScore1)
        #print(bestScore2)
        #print(inertia_weight - ((w_max-w_min)/num_iterations)*i*(bestScore1-bestScore2))
    
    #constriction factor
    psi = 2/abs(2-(c2+c3)**2-5*(c2+c3))
    
    particle.velocity = psi * (inertia_weight * particle.velocity 
                               + cognitive_weight * r1 * cognitive_component 
                               + social_weight * r2 * c2 * social_component 
                               + randPart_weight * r3 * c3 * randPart_component 
                               + rand_velocity_weight * rand_velocity)
    print("Inertia",inertia_weight * particle.velocity)
    print("Cognitve",cognitive_weight * r1 *cognitive_component)
    print("Social", social_weight * r2 * social_component)
    print("RandParticle", randPart_weight * r3 * c3 * randPart_component)
    print("RandVelocity", rand_velocity_weight * rand_velocity)
    print("Velocity:"+str(particle.velocity))
    particle.position += particle.velocity
    particle.position = checkBound(particle.position, dimension)
    
#return values from runParticle are stored here   
result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)

#main optimization method
#dimension: number of ion channels to be optimized
#num_particles: number of particles
#num_iterations: number of iterations, how many times the particles are updated
#inertia_weight: weight for velocity
#cognitive_weight: weight for best position of particle
#social_weight: weight for best global position of all particles
#name: name for saving particle positions and scores, to differentiate between different runs of the optimization
#prot: number of protocol to be run and evaluated
#runs the model for each particle in swarm in parallel, evaluates score of each particle, saves all particle positions and scores in file, updates particles
def particle_swarm_optimization(dimension, num_particles, num_iterations, inertia_weight, cognitive_weight, social_weight, c2, randPart_weight, c3, rand_velocity_weight, name, prot=-1):
    #start timer
    tic = time.perf_counter()
    
    #initialize particles
    swarm = [Particle(dimension) for _ in range(num_particles)]
    global_best_position = None
    global_best_score = float('inf')
    global_best_score_prev = float('inf')
    
    #create folder
    if not os.path.exists('Results/'+str(name)):
        os.mkdir('Results/'+str(name))
    
    #save particle position and score
    for i in range(num_particles):
        filename = 'Results/'+str(name)+'/'+str(name)+'particle'+str(i)+'.csv'
        #creates file, deletes content, if file already exists
        with open(filename,'w', newline='') as f:
            csv.writer(f).writerow(["Position", "Score"])
    #save best scores
    filenameBest = 'Results/'+str(name)+'/'+str(name)+'_bestParticles.csv'
    with open(filenameBest,'w', newline='') as f:
            csv.writer(f).writerow(["Position", "Score"])

    with Pool() as pool:
        for j in range(num_iterations):
            print('***************************************************** iteration:'+str(j))
            tic2 = time.perf_counter()
            results = []
            #global result_list
            #result_list=[]
            nr=0#identifier for particles
            for particle in swarm:
                print("Pos:"+str(particle.position))
                #run model for each particle in parallel
                args=copy.copy(particle.position)
                args.append(name)
                args.append(nr)
                nr=nr+1
                args.append(prot)
                #print(args)
                result=pool.apply_async(evaluate.runParticle, args=args, callback = log_result)
                #result=pool.apply_async(evaluate.runParticle, args=args)
                #result=pool.apply(evaluate.runParticle, args=args)#blocks until ready, not what i need
                results.append(result)
             
            
            #for result in results:
                #print(result)
                #result.get()
            
            while True:
                time.sleep(1)
                # catch exception if results are not ready yet
                try:
                    ready = [result.ready() for result in results]
                    successful = [result.successful() for result in results]
                except Exception:
                    continue
                # exit loop if all tasks returned success
                if all(successful):
                    break
                # raise exception reporting exceptions received from workers
                if all(ready) and not all(successful):
                    print(result_list)
                    for result in results:
                        print(result)
                    raise Exception(f'Workers raised following exceptions {[result._value for result in results if not result.successful()]}')
            
            #print("Result:"+str(result_list))
            toc2 = time.perf_counter()
            print(f"Simulation time Run: {(toc2 - tic2)/60:0.4f} min")
            tic3 = time.perf_counter()
            
            i=0
            nr=0
            for particle in swarm:
                #evaluate score of particle
                args=copy.copy(particle.position)
                args.append(name)
                args.append(nr)
                nr=nr+1
                args.append(prot)
                score=evaluate.evalParticle(args)
                
                #use score that is stored in results array
                '''
                for res in result_list:
                    if i==res[0]:#particle corresponding to result list
                        score=res[1]
                '''
                print("Score="+str(score))
                
                #save position and score
                filename = 'Results/'+str(name)+'/'+str(name)+'particle'+str(i)+'.csv'
                i=i+1
                with open(filename,'a', newline='') as f:
                    csv.writer(f).writerow([particle.position, score])
                
                #update scores
                if score < particle.best_score:
                    print("particle is better")
                    particle.best_position = particle.position.copy()
                    particle.best_score = score
                if score < global_best_score:
                    print("particle is best")
                    global_best_position = particle.position.copy()
                    global_best_score_prev = global_best_score
                    global_best_score = score
                    
            #save best scores 
            with open(filenameBest,'a', newline='') as f:
                csv.writer(f).writerow([global_best_position, global_best_score])
                    
            #stop optimization if score is good enough
            if global_best_score < 0.001:
                print("best score < 0.001: break")
                break
                
            #calculate next position of particle
            for particle in swarm:
                r1 = int(random.random()*len(swarm))
                update_particle(particle, global_best_position, swarm[r1], inertia_weight, cognitive_weight, social_weight, c2, randPart_weight, c3, rand_velocity_weight, dimension, num_iterations, j, global_best_score, global_best_score_prev)
            
            toc3 = time.perf_counter()
            print(f"Simulation time Score and update: {(toc3 - tic3)/60:0.4f} min")
            
    toc = time.perf_counter()
    print(f"Simulation time Full: {(toc - tic)/60:0.4f} min")
    
    return global_best_position, global_best_score

if __name__ == "__main__":
    # Beispielaufruf
    dimension = 9 #number of ion channels
    num_particles = 10
    num_iterations = 10
    lower_bound = 0
    upper_bound = 1
    inertia_weight = 0.7
    cognitive_weight = 1.5
    social_weight = 1.5

    best_position, best_score = particle_swarm_optimization(dimension, num_particles, num_iterations, lower_bound, upper_bound, inertia_weight, cognitive_weight, social_weight)

    print("Beste Position gefunden:", best_position)
    print("Beste Punktzahl:", best_score)
    