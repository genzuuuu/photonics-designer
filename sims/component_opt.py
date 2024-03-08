import time
import scipy
import np
import pyswarms as ps

def ScipyOptMin(component):
    t0 = time.time()
    res = scipy.optimize.minimize(component.fitness_function, (5.5,0.25),method='COBYLA', options={"maxiter": 100})
    component.insert_into_database()
    print(res)
    print(time.time() - t0)


def particleswarm(component, n_processes):
    max_bound = np.array([10, 10])
    min_bound = np.array([0.1,0.1])
    bounds = (min_bound, max_bound)

    # Set options for the PSO optimizer
    options = {"c1": 0.5, "c2": 0.3, "w": 0.9} #?

    # Create an instance of the PSO optimizer
    t0 = time.time()
    optimizer = ps.single.GlobalBestPSO(
        n_particles=10, dimensions=2, options=options, bounds=bounds
    )

    cost, pos = optimizer.optimize(component.fitness_function_swarm, iters=100, n_processes=n_processes)
    print(time.time() - t0) #print completion time

    print(cost)
    print(pos)

