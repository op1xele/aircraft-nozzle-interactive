"""
Optimization Module
Automated nozzle parameter optimization for target performance.
Algorithms: Gradient descent, Nelder-Mead simplex, Genetic algorithm (simple)
"""

import numpy as np


class NozzleOptimizer:
    """
    Optimize nozzle parameters to meet target Isp, thrust, or efficiency.
    Works by calling the MOC engine repeatedly with different parameters.
    """

    def __init__(self, moc_engine, performance_calc, gamma=1.4):
        self.moc = moc_engine
        self.perf = performance_calc
        self.gamma = gamma
        self.history = []  # iteration history

    # ------------------------------------------------------------------
    def _evaluate(self, params, target, fixed):
        """
        Evaluate objective function for a given parameter set.
        params: [mach_exit, num_rays] or [mach_exit] depending on what's free
        target: dict with 'Isp', 'thrust', 'efficiency' targets
        fixed: dict with fixed parameters (p0, T0, throat_height, etc.)
        """
        try:
            mach_exit  = float(np.clip(params[0], 1.1, 8.0))
            num_rays   = int(np.clip(params[1], 10, 60)) if len(params) > 1 else fixed.get('num_rays', 30)
            throat_h   = fixed.get('throat_height', 1.0)

            result = self.moc.calculate_2d_nozzle(mach_exit, num_rays, throat_h)
            if result is None:
                return 1e9

            exit_h   = result['wall_y'][-1]
            throat_a = throat_h / 12.0
            exit_a   = exit_h / 12.0
            p0_psf   = fixed['p0'] * 144.0
            p_amb    = fixed.get('p_ambient', 14.7) * 144.0
            T0       = fixed['T0']

            perf = self.perf.calculate_all_performance(
                result, p0_psf, T0, p_amb, throat_a, exit_a)

            # Objective: weighted sum of target deviations
            obj = 0.0
            if 'Isp' in target:
                obj += ((perf['Isp'] - target['Isp']) / target['Isp'])**2 * target.get('w_Isp', 1.0)
            if 'thrust' in target:
                obj += ((perf['thrust'] - target['thrust']) / max(target['thrust'],1))**2 * target.get('w_thrust', 1.0)
            if 'efficiency' in target:
                obj += ((perf['efficiency'] - target['efficiency']) / 100.0)**2 * target.get('w_eff', 1.0)

            self.history.append({
                'mach_exit': mach_exit, 'num_rays': num_rays,
                'Isp': perf['Isp'], 'thrust': perf['thrust'],
                'efficiency': perf['efficiency'], 'obj': obj
            })
            return obj

        except Exception:
            return 1e9

    # ------------------------------------------------------------------
    def nelder_mead(self, x0, target, fixed, max_iter=50, tol=1e-4):
        """
        Nelder-Mead simplex optimization (no scipy needed).
        x0: initial guess [mach_exit, num_rays]
        """
        n = len(x0)
        # Build initial simplex
        simplex = [np.array(x0, dtype=float)]
        for i in range(n):
            x = np.array(x0, dtype=float)
            x[i] *= 1.1 if x[i] != 0 else 0.1
            simplex.append(x)

        f = [self._evaluate(x, target, fixed) for x in simplex]
        self.history = []

        alpha, gamma_nm, rho, sigma = 1.0, 2.0, 0.5, 0.5

        for iteration in range(max_iter):
            # Sort
            order = np.argsort(f)
            simplex = [simplex[i] for i in order]
            f = [f[i] for i in order]

            if f[0] < tol:
                break

            # Centroid (excluding worst)
            xo = np.mean(simplex[:-1], axis=0)

            # Reflection
            xr = xo + alpha*(xo - simplex[-1])
            fr = self._evaluate(xr, target, fixed)

            if f[0] <= fr < f[-2]:
                simplex[-1] = xr; f[-1] = fr
            elif fr < f[0]:
                # Expansion
                xe = xo + gamma_nm*(xr - xo)
                fe = self._evaluate(xe, target, fixed)
                if fe < fr:
                    simplex[-1] = xe; f[-1] = fe
                else:
                    simplex[-1] = xr; f[-1] = fr
            else:
                # Contraction
                xc = xo + rho*(simplex[-1] - xo)
                fc = self._evaluate(xc, target, fixed)
                if fc < f[-1]:
                    simplex[-1] = xc; f[-1] = fc
                else:
                    # Shrink
                    for i in range(1, n+1):
                        simplex[i] = simplex[0] + sigma*(simplex[i]-simplex[0])
                        f[i] = self._evaluate(simplex[i], target, fixed)

        best = simplex[0]
        return {
            'mach_exit':  float(best[0]),
            'num_rays':   int(best[1]) if len(best) > 1 else fixed.get('num_rays', 30),
            'obj_final':  float(f[0]),
            'iterations': iteration+1,
            'converged':  f[0] < tol,
            'history':    self.history,
        }

    # ------------------------------------------------------------------
    def genetic_algorithm(self, target, fixed, pop_size=20, generations=30):
        """
        Simple genetic algorithm for global search.
        """
        # Random initial population
        pop = np.column_stack([
            np.random.uniform(1.5, 6.0, pop_size),   # mach_exit
            np.random.randint(15, 50, pop_size).astype(float)  # num_rays
        ])
        self.history = []

        for gen in range(generations):
            fitness = np.array([self._evaluate(ind, target, fixed) for ind in pop])
            order   = np.argsort(fitness)
            pop     = pop[order]
            fitness = fitness[order]

            # Elitism: keep top 25%
            n_elite = max(2, pop_size//4)
            new_pop = list(pop[:n_elite])

            # Crossover + mutation to fill rest
            while len(new_pop) < pop_size:
                i, j = np.random.choice(n_elite, 2, replace=False)
                alpha_cross = np.random.random()
                child = alpha_cross*pop[i] + (1-alpha_cross)*pop[j]
                # Mutation
                if np.random.random() < 0.2:
                    child += np.random.normal(0, 0.1, len(child))
                new_pop.append(child)

            pop = np.array(new_pop)

        best = pop[0]
        best_fitness = self._evaluate(best, target, fixed)
        return {
            'mach_exit':  float(np.clip(best[0], 1.1, 8.0)),
            'num_rays':   int(np.clip(best[1], 10, 60)),
            'obj_final':  float(best_fitness),
            'generations': generations,
            'history':    self.history,
        }