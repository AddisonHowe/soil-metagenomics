"""Modeling metabolic dynamics

"""

import numpy as np


def timestep(rhs, y_n, params, dt):
    """RK4 timestep
    
    Args:
        rhs (callable): Right hand side derivative function dy/dt. Should return
            an array of length d.
        y_n (np.ndarray[float]): Current state. Length d.
        params (list): variable length list of parameters.
        dt (float): stepsize.
    
    Returns:
        (np.ndarray) array of length d, the updated state, using RK4.
    """
    k_1 = rhs(y_n, params)
    k_2 = rhs(y_n + dt*k_1/2, params)
    k_3 = rhs(y_n + dt*k_2/2, params)
    k_4 = rhs(y_n + dt*k_3, params)
    y_new = y_n + (dt/6) * (k_1 + 2*k_2 + 2*k_3 + k_4)
    return y_new


def simulate(rhs, y0, params, dt, T):
    """Simulate an ODE from an initial condition over a specified interval.
    
    Args:
        rhs (callable): Right hand side derivative function dy/dt. Should return
            an array of length d.
        y0 (np.ndarray[float]): Initial condition. Length d.
        params (list): variable length list of parameters.
        dt (float): stepsize.
        T (float): final time.
    
    Returns:
        (np.ndarray) Array of shape (d, N). The state history over the interval.
    """
    N = int(T/dt)
    d = len(y0)
    ys = np.zeros([d, N])
    ys[:,0] = y0
    y = ys[:,0]
    for i in range(N-1):
        y_new = timestep(rhs, y, params, dt)
        for x in y_new:
            x = max(x, 1e-5)
        ys[:, i + 1] = y_new
        y = y_new
    return ys


class Model:
    """Model class.

    Attributes:
        rhs: ODE defining the model.
        param_names: List of parameter names.
        infer_screen: boolean screen specifying which parameters are inferred.
        fixed_params: dictionary specifying values of fixed parameters.
    """

    def __init__(
            self, 
            rhs, 
            param_names,  
            fixed_params={},
    ):
        self.rhs = rhs
        self.param_names = np.array(param_names, dtype=str)
        self.set_fixed_parameters(fixed_params, quiet=True)

    def __str__(self):
        s = f"Model [num_params={len(self.param_names)}, num_inferred={self.num_inferred}]"
        s += f"\n\tFixed params: {self.fixed_params}"
        return s

    def simulate(self, y0, params, dt, tfin):
        p = self.fixed_params | params
        return simulate(self.rhs, y0, p, dt, tfin)
    
    def set_fixed_parameters(self, fixed_params, quiet=False):
        infer_params = {k for k in self.param_names if k not in fixed_params}
        self._set_infer_screen(infer_params, quiet=quiet)
        self._set_fixed_params(fixed_params, quiet=quiet)

    def get_inferred_param_names(self) -> list[str]:
        return list(self.param_names[self.infer_screen])

    def get_fixed_param_names(self) -> list[str]:
        return list(self.param_names[~self.infer_screen])
    
    def _set_infer_screen(self, infer_screen, quiet=False):
        # Initialize the screen to infer no parameters.
        self.infer_screen = np.zeros(len(self.param_names), dtype=bool)
        if infer_screen is None:
            # If None, infer all parameters.
            self.infer_screen[:] = True
        elif isinstance(infer_screen, dict):
            # If a dictionary, use value associated with parameter name.
            for i, k in enumerate(self.param_names):
                self.infer_screen[i] = infer_screen.get(k, 0)
        elif isinstance(infer_screen, set):
            # If a set, infer if parameter name is included.
            for i, k in enumerate(self.param_names):
                if k in infer_screen:
                    self.infer_screen[i] = True
        elif isinstance(infer_screen, (list, np.ndarray)):
            # If list or array, treat as a boolean screen.
            self.infer_screen = np.array(infer_screen, dtype=bool)
            if len(self.infer_screen) != len(self.param_names):
                msg = "List input for infer_screen must match length of param_names"
                raise RuntimeError(msg)
        self.num_inferred = np.sum(self.infer_screen)

    def _set_fixed_params(self, fixed_params: dict, quiet=False):
        self.fixed_params = {}
        if fixed_params is None:
            return
        elif isinstance(fixed_params, dict):
            # If a dictionary, store values
            for k, v in fixed_params.items():
                self.fixed_params[k] = v
        elif isinstance(fixed_params, (list, np.ndarray)):
            # If list or array, treat as corresponding to the order of 
            # NON-inferred parameters.
            if len(fixed_params) != np.sum(~self.infer_screen):
                msg = "Length of fixed_params must be equal to number of "
                msg += "NON-inferred parameters."
                raise RuntimeError(msg)
            counter = 0
            for i, p in enumerate(self.param_names):
                if not self.infer_screen[i]:
                    self.fixed_params[p] = fixed_params[counter]
        if not quiet:
            print("Fixed parameters:", self.fixed_params)
        