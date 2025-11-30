import numpy as np
import matplotlib.pyplot as plt

# parameters / initial values

beta   = 24.0
delta_E = 0.6
mu_E   = 0.15
delta_J = 0.08
mu_J   = 0.05
alpha  = 0.003
omega  = 0.5
mu_A   = 0.1

t0    = 0.0
t_end = 365.0
h     = 0.01
y0    = np.array([10.0, 0.0, 0.0], dtype=float)

#algorithm 2 

def f_algorithm2(t, y, beta, delta_E, mu_E, delta_J, mu_J, alpha, omega, mu_A):

    """
    Implements Algorithm 2 as follows:
    function f(t_n, y_n, β, δ_E, μ_E, δ_J, μ_J, α, ω, μ_A)
    """
    E = y[0]  
    J = y[1]  
    A = y[2]  

    dE = beta * A - delta_E * E - mu_E * E
    dJ = delta_E * E - delta_J * J - alpha * J**2 - mu_J * J
    dA = omega * delta_J * J - mu_A * A

    dy = np.array([dE, dJ, dA], dtype=float)
    return dy

#algorithm 1 

def RK4(f, y_n, t_n, h):

    """
    Implements Algorithm 1 as follows:

    function RK4(f, y_n, t_n, h)
    """
    k1 = f(t_n,y_n)
    k2 = f(t_n + 0.5 * h,y_n + 0.5 * h * k1)
    k3 = f(t_n + 0.5 * h,y_n + 0.5 * h * k2)
    k4 = f(t_n + h,y_n + h * k3)

    y_next = y_n + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return y_next

#Time vector T as : T_i = t0 + h (i-1), i = 1,...,N+1

N_steps = int((t_end - t0) / h)            # = tend/h and length N_steps+1 from t0 to t_end
T = t0 + h * np.arange(N_steps + 1)       

#wrapping f_algorithm2 with fixed parameters so RK4 sees function f(t, y)

def f(t, y):
    return f_algorithm2(t, y,
                        beta, delta_E, mu_E,
                        delta_J, mu_J,
                        alpha, omega, mu_A)

#RK4 time-stepping loop and Y

Y = np.zeros((N_steps + 1, 3), dtype=float)
Y[0] = y0

for i in range(1, N_steps + 1):
    t_n = T[i - 1]   
    y_n = Y[i - 1]  
    Y[i] = RK4(f, y_n, t_n, h) #Y contains y0 and all RK4 outputs now

#code for result plotting

fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

#T vs E subplot, T vs J subplot, and T vs A subplot
axes[0].plot(T, Y[:, 0])
axes[0].set_ylabel("E")
axes[0].set_xlabel("T")

axes[1].plot(T, Y[:, 1])
axes[1].set_ylabel("J")
axes[1].set_xlabel("T")

axes[2].plot(T, Y[:, 2])
axes[2].set_ylabel("A")
axes[2].set_xlabel("T")

fig.suptitle("RK4 solution of the system", y=1.02)
plt.tight_layout()
plt.show()
