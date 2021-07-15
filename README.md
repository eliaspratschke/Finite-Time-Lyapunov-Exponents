# Finite-Time-Lyapunov-Exponents

This is a short piece of code that computes the Finite-time Lypunov exponents (FTLE) for a system of inhomogenous, non-linear ODE's. The FTLE is often used to identify areas with different long-term behavious in the phase space of a system. The code features a Runge-Kutta scheme of order 4 and uses Eigens's Matrix class as the main workhorse. There is also an option for multithreading, since this task can be parallelized trivially. Included is also a short piece of Python code to plot the result. 
