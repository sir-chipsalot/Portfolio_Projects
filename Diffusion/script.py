import numpy as np
import matplotlib.pyplot as plt

# forcing term f(x,t)
def f_func(x, t):
    return (np.pi**2 - 1) * np.exp(-t) * np.sin(np.pi * x)

# gaussian points for 2 point rule
gp = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
gw = np.array([1, 1])

# number of nodes
N = 11

# mesh on [0,1]
x_nodes = np.linspace(0, 1, N)
num_elems = N - 1
h = x_nodes[1] - x_nodes[0]

# global matrices
M = np.zeros((N, N))
K = np.zeros((N, N))

# assemble M and K
for e in range(num_elems):
    xL = x_nodes[e]
    xR = x_nodes[e+1]

    detJ = (xR - xL) / 2.0

    Me = np.zeros((2,2))
    Ke = np.zeros((2,2))

    for q in range(2):
        xi = gp[q]
        w = gw[q]

        phi = np.array([(1 - xi)/2, (1 + xi)/2])
        dphi_dxi = np.array([-0.5, 0.5])
        dphi_dx = dphi_dxi / detJ

        Me += w * detJ * np.outer(phi, phi)
        Ke += w * detJ * np.outer(dphi_dx, dphi_dx)

    inds = [e, e+1]
    for i in range(2):
        for j in range(2):
            M[inds[i], inds[j]] += Me[i,j]
            K[inds[i], inds[j]] += Ke[i,j]

# build F vector
def build_F(t):
    F = np.zeros(N)
    for e in range(num_elems):
        xL = x_nodes[e]
        xR = x_nodes[e+1]
        detJ = (xR - xL) / 2.0

        Fe = np.zeros(2)

        for q in range(2):
            xi = gp[q]
            w = gw[q]
            phi = np.array([(1 - xi)/2, (1 + xi)/2])
            x_q = detJ * xi + (xL + xR) / 2.0
            Fe += w * detJ * phi * f_func(x_q, t)

        inds = [e, e+1]
        for i in range(2):
            F[inds[i]] += Fe[i]

    return F

# initial condition
u0 = np.sin(np.pi * x_nodes)

# time settings
dt = 1.0 / 551.0
t_final = 1.0
steps = int(t_final / dt)

# explicit forward euler

# lumped mass matrix
M_lumped = np.sum(M, axis=1)
M_inv = 1.0 / M_lumped

u_exp = u0.copy()

for step in range(steps):
    t = step * dt
    F = build_F(t)
    rhs = F - K @ u_exp
    u_exp = u_exp + dt * (rhs * M_inv)

    # enforce BC
    u_exp[0] = 0
    u_exp[-1] = 0

# implicit backward euler

u_imp = u0.copy()

A = M + dt * K
A[0,:] = 0
A[-1,:] = 0
A[:,0] = 0
A[:,-1] = 0
A[0,0] = 1
A[-1,-1] = 1

A_inv = np.linalg.inv(A)

for step in range(steps):
    t = step * dt
    F = build_F(t)
    rhs = M @ u_imp + dt * F
    rhs[0] = 0
    rhs[-1] = 0
    u_imp = A_inv @ rhs

# exact solution
u_exact = np.exp(-1.0) * np.sin(np.pi * x_nodes)

# plot main comparison
# 1) Explicit Forward Euler only
plt.figure()
plt.plot(x_nodes, u_exp, label="explicit/Forward Euler (dt=1/551)")
plt.xlabel("x")
plt.ylabel("u(x,1)")
plt.title("Forward Euler Solution (N = 11)")
plt.grid()
plt.legend()
plt.show()

# 2) Implicit Backward Euler only
plt.figure()
plt.plot(x_nodes, u_imp, label="implicit/Backward Euler (dt=1/551)")
plt.xlabel("x")
plt.ylabel("u(x,1)")
plt.title("Backward Euler Solution (N = 11)")
plt.grid()
plt.legend()
plt.show()

# 3) Exact only
plt.figure()
plt.plot(x_nodes, u_exact, "--k", label="exact")
plt.xlabel("x")
plt.ylabel("u(x,1)")
plt.title("Exact Solution ( N = 11)")
plt.grid()
plt.legend()
plt.show()

print("As we can see the Forward Euler and Backward Euler all work and match the exact solution very well. I am very proud of myself!.\n")
print("Now onto analyzing.\n")

# explicit euler instability test
dt_list = [1/551, 0.001, 0.003, 0.005,0.0057]

plt.figure()

for dt_test in dt_list:
    u_test = u0.copy()
    steps_test = int(t_final / dt_test)
    M_lumped = np.sum(M, axis=1)
    M_inv = 1.0 / M_lumped

    for step in range(steps_test):
        t = step * dt_test
        F = build_F(t)
        rhs = F - K @ u_test
        u_test = u_test + dt_test * (rhs * M_inv)
        u_test[0] = 0
        u_test[-1] = 0

        if np.isnan(np.max(u_test)):
            break

    plt.plot(x_nodes, u_test, label="dt = " + str(dt_test))

plt.title("Explicit Euler Instability Test")
plt.xlabel("x")
plt.ylabel("u(x,1)")
plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

# explicit euler instability test pt2
dt_list = [1/551, 0.005,0.00535,0.0057,0.00571, 0.00575]

plt.figure()

for dt_test in dt_list:
    u_test = u0.copy()
    steps_test = int(t_final / dt_test)
    M_lumped = np.sum(M, axis=1)
    M_inv = 1.0 / M_lumped

    for step in range(steps_test):
        t = step * dt_test
        F = build_F(t)
        rhs = F - K @ u_test
        u_test = u_test + dt_test * (rhs * M_inv)
        u_test[0] = 0
        u_test[-1] = 0

        if np.isnan(np.max(u_test)):
            break

    plt.plot(x_nodes, u_test, label="dt = " + str(dt_test))

plt.title("Explicit Euler Instability Test (pt.2)")
plt.xlabel("x")
plt.ylabel("u(x,1)")
plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

print("CONTEXT: These are the graphs for the explicit Euler instability Test. The first figure is the Explicit Euler instability test while everything does not blow up, the second graph shows everything blowing up. I had to make different graphs because even having a difference of .0001, caused the graph to explode such that you can barely see the heat equation anymore.\n")

print("As we can see from the first graph, we have a smooth continous stable graph from dt's up until 0.0057, so much so that the lines overlap on eachother. From there we can look at the second graph for values above 0.0057, in these lines the start and end is the same ('0') which is enforced by our boundary conditions but line of the dt's greater than 0.0057, start oscillatiing indicating instability. From this, we can logically conclude that the critical dt value for this heat diffusion function is 0.0057, and for all dt's greater than 0.0057 there graphs will be unstable and thus unreliable.\n")

# forward euler solution for different N values
node_list = [3, 5, 8, 11,13]

plt.figure()

for Node in node_list:
    # build mesh
    x_nodes2 = np.linspace(0, 1, Node)
    num_elems2 = Node - 1
    h2 = x_nodes2[1] - x_nodes2[0]

    # new matrices for this N
    M2 = np.zeros((Node, Node))
    K2 = np.zeros((Node, Node))

    # assemble for this N
    for e in range(num_elems2):
        xL = x_nodes2[e]
        xR = x_nodes2[e+1]
        detJ = (xR - xL) / 2.0

        Me2 = np.zeros((2,2))
        Ke2 = np.zeros((2,2))

        for q in range(2):
            xi = gp[q]
            w = gw[q]
            phi = np.array([(1 - xi)/2, (1 + xi)/2])
            dphi_dxi = np.array([-0.5, 0.5])
            dphi_dx = dphi_dxi / detJ

            Me2 += w * detJ * np.outer(phi, phi)
            Ke2 += w * detJ * np.outer(dphi_dx, dphi_dx)

        inds = [e, e+1]
        for i in range(2):
            for j in range(2):
                M2[inds[i], inds[j]] += Me2[i,j]
                K2[inds[i], inds[j]] += Ke2[i,j]

    # initial condition for this N
    u2 = np.sin(np.pi * x_nodes2)

    # lumped mass for forward euler
    M2_lumped = np.sum(M2, axis=1)
    M2_inv = 1.0 / M2_lumped

    # explicit time stepping
    steps2 = int(t_final / dt)
    for step in range(steps2):
        t = step * dt

        # build F2
        F2 = np.zeros(Node)
        for e in range(num_elems2):
            xL = x_nodes2[e]
            xR = x_nodes2[e+1]
            detJ = (xR - xL) / 2.0

            Fe2 = np.zeros(2)
            for q in range(2):
                xi = gp[q]
                w = gw[q]
                phi = np.array([(1 - xi)/2, (1 + xi)/2])
                x_q = detJ * xi + (xL + xR) / 2.0
                Fe2 += w * detJ * phi * f_func(x_q, t)

            inds = [e, e+1]
            for i in range(2):
                F2[inds[i]] += Fe2[i]

        rhs2 = F2 - K2 @ u2
        u2 = u2 + dt * (rhs2 * M2_inv)

        # enforce BCs
        u2[0] = 0
        u2[-1] = 0

    plt.plot(x_nodes2, u2, label="N = " + str(Node))

plt.title("Forward Euler for different Numbers of Nodes")
plt.xlabel("x")
plt.ylabel("u(x,1)")
plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

print("In the plot above, I compared different values of N (nodes) to see how the graph of my heat equation changes with regards to different nodes. When I decrease N (so the mesh is coarser and each element is bigger), the solution becomes more jagged or 'stiff' with steeper and more abrupt turns.The peak near x=0.5 is also not captured well. Both the peak and the curve look too triangular, or just not match the smooth sinusoidal exact shape for small number of Nodes. As the amount of Nodes increase however, we see a 'smoothening' effect where the peaks are no longer as point and it follows the smooth heat diffusion we expect. From this we can see that The Number of nodes is sort of like the resolution of a camera, where each increase of nodes exponentially increases the resolution of our 'image' or graph of the function that we want to capture.\n")

dt_list_backward = [1/551, 0.001, 0.01, 0.05, 0.1,0.15, 0.2, 0.22, 0.3]

plt.figure()

for dt_tester in dt_list_backward:
    u_be = u0.copy()

    A = M + dt_tester * K
    A[0,:] = 0
    A[-1,:] = 0
    A[:,0] = 0
    A[:,-1] = 0
    A[0,0] = 1
    A[-1,-1] = 1

    A_inv = np.linalg.inv(A)

    steps_test = int(t_final / dt_tester)

    for step in range(steps_test):
        t = step * dt_tester
        F = build_F(t)
        rhs = M @ u_be + dt_tester * F
        rhs[0] = 0
        rhs[-1] = 0
        u_be = A_inv @ rhs

    plt.plot(x_nodes, u_be, label="dt = " + str(dt_tester))

plt.title("Backward Euler for Different delta t's")
plt.xlabel("x")
plt.ylabel("u(x,1)")
plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

print("In this graph, I graphed the behavior of Backward Euler at different dt's. From the graph that regardless of the dt I provided, the graph stays smooth, continous, and very much so models, at least in apperance, the correct graph. However, we can also see that as the dt increases, even though the curves are smooth, the larger dt's begin to overshoot the correct curve. In other words, the dt stays stable but the error gets larger the greater the dt.\n")

print("This result appears due to a combination of things, The graph stays stable because of the dampening effect of backward Euler. Because it is an energy sink, the graph will never blow up but this can be dangerous because the graph could return a stable but incorrect curve. Now onto the error. The error comes from the fact that as explained in notes earlier in class, backward euler is first order accurate in time which means that the bigger the jumps in dt between the less reliable data the backward euler function will use thus meaning the error will also increase. This can be seen especially when the dt is greater than the special step size. The spatial step size for the backward Euler heat equation is 0.1, derived from the 1/(N-1) equation where N is 11. After the critical dt of 0.1, derived from dt>=h, the curve remains stable but the error starts to increase exponentially.\n")

# test how the solution changes as the number of nodes changes for backward euler
node_list = [3, 5, 8, 11,13]

plt.figure()

for Node in node_list:

    # remake the mesh
    x_nodes2 = np.linspace(0, 1, Node)
    num_elems2 = Node - 1
    h2 = x_nodes2[1] - x_nodes2[0]

    # make new matrices
    M2 = np.zeros((Node, Node))
    K2 = np.zeros((Node, Node))

    # assemble for this N
    for e in range(num_elems2):
        xL = x_nodes2[e]
        xR = x_nodes2[e+1]
        detJ = (xR - xL) / 2.0

        Me2 = np.zeros((2,2))
        Ke2 = np.zeros((2,2))

        for q in range(2):
            xi = gp[q]
            w = gw[q]
            phi = np.array([(1 - xi)/2, (1 + xi)/2])
            dphi_dxi = np.array([-0.5, 0.5])
            dphi_dx = dphi_dxi / detJ
            Me2 += w * detJ * np.outer(phi, phi)
            Ke2 += w * detJ * np.outer(dphi_dx, dphi_dx)

        inds = [e, e+1]
        for i in range(2):
            for j in range(2):
                M2[inds[i], inds[j]] += Me2[i,j]
                K2[inds[i], inds[j]] += Ke2[i,j]

    # initial condition for this N
    u2 = np.sin(np.pi * x_nodes2)

    # implicit since it's always stable
    A2 = M2 + dt * K2
    A2[0,:] = 0
    A2[-1,:] = 0
    A2[:,0] = 0
    A2[:,-1] = 0
    A2[0,0] = 1
    A2[-1,-1] = 1

    A2_inv = np.linalg.inv(A2)

    # time stepping for this N
    steps2 = int(t_final / dt)
    for step in range(steps2):
        t = step * dt
        F2 = np.zeros(Node)
        # build F for this N
        for e in range(num_elems2):
            xL = x_nodes2[e]
            xR = x_nodes2[e+1]
            detJ = (xR - xL) / 2.0
            Fe2 = np.zeros(2)
            for q in range(2):
                xi = gp[q]
                w = gw[q]
                phi = np.array([(1 - xi)/2, (1 + xi)/2])
                x_q = detJ * xi + (xL + xR) / 2.0
                Fe2 += w * detJ * phi * f_func(x_q, t)
            inds = [e, e+1]
            for i in range(2):
                F2[inds[i]] += Fe2[i]

        rhs2 = M2 @ u2 + dt * F2
        rhs2[0] = 0
        rhs2[-1] = 0
        u2 = A2_inv @ rhs2

    plt.plot(x_nodes2, u2, label="N = " + str(Node))

plt.title("Backward Euler for Different Numbers of Nodes")
plt.xlabel("x")
plt.ylabel("u(x,1)")
plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

print("In the plot above, very similalry to the Forward Euler Node test, we see that as the amount of nodes increase, the 'resolution' of the graph becomes smoother, less jagged, and overall more correct. That being said, I would like to highlight an interesting difference I have seen in the graphs, for Forward Euler, at peaks of less 'resolution' nodes i.e lower nodes, we would see an overshoot at the peaks but in this graph we see an undershoot. This further demonstrates how the Forward Euler is an energy adder whereas the Backward Euler is an energy sinker and dampner. This was a very cool thing that I just happened to see.\n")
