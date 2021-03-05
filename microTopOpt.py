from fenics import *
from fenics_adjoint import *
import matplotlib.pyplot as plt

class Periodic2DBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool(near(x[0], 0) or near(x[1], 0)) and \
            not (bool(near(x[0], 1)) and near(x[1], 0)) and \
            not (bool(near(x[0], 0)) and near(x[1], 1)) and \
            not (bool(near(x[0], 1)) and near(x[1], 1)) and \
            on_boundary

    def map(self, x, y):
        if near(x[0], 1):
            y[0] = x[0] - 1
        else:
            y[0] = x[0]
        if near(x[1], 1):
            y[1] = x[1] - 1
        else:
            y[1] = x[1]

def eps(v):
    return sym(grad(v))

def sigma(v, Eps, phi):
    E0, E1, nu = 0.91, 0.0001, 0.3  # E0 is material, E1 void
    mu0 = E0/(2*(1+nu))
    mu1 = E1/(2*(1+nu))
    lmbda0 = E0*nu/((1+nu)*(1-2*nu))
    lmbda1 = E1*nu/((1+nu)*(1-2*nu))
    cmat0 = lmbda0 * tr(eps(v) + Eps) * Identity(2) + 2 * mu0 * (eps(v) + Eps)
    cmat1 = lmbda1 * tr(eps(v) + Eps) * Identity(2) + 2 * mu1 * (eps(v) + Eps)
    return phi**4*cmat0+(1-phi**4)*cmat1

def micro_elast(ii, phi):

    Ve = VectorElement("CG", mesh.ufl_cell(), 2)
    Re = VectorElement("R", mesh.ufl_cell(), 0)
    W = FunctionSpace(mesh, MixedElement([Ve, Re]), constrained_domain=Periodic2DBoundary())

    dv, dlamb = TrialFunctions(W)
    v_, lamb_ = TestFunctions(W)
    w = Function(W)

    if(ii == 11):  
        Eij = E11
    else:          
        Eij = E22

    F = inner(sigma(dv, Eij, phi), eps(v_)) * dx
    a, L = lhs(F), rhs(F)
    a += dot(lamb_, dv) * dx + dot(dlamb, v_) * dx

    solve(a == L, w)
    (u, lamb) = split(w)
    return u

def microTopOpt():

    # Target values and weights
    w1111, w1122, w2222    =   1,   30,   1
    AT1111, AT1122, AT2222 = 0.2, -0.1, 0.2

    # Parameters
    vc = 0.6
    GL_gamma = 0.00001
    GL_eps = 1
    niter = 50

    C = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=Periodic2DBoundary())

    initial_guess = Expression("(sin(2*pi*3*x[0]-0.5*pi)*sin(2*pi*3*x[1]-0.5*pi)+1)/2", degree=2)

    phi = interpolate(initial_guess, C)

    plt.figure(1)
    c = plot(phi, mode='color', vmin=0, vmax=1, cmap="coolwarm")
    plt.colorbar(c)

    u11 = micro_elast(11, phi)
    u22 = micro_elast(22, phi)

    allctrls = File("homogenized/allcontrols.pvd")
    rho_viz = Function(C)

    def eval_cb(j, phi):
        plt.figure(1)
        plot(phi, mode='color', vmin=0, vmax=1,
             cmap="coolwarm", title='Phasefield')
        plt.pause(0.01)
        rho_viz.assign(phi)
        allctrls << rho_viz

    J = 0.5*w1111*(assemble(inner(sigma(u11, E11, phi), eps(u11) + E11) * dx)-AT1111)**2  \
      + 0.5*w1122*(assemble(inner(sigma(u11, E11, phi), eps(u22) + E22) * dx)-AT1122)**2  \
      + 0.5*w2222*(assemble(inner(sigma(u22, E22, phi), eps(u22) + E22) * dx)-AT2222)**2  \
      + GL_gamma*assemble(GL_eps*dot(grad(phi), grad(phi)) * dx+0.25/GL_eps*(phi*phi-phi)**2*dx)

    m = Control(phi)              

    Jhat = ReducedFunctional(J, m, eval_cb_post=eval_cb)  

    lb = 0.0 
    ub = 1.0 

    volume_constraint = UFLInequalityConstraint((Constant(vc) - phi)*dx, m)

    problem = MinimizationProblem(Jhat, bounds=(lb, ub), constraints=volume_constraint)
    parameters = {"acceptable_tol": 1.0e-16,"maximum_iterations": niter, "print_level": 6}
    solver = IPOPTSolver(problem, parameters=parameters)
    phi_opt = solver.solve()

    u11_opt = micro_elast(11, phi_opt)
    u22_opt = micro_elast(22, phi_opt)

    A1111_opt = assemble(inner(sigma(u11_opt, E11, phi_opt), eps(u11_opt) + E11) * dx)
    A1122_opt = assemble(inner(sigma(u11_opt, E11, phi_opt), eps(u22_opt) + E22) * dx)
    A2222_opt = assemble(inner(sigma(u22_opt, E22, phi_opt), eps(u22_opt) + E22) * dx)

    print(f"Target tensor value: {AT1111} Final value: {A1111_opt:.3f}")
    print(f"Target tensor value: {AT1122} Final value: {A1122_opt:.3f}")
    print(f"Target tensor value: {AT2222} Final value: {A2222_opt:.3f}")

    print(f"First  target poisson ratio: {AT1122/AT2222} Final value: {A1122_opt/A2222_opt:.3f}")
    print(f"Second target poisson ratio: {AT1122/AT1111} Final value: {A1122_opt/A1111_opt:.3f}")

E11 = Constant(((1, 0), (0, 0))) 
E22 = Constant(((0, 0), (0, 1))) 

ndof = 50
mesh = UnitSquareMesh(ndof, ndof, "crossed")

microTopOpt()
# plt.show()