include("kernel.jl")

# experiment name (determines folder in runs/)
name = "test2"

# aspect ratio (alpha)
a = .2

# friction coefficient (rho)
r = 1.0

# domain boundaries
x0 = 0
x1 = 0.5
y0 = 0.1
y1 = 0.6

# depth and its derivatives
d = .02
g(x) = 1-exp.(-x.^2./(2*d^2))
gp(x) = x/d^2.*exp.(-x.^2./(2*d^2))
h(x,y) = g(x-x0).*g(x1-x).*g(y-y0).*g(y1-y)
hx(x,y) = (gp(x-x0).*g(x1-x)-g(x-x0).*gp(x1-x)).*g(y-y0).*g(y1-y)
hy(x,y) = g(x-x0).*g(x1-x).*(gp(y-y0).*g(y1-y)-g(y-y0).*gp(y1-y))

# diffusivity map
#   uniform diffusivity:
#     k(x,y,s) = 1.00e-2*ones(s)
#   bottom-intensified diffusivity:
#     k(x,y,s) = 1.00e-1*exp(-(s+1).*h(x,y)/.1)
k(x,y,s) = 1e-2*ones(s)

# meridional restoring profile
#   no restoring:
#     c(y) = zeros(y)
#   restoring in "Southern Ocean":
#     c(y) = (1-tanh((y+.5)/.1))/2
#   pick restoring constant 100x the average diffusivity
c(y) = zeros(y)

# surface BC buoyancy
b0(x,y) = cos.(pi*y)
#b0(x,y) = zeros(y)

# simulation length
T = 2.5

# number of grid points
nx = 100
ny = 100
ns = 20

# time step
dt = 2.5e-5

# set up model
m = ModelSetup(a, r, x0, x1, y0, y1, T, h, k, c, b0, nx, ny, ns, dt)

# initialize model state
s = ModelState(m)

# save model setup and initial state
path = @sprintf("runs/%s", name)
mkpath(path)
save(m, path)
save(s, path)
