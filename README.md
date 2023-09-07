General information about the CTL course available at https://ctl.polyphys.mat.ethz.ch/

# :wave: PROJECT SoftSolid

In this project you implement a known 2D model for a soft solid, the so-called elastic Lennard-Jones model [1], and explore its features. The model consists of a number of particles, whose positions are initially coinciding with the node positions of a regular *N* $\times$ *N* grid with grid constant *g*. The size of the simulation box is therefore *L* $\times$ *L* with *L=Ng*, and periodic boundary conditions apply. Each of the particles is permanently spring-connected to those four particles that are nearest neighbors in the initial configuration. If **r** is the relative vector between two spring-connected particles, and *k* the spring coefficient, the potential energy of the connecting harmonic spring is 

*E<sub>spring</sub>(r)=* $\frac{1}{2}$*kr<sup>2</sup>* 

where *r* is the length of **r**, i.e., *r=|**r**|*. Moreover, each particle interacts with all those particles whose centers are located within a circle of radius *r<sub>c</sub>* around its own position (note that *L>2r<sub>c</sub>* must hold). This interaction is captured by the radially symmetric Lennard-Jones potential

*E<sub>pair</sub>(r) = 4( r<sup>-12</sup> - r<sup>-6</sup> - r<sub>c</sub><sup>-12</sup> + r<sub>c</sub><sup>-6</sup> )*

which reaches zero at *r=r<sub>c</sub>* by construction. The contribution to the total Lennard-Jones energy from a single pair of particles separated by distance vector **r** is *E<sub>pair</sub>(r)*. As obvious from these expressions for energies, all quantities are considered dimensionless, while dimensions may also be introduced. 

Qualitatively, particles experience a long-range repulsion caused by the elastic springs to the original neighbors, and they attract each other due to the Lennard-Jones potential at small separation distances. This competition gives rise to a phase transition, if parameters are suitably chosen. The four model parameters are *g* (grid constant), *k* (spring coefficient), *T* (temperature), and integer-valued *N=L/g*. Typical parameters are: *N*=100 or smaller, *k=0.1*, *g=3.5*, and temperatues *T* in the range from 0.1 to 1.5.

The model just described is solved using so-called (overdamped) Brownian dynamics, which ignores the usual accelerations provided by Newton's equation, but provides displacements of all particles at each simulation (time) step, and at the same time ensures that the system operates at a constant temperature (which is not the case if Newton's equations were solved). If the time step is denoted by $\Delta$*t*, the new coordinates of particles are given in terms of the existing coordinates via

**x** = **x** +  **F**$\Delta t$ + $\sqrt{2 T \Delta t}$ **W**

In case you are browsing through publications, you will notice that we absorbed the friction coefficient into the time step, and that some people call it Langevin equation, others include acceleration terms. In every case the application of this displacement rule requires knowledge of all forces **F** on all particles, which is the main part of the algorithm. 


To this end one chooses a sufficiently small time step $\Delta$*t*=0.005 and discretizes Newton's equation of motion: acceleration = force to obtain new coordinates and velocities for each particle in the course of discretized time. The main ingredient for such a scheme is the calculation of the force on a given particle, for a given configuration (coordinates of all particles). Using the above expressions for energies, one has to first write down the corresponding forces analytically via force = negative gradient of the potential energy. Once all forces are known at time *t*, one can apply the so-called velocity Verlet scheme to obtain new coordinates and new velocities at the later time *t*+ $\Delta$*t*. This scheme does not only change the velocities of particles, but also their coordinates, and one has to make sure that periodic boundary conditions are respected. Afterwards, all velocities are rescaled by a common factor to ensure $\langle$*v*<sup>2</sup>$\rangle$/*2=T* holds at the new time *t*+ $\Delta$*t*, and one proceeds with the next time step.

The molecular dynamics methods thus generates particle trajectories (coordinates + velocities) for all particles, while all particles effectively stay within the central simulation box, or have at least a periodic image inside the central simulation box. 



<img src="https://www.complexfluids.ethz.ch/images/PROJECT-soft-solid.PNG" width="50%">

Movie 1: https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/Media1.mp4

Movie 2: https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/Media2.mp4


<img src="https://www.complexfluids.ethz.ch/images/soft-solid/Media1.gif" width="50%">

The following tasks are meant to become python-functions, with some input like the set of *g,k,T,N* and some output like *X,Y,L". The reason for creating such functions becomes most evident lateron in the molecular_dynamics function, where you may start to read. 

## Task 1: create_initial_configuration(*g,k,T,N*)

1. Place *N* $\times$ *N* nodes on a regular two-dimensional grid with lattice spacing *g* in a square box of size *L* $\times$ *L* where *L = gN*. The box is centered at the origin, i.e., the coordinates are all in [-*L*/2,*L*/2]. 
2. Collect the <em>x</em>- and <em>y</em>-components of the node positions in *N* $\times$ *N* arrays *X* and *Y*. The component *X*[*n*,*m*] then contains the *x*-component of the position of the node (*n,m*). In python it is most convenient to use *n,m* $\in$ {0,..,*N*-1}.
3. Return *X* and *Y* and *L*

## Task 2: create_initial_velocities(*N,T*)

This function creates and returns two *N* $\times$ *N* random velocity arrays (VX and VY). The element (*n,m*) of VX carries the *x*-component of the velocity of the node (*n,m*). Make sure that the mean(VX)=0 and mean(VY)=0, so that there is no drift of your systems center of mass. 

## Task 3: thermostat(*N,T*,VX,VY)

This function 
1. multiplies both VX and VY with a unique factor so that the kinetic energy is *N* $\times$ *T*. The kinetic energy is the sum of all squared velocity components.
2. Calculate and print the new kinetic energy and check that it now indeed equals *N* $\times$ *T*. This command can be removed once the test has been successfully passed. 
3. Calculates the mean of all velocity vectors meanVX, meanVY and afterwards sets VX -= meanVX and VY -= meanVY to make sure the mean velocity is zero.
2. returns the modified VX and VY.

## Task 4: fold(*L,x*)

Because coordinates may leave the central simulation box, one has to fold them back to the central simulation box 

     def fold(L,x):
        x -= - L*round(x/L)
        return x
        
## Task 5: get_dist(L,x1,y1,x2,y2)

This function calculates the true vector between two nodes at positions (x1,y1) and (x2,y2). Because the simulation box has periodic boundaries, the true distance vector and its length are 

    def get_dist(L,x1,y1,x2,y2):
        distance_vector = fold(L,[x1,y1] - [x2,y2])
        distance = norm(distance_vector)
        return distance_vector,distance
    
## Task 6: get_forces(N,L,X,Y)

This function calculates force components FX and FY for all nodes. The force on node (*n,m*) has two qualitatively different contributions
1. an elastic force from the nodes that have been direct neighbors of (*n,m*) directly after create_initial_configuration(). This means, that node (*n,m*) is spring-connected (permanently bonded) to nodes (*n*+1,*m*), (*n*-1,*m*), (*n*,*m*+1), and (*n*,*m*-1). Because of periodic boundary conditions, *n*+1 stands for 1, if *n=N*, and *n*-1 stands for *N*-1, if *n*=0. 
2. a Lennard-Jones-type repulsive force from all those nodes that are currently closer to (*n,m*) than distance 3.5. 

This is the difficult part. To calculate distances and distance vectors, always use get_dist. If the distance vector **d** between two bonded particles is known, the elastic force is **F** = -*k***d** on one of the two particles (and -**F** on the other), where *k* is the spring coefficient. The Lennard-Jones force between any pair of particles (two particles form a pair if their distance is &le; 3.5) is **F** = -&nabla;*V* where *V*=*V*(*r*) is the radially symmetric Lennard-Jones interaction potential *V*(*r*)=4(*r*<sup>-12</sup>-*r*<sup>-6</sup>). 

## Task 7: molecular_dynamics(*N,L,g,k,T*,MDsteps,*dt*)

Now we are ready to implement a full molecular dynamics with MDsteps steps, each of time duration $dt$. The overall structure looks like this, using the above functions. The integration scheme implemented here is known as the velocity-Verlet algorithm. 

    create_initial_configuration
    create_initial_velocities
    thermostat
    get_forces
    t = 0   # time

    # Loop over the following for MDsteps steps: 

    t += dt
    VX += dt*FX/2;  VY += dt*FY/2
    X += dt*VX;     Y += dt*VY
    get_forces
    VX = dt*FX/2;   VY += dt*FY/2
    thermostat
    
return the final t,coordinates, velocities, and forces. 

## Task 8: order_parameter(*N,L,X,Y*)

For each of the nodes check if it has at least one other node in its neighborhood, i.e., at a distance smaller 1.5. Count the number of nodes belonging to this species and divide by N<sup>2</sup> (the total number of nodes). This defines the so-called high density order parameter &Phi;. Having defined the order parameter
1. plot &Phi; versus time to check if &Phi; reaches a stationary value,
2. plot steady-state values of &Phi; versus temperature *T*,

using g=3.5, k=0.1, dt=0.005, and some N>10. 

[1] http://doi.org/10.1209/0295-5075/77/58007

