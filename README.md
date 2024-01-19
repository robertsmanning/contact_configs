# contact_configs
MATLAB code for computing energy-minimizing self-contact configurations of a Cosserat rod

Follows the model in "Energy-minimizing configurations for an elastic rod with self-contact
energy close to the inextensible-unshearable and hard-contact limits", by 
Robert Manning, Kathleen Hoffman, Michael Merkle, Li Fan, and Anubhav Sharma (under review
at CMAME as of January 2024)

Breakdown of files:

* Main programs are called twisted*.m.  Each solves a different problem, aligned with examples in paper cited above.  Otherwise .m files are called by these scripts

* find_minimum_homegrown.m is the energy-minimization code

* add_twist_to_z.m is used to generate initial guess from previous solution, if adding some imposed twist
* interpolate_z.m is used to generate initial guess from previous solution, if changing number of rod segments

* quatmult.m is quaternion multiplication (used in setting some boundary conditions)
 
* Four files named compute*.m compute terms that appear in the energy
* Three files named discrete*.m compute energy and/or its gradient and/or its Hessian
* five_energies_plus.m computes components of energy and some key features of a configuration

* make_H_posdef.m takes the approximate Hessian and finds a nearby positive definite matrix
