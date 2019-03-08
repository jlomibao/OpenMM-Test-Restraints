from simtk.openmm.app import *
from simtk.openmm import *
import simtk.openmm.openmm as mm
from simtk.unit import *
from sys import stdout
import parmed as pmd
import pytraj as pt
from mdtraj.reporters import NetCDFReporter

# Setup simulation
struct_dir = "cb6-but-dum"
old_topology = "cb6-but-dum.prmtop"
old_coordinates = "cb6-but-dum.rst7"

host = ":CB6"
guest = ":BUT"

H = [host+"@C7", host+"@C31", host+"@C19"]
G = [guest+"@C", guest+"@C3"]
C1 = [host+"@O", host+"@O6"]
C2 = [host+"@O2", host+"@O8"]
C3 = [host+"@O4", host+"@O10"]
D = [":DM1", ":DM2", ":DM3"]

# To store index values for atoms
H_i = [0,0,0]
G_i = [0,0]
C1_i = [0,0]
C2_i = [0,0]
C3_i = [0,0]
D_i = [0,0,0]
structure = pmd.load_file(
    struct_dir+'/'+old_topology,
    struct_dir+'/'+old_coordinates,
    structure=True
)

# There's probably a shorter way to get index of each atom
for atom in structure.atoms:
    if atom.residue.name == "CB6":
        if atom.name == "C7":
            H_i[0] = atom.idx
        elif atom.name == "C31":
            H_i[1] = atom.idx
        elif atom.name == "C19":
            H_i[2] = atom.idx
        elif atom.name == "O":
            C1_i[0] = atom.idx
        elif atom.name == "O6":
            C1_i[1] = atom.idx
        elif atom.name == "02":
            C2_i[0] = atom.idx
        elif atom.name == "O8":
            C2_i[1] = atom.idx
        elif atom.name == "O4":
            C3_i[0] = atom.idx
        elif atom.name == "O10":
            C3_i[1] = atom.idx
    elif atom.residue.name == "BUT":
        if atom.name == "C":
            G_i[0] = atom.idx
        elif atom.name == "C3":
            G_i[1] = atom.idx
    elif atom.residue.name == "DM1":
        D_i[0] = atom.idx
    elif atom.residue.name == "DM2":
        D_i[1] = atom.idx
    elif atom.residue.name == "DM3":
        D_i[2] = atom.idx

# Set mass of Dummy atoms to 0 so they are non-interacting
for i, atom in enumerate(structure.atoms):
    if atom.name == 'DUM':
        atom.mass=0.0

topology = "cb6-but-0m_dum.prmtop"
coordinates = "cb6-but-0m_dum.rst7"
structure.save(struct_dir+'/'+topology, overwrite=True)
structure.save(struct_dir+'/'+coordinates, overwrite=True)

prmtop = AmberPrmtopFile(struct_dir+'/'+topology)
inpcrd = AmberInpcrdFile(struct_dir+'/'+coordinates)

system = prmtop.createSystem(
    nonbondedMethod=NoCutoff,
    constraints=HBonds,
    implicitSolventSaltConc=0.1*moles/liter
)

integrator = LangevinIntegrator(
    300*kelvin,
    1/picosecond,
    0.002*picoseconds
)

# May need to add positional restraints to dummy atoms
'''
# Create Positional Restraints
pos_rest = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
pos_rest.addGlobalParameter("k",50.0*kilocalories_per_mole/angstroms**2) 
pos_rest.addPerParticleParameter("x0")
pos_rest.addPerParticleParameter("y0")
pos_rest.addPerParticleParameter("z0")
for i, atom in enumerate(structure.positions):
    if structure.atoms[i].name == 'DUM':
        pos_rest.addParticle(i, atom.value_in_unit(nanometers))
'''
# Distance restraint

# Use original prmtop and coordinates. Pytraj does not give correct distance to dummy atoms otherwise
traj = pt.load(struct_dir+'/'+old_coordinates,struct_dir+'/'+old_topology)
static_distance_rest = [D[0], H[0]]
static_init_dist = pt.distance(traj,D[0]+' '+H[0])[0]

dist_restraint = mm.CustomBondForce('k * (r - r_0)^2')
dist_restraint.addPerBondParameter('k')
dist_restraint.addPerBondParameter('r_0')

r_0 = static_init_dist * angstroms
k = 5 * kilojoules_per_mole / angstroms**2

dist_restraint.addBond(D_i[0],H_i[0], [k, r_0])

# Angle restraint 1
static_angle_rest1 = [D[1], D[0], H[0]]
static_init_angle1 = pt.angle(traj,D[1]+' '+D[0]+' '+H[0])[0]

angle_restraint1 = mm.CustomAngleForce('0.5*k*(theta-theta_0)^2')
angle_restraint1.addPerAngleParameter('k')
angle_restraint1.addPerAngleParameter('theta_0')

theta_0 = static_init_angle1 * degrees
k = 100 * kilojoules_per_mole / angstroms**2

angle_restraint1.addAngle(D_i[1],D_i[0],H_i[0], [k, theta_0])

# Dihedral restraint 1
static_dihedral_rest1 = [D[2], D[1], D[0], H[0]]
static_init_dihedral1 = pt.dihedral(traj,D[2]+' '+D[1]+' '+D[0]+' '+H[0])[0]

dihedral_restraint1 = mm.CustomTorsionForce('0.5*k*(theta-theta_0)^2')
dihedral_restraint1.addPerTorsionParameter('k')
dihedral_restraint1.addPerTorsionParameter('theta_0')

theta_0 = static_init_dihedral1 * degrees
k = 100 * kilojoules_per_mole / angstroms**2

dihedral_restraint1.addTorsion(D_i[2],D_i[1],D_i[0],H_i[0], [k, theta_0])

# Angle restraint 2
static_angle_rest2 = [D[0], H[0], H[1]]
static_init_angle2 = pt.angle(traj,D[0]+' '+H[0]+' '+H[1])[0]

angle_restraint2 = mm.CustomAngleForce('0.5*k*(theta-theta_0)^2')
angle_restraint2.addPerAngleParameter('k')
angle_restraint2.addPerAngleParameter('theta_0')

theta_0 = static_init_angle2 * degrees
k = 100 * kilojoules_per_mole / angstroms**2

angle_restraint2.addAngle(D_i[0],H_i[0],H_i[1], [k, theta_0])

# Dihedral restraint 2
static_dihedral_rest2 = [D[1], D[0], H[0], H[1]]
static_init_dihedral2 = pt.dihedral(traj,D[1]+' '+D[0]+' '+H[0]+' '+H[1])[0]

dihedral_restraint2 = mm.CustomTorsionForce('0.5*k*(theta-theta_0)^2')
dihedral_restraint2.addPerTorsionParameter('k')
dihedral_restraint2.addPerTorsionParameter('theta_0')

theta_0 = static_init_dihedral2 * degrees
k = 100 * kilojoules_per_mole / angstroms**2

dihedral_restraint2.addTorsion(D_i[1],D_i[0],H_i[0],H_i[1], [k, theta_0])

# Dihedral restraint 3
static_dihedral_rest3 = [D[0], H[0], H[1], H[2]]
static_init_dihedral3 = pt.dihedral(traj,D[0]+' '+H[0]+' '+H[1]+' '+H[2])[0]

dihedral_restraint3 = mm.CustomTorsionForce('0.5*k*(theta-theta_0)^2')
dihedral_restraint3.addPerTorsionParameter('k')
dihedral_restraint3.addPerTorsionParameter('theta_0')

theta_0 = static_init_dihedral3 * degrees
k = 100 * kilojoules_per_mole / angstroms**2

dihedral_restraint3.addTorsion(D_i[0],H_i[0],H_i[1],H_i[2], [k, theta_0])

# Add Forces to system
system.addForce(dist_restraint)
system.addForce(angle_restraint1)
system.addForce(dihedral_restraint1)
system.addForce(angle_restraint2)
system.addForce(dihedral_restraint2)
system.addForce(dihedral_restraint3)

simulation = Simulation(prmtop.topology, system, integrator, mm.Platform.getPlatformByName('Reference'))
simulation.context.setPositions(inpcrd.positions)

if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

simulation.minimizeEnergy()
simulation.reporters.append(NetCDFReporter('cb6-but-dum_output_6rest.nc', 500))
simulation.reporters.append(
    StateDataReporter(
        stdout, 500,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True
    )
)

simulation.step(100000)
