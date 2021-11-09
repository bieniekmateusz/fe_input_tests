from simtk import openmm, unit
from openmmtools import testsystems
from openmmtools import alchemy
from openmmtools import integrators
from openmmtools import cache
from openmmtools import mcmc
from openmmtools import multistate
from simtk.openmm.app import GromacsGroFile, GromacsTopFile, PME, HBonds, Simulation, PDBReporter, StateDataReporter
import simtk.openmm as omm
# from simtk.openmm.app import *
# from simtk.openmm import *
# from simtk import *
from simtk import unit
from sys import stdout
import sys
import mdtraj
from openmmtools import states
import tempfile

pdb = omm.app.PDBFile('mol.pdb')
modeller = omm.app.Modeller(pdb.topology, pdb.positions)

# add virtual sites
forcefield = omm.app.ForceField('mol.xml')
modeller.addExtraParticles(forcefield)

system = forcefield.createSystem(modeller.topology, constraints=HBonds)

omm.app.PDBFile.writeFile(modeller.topology, modeller.positions, open("sys.pdb", "w"))

# Define the region of the System to be alchemically modified
mdtraj_top = mdtraj.Topology.from_openmm(modeller.topology)

# create the alchemical region
factory = alchemy.AbsoluteAlchemicalFactory(alchemical_pme_treatment='exact')
alchemical_region = alchemy.AlchemicalRegion(alchemical_atoms=[a.index for a in mdtraj_top.atoms])
print(f'Alchemical region: {alchemical_region}')
alchemical_system = factory.create_alchemical_system(system, alchemical_region)
alchemical_state = alchemy.AlchemicalState.from_system(alchemical_system)

# define the lambdas and the protocol
lambda_electrostatics = [1.00, 0.75, 0.50, 0.25] + [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
lambda_sterics  =  [1.00, 1.00, 1.00, 1.00, 1.00] + [0.95, 0.90, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2, 0.20, 0.1, 0.05, 0.00]
protocol = {'lambda_electrostatics': lambda_electrostatics, 
			'lambda_sterics': lambda_sterics}

# create the states for the replica exchange
compound_states = states.create_thermodynamic_state_protocol(
	alchemical_system, 
	protocol=protocol, 
	composable_states=[alchemical_state],
	constants={'temperature': 298.0 *unit.kelvin}
	)

# we define here 3k iterations for the entire RE sampler
# use the default Langevin for every step (1 ps, 2fm step, 5.0/ps collision rate)
simulation = multistate.ReplicaExchangeSampler(number_of_iterations=3000)

reporter = multistate.MultiStateReporter('vac.nc', checkpoint_interval=1)
# TODO add structure reporter
simulation.create(
	thermodynamic_states=compound_states, 
	sampler_states=states.SamplerState(modeller.positions, box_vectors=None),
	storage=reporter
	)

simulation.run()
print(simulation.iteration)


