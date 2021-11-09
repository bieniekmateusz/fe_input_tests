from sys import stdout
import sys
import os

from simtk import unit
from simtk import openmm as omm
import openmmtools as ommt
import mdtraj

#cache.global_context_cache.platform = openmm.Platform.getPlatformByName('CUDA')

pdb = omm.app.PDBFile('mol.pdb')
modeller = omm.app.Modeller(pdb.topology, pdb.positions)
forcefield = omm.app.ForceField('mol.xml', 'mytip3p.xml') # tip3p tip4pew tip4pfb

# add virtual sites
modeller.addExtraParticles(forcefield)

modeller.addSolvent(forcefield, padding=1.6, model='tip3p')

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=omm.app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=omm.app.HBonds
)

omm.app.PDBFile.writeFile(modeller.topology, modeller.positions, open("sol.pdb", "w"))

# Define the region of the System to be alchemically modified
mdtraj_top = mdtraj.Topology.from_openmm(modeller.topology)
lig_atoms = mdtraj_top.select('resname mol')
print(f'Selected ligand: {lig_atoms}')

# create the alchemical region
factory = ommt.alchemy.AbsoluteAlchemicalFactory(alchemical_pme_treatment='exact')
alchemical_region = ommt.alchemy.AlchemicalRegion(alchemical_atoms=lig_atoms)
print(f'Alchemical region: {alchemical_region}')
alchemical_system = factory.create_alchemical_system(system, alchemical_region)
alchemical_state = ommt.alchemy.AlchemicalState.from_system(alchemical_system)

# define the lambdas and the protocol
lambda_electrostatics = [1.00, 0.75, 0.50, 0.25] + [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
lambda_sterics  =  		[1.00, 1.00, 1.00, 1.00] + [1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.35, 0.30, 0.25, 0.20, 0.20, 0.10, 0.05, 0.00]
protocol = {'lambda_electrostatics': lambda_electrostatics, 
			'lambda_sterics': lambda_sterics}

# create the states for the replica exchange
compound_states = ommt.states.create_thermodynamic_state_protocol(
	alchemical_system, 
	protocol=protocol, 
	composable_states=[alchemical_state],
	constants={'temperature': 298.0*unit.kelvin, 'pressure': 1.0*unit.atmosphere}
	)

# we define here 3k iterations for the entire RE sampler
# use the default Langevin for every step (1 ps, 2fm step, 5.0/ps collision rate)
simulation = ommt.multistate.ReplicaExchangeSampler(number_of_iterations=3000)

reporter = ommt.multistate.MultiStateReporter('sol.nc', checkpoint_interval=1)
# TODO add structure reporter
simulation.create(
	thermodynamic_states=compound_states, 
	sampler_states=ommt.states.SamplerState(
							modeller.positions, 
							box_vectors=modeller.topology.getPeriodicBoxVectors()),
	storage=reporter
	)

simulation.run()
print(simulation.iteration)


