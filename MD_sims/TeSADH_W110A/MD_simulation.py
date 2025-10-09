from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import parmed as pmd
from pdbfixer import PDBFixer
from mdtraj.reporters import HDF5Reporter
import mdtraj as mdt
import os
import numpy as np

boxtype = 'rectangular'
box_padding = 1.78
traj_folder = 'TeSADH_W110A'
production_steps = 50000000
gpu_index = '0'
sim_force_field = 'amber14-all.xml'
sim_water_model = 'amber14/tip3p.xml'
sim_gaff = 'gaff.xml'
prot = 'TeSADH_W110A.pdb'
ligand_names = ["ISO","DFX"]
ligand_xml_files = ["dfx.xml","isopropanol.xml"]
ligand_pdb_files = ["dfx_iso.pdb"]
additional_residue_definitions_file = "add_residue_def.xml"
sim_ph = 8
sim_temperature = 323.15
trajectory_out_atoms = 'protein or resn Zn or resname DFX or resn ISO'
trajectory_out_interval = 1000
restrained_eq_atoms = "chainid 0 and name CA or chainid 1 and name CA or chainid 2 and name CA or chainid 3 and name CA"

os.system(f'mkdir {traj_folder}')

xml_list = [sim_force_field, sim_water_model, sim_gaff] + ligand_xml_files
forcefield = app.ForceField(*xml_list)

template = PDBFixer(filename=prot)
template.findNonstandardResidues()
template.findMissingResidues()
template.findMissingAtoms()
template.addMissingAtoms()
app.PDBFile.writeFile(template.topology, template.positions, open(f"{traj_folder}/variant.pdb", 'w'), keepIds=True)

protein_pdb = app.PDBFile(f"{traj_folder}/variant.pdb")
if additional_residue_definitions_file:
    protein_pdb.topology.loadBondDefinitions(additional_residue_definitions_file)
protein_pdb.topology.createStandardBonds()

protein_mod = app.Modeller(protein_pdb.topology, protein_pdb.positions)
protein_mod.addHydrogens(forcefield, pH=sim_ph)

x_list = [pos[0]._value for pos in protein_mod.positions]
y_list = [pos[1]._value for pos in protein_mod.positions]
z_list = [pos[2]._value for pos in protein_mod.positions]

for lig_pdb_file in ligand_pdb_files:
    ligand_pdb = app.PDBFile(lig_pdb_file)
    protein_mod.add(ligand_pdb.topology, ligand_pdb.positions)

x_span = max(x_list) - min(x_list)
y_span = max(y_list) - min(y_list)
z_span = max(z_list) - min(z_list)

d = max(x_span, y_span, z_span) + (2 * box_padding)
d_x = x_span + (2 * box_padding)
d_y = y_span + (2 * box_padding)
d_z = z_span + (2 * box_padding)

prot_x_mid = min(x_list) + (0.5 * x_span)
prot_y_mid = min(y_list) + (0.5 * y_span)
prot_z_mid = min(z_list) + (0.5 * z_span)

box_x_mid = d_x * 0.5
box_y_mid = d_y * 0.5
box_z_mid = d_z * 0.5

shift_x = box_x_mid - prot_x_mid
shift_y = box_y_mid - prot_y_mid
shift_z = box_z_mid - prot_z_mid

solvated_model = app.Modeller(protein_mod.topology, protein_mod.positions)
for index in range(len(solvated_model.positions)):
    solvated_model.positions[index] = (solvated_model.positions[index][0]._value + shift_x, 
                                     solvated_model.positions[index][1]._value + shift_y, 
                                     solvated_model.positions[index][2]._value + shift_z) * unit.nanometers

if boxtype == 'cubic':
    solvated_model.addSolvent(forcefield, model='tip3p', neutralize=True, ionicStrength=0*unit.molar, 
                            boxVectors=(mm.Vec3(d, 0., 0.), mm.Vec3(0., d, 0.), mm.Vec3(0, 0, d)))
elif boxtype == 'rectangular':
    solvated_model.addSolvent(forcefield, model='tip3p', neutralize=True, ionicStrength=0*unit.molar, 
                            boxVectors=(mm.Vec3(d_x, 0., 0.), mm.Vec3(0., d_y, 0.), mm.Vec3(0, 0, d_z)))

selection_reference_topology = mdt.Topology().from_openmm(solvated_model.topology)
trajectory_out_indices = selection_reference_topology.select(trajectory_out_atoms)
restrained_eq_indices = selection_reference_topology.select(restrained_eq_atoms)

system = forcefield.createSystem(solvated_model.topology, nonbondedMethod=app.PME, 
                               nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, 
                               ewaldErrorTolerance=0.0005, rigidWater=True)

integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'Precision': 'single', 'DeviceIndex': gpu_index}
simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
simulation.context.setPositions(solvated_model.positions)
simulation.minimizeEnergy()

min_pos = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, min_pos, open(f'{traj_folder}/MIN.pdb', 'w'))
del simulation

force = mm.CustomExternalForce("(k/2)*periodicdistance(x, y, z, x0, y0, z0)^2")
force.addGlobalParameter("k", 50.0*unit.kilojoules_per_mole/unit.angstroms**2)
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")

for res_atom_index in restrained_eq_indices:
    force.addParticle(int(res_atom_index), min_pos[int(res_atom_index)].value_in_unit(unit.nanometers))
system.addForce(force)

integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
simulation.context.setPositions(min_pos)
simulation.context.setVelocitiesToTemperature(sim_temperature*unit.kelvin)

simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, 
                                                temperature=True, progress=True, remainingTime=True, 
                                                speed=True, totalSteps=7500000, separator='\t'))
simulation.reporters.append(HDF5Reporter(f'{traj_folder}/EQ_NVT.h5', 1000, atomSubset=trajectory_out_indices))
simulation.step(25000000)

state_npt_EQ = simulation.context.getState(getPositions=True, getVelocities=True)
simulation.saveCheckpoint(f"{traj_folder}/EQ1.chk")
positions = state_npt_EQ.getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(f'{traj_folder}/NVT_EQ.pdb', 'w'), keepIds=True)
del simulation

system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, sim_temperature*unit.kelvin, 25))
integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
simulation.context.setState(state_npt_EQ)

simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, 
                                                temperature=True, progress=True, remainingTime=True, 
                                                speed=True, totalSteps=7500000, separator='\t'))
simulation.reporters.append(HDF5Reporter(f'{traj_folder}/EQ_NPT.h5', 1000, atomSubset=trajectory_out_indices))
simulation.step(2500000)

state_npt_EQ = simulation.context.getState(getPositions=True, getVelocities=True)
simulation.saveCheckpoint(f"{traj_folder}/EQ2.chk")
positions = state_npt_EQ.getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(f'{traj_folder}/NPT_EQ.pdb', 'w'), keepIds=True)
del simulation

n_forces = len(system.getForces())
system.removeForce(n_forces-2)
integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
simulation.context.setState(state_npt_EQ)

simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, 
                                                temperature=True, progress=True, remainingTime=True, 
                                                speed=True, totalSteps=5000000, separator='\t'))
simulation.reporters.append(HDF5Reporter(f'{traj_folder}/EQ_NPT_free.h5', 1000, atomSubset=trajectory_out_indices))
simulation.step(2500000)

state_free_EQ = simulation.context.getState(getPositions=True, getVelocities=True)
simulation.saveCheckpoint(f"{traj_folder}/EQ3.chk")
positions = state_free_EQ.getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(f'{traj_folder}/free_NPT_EQ.pdb', 'w'), keepIds=True)
del simulation

for i in range(5):
    integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
    simulation.context.setState(state_free_EQ)
    simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, 
                                                    temperature=True, progress=True, remainingTime=True, 
                                                    speed=True, totalSteps=production_steps, separator='\t'))
    simulation.reporters.append(HDF5Reporter(f'{traj_folder}/production_{i}.h5', trajectory_out_interval, 
                                           atomSubset=trajectory_out_indices))
    simulation.step(production_steps)
    state_prod = simulation.context.getState(getPositions=True, getVelocities=True)
    positions = state_prod.getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(f'{traj_folder}/production_{i}.pdb', 'w'), keepIds=True)
    simulation.saveCheckpoint(f'{traj_folder}/production_{i}.chk')
    del simulation