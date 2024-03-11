/**
 * @file 	particle_relaxation_single_resolution.cpp
 * @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Setting for the second geometry.
//	To use this, please commenting the setting for the first geometry.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/blade.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.06;
Real BW = 4 * resolution_ref;
Vec3d domain_lower_bound(2.0, -2.0, 4.5);
Vec3d domain_upper_bound(8.0, 2.0, 8.0);
Vec3d translation(-30, 0.0, 0.0);
Real scaling = 0.1;
//----------------------------------------------------------------------
//	Below are common parts for the two test geometries.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	define the imported model.
//----------------------------------------------------------------------
class SolidBodyFromMesh : public ComplexShape
{
  public:
    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
        //add<ExtrudeShape<TriangleMeshShapeSTL>>(2.0 * resolution_ref, full_path_to_file, translation, scaling);
        //subtract<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
    }
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody imported_model(sph_system, makeShared<SolidBodyFromMesh>("symbol"));
    // level set shape is used for particle relaxation
    imported_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(io_environment);
    imported_model.defineAdaptationRatios(1.3, 1.0);
    imported_model.defineParticlesAndMaterial();
    imported_model.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToPlt write_imported_model_to_plt(io_environment, { imported_model });
    MeshRecordingToPlt write_cell_linked_list(io_environment, imported_model.getCellLinkedList());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation imported_model_inner(imported_model);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    /** Calculate the correction matrix. */
    InteractionWithUpdate <KernelCorrectionMatrixInnerWithLevelSet> correction_matrix(imported_model_inner);
    imported_model.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
    SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(imported_model);
    /** A Physics relaxation step . */
    //relax_dynamics::RelaxationStepInner relaxation_inner(imported_model_inner, true);
    //relax_dynamics::RelaxationStepInnerImplicit relaxation_inner(imported_model_inner, true);
    relax_dynamics::RelaxationStepInner<CorrectionMatrixRelaxation> relaxation_inner(imported_model_inner, true);
    //relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_inner(imported_model_inner, true);
    /** Update the kinetic energy for the relaxation. */
    SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_body_kinetic_energy(imported_model_inner);
    ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_kinetic_energy(imported_model, "ParticleKineticEnergy");
    ReduceDynamics<QuantityMaximum<Real>> calculate_body_maximum_kinetic_energy(imported_model, "ParticleKineticEnergy");
    std::string filefullpath_error_analysis = io_environment.output_folder_ + "/" + "error_analysis.dat";
    std::ofstream out_file_error_analysis(filefullpath_error_analysis.c_str(), std::ios::app);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_imported_model_particles.exec(0.25);
    relaxation_inner.SurfaceBounding().exec();
    write_imported_model_to_plt.writeToFile(0.0);
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    write_cell_linked_list.writeToFile(0.0);
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    int ite = 0;

    Real body_average_kinetic_energy = 0.0;
    Real body_maximum_kinetic_energy = 0.0;

    GlobalStaticVariables::physical_time_ = ite;

    while (ite < 60000)
    {
        correction_matrix.exec();
        relaxation_inner.exec();
        ite += 1;
        if (ite % 500 == 0)
        {
            update_body_kinetic_energy.exec();
            body_average_kinetic_energy = calculate_body_average_kinetic_energy.exec();
            body_maximum_kinetic_energy = calculate_body_maximum_kinetic_energy.exec();
            std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
            out_file_error_analysis << std::fixed << std::setprecision(12) << ite << "   " << body_average_kinetic_energy <<
                "   " << body_maximum_kinetic_energy << "\n";
            write_imported_model_to_plt.writeToFile(ite);
        }
    }
    ite++;
    write_imported_model_to_plt.writeToFile(ite);
    TickCount t2 = TickCount::now();
    TickCount::interval_t tt;
    tt = t2 - t1;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    return 0;
}
