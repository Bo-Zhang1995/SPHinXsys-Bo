/**
 * @file 	particle_relaxation.cpp
 * @brief 	This is the test of using levelset to generate body fitted particles (3D).
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Bo Zhang and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real particle_spacing_ref = 1.0 / 20.0;
Real BW = particle_spacing_ref * 10.0;
Real inner_circle_radius = 1.0;

Vec3d domain_lower_bound(-1.0 - BW, -1.0 - BW, -1.0 - BW);
Vec3d domain_upper_bound(1.0 + BW, 1.0 + BW, 1.0 + BW);

BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(5);

//----------------------------------------------------------------------
//	Define the body.
//----------------------------------------------------------------------
class SphereShape : public ComplexShape
{
public:
    explicit SphereShape(const std::string& shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_sphere(0.0, 0.0, 0.0);
        add<TriangleMeshShapeSphere>(inner_circle_radius, resolution, translation_sphere);
    }
};

int main(int ac, char* av[])
{
    /** Setup the system. */
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(sph_system);

    /** create a body with corresponding material, particles and reaction model. */
    SolidBody sphereshape(sph_system, makeShared<SphereShape>("SphereShape"));
    sphereshape.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    sphereshape.defineAdaptationRatios(0.8, 1.0);
    sphereshape.defineParticlesAndMaterial<>();
    sphereshape.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    /** body relation topology */
    InnerRelation sphere_inner(sphereshape);
    /**
     * @brief 	Methods used for particle relaxation.
     */
    /** Random reset the insert body particle position. */
    SimpleDynamics<RandomizeParticlePosition> random_sphere_particles(sphereshape);
    /** Write the body state to Vtp file. */
    BodyStatesRecordingToPlt write_sphere_to_plt(io_environment, sphereshape);
    /** Write the particle reload files. */
    ReloadParticleIO write_particle_reload_files(io_environment, sphereshape);
    /** Calculate the correction matrix. */
    InteractionWithUpdate<KernelCorrectionMatrixInnerWithLevelSet> correction_matrix(sphere_inner);
    sphereshape.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
    /** A Physics relaxation step . */
    //relax_dynamics::RelaxationStepInner relaxation_inner(sphere_inner, true);
    //relax_dynamics::RelaxationStepInnerImplicit relaxation_inner(sphere_inner, true);
    //relax_dynamics::RelaxationStepInner<CorrectionMatrixRelaxation> relaxation_inner(sphere_inner, true);
    relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_inner(sphere_inner, true);
    /** Update the kinetic energy for the relaxation. */
    SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_body_kinetic_energy(sphere_inner);
    ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_kinetic_energy(sphereshape, "ParticleKineticEnergy");
    ReduceDynamics<QuantityMaximum<Real>> calculate_body_maximum_kinetic_energy(sphereshape, "ParticleKineticEnergy");
    std::string filefullpath_error_analysis = io_environment.output_folder_ + "/" + "error_analysis.dat";
    std::ofstream out_file_error_analysis(filefullpath_error_analysis.c_str(), std::ios::app);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_sphere_particles.exec(0.25);
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    write_sphere_to_plt.writeToFile(0);
    //----------------------------------------------------------------------
    //	Relax particles of the insert body.
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    int ite = 0;

    Real body_average_kinetic_energy = 0.0;
    Real body_maximum_kinetic_energy = 0.0;

    GlobalStaticVariables::physical_time_ = ite;
    while (ite < 100000)
    {
        correction_matrix.exec();
        relaxation_inner.exec();
        ite++;

        if (ite % 500 == 0)
        {
            update_body_kinetic_energy.exec();
            body_average_kinetic_energy = calculate_body_average_kinetic_energy.exec();
            body_maximum_kinetic_energy = calculate_body_maximum_kinetic_energy.exec();
            std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
            out_file_error_analysis << std::fixed << std::setprecision(12) << ite << "   " << body_average_kinetic_energy <<
                "   " << body_maximum_kinetic_energy << "\n";
            write_sphere_to_plt.writeToFile(ite);
        }
    }
    ite++;
    write_sphere_to_plt.writeToFile(ite);
    write_particle_reload_files.writeToFile(0);
    TickCount t2 = TickCount::now();
    TickCount::interval_t tt;
    tt = t2 - t1;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    return 0;
}