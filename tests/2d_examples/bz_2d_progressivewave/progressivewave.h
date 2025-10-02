/**
 * @file 	owsc.h
 * @brief 	This is the case file for the test of Oscillating Wave Surge Converter (OWSC).
 * @author   Chi Zhang and Xiangyu Hu
 */
#ifndef TEST_2D_PROGRESSIVEWAVE_H
#define TEST_2D_PROGRESSIVEWAVE_H

#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real total_physical_time = 100.0;

Real DL = 60; // tank length
Real DL_FLAT = 50; // flat tank length
Real DL_SLOPE = 10; // slope tank length
Real DH = 1.5; // tank height

Real Water_H = 1.0; // water height
Real DL_extra = 1.0; // for wave maker
Real particle_spacing_ref = Water_H / 128; // particle spacing
Real BW = particle_spacing_ref * 4.0; //boundary width
BoundingBox system_domain_bounds(Vec2d(-DL_extra - BW, -BW), Vec2d(DL + BW, DH + BW));
Real gravity_g = 9.81; // gravity
// for material properties of the fluid
Real rho0_f = 1000.0;
Real U_f = 2.0 * sqrt(0.79 * gravity_g);
Real c_f = 10.0 * U_f;
Real mu_f = 1.0e-6;
//------------------------------------------------------------------------------
// geometric shape elements used in the case
//------------------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> waterblock;
    waterblock.push_back(Vec2d(0.0, 0.0));
    waterblock.push_back(Vec2d(0.0, Water_H));
    waterblock.push_back(Vec2d(DL, Water_H));
    waterblock.push_back(Vec2d(DL_FLAT, 0.0));
    waterblock.push_back(Vec2d(0.0, 0.0));

    return waterblock;
}

std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outerwallshape;
    outerwallshape.push_back(Vec2d(-DL_extra - BW, -BW));
    outerwallshape.push_back(Vec2d(-DL_extra - BW, DH));
    outerwallshape.push_back(Vec2d(DL + BW, DH));
    outerwallshape.push_back(Vec2d(DL + BW, Water_H-BW));
    outerwallshape.push_back(Vec2d(DL_FLAT, -BW));
    outerwallshape.push_back(Vec2d(-DL_extra - BW, -BW));

    return outerwallshape;
}

std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> innerwallshape;
    innerwallshape.push_back(Vec2d(-DL_extra, 0.0));
    innerwallshape.push_back(Vec2d(-DL_extra, DH + BW));
    innerwallshape.push_back(Vec2d(-BW, DH + BW));
    innerwallshape.push_back(Vec2d(-BW, 0.0));
    innerwallshape.push_back(Vec2d(-DL_extra, 0.0));

    return innerwallshape;
}

std::vector<Vecd> createInnerWallShape2()
{
    std::vector<Vecd> innerwallshape2;
    innerwallshape2.push_back(Vec2d(0.0, 0.0));
    innerwallshape2.push_back(Vec2d(0.0, DH));
    innerwallshape2.push_back(Vec2d(DL, DH));
    innerwallshape2.push_back(Vec2d(DL, Water_H));
    innerwallshape2.push_back(Vec2d(DL_FLAT, 0.0));
    innerwallshape2.push_back(Vec2d(0.0, 0.0));

    return innerwallshape2;
}

MultiPolygon createWaveMakerShape()
{
    std::vector<Vecd> wave_make_shape;
    wave_make_shape.push_back(Vec2d(-BW, 0.0));
    wave_make_shape.push_back(Vec2d(-BW, DH + BW));
    wave_make_shape.push_back(Vec2d(0.0, DH + BW));
    wave_make_shape.push_back(Vec2d(0.0, 0.0));
    wave_make_shape.push_back(Vec2d(-BW, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(wave_make_shape, ShapeBooleanOps::add);
    return multi_polygon;
}
//------------------------------------------------------------------------------
// geometric shapes for the bodies used in the case
//------------------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape2(), ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
//------------------------------------------------------------------------------
// Body parts used in the case
//------------------------------------------------------------------------------
class WaveMaking : public solid_dynamics::BaseMotionConstraint<BodyPartByParticle>
{
    Real model_scale_;
    Real gravity_;
    Real water_depth_;
    Real wave_height_;
    Real wave_period_;
    Real wave_freq_;
    Real wave_stroke_;

    Vecd getDisplacement(const Real& time)
    {
        Vecd displacement{ Vecd::Zero() };
        displacement[0] = 0.5 * wave_stroke_ * sin(wave_freq_ * time);
        return displacement;
    }

    Vec2d getVelocity(const Real& time)
    {
        Vec2d velocity{ Vecd::Zero() };
        velocity[0] = 0.5 * wave_stroke_ * wave_freq_ * cos(wave_freq_ * time);
        return velocity;
    }

    Vec2d getAcceleration(const Real& time)
    {
        Vec2d acceleration{ Vecd::Zero() };
        acceleration[0] = -0.5 * wave_stroke_ * wave_freq_ * wave_freq_ * sin(wave_freq_ * time);
        return acceleration;
    }

    void computeWaveStrokeAndFrequency()
    {
        Real scaled_wave_height = wave_height_ / model_scale_;
        Real scaled_wave_period = wave_period_ / sqrt(model_scale_);
        Real scaled_wave_freq = 2.0 * PI / scaled_wave_period;
        Real scaled_wave_amp = 0.5 * scaled_wave_height;

        int iterator = 20;
        Real Tol = 1.0e-6;
        Real scaled_wave_number = 1.0;
        for (int i = 1; i < iterator; i++)
        {
            Real term1 = tanh(scaled_wave_number * water_depth_);
            Real term2 = scaled_wave_freq * scaled_wave_freq / gravity_;
            Real term3 = scaled_wave_number * term1 - term2;
            Real term4 = term1 + scaled_wave_number * water_depth_ * (1.0 - term1 * term1);
            Real wave_number_old = scaled_wave_number;
            scaled_wave_number = wave_number_old - term3 / term4;
            Real error = abs(scaled_wave_number - wave_number_old) / abs(scaled_wave_number);
            if (error <= Tol)
                break;
        }

        Real term_1 = gravity_ / scaled_wave_freq / scaled_wave_freq;
        Real term_2 = 2.0 * scaled_wave_number * water_depth_;
        Real term_3 = scaled_wave_number * water_depth_;
        Real scaled_wave_stroke = 0.5 * scaled_wave_amp * scaled_wave_number * term_1 *
            (term_2 + sinh(term_2)) / (cosh(term_3) * sinh(term_3));

        wave_stroke_ = scaled_wave_stroke;
        wave_freq_ = scaled_wave_freq;
        std::cout << "Wave stroke: " << wave_stroke_ << " Wave frequency: " << wave_freq_ << std::endl;
    }

public:
    WaveMaking(BodyPartByParticle& body_part)
        : solid_dynamics::BaseMotionConstraint<BodyPartByParticle>(body_part),
        model_scale_(1.0), gravity_(gravity_g), water_depth_(Water_H), wave_height_(0.08),
        wave_period_(0.98)
    {
        computeWaveStrokeAndFrequency();
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        Real time = GlobalStaticVariables::physical_time_;
        pos_[index_i] = pos0_[index_i] + getDisplacement(time);
        vel_[index_i] = getVelocity(time);
        acc_[index_i] = getAcceleration(time);
    };
};

Real h = 1.3 * particle_spacing_ref;
MultiPolygon createWaveProbeShape(size_t index)
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(0.1 * index - h, 0.0));
    pnts.push_back(Vecd(0.1 * index - h, 1.5));
    pnts.push_back(Vecd(0.1 * index + h, 1.5));
    pnts.push_back(Vecd(0.1 * index + h, 0.0));
    pnts.push_back(Vecd(0.1 * index - h, 0.0));
    
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}
#endif // TEST_2D_OWSC_CASE_H
