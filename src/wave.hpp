//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
//! @date 2023-10-31
//! @copyright 2023 Feel++ Consortium
//! @copyright 2023 Université de Strasbourg
//!
#pragma once
#include <iostream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/minmax.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/measure.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/sensors.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include "wavelet.hpp"

namespace Feel
{
inline const int FEELPP_DIM=2;
inline const int FEELPP_ORDER=1;

static inline const bool do_print = true;
static inline const bool dont_print = false;

/**
 * @brief compute the summary of a container
 *
 * @tparam Container type of the container
 * @param c container
 * @param print boolean, true print the summary, false otherwise
 * @return nl::json json object containing the summary
 */
template<typename Container>
nl::json summary( Container const& c, bool print = do_print )
{
    using namespace Feel;
    using namespace Feel::vf;
    nl::json j;
    j["size"] = c.size();
    auto r = minmaxelt(_range = elements(support(c.functionSpace())), _element = c);
    j["min"] = r[0];
    j["max"] = r[1];
    j["mean"] = mean( _range = elements( c.mesh() ), _expr = idv( c ) );

    if (print)
    {
        if (Environment::isMasterRank())
            std::cout << j.dump(2) << std::endl;
    }
    return j;
}
inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "wave options" );
    options.add_options()

        // mesh parameters
        ( "specs", Feel::po::value<std::string>(),
          "json spec file for rht" )

        ( "steady", Feel::po::value<bool>()->default_value( 1 ),
          "if 1: steady else unsteady" );

    return options.add( Feel::feel_options() );
}

template<typename T>
T get_value(const nl::json& specs, const std::string& path, const T& default_value)
{
    auto json_pointer = nl::json::json_pointer(path);
    return specs.contains(json_pointer) ? specs[json_pointer].get<T>() : default_value;
}

template <int Dim, int Order, int OrderGeo=1>
class Wave
{
public:
    using mesh_t = Mesh<Simplex<Dim,Order,OrderGeo>>;
    using mesh_P1_t = Mesh<Simplex<Dim,1,1>>;
    using space_t = Pch_type<mesh_t, Order>;
    using space_ptr_t = Pch_ptrtype<mesh_t, Order>; // Define the type for Pch_ptrtype
    using space_P1_t = Pch_type<mesh_P1_t, 1>;
    using space_P1_ptr_t = Pch_ptrtype<mesh_P1_t, 1>;
    using element_ = typename space_P1_t::element_type;
    using form2_type = form2_t<space_t,space_t>; // Define the type for form2
    using form1_type = form1_t<space_t>; // Define the type for form1
    using bdf_ptrtype = std::shared_ptr<Bdf<space_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_P1_t>>; // Define the type for exporter_ptrtype
    using expr_type = typename Feel::vf::Expr<Feel::vf::GinacEx<2>>; // Define the type for expr_type
    using opI_type = I_t<space_t,space_P1_t>;

    Wave() = default;
    Wave(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    space_ptr_t const& Xh() const { return Xh_; }
    element_ const& u() const { return u_; }
    element_ const& v() const { return v_; }
    form2_type const& a() const { return a_; }
    form2_type const& at() const { return at_; }
    form1_type const& l() const { return l_; }
    form1_type const& lt() const { return lt_; }
    bdf_ptrtype const& bdf() const { return bdf_; }
    exporter_ptrtype const& exporter() const { return e_; }
    nl::json measures() const { return meas_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(element_ const& u) { u_ = u; }

    void initialize();
    void processMaterials();
    void processBoundaryConditions();
    void processWavelet(form1_type& l, form2_type& a,double t, int it);
    void run();
    void timeLoop();
    void exportResults();
    void summary(/*arguments*/);
    void writeResultsToFile(const std::string& filename) const;
    
    // Accessors and mutators for members
    /* ... */

private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
    std::shared_ptr<mesh_P1_t> mesh_P1_;
    space_ptr_t Xh_;
    space_P1_ptr_t Xh_P1_;
    element_ u_, v_;
    form2_type a_, at_;
    form1_type l_, lt_;
    bdf_ptrtype bdf_;
    exporter_ptrtype e_;
    nl::json meas_;
    expr_type mu, rho, s, g;
    node_type n;
    form1_type l_dirac_;
    bool dirac;
    bool gamma;
    bool abc;
    bool bat;
    bool circle;
    std::shared_ptr<opI_type> Ih_;
};

// Constructor
template <int Dim, int Order, int OrderGeo>
Wave<Dim,Order, OrderGeo>::Wave(nl::json const& specs) : specs_(specs)
{
    initialize();
}

// Initialization
template <int Dim, int Order,int OrderGeo>
void Wave<Dim, Order,OrderGeo>::initialize()
{
    double H = specs_["/Meshes/wave/Import/h"_json_pointer].get<double>();
    // Load mesh and initialize Xh, a, l, etc.
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/wave/Import/filename"_json_pointer].get<std::string>(), _h = H);
    if constexpr (OrderGeo==1)
    {
        mesh_P1_ = mesh_;
    }
    else
    {
        mesh_P1_ = loadMesh( _mesh = new mesh_P1_t, _filename = specs_["/Meshes/wave/Import/filename"_json_pointer].get<std::string>(), _h = H);
    }
    // define Xh on a marked region
    if ( specs_["/Spaces/wave/Domain"_json_pointer].contains("marker") )
    {
        Xh_ = Pch<Order>(mesh_, markedelements(mesh_, specs_["/Spaces/wave/Domain/marker"_json_pointer].get<std::vector<std::string>>()));
        Xh_P1_ = Pch<1>(mesh_P1_, markedelements(mesh_P1_, specs_["/Spaces/wave/Domain/marker"_json_pointer].get<std::vector<std::string>>()));
    }
    // define Xh via a levelset phi where phi < 0 defines the Domain and phi = 0 the boundary
    else if (specs_["/Spaces/wave/Domain"_json_pointer].contains("levelset"))
    {
        Xh_ = Pch<Order>(mesh_, elements(mesh_, expr(specs_["/Spaces/wave/Domain/levelset"_json_pointer].get<std::string>())));
        Xh_P1_ = Pch<1>(mesh_P1_, elements(mesh_P1_, expr(specs_["/Spaces/wave/Domain/levelset"_json_pointer].get<std::string>())));
    }
    // define Xh on the whole mesh
    else
    {
        Xh_ = Pch<Order>(mesh_);
        Xh_P1_ = Pch<1>(mesh_P1_);
    }

    Ih_ = I( _domain = Xh_, _image = Xh_P1_ );

    u_ = Xh_->element();
    v_ = Xh_->element();

    a_ = form2( _test = Xh_, _trial = Xh_ );
    at_ = form2( _test = Xh_, _trial = Xh_ );
    l_ = form1( _test = Xh_ );
    lt_ = form1( _test = Xh_ );

    bool steady = get_value(specs_, "/TimeStepping/wave/steady", true);
    int time_order = get_value(specs_, "/TimeStepping/wave/order", 2);
    double initial_time = get_value(specs_, "/TimeStepping/wave/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/wave/end", 4.0);
    double time_step = get_value(specs_, "/TimeStepping/wave/step", 0.002);

    // CFL condition
    ////////////////////////////////
    // C = max(C(x,y))
    // Later find a way to setup c expression in json file and compute C
    ////////////////////////////////
    double C = specs_["/Parameters/wave/c"_json_pointer].get<double>();
    time_step = std::min(time_step, H/C);

    bdf_ = Feel::bdf( _space = Xh_, _steady=steady, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _order=time_order );

    bdf_->start();
    if ( steady )
        bdf_->setSteady();

    // Initialize u0_ and u1_ with initial conditions
    auto u0_ = Xh_->element();
    u0_.on(_range = elements(mesh_), _expr = expr( specs_["/InitialConditions/wave/pressure/Expression/Omega/expr"_json_pointer].get<std::string>() ));
    auto w0_ = Xh_->element();
    w0_.on(_range = elements(mesh_), _expr = expr( specs_["/InitialConditions/wave/velocity/Expression/Omega/expr"_json_pointer].get<std::string>() ));
    // Parameters
    auto Mu = specs_["/Parameters/wave/mu"_json_pointer].get<std::string>();
    auto Rho = specs_["/Parameters/wave/rho"_json_pointer].get<std::string>();
    auto S = specs_["/Parameters/wave/s"_json_pointer].get<std::string>();
    auto G = specs_["/BoundaryConditions/wave/flux/Gamma/expr"_json_pointer].get<std::string>();
    mu = expr(Mu);
    rho = expr(Rho);
    s = expr(S);
    g = expr(G);

    std::cout << mu << std::endl;
    // Dirac
    if ( specs_["/Parameters/wave"_json_pointer].contains("dirac") )
    {
        auto coords = specs_["/Parameters/wave/dirac"_json_pointer].get<std::vector<double>>();
        n = node_type(coords.size());
        dirac = true;
        for (int i=0; i<coords.size(); i++)
            n(i) = coords[i];
        auto s_dirac = std::make_shared<SensorPointwise<space_t>>(Xh_, n, "S"); // "S" is the name of the sensor
        l_dirac_ = form1( _test = Xh_, _vector = s_dirac->containerPtr() );
    }
    else
        dirac = false;
    if (specs_["/Parameters/wave"_json_pointer].contains("abc"))
    {
        auto variety = specs_["/Parameters/wave/abc"_json_pointer].get<std::string>();
        if (variety == "circle")
        {
            circle = true;
        }
        else
            circle = false;
        abc = true;
    }
    else
        abc = false;

    if (specs_["/Parameters/wave"_json_pointer].contains("gamma"))
    {
        gamma = true;
    }
    else
        gamma = false;

    if (specs_["/Parameters/wave"_json_pointer].contains("bat"))
    {
        bat = true;
    }
    else
        bat = false;

    // Compute u1_
    a_.zero();
    l_.zero();
    a_ += integrate( _range = elements(mesh_), _expr = 1/mu * idt(u_) * id(v_) );
    l_ += integrate( _range = elements(mesh_),
            _expr = 1/mu * idv(u0_) * id(v_)
            + expr(bdf_->timeStep()) * 1/mu * idv(w0_) * id(v_)
            + expr(bdf_->timeStep()) * expr(bdf_->timeStep()) * s * id(v_) / 2
            + expr(bdf_->timeStep()) * expr(bdf_->timeStep()) * -1/mu * inner(gradv(u0_),gradv(v_)) /2);
    // add dirac
    if ( dirac )
    {
        auto tmp_l = l_dirac_;
        tmp_l.scale(bdf_->timeStep() * bdf_->timeStep() / 2);
        l_ += tmp_l;
    }
    l_ += integrate( _range = markedfaces(mesh_, "Gamma"), _expr = expr(bdf_->timeStep()) * expr(bdf_->timeStep()) * 1/rho * g * id(v_) / 2);
    a_.solve( _rhs = l_, _solution = u_ );

    // Initialize bdf
    bdf_->initialize( u0_ );
    bdf_->shiftRight( u_ );

    if ( steady )
        std::cout << "\n***** Compute Steady state *****" << std::endl;
    else
    {
        std::cout << "\n***** Compute Transient state *****" << std::endl;
        std::cout << "The step is  " << bdf_->timeStep() << "\n"
                  << "The initial time is " << bdf_->timeInitial() << "\n"
                  << "The final time is " << bdf_->timeFinal() << "\n"
                  << "BDF order :  " << bdf_->timeOrder() << "\n" << std::endl
                  << "BDF coeff :  " << bdf_->polyDerivCoefficient( 0 ) << "\n" << std::endl;
    }

    a_.zero();
    at_.zero();
    l_.zero();
    lt_.zero();

    e_ = Feel::exporter(_mesh = mesh_P1_);
}

// Process materials
template <int Dim, int Order, int OrderGeo>
void Wave<Dim,Order, OrderGeo>::processMaterials()
{
    for ( auto [key, material] : specs_["/Models/wave/Materials"_json_pointer].items() )
    {
        LOG( INFO ) << fmt::format( "Material {} found", material );
        std::string mat = fmt::format( "/Materials/{}/k", material.get<std::string>() );
        auto k = specs_[nl::json::json_pointer( mat )].get<std::string>();
        std::string matRho = fmt::format( "/Materials/{}/rho", material.get<std::string>() );
        auto Rho = specs_[nl::json::json_pointer( matRho )].get<std::string>();
        std::string matCp = fmt::format( "/Materials/{}/Cp", material.get<std::string>() );
        auto Cp = specs_[nl::json::json_pointer( matCp )].get<std::string>();

        a_ += integrate( _range = markedelements( support( Xh_ ), material.get<std::string>() ),
                _expr = bdf_->polyDerivCoefficient( 0 ) * expr( Rho ) * expr( Cp ) * idt( u_ ) * id( v_ ) + expr( k ) * gradt( u_ ) * trans( grad( v_ ) ) );
    }
}

// Process boundary conditions
template <int Dim, int Order, int OrderGeo>
void Wave<Dim,Order, OrderGeo>::processBoundaryConditions()
{
    // BC Neumann
    if ( specs_["/BoundaryConditions/wave"_json_pointer].contains( "flux" ) )
    {
        for ( auto& [bc, value] : specs_["/BoundaryConditions/wave/flux"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "flux {}: {}", bc, value.dump() );
            auto flux = value["expr"].get<std::string>();

            l_ += integrate( _range = markedfaces( support( Xh_ ), bc ),
                    _expr = expr( flux ) * id( v_ ) );
            if ( dirac )
                l_ += l_dirac_;
        }
    }

    // BC Robin
    if ( specs_["/BoundaryConditions/wave"_json_pointer].contains( "convective_wave_flux" ) )
    {
        for ( auto& [bc, value] : specs_["/BoundaryConditions/wave/convective_wave_flux"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "convective_wave_flux {}: {}", bc, value.dump() );
            auto h = value["h"].get<std::string>();
            auto Text = value["Text"].get<std::string>();

            a_ += integrate( _range = markedfaces( support( Xh_ ), bc ),
                    _expr = expr( h ) * id( v_ ) * idt( u_ ) );
            l_ += integrate( _range = markedfaces( support( Xh_ ), bc ),
                    _expr = expr( h ) * expr( Text ) * id( v_ ) );
        }
    }
}

// Run method (main method to run Wave process)
template <int Dim, int Order, int OrderGeo>
void Wave<Dim,Order, OrderGeo>::run()
{
    std::cout << "\n***** Initialize *****" << std::endl;
    initialize();
    std::cout << "\n***** Process materials *****" << std::endl;
    processMaterials();
    std::cout << "\n***** Process boundary conditions *****" << std::endl;
    processBoundaryConditions();
    std::cout << "\n***** Time loop *****" << std::endl;
    timeLoop();
    std::cout << "\n***** Export results *****" << std::endl;
    exportResults();
}

template <int Dim, int Order, int OrderGeo>
void Wave<Dim,Order,OrderGeo>::processWavelet(form1_type& l, form2_type& a, double t, int it)
{
    if ( specs_["/Models/wave"_json_pointer].contains("loading"))
    {
        for ( auto [key, loading] : specs_["/Models/wave/loading"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Loading {} found", key );
            std::string loadtype = fmt::format( "/Models/wave/loading/{}/type", key );
            if (specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Wavelet")
            {
                LOG( INFO ) << fmt::format( "Dirichlet conditions found" );
                std::string loadexpr = fmt::format( "/Models/wave/loading/{}/parameters/expr", key );
                // auto e = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
                double force = wavelet(t);
                std::cout << "force: " << force << std::endl;
                std::string e = "{";
                // for (int i = 0; i < FEELPP_DIM-1;i++)
                // {
                //     e.append("0,");
                // }
                e.append(std::to_string(wavelet(t)));
                e.append("}");
                std::cout << e << std::endl;
                // TODO: invert e and p, e is the expression of the force and p is the position of the force
                auto loadpos = fmt::format( "/Models/wave/loading/{}/parameters/location", key );
                std::vector<double> p = specs_[nl::json::json_pointer( loadpos )].get<std::vector<double>>();
                node_type n(p.size());
                for (int i = 0; i < p.size(); i++)
                    n(i) = p[i];
                auto dirac = std::make_shared<SensorPointwise<space_t>>(Xh_, n, key, e);
                auto f = form1( _test = Xh_, _vector = dirac->containerPtr() );
                l += f;
            }
        }
    }
}

// Time loop
template <int Dim, int Order, int OrderGeo>
void Wave<Dim,Order, OrderGeo>::timeLoop()
{
    int it = 0;
    // time loop
    for ( bdf_->start(); bdf_->isFinished()==false; bdf_->next(u_) )
    {
        // std::cout << "Time " << bdf_->time() << std::endl;

        at_ += integrate( _range = elements(mesh_), _expr = (1/mu) * idt(u_) * id(v_) );
        if (abc)
        {
            at_ += integrate( _range = markedfaces( mesh_, "ABC" ), _expr = expr(bdf_->timeStep()) * idt(u_) * id(v_) / 2.);
        }
        auto un = bdf_->unknown(0);
        auto un_1 = bdf_->unknown(1);
        lt_ += integrate( _range = elements(mesh_),
                          _expr = (1/mu) * (2 * idv(un) - idv(un_1) ) * id(v_)
                          + expr(bdf_->timeStep()) * expr(bdf_->timeStep()) * ((-1)/mu) * inner(gradv(un), grad(v_))
                          + expr(bdf_->timeStep()) * expr(bdf_->timeStep()) * s * id(v_));
        if (bat)
        {
            lt_ += integrate( _range = elements(mesh_,"Bat"),
                          _expr = expr(bdf_->timeStep()) * expr(bdf_->timeStep()) * dnv(un) * id(v_));
        }

        if (abc)
        {
            if (circle)
            {
                lt_ += integrate( _range = markedfaces( mesh_, "ABC" ),
                                  _expr = - expr(bdf_->timeStep()) * expr(bdf_->timeStep()) * idv(un) * id(v_) / 2.
                                  + expr(bdf_->timeStep()) * idv(un_1)*id(v_)/2.);
            }
            else
                lt_ += integrate( _range = markedfaces( mesh_, "ABC" ),
                                  _expr = - expr(bdf_->timeStep()) * idv(un_1) * id(v_) / 2.);
        }
        // add dirac
        if ( dirac )
        {
            auto tmp_l = l_dirac_;
            tmp_l.scale(bdf_->timeStep() * bdf_->timeStep() / 2);
            lt_ += tmp_l;
        }
        if (gamma)
        {
            lt_ += integrate( _range = markedfaces(mesh_, "Gamma"), _expr = expr(bdf_->timeStep()) * expr(bdf_->timeStep()) * (1/rho) * g * id(v_));
        }

        processWavelet(lt_, at_, bdf_->time(), it);
        it++;

        at_.solve( _rhs = lt_, _solution = u_ );

        this->exportResults();

        at_.zero();
        lt_.zero();
    }
}

// Export results
template <int Dim, int Order, int OrderGeo>
void Wave<Dim,Order, OrderGeo>::exportResults()
{
    // auto u_thin = Ih_->interpolate(u_);
    auto u_thin = Ih_(u_);
    e_->step(bdf_->time())->addRegions();
    e_->step(bdf_->time())->add("u", u_thin);
    e_->save();


    auto totalQuantity = integrate(_range=elements(mesh_), _expr=idv(u_)).evaluate()(0,0);
    auto totalFlux = integrate(_range=boundaryfaces(mesh_), _expr=gradv(u_)*N()).evaluate()(0,0);
    double meas=measure(_range=elements(mesh_), _expr=cst(1.0));
    meas_["time"].push_back(bdf_->time());
    meas_["totalQuantity"].push_back(totalQuantity);
    meas_["totalFlux"].push_back(totalFlux);
    meas_["mean"].push_back(totalQuantity/meas);
    meas_["min"].push_back(u_.min());
    meas_["max"].push_back(u_.max());
    for( auto [key,values] : mesh_->markerNames())
    {
        if ( values[1] == Dim )
        {
            double meas=measure(_range=markedelements(mesh_,key), _expr=cst(1.0));
            auto quantity = integrate(_range=markedelements(mesh_,key), _expr=idv(u_)).evaluate()(0,0);
            meas_[fmt::format("quantity_{}",key)].push_back(quantity);
            meas_[fmt::format("mean_{}",key)].push_back(quantity/meas);
        }
        else if ( values[1] == Dim-1 )
        {
            double meas=measure(_range=markedfaces(mesh_,key), _expr=cst(1.0));
            auto quantity = integrate(_range=markedfaces(mesh_,key), _expr=idv(u_)).evaluate()(0,0);
            meas_[fmt::format("quantity_{}",key)].push_back(quantity);
            meas_[fmt::format("mean_{}",key)].push_back(quantity/meas);
            auto flux = integrate(_range=markedfaces(mesh_,key), _expr=gradv(u_)*N()).evaluate()(0,0);
            meas_[fmt::format("flux_{}",key)].push_back(flux);
        }
        
    }
    
}
template <int Dim, int Order, int OrderGeo>
void Wave<Dim,Order, OrderGeo>::writeResultsToFile(const std::string& filename) const
{
    std::ofstream file(filename);
    if (file.is_open()) {
        file << meas_.dump(4);  // Indent of 4 spaces for readability
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

// Summary method
template <int Dim, int Order, int OrderGeo>
void Wave<Dim,Order, OrderGeo>::summary(/*arguments*/) {
    /* ... summary code ... */
}

// Accessors and Mutators
/* ... */

} // namespace Feel
