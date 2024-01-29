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
#include <feel/feeldiscr/sensors.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/measure.hpp>
#include <feel/feelts/bdf.hpp>


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
    Feel::po::options_description options( "laplacian options" );
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

template <int Dim, int Order>
class Elastic
{
public:
    using mesh_t = Mesh<Simplex<Dim>>;
    using space_t = Pchv_type<mesh_t, Order>;
    using space_ptr_t = Pchv_ptrtype<mesh_t, Order>; // Define the type for Pchv_ptrtype
    using element_ = typename space_t::element_type;
    using form2_type = form2_t<space_t,space_t>; // Define the type for form2
    using form1_type = form1_t<space_t>; // Define the type for form1
    using bdf_ptrtype = std::shared_ptr<Bdf<space_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; // Define the type for exporter_ptrtype

    Elastic() = default;
    Elastic(nl::json const& specs);

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
    space_ptr_t Xh_;
    element_ u_, v_;
    form2_type a_, at_;
    form1_type l_, lt_;
    bdf_ptrtype bdf_;
    exporter_ptrtype e_;
    nl::json meas_;
    double E, nu, lambda, mu, rho;
    std::string F, G;
};

// Constructor
template <int Dim, int Order>
Elastic<Dim, Order>::Elastic(nl::json const& specs) : specs_(specs)
{
    initialize();
}

// Initialization
template <int Dim, int Order>
void Elastic<Dim, Order>::initialize()
{
    double H = specs_["/Meshes/elastic/Import/h"_json_pointer].get<double>();
    // Load mesh and initialize Xh, a, l, etc.
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/elastic/Import/filename"_json_pointer].get<std::string>(), _h = H);
    // define Xh on a marked region
    if ( specs_["/Spaces/elastic/Domain"_json_pointer].contains("marker") )
        Xh_ = Pchv<Order>(mesh_, markedelements(mesh_, specs_["/Spaces/elastic/Domain/marker"_json_pointer].get<std::vector<std::string>>()));
    // define Xh via a levelset phi where phi < 0 defines the Domain and phi = 0 the boundary
    else if (specs_["/Spaces/elastic/Domain"_json_pointer].contains("levelset"))
        Xh_ = Pchv<Order>(mesh_, elements(mesh_, expr(specs_["/Spaces/elastic/Domain/levelset"_json_pointer].get<std::string>())));
    // define Xh on the whole mesh
    else
        Xh_ = Pchv<Order>(mesh_); // Pchv : vectoriel

    u_ = Xh_->element();
    v_ = Xh_->element();

    a_ = form2( _test = Xh_, _trial = Xh_ );
    at_ = form2( _test = Xh_, _trial = Xh_ );
    l_ = form1( _test = Xh_ );
    lt_ = form1( _test = Xh_ );

    bool steady = get_value(specs_, "/TimeStepping/wave/steady", true);
    int time_order = get_value(specs_, "/TimeStepping/wave/order", 2);
    double initial_time = get_value(specs_, "/TimeStepping/wave/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/wave/end", 0.3);
    double time_step = get_value(specs_, "/TimeStepping/wave/step", 0.002);

    /////////////////////////////////////////////////////
    //  Parameters and Initializations for Elasticity  //
    /////////////////////////////////////////////////////
    E = get_value(specs_, "/Parameters/elastic/E/expr", 1.0e9); // Young modulus
    nu = get_value(specs_, "/Parameters/elastic/nu/expr", 0.3); // Poisson ratio
    rho = get_value(specs_, "/Parameters/elastic/rho/expr", 7800.0); // Density in kg.m^-3 (default : steel)
    lambda = E*nu/( (1+nu)*(1-2*nu) );
    mu = E/(2*(1+nu));
    F = specs_["/InitialConditions/elastic/externalF/Expression/Omega/expr"_json_pointer].get<std::string>();
    G = specs_["/BoundaryConditions/elastic/Gamma/g/expr"_json_pointer].get<std::string>();
    std::string DTU0 = specs_["/InitialConditions/elastic/displacement/Expression/Omega/expr"_json_pointer].get<std::string>();
    auto f = expr<FEELPP_DIM,1>(F);
    auto g = expr<FEELPP_DIM,1>(G);
    auto dtu0 = expr<FEELPP_DIM,1>(DTU0);
    auto u0_ = Xh_->element();
    u0_.zero();
    std::cout << "**** elasticity parameters initialized **** \n" << std::endl;


    /////////////////////////////////////////////////////
    //               Initial Condition                 //
    /////////////////////////////////////////////////////
    std::cout << "**** Initialize Dirac **** \n" << std::endl;
    node_type n(FEELPP_DIM);
    for (int i = 0 ; i < FEELPP_DIM ; i++){
        n(i) = 0.1;
    }
    auto s = std::make_shared<SensorPointwise<space_t>>(Xh_, n, "S");
    auto f_0 = form1( _test = Xh_, _vector = s->containerPtr() ); // contient la contribution du dirac

    /////////////////////////////////////////////////
    //         Test using different dirac          //
    /////////////////////////////////////////////////
    // // Initialize the Dirac as a smoothstep
    // std::cout << "**** Initialization of the Dirac **** \n" << std::endl;
    // auto f_0 = form1(_test = Xh_);
    // double diracMagnitude = get_value(specs_, "/Parameters/elastic/dirac_magintude/expr", 0.1);
    // auto dirac = vec(
    //     diracMagnitude * expr("smoothstep(x,0.01,0.02):x"),
    //     diracMagnitude * expr("smoothstep(y,0.01,0.02):y"));
    // std::cout << "**** Dirac initialized **** \n" << std::endl;
    // u0_.on(_range=elements(mesh_), _expr=dirac);

    double rho = 7800.; // kg.m^-3

    std::cout << "Dim : " << FEELPP_DIM << std::endl;
    std::cout << "Order : " << FEELPP_ORDER << std::endl;


    l_.zero();
    a_.zero();
    a_ += integrate( _range = elements(mesh_), _expr=trans(idt(u_))*id(v_)); //1
    // f_0 += integrate(_range = elements(mesh_), _expr = time_step * time_step / rho / 2 * id(v_)); // contribution du dirac évalué sur v
    f_0 += integrate(_range = elements(mesh_), _expr = -time_step * time_step / rho / 2 * lambda * inner(gradv(u0_),grad(v_))); //3
    f_0 += integrate(_range = elements(mesh_), _expr = -time_step * time_step / rho * mu * trace(sym(gradv(u0_)) * trans(sym(grad(v_))))); //4
    f_0 += integrate(_range = markedfaces(mesh_, "Gamma"), _expr = time_step * time_step / rho / 2 * inner(g, id(v_))); //5
    a_.solve( _rhs = f_0, _solution = u_ );


    ////////////////////////////////////////////////////
    //                Initialize BDF                  //
    ////////////////////////////////////////////////////

    bdf_ = Feel::bdf( _space = Xh_, _steady=steady, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _order=time_order );

    bdf_->start();
    if ( steady )
        bdf_->setSteady();

    bdf_->initialize( u0_ ); // set u0_
    bdf_->shiftRight( u_ ); // set u1_


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

    e_ = Feel::exporter(_mesh = mesh_);
}

// Process materials
template <int Dim, int Order>
void Elastic<Dim, Order>::processMaterials()
{
    auto def = sym(grad(u_));
    auto Id = eye<FEELPP_DIM,FEELPP_DIM>();
    for ( auto [key, material] : specs_["/Models/elastic/Materials"_json_pointer].items() )
    {
        LOG( INFO ) << fmt::format( "Material {} found", material );

        std::string matMu = fmt::format( "/Materials/{}/mu", material.get<std::string>() );
        auto MU = specs_[nl::json::json_pointer( matMu )].get<double>();
        std::string matLambda = fmt::format( "/Materials/{}/lambda", material.get<std::string>() );
        auto Lambda = specs_[nl::json::json_pointer( matLambda )].get<double>();

        auto sigma = Lambda*trace(def)*Id + 2*MU*def;
        a_ += integrate( _range = elements(mesh_), _expr = inner( sigma, grad(v_) ) );
        // a_ += integrate( _range=markedfaces(mesh_,"Gamma"), _expr=-inner(sigma*N(),id(v_)));
        a_ += integrate( _range=markedfaces(mesh_,"Gamma"), _expr=-inner(expr<FEELPP_DIM,1>(G),id(v_)));
    }
}

// Process boundary conditions
template <int Dim, int Order>
void Elastic<Dim, Order>::processBoundaryConditions()
{
    // BC Dirichlet (homogenous if g = 0)
    if ( specs_["/BoundaryConditions/elastic/Gamma/"_json_pointer].contains( "g" ) )
    {
        for ( auto& [bc, value] : specs_["/BoundaryConditions/elastic/Gamma/g"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "g {}: {}", bc, value.dump() );
            // auto flux = value["expr"].get<std::string>();

            // l_ += integrate( _range = markedfaces( support( Xh_ ), bc ),
            //         _expr = inner( expr<FEELPP_DIM,1>( G ), id( v_ )) );
            a_+=on(_range=markedfaces(mesh_,"Gamma"), _rhs=l_, _element=u_, _expr=expr<FEELPP_DIM,1>( G ) );
        }
    }

    // BC Neuman (Pure Traction Problem)
    if ( specs_["/BoundaryConditions/elastic/Gamma"_json_pointer].contains( "pure_traction" ) )
    {
        for ( auto& [bc, value] : specs_["/BoundaryConditions/elsatic/Gamma/pure_traction/g"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "pure traction {}: {}", bc, value.dump() );
            auto g = value["g"].get<std::string>();

            l_ += integrate( _range = markedfaces( support( Xh_ ), bc ),
                    _expr = inner( expr<FEELPP_DIM,1>( g ), id( v_ )) );
        }
    }
}

// Run method (main method to run Laplacian process)
template <int Dim, int Order>
void Elastic<Dim, Order>::run()
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

// Time loop
template <int Dim, int Order>
void Elastic<Dim, Order>::timeLoop()
{
    ////////////////////////////////////////////////////
    //          Linear Elasticity Time Loop           //
    ////////////////////////////////////////////////////
    auto f = expr<FEELPP_DIM,1>(F);
    auto g = expr<FEELPP_DIM,1>(G);
    double rho = 1.0;
    for ( bdf_->start(); bdf_->isFinished()==false; bdf_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << "time " << bdf_->time() << std::endl;
        // u_ = u_n+1
        auto un = bdf_->unknown(0); // un = u_n
        auto un_1 = bdf_->unknown(1);  // un_1 = u_{n-1}

        auto dt = expr(bdf_->timeStep());

        at_ += integrate( _range = elements(mesh_), _expr = trans(idt(u_))*id(v_) );

        lt_ += integrate( _range = elements(mesh_), _expr = dt*dt/rho * inner(f,id(v_))); //1
        lt_ += integrate( _range = elements(mesh_), _expr = -dt*dt/rho * lambda * inner(gradv(un),grad(v_))); //2
        lt_ += integrate( _range = elements(mesh_), _expr = -dt*dt/rho * 2 * mu * trace(sym(gradv(un)) * trans(sym(grad(v_)) ))); //3
        lt_ += integrate( _range = markedfaces(mesh_, "Gamma"), _expr = dt*dt/rho * inner(g, id(v_))); //4
        lt_ += integrate( _range = elements(mesh_), _expr = 2*inner(id(un),id(v_)) - inner(id(un_1),id(v_))); //5

        at_.solve( _rhs = lt_, _solution = u_ );

        this->exportResults();

        at_.zero();
        lt_.zero();
    }
}


// Export results
template <int Dim, int Order>
void Elastic<Dim, Order>::exportResults()
{
    e_->step(bdf_->time())->addRegions();
    e_->step(bdf_->time())->add("u", u_);
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

template <int Dim, int Order>
void Elastic<Dim, Order>::writeResultsToFile(const std::string& filename) const
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
template <int Dim, int Order>
void Elastic<Dim, Order>::summary(/*arguments*/) {
    /* ... summary code ... */
}

// Accessors and Mutators
/* ... */

} // namespace Feel
