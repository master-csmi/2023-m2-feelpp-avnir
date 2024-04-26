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
#include <feel/feelts/newmark.hpp>
#include "wavelet.hpp"
#include <typeinfo>


namespace Feel
{
inline const int FEELPP_DIM=2;
inline const int FEELPP_ORDER=3;

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

template <typename ... Ts>
auto hoexporter( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && mesh = args.get(_mesh);
    bool fileset = args.get_else_invocable(_fileset,[](){ return boption(_name="exporter.fileset"); } );
    using mesh_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(mesh)>>>;
    using mesh_fragmentation_type = MeshFragmentation<mesh_type>;
    auto && meshFragmentation = args.get_else_invocable(_byparts,[](){ return mesh_fragmentation_type(boption(_name="exporter.byparts")? mesh_fragmentation_type::Strategy::AllMarkedElements : mesh_fragmentation_type::Strategy::None ); } );

    std::string const& name = args.get_else_invocable(_name,[](){ return Environment::about().appName(); } );
    std::string const& geo = args.get_else_invocable(_geo, [](){ return soption(_name="exporter.geometry"); } );
    auto && path = args.get_else_invocable(_path, [&name](){ return std::string((fs::path(Environment::exportsRepository())/fs::path(soption("exporter.format"))/name).string()); } );

    using exporter_type = Exporter<mesh_type,1>;

    auto e =  exporter_type::New( name,mesh->worldCommPtr() );
    e->setPrefix( name );
    e->setUseSingleTransientFile( fileset );
    e->setMeshFragmentation( meshFragmentation );
    if ( std::string(geo).compare("change_coords_only") == 0 )
        e->setMesh( mesh, EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY );
    else if ( std::string(geo).compare("change") == 0 )
        e->setMesh( mesh, EXPORTER_GEOMETRY_CHANGE );
    else // default
        e->setMesh( mesh, EXPORTER_GEOMETRY_STATIC );
    e->setPath( path );
    // addRegions not work with transient simulation!
    //e->addRegions();
    return e;
    //return Exporter<Mesh<Simplex<2> >,1>::New();
}


template <int Dim, int Order>
class Elastic
{
public:
    // using mesh_t = Mesh<Simplex<Dim,Order>>;
    using mesh_t = Mesh<Simplex<Dim>>;
    using space_t = Pchv_type<mesh_t, Order>;
    using space_ptr_t = Pchv_ptrtype<mesh_t, Order>; // Define the type for Pchv_ptrtype
    using element_t = typename space_t::element_type;
    using form2_type = form2_t<space_t,space_t>; // Define the type for form2
    using form1_type = form1_t<space_t>; // Define the type for form1
    using ts_ptrtype = std::shared_ptr<Newmark<space_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; // Define the type for exporter_ptrtype

    Elastic() = default;
    Elastic(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    std::shared_ptr<mesh_t> const& meshThin() const { return mesh_thin_; }
    space_ptr_t const& Xh() const { return Xh_; }
    element_t const& u() const { return u_; }
    element_t const& v() const { return v_; }
    form2_type const& a() const { return a_; }
    form2_type const& at() const { return at_; }
    form1_type const& l() const { return l_; }
    form1_type const& lt() const { return lt_; }
    ts_ptrtype const& bdf() const { return ts_; }
    exporter_ptrtype const& exporter() const { return e_; }
    nl::json measures() const { return meas_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setMeshThin(std::shared_ptr<mesh_t> const& mesh) { mesh_thin_ = mesh; }
    void setU(element_t const& u) { u_ = u; }

    void initialize();
    void processLoading(form1_type& l);
    void processMaterials(form2_type &a);
    void processBoundaryConditions(form1_type& l, form2_type& a,double t, int it);
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
    std::shared_ptr<mesh_t> mesh_thin_;
    space_ptr_t Xh_;
    element_t u_, v_, dtun;
    form2_type a_, at_;
    form1_type l_, lt_;
    ts_ptrtype ts_;
    exporter_ptrtype e_;
    nl::json meas_;
    double E, nu, lambda, mu, rho, beta, gamma;
    std::string F, G;
};

// Constructor
template <int Dim, int Order>
Elastic<Dim, Order>::Elastic(nl::json const& specs) : specs_(specs)
{
}

// Initialization
template <int Dim, int Order>
void Elastic<Dim, Order>::initialize()
{
    double H = specs_["/Meshes/elastic/Import/h"_json_pointer].get<double>();
    // Load mesh and initialize Xh, a, l, etc.
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/elastic/Import/filename"_json_pointer].get<std::string>(), _h = H);
    mesh_thin_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/elastic/Import/filename"_json_pointer].get<std::string>(), _h = H);
    // define Xh on a marked region
    if ( specs_["/Spaces/elastic/Domain"_json_pointer].contains("marker") )
        Xh_ = Pchv<Order>(mesh_, markedelements(mesh_, specs_["/Spaces/elastic/Domain/marker"_json_pointer].get<std::vector<std::string>>()));
    // define Xh via a levelset phi where phi < 0 defines the Domain and phi = 0 the boundary
    else if (specs_["/Spaces/elastic/Domain"_json_pointer].contains("levelset"))
        Xh_ = Pchv<Order>(mesh_, elements(mesh_, expr(specs_["/Spaces/elastic/Domain/levelset"_json_pointer].get<std::string>())));
    // define Xh on the whole mesh
    else
        Xh_ = Pchv<Order>(mesh_); // Pchv : vectoriel

    std::cout << "**** About the mesh ****" << std::endl;
    std::cout << "degree of freedom = " << Xh_->nDof() << std::endl;
    std::cout << "Nel = " << mesh_->numGlobalElements() << std::endl;
    std::cout << "Number of points = " << mesh_->numGlobalPoints() << std::endl;
    std::cout << "degree of freedom per points = " << Xh_->nDof()/(mesh_->numGlobalPoints()-1) << std::endl;

    u_ = Xh_->element();
    v_ = Xh_->element();
    dtun = Xh_->element();

    a_ = form2( _test = Xh_, _trial = Xh_ );
    at_ = form2( _test = Xh_, _trial = Xh_ );
    l_ = form1( _test = Xh_ );
    lt_ = form1( _test = Xh_ );

    bool steady = get_value(specs_, "/TimeStepping/elastic/steady", true);
    int time_order = get_value(specs_, "/TimeStepping/elastic/order", 2);
    double initial_time = get_value(specs_, "/TimeStepping/elastic/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/elastic/end", 0.3);
    double time_step = get_value(specs_, "/TimeStepping/elastic/step", 0.002);

    /////////////////////////////////////////////////////
    //  Parameters and Initializations for Elasticity  //
    /////////////////////////////////////////////////////
    E = get_value(specs_, "/Parameters/elastic/E/expr", 2.1e11); // Young's modulus
    nu = get_value(specs_, "/Parameters/elastic/nu/expr", 0.3); // Poisson ratio
    rho = get_value(specs_, "/Parameters/elastic/rho/expr", 7800.0); // Density in kg.m^-3 (default : steel)
    beta = get_value(specs_, "/Parameters/elastic/beta/expr", 0.25); // Newmark beta
    gamma = get_value(specs_, "/Parameters/elastic/gamma/expr", 0.5); // Newmark gamma
    lambda = E*nu/( (1+nu)*(1-2*nu) );
    mu = E/(2*(1+nu));
//    F = specs_["/InitialConditions/elastic/externalF/Expression/Omega/expr"_json_pointer].get<std::string>();
//    G = specs_["/InitialConditions/elastic/displacement/Expression/Omega/expr"_json_pointer].get<std::string>();
//    auto f = expr<FEELPP_DIM,1>(F);
//    auto g = expr<FEELPP_DIM,1>(G);
    auto u0_ = Xh_->element();
    u0_.zero();
    dtun.zero();
    std::cout << "**** elasticity parameters initialized **** \n" << std::endl;


    /////////////////////////////////////////////////////
    //               Initial Condition                 //
    /////////////////////////////////////////////////////
    std::cout << "**** Initialize Dirac **** \n" << std::endl;
    node_type n(FEELPP_DIM);
    n(0) = -1;
    for (int i = 1; i < FEELPP_DIM; i++)
        n(i) = 0.2312;


    std::cout << " n = " << n << "\n" << std::endl;
    // if (Dim == 2)
    // {
    //     auto s_dirac = std::make_shared<SensorPointwise<space_t>>(Xh_, n, "S", "{1e8,0}");
    //     auto f_0 = form1( _test = Xh_, _vector = s_dirac->containerPtr() ); // contient la contribution du dirac
    // }
    // else if (Dim == 3)
    // {
    //     auto s_dirac = std::make_shared<SensorPointwise<space_t>>(Xh_, n, "S", "{1e8,0,0}");
    //     auto f_0 = form1( _test = Xh_, _vector = s_dirac->containerPtr() ); // contient la contribution du dirac
    // }
    // else
    // {
    //     std::cerr << "Dimension not supported" << std::endl;
    //     return;
    // }

    /////////////////////////////////////////////////
    //         Test using different dirac          //
    /////////////////////////////////////////////////
    // Initialize the Dirac as a smoothstep

    // std::cout << "**** Initialization of the Dirac **** \n" << std::endl;
    // auto f_0 = form1(_test = Xh_);
    // double diracMagnitude = get_value(specs_, "/Parameters/elastic/dirac_magintude/expr", 0.1);
    // auto dirac = vec(
    //     diracMagnitude * expr("smoothstep(x,0.01,0.02):x") * expr("smoothstep(y,0.01,0.02):y"),
    //     diracMagnitude * expr("smoothstep(x,0.01,0.02):x") * expr("smoothstep(y,0.01,0.02):y"));
    // std::cout << "**** Dirac initialized **** \n" << std::endl;
    // u0_.on(_range=elements(mesh_), _expr=dirac);


    // auto localizedPerturbation = expr<FEELPP_DIM,1>("{1.,1.}*triangle(x,0.01,0.02)*triangle(y,0.01,0.02):x:y");
    // auto localizedPerturbation = expr("triangle(x,0.01,0.02)*triangle(y,0.01,0.02):x:y");
    // u0_.on(_range=elements(mesh_), _expr=localizedPerturbation);


    std::cout << "Dim : " << FEELPP_DIM << std::endl;
    std::cout << "Order : " << FEELPP_ORDER << std::endl;
    
    ts_ = newmark( _space = Xh_, _steady=steady, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _order=time_order );


    e_ = Feel::exporter(_mesh = mesh_, _name = "elasticity");

    ////////////////////////////////////////////////////
    //          Newmark beta-model for dttun          //
    ////////////////////////////////////////////////////
    l_.zero();
    a_.zero();
    
    ts_->start();
    if ( steady )
        ts_->setSteady();

    ts_->initialize( u0_ ); // set u0_
    u_ = u0_;
    //ts_->shiftRight( u_ ); // set u1_


    if ( steady )
        std::cout << "\n***** Compute Steady state *****" << std::endl;
    else
    {
        std::cout << "\n***** Compute Transient state *****" << std::endl;
        std::cout << "The step is  " << ts_->timeStep() << "\n"
                  << "The initial time is " << ts_->timeInitial() << "\n"
                  << "The final time is " << ts_->timeFinal() << "\n";
                  //<< "BDF order :  " << ts_->timeOrder() << "\n" << std::endl
                  //<< "BDF coeff :  " << ts_->polyDerivCoefficient( 0 ) << "\n" << std::endl;
    }

    at_.zero();
    lt_.zero();

    ts_->updateFromDisp(u_);
    auto u_proj = project(_space=Xh_,_range=elements(mesh_),_expr=idv(u_));
    e_->step(0)->add( "displacement", u_proj );
    e_->step(0)->add( "velocity", ts_->currentVelocity() );
    e_->step(0)->add( "acceleration", ts_->currentAcceleration() );
    e_->save();
    for (auto [key, sensor] : specs_["/Sensors"_json_pointer].items())
    {
        // Create a csv file for each sensor
        std::ofstream file;
        file.open(fmt::format( "{}.csv", key ));
        if (!file.is_open())
        {
            std::cerr << "Unable to open file : " << fmt::format( "{}.csv", key ) << std::endl;
            return;
        }
        file << "time,";
        for (int i = 0; i < FEELPP_DIM-1; i++)
        {
            // TODO: Use scientific notation and add the precision
            file << "u(" << i << "),";
        }
        file << "u(" << FEELPP_DIM-1 << ")";
        file << "\n";
        auto loadpos = fmt::format( "/Sensors/{}/location", key );
        std::vector<double> p = specs_[nl::json::json_pointer( loadpos )].get<std::vector<double>>();
        // std::cout << "p(0):" << p[0] << " p(1):" << p[1] << std::endl;
        Eigen::Matrix<double,Dim,1> __x;
        for (int i = 0; i < Dim; i++)
            __x(i) = p[i];
        file << 0 << ",";
        for (int i = 0; i < FEELPP_DIM-1; i++)
        {
            file << u_(__x)(i,0,0) << ",";
        }
        file << u_(__x)(FEELPP_DIM-1,0,0);
        file << "\n";
        file.close();
    }

}

template <int Dim, int Order>
void Elastic<Dim, Order>::processLoading(form1_type& l)
{
    if ( specs_["/Models/elastic"_json_pointer].contains("loading") )
    {
        for ( auto [key, loading] : specs_["/Models/elastic/loading"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Loading {} found", key );
            std::string loadtype = fmt::format( "/Models/elastic/loading/{}/type", key );
            
            if ( specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Dirac" )
            {
                LOG( INFO ) << fmt::format( "Loading {}: Dirac found", key );
                std::string loadexpr = fmt::format( "/Models/elastic/loading/{}/parameters/expr", key );
                auto e = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
                std::cout << fmt::format("Loading expr : {} ", e) << std::endl;
                auto loadpos = fmt::format( "/Models/elastic/loading/{}/parameters/location", key );
                std::vector<double> p = specs_[nl::json::json_pointer( loadpos )].get<std::vector<double>>();
                std::cout << fmt::format("Loading position : {} ", p) << std::endl;
                node_type n(p.size());
                for (int i = 0; i < p.size(); i++)
                    n(i) = p[i];
                auto dirac = std::make_shared<SensorPointwise<space_t>>(Xh_, n, key, e);
                auto f = form1( _test = Xh_, _vector = dirac->containerPtr() );
                l += f;
            }
            else if(specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Curve")
            {
                LOG(INFO) << fmt::format("Loading {}: Curve found", key);
                std::string loadexpr = fmt::format("/Models/elastic/loading/{}/parameters/expr", key);
                auto e = specs_[nl::json::json_pointer(loadexpr)].get<std::string>();
                std::cout << fmt::format("Loading expr : {} ", e) << std::endl;
                auto loadpos = fmt::format("/Models/elastic/loading/{}/parameters/location", key);
                // l+=on(_range=markedfaces(mesh_,loadpos),rhs=l,element=u_,expr=expr<FEELPP_DIM,1>(e));
                l_ += integrate(_range = markedfaces( mesh_ , loadpos), _expr = inner( expr<FEELPP_DIM,1>(e),id(v_)));
                // for ( auto [key, bc] : specs_["/BoundaryConditions/elastic/Dirichlet"_json_pointer].items() )
                // {
                //     LOG( INFO ) << fmt::format( "Dirichlet conditions found: {}", key );
                //     std::string e = fmt::format("/BoundaryConditions/elastic/Dirichlet/{}/g/expr",key);
                //     auto bc_dir = specs_[nl::json::json_pointer( e )].get<std::string>();
                //     std::cout << "BoundaryCondition Dirichlet : " << bc_dir << std::endl;
                //     a+=on(_range=markedfaces(mesh_,key), _rhs=l, _element=u_, _expr=expr<FEELPP_DIM,1>( bc_dir ) );
                // }
            }
            else if (specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Wavelet")
            {
                LOG( INFO ) << fmt::format( "Dirichlet conditions found" );
                std::string loadexpr = fmt::format( "/Models/elastic/loading/{}/parameters/expr", key );
                // auto e = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
                double force = wavelet(0);
                std::string e = "{";
                for (int i = 0; i < FEELPP_DIM-1;i++)
                {
                    e.append("0,");
                }
                e.append(std::to_string(wavelet(0)));
                e.append("}");
                // std::cout << e << std::endl;
                // TODO: invert e and p, e is the expression of the force and p is the position of the force
                auto loadpos = fmt::format( "/Models/elastic/loading/{}/parameters/location", key );
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

// Process materials
template <int Dim, int Order>
void Elastic<Dim, Order>::processMaterials( form2_type &a )
{
    for ( auto [key, material] : specs_["/Models/elastic/Materials"_json_pointer].items() )
    {
        LOG( INFO ) << fmt::format( "Material {} found", material.get<std::string>() );

        std::string matE = fmt::format( "/Materials/{}/parameters/E/value", material.get<std::string>() );
        auto E = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());
        std::string matNu = fmt::format( "/Materials/{}/parameters/nu/value", material.get<std::string>() );
        auto nu = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());
        std::string matRho = fmt::format( "/Materials/{}/parameters/rho/value", material.get<std::string>() );
        auto rho = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
        double lambda = E*nu/( (1+nu)*(1-2*nu) );
        double mu = E/(2*(1+nu));

        a += integrate( _range = markedelements(mesh_, material.get<std::string>()), 
                         _expr= rho*inner( ts_->polyDerivCoefficient()*idt(u_),id( v_ ) )
                              + lambda * divt(u_)*div(v_) + 2 * mu * trace(sym(gradt(u_)*trans(sym(grad(v_))))));
    }
}

// Process boundary conditions
template <int Dim, int Order>
void Elastic<Dim, Order>::processBoundaryConditions(form1_type& l, form2_type& a,double t, int it)
{
    // Boundary Condition Dirichlet
    if ( specs_["/BoundaryConditions/elastic"_json_pointer].contains("Dirichlet") )
    {
        LOG( INFO ) << fmt::format( "Dirichlet conditions found" );
        for ( auto [key, bc] : specs_["/BoundaryConditions/elastic/Dirichlet"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Dirichlet conditions found: {}", key );
            std::string e = fmt::format("/BoundaryConditions/elastic/Dirichlet/{}/g/expr",key);
            auto bc_dir = specs_[nl::json::json_pointer( e )].get<std::string>();
            // std::cout << "BoundaryCondition Dirichlet : " << bc_dir << std::endl;
            a+=on(_range=markedfaces(mesh_,key), _rhs=l, _element=u_, _expr=expr<FEELPP_DIM,1>( bc_dir ) );
        }
    }

    if ( specs_["/Models/elastic"_json_pointer].contains("loading") )
    {
        for ( auto [key, loading] : specs_["/Models/elastic/loading"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Loading {} found", key );
            std::string loadtype = fmt::format( "/Models/elastic/loading/{}/type", key );
            if (specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Dirac")
            {
                LOG( INFO ) << fmt::format( "Dirichlet conditions found" );
                std::string loadexpr = fmt::format( "/Models/elastic/loading/{}/parameters/expr", key );
                auto e = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
                auto loadpos = fmt::format( "/Models/elastic/loading/{}/parameters/location", key );
                std::vector<double> p = specs_[nl::json::json_pointer( loadpos )].get<std::vector<double>>();
                node_type n(p.size());
                for (int i = 0; i < p.size(); i++)
                    n(i) = p[i];
                auto dirac = std::make_shared<SensorPointwise<space_t>>(Xh_, n, key, e);
                auto f = form1( _test = Xh_, _vector = dirac->containerPtr() );
                l += f;
            }
            else if (specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Wavelet")
            {
                LOG( INFO ) << fmt::format( "Dirichlet conditions found" );
                std::string loadexpr = fmt::format( "/Models/elastic/loading/{}/parameters/expr", key );
                // auto e = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
                double force = wavelet(t);
                std::string e = "{";
                for (int i = 0; i < FEELPP_DIM-1;i++)
                {
                    e.append("0,");
                }
                e.append(std::to_string(wavelet(t)));
                e.append("}");
                // std::cout << e << std::endl;
                // TODO: invert e and p, e is the expression of the force and p is the position of the force
                auto loadpos = fmt::format( "/Models/elastic/loading/{}/parameters/location", key );
                std::vector<double> p = specs_[nl::json::json_pointer( loadpos )].get<std::vector<double>>();
                node_type n(p.size());
                for (int i = 0; i < p.size(); i++)
                    n(i) = p[i];
                auto dirac = std::make_shared<SensorPointwise<space_t>>(Xh_, n, key, e);
                auto f = form1( _test = Xh_, _vector = dirac->containerPtr() );
                l += f;
            }
            else if (specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Sensor")
            {
                LOG( INFO ) << fmt::format( "Dirichlet conditions found" );
                std::string loadexpr = fmt::format( "/Models/elastic/loading/{}/parameters/expr", key );
                // auto e = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
                std::string loadcsv = fmt::format( "/Models/elastic/loading/{}/parameters/csv", key );
                // std::cout << "loadcsv : " << specs_[nl::json::json_pointer( loadcsv )].get<std::string>() << std::endl;
                // Read the it line of the csv file
                std::ifstream file(specs_[nl::json::json_pointer( loadcsv )].get<std::string>());
                std::string line;
                // Check if the file is open
                if (!file.is_open())
                {
                    std::cerr << "Unable to open file : " << specs_[nl::json::json_pointer( loadcsv )].get<std::string>() << std::endl;
                    return;
                }

                std::vector<std::string> values;
                int i = 0;
                std::string e = "{";
                while (std::getline(file, line))
                {
                    if (i == it+1)
                    {
                        std::stringstream ss(line);
                        std::string value;
                        int j = 0;
                        while (std::getline(ss, value, ','))
                        {
                            if ((j>0) && (j<Dim))
                            {
                                e.append(value);
                                e.append(",");
                            }
                            else if (j==Dim)
                            {
                                e.append(value);
                                e.append("}");
                            }
                            values.push_back(value);
                            j++;
                        }
                        break;
                    }
                    i++;
                }
                file.close();
                // std::cout << "e : " << e << std::endl;
                // e.append(std::to_string(force));
                // e.append(",");
                // for (int i = 1; i < FEELPP_DIM-1;i++)
                // {
                //     e.append("0,");
                // }
                // e.append("0}");
                auto loadpos = fmt::format( "/Models/elastic/loading/{}/parameters/location", key );
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


    //// Boundary Condition Neumann
    //auto bc_neu = specs_["/BoundaryConditions/elastic/Neumann/Gamma/pure_traction/g/expr"_json_pointer].get<std::string>();
    //std::cout << "Boundary Condition Neumann : " <<  bc_neu << std::endl;
    //l_ += integrate( _range = markedfaces( mesh_ , "Gamma"), _expr = inner( expr<FEELPP_DIM,1>( bc_neu ), id( v_ )) );
}

// Run method (main method to run Laplacian process)
template <int Dim, int Order>
void Elastic<Dim, Order>::run()
{
    std::cout << "\n***** Initialize *****" << std::endl;
    initialize();
    //std::cout << "\n***** Process materials *****" << std::endl;
    //processMaterials();
    std::cout << "\n***** Process boundary conditions *****" << std::endl;
    //processBoundaryConditions();
    std::cout << "\n***** Time loop *****" << std::endl;
    timeLoop();
    std::cout << "\n***** Export results *****" << std::endl;
    writeResultsToFile("results.json");
    //exportResults();
}

// Time loop
template <int Dim, int Order>
void Elastic<Dim, Order>::timeLoop()
{
    int it = 0; // initialization of the counter for time backtracing when evaluating the disturbance on the point of interest
    processLoading(lt_);
    processMaterials(a_);

    for ( ts_->start(); ts_->isFinished()==false; ts_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << "time " << ts_->time() << std::endl;

        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////
        at_ = a_;

        for ( auto [key, material] : specs_["/Models/elastic/Materials"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Material {} found", material );

            std::string matRho = fmt::format( "/Materials/{}/parameters/rho/value", material.get<std::string>() );
            auto rho = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
            lt_ +=  integrate( _range=markedelements( mesh_, material.get<std::string>() ), _expr= rho*inner( idv(ts_->polyDeriv()),id( v_ ) ) );
        }

        processBoundaryConditions(lt_, at_,ts_->time(),it);

        at_.solve( _rhs = lt_, _solution = u_ );

        // Uncomment and export ts_
        // this->exportResults();
        ts_->updateFromDisp(u_);
        auto u_proj = project(_space=Xh_,_range=elements(mesh_),_expr=idv(u_));
        e_->step(ts_->time())->add( "displacement", u_proj);
        e_->step(ts_->time())->add( "velocity", ts_->currentVelocity() );
        e_->step(ts_->time())->add( "acceleration", ts_->currentAcceleration() );
        e_->save();

        ////////////////////////////////////////////////////
        //         Newmark gamma-model for dtun           //
        ////////////////////////////////////////////////////
        at_.zero();
        lt_.zero();

        it += 1;

        for (auto [key, sensor] : specs_["/Sensors"_json_pointer].items())
        {
            auto loadpos = fmt::format( "/Sensors/{}/location", key );
            std::vector<double> p = specs_[nl::json::json_pointer( loadpos )].get<std::vector<double>>();
            // std::cout << "p(0):" << p[0] << " p(1):" << p[1] << std::endl;
            Eigen::Matrix<double,FEELPP_DIM,1> __x;
            for (int i = 0; i < FEELPP_DIM; i++)
                __x(i) = p[i];
            // std::cout << "u(" << p[0] << "," << p[1] <<"): " << u_(__x) << std::endl;
            std::ofstream file;
            file.open(fmt::format( "{}.csv", key ), std::ios::app);
            // Check if the file is open
            if (!file.is_open())
            {
                std::cerr << "Unable to open file : " << fmt::format( "{}.csv", key ) << std::endl;
                return;
            }
            file << ts_->time() << ",";
            for (int i = 0; i < FEELPP_DIM-1; i++)
            {
                file << u_(__x)(i,0,0) << ",";
            }
            file << u_(__x)(FEELPP_DIM-1,0,0);
            file << "\n";
            file.close();
        }
    }
}


// Export results
template <int Dim, int Order>
void Elastic<Dim, Order>::exportResults()
{
    // e_->step(ts_->time())->addRegions();
    // e_->step(ts_->time())->add("u", u_);
    // e_->save();


    // auto totalQuantity = integrate(_range=elements(mesh_), _expr=idv(u_)).evaluate()(0,0);
    // auto totalFlux = integrate(_range=boundaryfaces(mesh_), _expr=gradv(u_)*N()).evaluate()(0,0);
    // double meas=measure(_range=elements(mesh_), _expr=cst(1.0));
    // meas_["time"].push_back(ts_->time());
    // meas_["totalQuantity"].push_back(totalQuantity);
    // meas_["totalFlux"].push_back(totalFlux);
    // meas_["mean"].push_back(totalQuantity/meas);
    // meas_["min"].push_back(u_.min());
    // meas_["max"].push_back(u_.max());
    for( auto [key,values] : mesh_->markerNames())
    {
        // std::cout << "key : " << key << std::endl;
        if (key=="S1")
        {
            double meas = measure(_range=markedelements(mesh_,key), _expr=cst(1.0));
            std::cout << "meas : " << meas << std::endl;
            meas_[fmt::format("quantity_{}",key)].push_back(meas);
        }
        // if ( values[1] == Dim )
        // {
        //     double meas=measure(_range=markedelements(mesh_,key), _expr=cst(1.0));
        //     auto quantity = integrate(_range=markedelements(mesh_,key), _expr=idv(u_)).evaluate()(0,0);
        //     meas_[fmt::format("quantity_{}",key)].push_back(quantity);
        //     meas_[fmt::format("mean_{}",key)].push_back(quantity/meas);
        // }
        // else if ( values[1] == Dim-1 )
        // {
        //     double meas=measure(_range=markedfaces(mesh_,key), _expr=cst(1.0));
        //     auto quantity = integrate(_range=markedfaces(mesh_,key), _expr=idv(u_)).evaluate()(0,0);
        //     meas_[fmt::format("quantity_{}",key)].push_back(quantity);
        //     meas_[fmt::format("mean_{}",key)].push_back(quantity/meas);
        //     auto flux = integrate(_range=markedfaces(mesh_,key), _expr=gradv(u_)*N()).evaluate()(0,0);
        //     meas_[fmt::format("flux_{}",key)].push_back(flux);
        // }
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
