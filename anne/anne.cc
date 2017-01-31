//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Driver for Anne's MSc project

// The oomphlib headers
#include "generic.h"
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;
using namespace oomph;


//===start_of_namespace=================================================
/// Namespace for global parameters
//======================================================================
namespace Global_Parameters
{
 /// Reynolds number
 double Re=100.0;

 /// x position of vortex
 double X_vortex=0.0;

 /// x position of vortex
 double Y_vortex=1.0;

 /// Initial condition for velocity
 void initial_condition(const Vector<double>& x, Vector<double>& u)
 {

  // hierher Anne: Update this!
  //               (1) All these numbers must have some relation to
  //                   the Reynolds number (I just set them to 1)
  //               (2) Add uniform background flow
  double a=1.0;
  double omega_0=1.0;
  double h=1.0;

  // Top vortex
  double r1=sqrt(pow(x[0]-X_vortex,2)+
                pow(x[1]-Y_vortex,2));
  double theta1=atan2(x[1]-Y_vortex,x[0]-X_vortex);
  double u_theta1=omega_0*h*h/(2.0*r1)*(1.0-exp(-r1*r1/(a*a)));

  // Bottom vortex
  double r2=sqrt(pow(x[0]-X_vortex,2)+
                 pow(x[1]+Y_vortex,2));
  double theta2=atan2(x[1]+Y_vortex,x[0]-X_vortex);
  double u_theta2=omega_0*h*h/(2.0*r2)*(1.0-exp(-r2*r2/(a*a)));

  u[0]=-u_theta1*sin(theta1)+u_theta2*sin(theta2);
  u[1]= u_theta1*cos(theta1)-u_theta2*cos(theta2);
 }

} // end of namespace


//===start_of_problem_class=============================================
/// Problem class for Anne's MSc problem
//======================================================================
template<class ELEMENT>
class AnneProblem : public Problem
{
public:

 /// Constructor:
 AnneProblem();

 //Update before solve is empty
 void actions_before_newton_solve() {}

 /// \short Update after solve is empty
 void actions_after_newton_solve() {}
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Impose no slip and re-assign eqn numbers
 void impose_no_slip_on_bottom_boundary();

private:

 /// oomph-lib iterative linear solver
 IterativeLinearSolver* Solver_pt;
 
 /// Preconditioner
 NavierStokesSchurComplementPreconditioner* Prec_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;
}; // end of problem_class


//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT>
AnneProblem<ELEMENT>::AnneProblem()
{

 //Allocate the timestepper
 add_time_stepper_pt(new BDF<2>); 

 // Number of elements in x direction
 unsigned nx=10;

 // Number of elements in y direction
 unsigned ny=10;

 // Left end of computational domain
 double x_min=-3.36;
 
 // Right end of computational domain
 double x_max=3.36;
 
 // Botton of computational domain
 double y_min=0.0;

 // Height of computational domain
 double y_max=5.13;
 
 //Now create the mesh 
 mesh_pt() = 
  new RectangularQuadMesh<ELEMENT>(nx,ny,x_min,x_max,y_min,y_max,
                                   time_stepper_pt());

 // Check enumeration of boundaries
 mesh_pt()->output_boundaries("boundaries.dat");

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here
 unsigned num_bound=mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Imposed velocity on top (2) and left (3) 
     if ((ibound==2)||(ibound==3))
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
     // Horizontal outflow on the left (1) and no penetration
     // at bottom (0)
     else if ((ibound==0)||(ibound==1))
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
    }
  } // end loop over boundaries

 //Complete the problem setup to make the elements fully functional

 //Loop over the elements
 unsigned n_el = mesh_pt()->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Parameters::Re;   
   el_pt->re_st_pt() = &Global_Parameters::Re;

  }

 //Assign equation numbers
 assign_eqn_numbers();

 // Setup initial condition
 //------------------------
 
 // Loop over nodes
 Vector<double> x(2);
 Vector<double> u(2);
 unsigned num_nod = mesh_pt()->nnode();
 for (unsigned n=0;n<num_nod;n++)
  {
   // Get nodal coordinates
   x[0]=mesh_pt()->node_pt(n)->x(0);
   x[1]=mesh_pt()->node_pt(n)->x(1);
   
   // Get initial velocity field
   Global_Parameters::initial_condition(x,u);
   
   // Assign solution
   mesh_pt()->node_pt(n)->set_value(0,u[0]);
   mesh_pt()->node_pt(n)->set_value(1,u[1]);
  }


 // Linear solver
 //--------------
 if (CommandLineArgs::command_line_flag_has_been_set("--use_oomph_gmres"))
  {
   // Use GMRES
   Solver_pt=new GMRES<CRDoubleMatrix>;   
   linear_solver_pt()=Solver_pt;
   
   // Set preconditioner
   Prec_pt=new NavierStokesSchurComplementPreconditioner(this);
   Prec_pt->set_navier_stokes_mesh(this->mesh_pt());  
   Solver_pt->preconditioner_pt()=Prec_pt;
  }
 
 
 
} // end of constructor




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void AnneProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);

} // end_of_doc_solution   



//========================================================================
/// Impose no slip and re-assign eqn numbers
//========================================================================
template<class ELEMENT>
void AnneProblem<ELEMENT>::impose_no_slip_on_bottom_boundary()
{


 // Pin horizontal velocity at bottom boundary
 unsigned ibound=0;
 unsigned num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
   mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
  }

 // Re-assign equation numbers
 oomph_info << std::endl 
            << "ndofs before applying no slip at bottom boundary: "
            << ndof() << std::endl;
 assign_eqn_numbers();
 oomph_info << "ndofs after applying no slip at bottom boundary : "
            << ndof() << std::endl << std::endl;


}



//===start_of_main======================================================
/// Driver code for Anne channel problem
//======================================================================
int main(int argc, char* argv[]) 
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Reynolds number
 CommandLineArgs::specify_command_line_flag("--re",
                                            &Global_Parameters::Re);

 // Use gmres?
 CommandLineArgs::specify_command_line_flag("--use_oomph_gmres");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Set up doc info
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 //Set up problem
 AnneProblem<QTaylorHoodElement<2> > problem;
  
// Initialise all history values for an impulsive start
 double dt=0.005; 
 problem.initialise_dt(dt);
 problem.assign_initial_values_impulsive();

 // Number of timesteps until switch-over to no slip
 unsigned ntsteps=10;

 // Doc initial condition
 problem.doc_solution(doc_info);
 
 // increment counter
 doc_info.number()++;

 //Loop over the timesteps
 for(unsigned t=1;t<=ntsteps;t++)
  {
   oomph_info << "TIMESTEP " << t << std::endl;
   
   //Take one fixed timestep
   problem.unsteady_newton_solve(dt);

   //Output the time
   oomph_info << "Time is now " << problem.time_pt()->time() << std::endl;

   // Doc solution
   problem.doc_solution(doc_info);

   // increment counter
   doc_info.number()++;
  }

 // Now do no slip
 problem.impose_no_slip_on_bottom_boundary();

 //Loop over the remaining timesteps
 ntsteps=1000;
 for(unsigned t=1;t<=ntsteps;t++)
  {
   oomph_info << "TIMESTEP " << t << std::endl;
   
   //Take one fixed timestep
   problem.unsteady_newton_solve(dt);

   //Output the time
   oomph_info << "Time is now " << problem.time_pt()->time() << std::endl;

   // Doc solution
   problem.doc_solution(doc_info);

   // increment counter
   doc_info.number()++;
  }



} // end of main
