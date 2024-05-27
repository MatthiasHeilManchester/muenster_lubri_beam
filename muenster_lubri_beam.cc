//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
//Driver function for a simple beam proble

//OOMPH-LIB includes
#include "generic.h"
#include "beam.h"
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;

using namespace oomph;


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


//=====================================================================
/// Upgraded Hermite Beam Element to incorporate surface tension driven
/// thin film flow
//=====================================================================
class HermiteLubriBeamElement : public virtual HermiteBeamElement
{
 
public:
 
 /// Constructor: 
 HermiteLubriBeamElement() 
  {
  }


 ///  Required  # of `values' (pinned or dofs)
 /// at node n
 inline unsigned required_nvalue(const unsigned& n) const
  {
   return 2; // Initial_Nvalue;
  }

 
 /// Film thickness dof (local node j, type k)
 double nodal_h_lubri(const unsigned& j, const unsigned& k) const
  {
   return this->node_pt(j)->value(k);
  }
 
 
 /// Film thickness at s
 double interpolated_h_lubri(const Vector<double>& s)
  {
   // Initialise storage h_lubri
   double h_lubri = 0.0;
   
   // Number of nodes
   unsigned n_node = this->nnode();
   
   // Number of degrees of freedom per node
   unsigned n_value = this->node_pt(0)->nvalue();
   
   // Get shape function
   Shape psi(n_node, n_value);
   this->shape(s, psi);
   
   // Sum up
   for (unsigned n = 0; n < n_node; n++)
    {
     for (unsigned k = 0; k < n_value; k++)
      {
       h_lubri += nodal_h_lubri(n, k) * psi(n, k);
      }
    }
   
   return h_lubri;
  }
 
 /// Overwrite otput function
 void output(std::ostream& outfile, const unsigned& n_plot)
  {
   // Local variables
   Vector<double> s(1);
   
   // Tecplot header info
   outfile << "ZONE I=" << n_plot << std::endl;
   
   // Set the number of lagrangian coordinates
   unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();
   
   // Set the dimension of the global coordinates
   unsigned n_dim = Undeformed_beam_pt->ndim();
   
   // Find out how many nodes there are
   unsigned n_node = nnode();
   
   // Find out how many positional dofs there are
   unsigned n_position_dofs = nnodal_position_type();
   Vector<double> posn(n_dim);
   
   // # of nodes, # of positional dofs
   Shape psi(n_node, n_position_dofs);


   // hierher

   // {
   //  DShape& dpsidx,
   //   DShape& d2psidx
   //  d2shape_eulerian(const Vector<double>& s,
   //                   Shape& psi,
   //                   DShape& dpsidx,
   //                   DShape& d2psidx) const;
    
    
   // }

   
   // Loop over element plot points
   for (unsigned l1 = 0; l1 < n_plot; l1++)
    {
     s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);
     
     // Get shape functions
     shape(s, psi);
     
     Vector<double> interpolated_xi(n_lagrangian);
     interpolated_xi[0] = 0.0;
     
     // Initialise
     double h_lubri=0.0;
     for (unsigned i = 0; i < n_dim; i++)
      {
       posn[i] = 0.0;
      }
     
     
     // Calculate values
     for (unsigned l = 0; l < n_node; l++)
      {
       // Loop over positional dofs
       for (unsigned k = 0; k < n_position_dofs; k++)
        {
         // Loop over Lagrangian coordinate directions [xi_gen[] are the
         // the *gen*eralised Lagrangian coordinates: node, type, direction]
         for (unsigned i = 0; i < n_lagrangian; i++)
          {
           interpolated_xi[i] +=
            raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }
         
         // Loop over components of the deformed position Vector
         for (unsigned i = 0; i < n_dim; i++)
          {
           posn[i] += raw_dnodal_position_gen_dt(0, l, k, i) * psi(l, k);
          }
         h_lubri+=nodal_h_lubri(l, k) * psi(l, k);
        }
      }
     

     // Beam
     for (unsigned i = 0; i < n_dim; i++)
      {
       outfile << posn[i] << " ";
      }
     // Surface of fluid and film thickness (measured
     // vertically upwards
     outfile << posn[1]+h_lubri << " "
             << h_lubri << " "
             << std::endl;
    }
  }
 
};



/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////



 
//========start_of_namespace========================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
 /// Non-dimensional thickness
 double H=0.05;

 /// 2nd Piola Kirchhoff pre-stress
 double Sigma0=0.0;

 /// Pressure load
 double P_ext=0.0;

 /// Square of timescale ratio (i.e. non-dimensional density)  
 /// -- 1.0 for default value of scaling factor
 double Lambda_sq=1.0;

 /// Load function: Apply a constant external pressure to the beam
 void load(const Vector<double>& xi, const Vector<double> &x,
           const Vector<double>& N, Vector<double>& load)
 {
  for(unsigned i=0;i<2;i++) {load[i] = -P_ext*N[i];}
 }

} // end of namespace

//======start_of_problem_class==========================================
/// Beam problem object
//======================================================================
class ElasticBeamProblem : public Problem
{
public:
 
 /// Constructor: The arguments are the number of elements, 
 /// the length of domain
 ElasticBeamProblem(const unsigned &n_elem, const double &length);
 
 /// Conduct a parameter study
 void parameter_study();
 
 /// Return pointer to the mesh
 OneDLagrangianMesh<HermiteLubriBeamElement>* mesh_pt() 
  {return dynamic_cast<OneDLagrangianMesh<HermiteLubriBeamElement>*>
    (Problem::mesh_pt());}

 /// No actions need to be performed after a solve
 void actions_after_newton_solve() {}

 /// No actions need to be performed before a solve
 void actions_before_newton_solve() {}

private:

 /// Pointer to the node whose displacement is documented
 Node* Doc_node_pt;

 /// Length of domain (in terms of the Lagrangian coordinates)
 double Length;

 /// Pointer to geometric object that represents the beam's undeformed shape
 GeomObject* Undef_beam_pt;

}; // end of problem class


//=============start_of_constructor=====================================
/// Constructor for elastic beam problem
//======================================================================
ElasticBeamProblem::ElasticBeamProblem(const unsigned &n_elem,
                                       const double &length) : Length(length)
{
 // Set the undeformed beam to be a straight line at y=0
 Undef_beam_pt=new StraightLine(0.0); 

 // Create the timestepper and add it to the Problem's collection of
 // timesteppers -- this creates the Problem's Time object.
 add_time_stepper_pt(new Newmark<3>()); // hierher why 3?

 
 // Create the (Lagrangian!) mesh, using the geometric object
 // Undef_beam_pt to specify the initial (Eulerian) position of the
 // nodes.
 Problem::mesh_pt() = 
  new OneDLagrangianMesh<HermiteLubriBeamElement>(n_elem,length,Undef_beam_pt,
                                                  Problem::time_stepper_pt());

 // Set the boundary conditions: Each end of the beam is fixed in space
 // Loop over the boundaries (ends of the beam)
 for(unsigned b=0;b<2;b++)
  {
   // Pin displacements in both x and y directions
   // [Note: The mesh_pt() function has been overloaded
   //  to return a pointer to the actual mesh, rather than
   //  a pointer to the Mesh base class. The current mesh is derived
   //  from the SolidMesh class. In such meshes, all access functions
   //  to the nodes, such as boundary_node_pt(...), are overloaded
   //  to return pointers to SolidNodes (whose position can be
   //  pinned) rather than "normal" Nodes.]
   mesh_pt()->boundary_node_pt(b,0)->pin_position(0); 
   mesh_pt()->boundary_node_pt(b,0)->pin_position(1);

   // Clamp
   mesh_pt()->boundary_node_pt(b,0)->pin_position(1,1);

  }


 // For now: Pin all lubri dofs and assign values
 // consistnt with constant film thickness
 double h_initial=0.1; // hierher global namespace
 unsigned nnod=Problem::mesh_pt()->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   // Value
   Problem::mesh_pt()->node_pt(j)->pin(0);
   Problem::mesh_pt()->node_pt(j)->set_value(0,h_initial);

   // Slope
   Problem::mesh_pt()->node_pt(j)->pin(1);
   Problem::mesh_pt()->node_pt(j)->set_value(1,0.0);
  }
 
 //Find number of elements in the mesh
 unsigned n_element = mesh_pt()->nelement();
 
 //Loop over the elements to set physical parameters etc.
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   HermiteLubriBeamElement *elem_pt = 
    dynamic_cast<HermiteLubriBeamElement*>(mesh_pt()->element_pt(e));
   
   // Set physical parameters for each element:
   elem_pt->sigma0_pt() = &Global_Physical_Variables::Sigma0;
   elem_pt->h_pt() = &Global_Physical_Variables::H;

   // Set the load Vector for each element
   elem_pt->load_vector_fct_pt() = &Global_Physical_Variables::load;

   // Pass pointer to square of timescale ratio (non-dimensional density)
   elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = Undef_beam_pt;
   
  } // end of loop over elements

 // Choose node at which displacement is documented (halfway along -- provided
 // we have an odd number of nodes; complain if this is not the
 // case because the comparison with the exact solution will be wrong 
 // otherwise!)
 unsigned n_nod=mesh_pt()->nnode();
 if (n_nod%2!=1)
  {
   cout << "Warning: Even number of nodes " << n_nod << std::endl;
   cout << "Comparison with exact solution will be misleading..." << std::endl;
  }
 Doc_node_pt=mesh_pt()->node_pt((n_nod+1)/2-1);
 
 // Assign the global and local equation numbers
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;

} // end of constructor


//=======start_of_parameter_study==========================================
/// Solver loop to perform parameter study
//=========================================================================
void ElasticBeamProblem::parameter_study()
{
 // Over-ride the default maximum value for the residuals
 Problem::Max_residuals = 1.0e10;
 
 // Set the increments in control parameters
 double pext_increment = -0.00001;
 
 // Set initial values for control parameters 
 Global_Physical_Variables::P_ext = 0.0 - pext_increment;
 
 // Set the 2nd Piola Kirchhoff prestress
 Global_Physical_Variables::Sigma0=0.01;
  
 
 // Create label for output
 DocInfo doc_info;
 
 // Set output directory -- this function checks if the output
 // directory exists and issues a warning if it doesn't.
 doc_info.set_directory("RESLT");
 
 // Open a trace file
 ofstream trace("RESLT/trace_beam.dat");
 
 // Write a header for the trace file
 trace <<  "VARIABLES=\"p_e_x_t\",\"sigma_0\""
       <<  ", \"d\"" << std::endl;
 
 // Output file stream used for writing results
 ofstream file;
 
 // String used for the filename
 char filename[100]; 

 // Counter for output
 unsigned counter=0;
 
 // STAGE 1: INFLATE, WITH GIVEN POSITIVE PRESTRESS
 //------------------------------------------------
 {
  
  // Set initial values for control parameters 
  Global_Physical_Variables::P_ext = 0.0 - pext_increment;
  
  // Loop over parameter increments
  unsigned nstep=10;
  for(unsigned i=1;i<=nstep;i++)
   {
    // Increment pressure
    Global_Physical_Variables::P_ext += pext_increment;
    
    oomph_info << "Solving for p_ext sigma_0 = "
               << Global_Physical_Variables::P_ext << " " 
               <<  Global_Physical_Variables::Sigma0
               << std::endl;

   
    // Solve the system
    steady_newton_solve();
    
    // Document the solution
    sprintf(filename,"RESLT/beam%i.dat",counter);
    counter++;
    file.open(filename);
    mesh_pt()->output(file,5);
    file.close();
    
    // Write trace file
    trace << Global_Physical_Variables::P_ext  << " "   
             <<  Global_Physical_Variables::Sigma0 << " "
            << abs(Doc_node_pt->x(1))
            << std::endl;
   }
 }
 
 
 // STAGE 2: CHANGE SIGN OF PRESTRESS
 //----------------------------------
 {
  unsigned nstep=10;
  double d_sigma=
   2.0*Global_Physical_Variables::Sigma0/double(nstep-1);
  
  // Loop over parameter increments
  for(unsigned i=1;i<=nstep;i++)
   {
    // Decrement pre-stress
    Global_Physical_Variables::Sigma0-=d_sigma;

    oomph_info << "Solving for p_ext sigma_0 = "
               << Global_Physical_Variables::P_ext << " " 
               <<  Global_Physical_Variables::Sigma0
               << std::endl;
    
    // Solve the system
    steady_newton_solve();
    
    // Document the solution
    sprintf(filename,"RESLT/beam%i.dat",counter);
    counter++;
    file.open(filename);
    mesh_pt()->output(file,5);
    file.close();
    
    // Write trace file
    trace << Global_Physical_Variables::P_ext  << " "   
          <<  Global_Physical_Variables::Sigma0 << " "
          << abs(Doc_node_pt->x(1))
          << std::endl;
    
   }
  
 }

 // STAGE 3: REDUCE PRESSURE TO ZERO
 //---------------------------------
 {

  
  // Loop over parameter increments
  unsigned nstep=10;
  double pext_increment=Global_Physical_Variables::P_ext/double(nstep);
  for(unsigned i=1;i<=nstep;i++)
   {
    // Increment pressure
    Global_Physical_Variables::P_ext -= pext_increment;
    
    oomph_info << "Solving for p_ext sigma_0 = "
               << Global_Physical_Variables::P_ext << " " 
               <<  Global_Physical_Variables::Sigma0
               << std::endl;
    
    // Solve the system
    steady_newton_solve();
    
    // Document the solution
    sprintf(filename,"RESLT/beam%i.dat",counter);
    counter++;
    file.open(filename);
    mesh_pt()->output(file,5);
    file.close();
    
    // Write trace file
    trace << Global_Physical_Variables::P_ext  << " "   
          <<  Global_Physical_Variables::Sigma0 << " "
          << abs(Doc_node_pt->x(1))
          << std::endl;
   }
 }


 // STAGE 4: TIMESTEP THE THING
 //----------------------------
 {
  // Set timestep
  double dt=1.0; 
  
  // Assign impulsive start
  assign_initial_values_impulsive(dt); // hierher try bypassing

  Global_Physical_Variables::P_ext -= pext_increment;


  unsigned nstep=1000;
  for (unsigned i=0;i<nstep;i++)
   {
    // Solve
    oomph_info << "Doing unsteady solve for t = "
               << time_stepper_pt()->time_pt()->time()
               << std::endl;
    unsteady_newton_solve(dt);
    
    // Document the solution
    sprintf(filename,"RESLT/beam%i.dat",counter);
    counter++;
    file.open(filename);
    mesh_pt()->output(file,5);
    file.close();
    
    // Write trace file
    trace << Global_Physical_Variables::P_ext  << " "   
          <<  Global_Physical_Variables::Sigma0 << " "
          << abs(Doc_node_pt->x(1))
          << std::endl;
    
   }
 }

 
} // end of parameter study

//========start_of_main================================================
/// Driver for beam (string under tension) test problem 
//=====================================================================
int main()
{

 // Set the non-dimensional thickness 
 Global_Physical_Variables::H=0.05; 
 
 // Set the length of domain
 double L = 10.0;

 // Number of elements (choose an even number if you want the control point 
 // to be located at the centre of the beam)
 unsigned n_element = 10;

 // Construst the problem
 ElasticBeamProblem problem(n_element,L);

 // Check that we're ready to go:
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Conduct parameter study
 problem.parameter_study();

} // end of main

