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
 HermiteLubriBeamElement() : Q_pt(0)
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
 
 /// Film thickness dof (local node j, type k) at previous time level t
 /// (t==: present)
 double nodal_h_lubri(const unsigned& t, const unsigned& j, 
                      const unsigned& k) const
  {
   return this->node_pt(j)->value(t,k);
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

   // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
   DShape dpsidxi(n_node, n_position_dofs, n_lagrangian);
   
   // # of nodes, # of positional dofs, # of derivs)
   DShape d2psidxi(n_node, n_position_dofs, n_lagrangian);
      
   // Loop over element plot points
   for (unsigned l1 = 0; l1 < n_plot; l1++)
    {
     s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);
     
     // Call the derivatives of the shape functions w.r.t. Lagrangian coords
     double J = d2shape_lagrangian(s, psi, dpsidxi, d2psidxi);
     
     // Initialise
     double h_lubri=0.0;
     double dh_dt=0.0;
     double dh_lubri_dxi=0.0;
     double curv=0.0;
     for (unsigned i = 0; i < n_dim; i++)
      {
       posn[i] = 0.0;
      }
     
     
     // Calculate values
     for (unsigned l = 0; l < n_node; l++)
      {
       // Loop over dof types
       for (unsigned k = 0; k < n_position_dofs; k++)
        {
         // Loop over components of the deformed position Vector
         for (unsigned i = 0; i < n_dim; i++)
          {
           posn[i] += raw_dnodal_position_gen_dt(0, l, k, i) * psi(l, k);
          }
         h_lubri+=nodal_h_lubri(l, k) * psi(l, k);
         dh_lubri_dxi+=nodal_h_lubri(l,k)*dpsidxi(l,k,0);

         // curvature of free surface, taking substrate (beam) curvature into
         // account hierher computed using (i) linearity (ii) ignoring
         // horizontal beam displacemen!
         curv+=(nodal_h_lubri(l, k)+raw_nodal_position_gen(l,k,1))*
          d2psidxi(l,k,0);

         // Number of timsteps (past & present)
         const unsigned n_time =node_pt(l)->time_stepper_pt()->ntstorage();
         
         // Add the contributions to the time derivative
         for (unsigned t = 0; t < n_time; t++)
          {
           dh_dt+=node_pt(l)->time_stepper_pt()->weight(1, t)*
            nodal_h_lubri(t,l,k);
          }
        }
      }
     

     // Beam
     for (unsigned i = 0; i < n_dim; i++)
      {
       outfile << posn[i] << " "; // 1,2
      }
     // Surface of fluid and film thickness (measured
     // vertically upwards
     outfile << posn[1]+h_lubri << " " // 3
             << h_lubri << " "         // 4
             << dh_lubri_dxi << " "    // 5
             << curv << " "            // 6
             << dh_dt << " "           // 7
             << J << " "               // 8
             << std::endl;
    }
  }


 /// Pointer to FSI parameter (const version)
 double* q_pt() const
  {
   return Q_pt;
  }
 
 /// Pointer to FSI parameter (read/write version)
 double*& q_pt()
  {
   return Q_pt;
  }
 
 /// FSI parameter 
 double q() const
  {
   if (Q_pt==0)
    {
     return 0.0;
    }
   
   return *Q_pt;
  }


 
 
 /// hierher
 void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
   fill_in_contribution_to_residuals_beam(residuals);
   fill_in_contribution_to_residuals_lubri(residuals);
  }

 /// Get the Jacobian and residuals. hierher overloaded yet againi
 virtual void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian)
  {
   // Call FD version in base class
   SolidFiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
  }

 
 /// Overloaded version of the load vector function to incorporate lubrication
 /// pressure. : Pass number of integration point,
 /// Lagr. and Eulerian coordinate and normal vector and return the load
 /// vector.
 void load_vector(const unsigned& ipt,
                  const Vector<double>& xi,
                  const Vector<double>& x,
                  const Vector<double>& N,
                  Vector<double>& load)
  {
   // Get external pressure
   Load_vector_fct_pt(xi, x, N, load);
   
   // Add lubri traction
   //-------------------
   
   // Set the number of lagrangian coordinates
   const unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();
   
   // Find out how many nodes there are
   const unsigned n_node = nnode();
   
   // Find out how many positional dofs there are
   const unsigned n_position_type = nnodal_position_type();
   
   // # of nodes, # of positional dofs
   Shape psi(n_node, n_position_type);
   
   // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
   DShape dpsidxi(n_node, n_position_type, n_lagrangian);
   
   // # of nodes, # of positional dofs, # of derivs)
   DShape d2psidxi(n_node, n_position_type, n_lagrangian);
   
   // Call the derivatives of the shape functions w.r.t. Lagrangian coords
   d2shape_lagrangian_at_knot(ipt, psi, dpsidxi, d2psidxi);
   
   // Get curvature
   double curv=0.0;
   for (unsigned l = 0; l < n_node; l++)
    {
     // Loop over dof types
     for (unsigned k = 0; k < n_position_type; k++)
      {
       // curvature of free surface, taking substrate (beam) curvature into
       // account hierher (i) linear (ii) ignore horizontal beam displacement!
       curv+=(nodal_h_lubri(l, k)+raw_nodal_position_gen(l,k,1))*d2psidxi(l,k,0);
      }
    }
   
   double q_fsi=q(); // 1.0e-4; // 1.0e-6; //1.0e-8;
   load[0]+=q_fsi*curv*N[0];
   load[1]+=q_fsi*curv*N[1];
   
  }

 /// Volume (per unit depgth) of fluid in the film
 void get_fluid_volume(double& v_fluid)
  {
   // Initialise
   v_fluid = 0.0;
      
   // Set the number of lagrangian coordinates
   unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();
   
   // Find out how many nodes there are
   unsigned n_node = nnode();
   
   // Find out how many positional dofs there are
   unsigned n_position_dofs = nnodal_position_type();
   
   // Set up memory for the shape functions:
   // # of nodes, # of positional dofs
   Shape psi(n_node, n_position_dofs);
   
   // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
   DShape dpsidxi(n_node, n_position_dofs, n_lagrangian);
   
   // # of nodes, # of positional dofs, # of derivs)
   DShape d2psidxi(n_node, n_position_dofs, n_lagrangian);
   
   // Set # of integration points
   unsigned n_intpt = integral_pt()->nweight();
   
   // Loop over the integration points
   for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
     // Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     // Call the derivatives of the shape functions w.r.t. Lagrangian coords
     double J = d2shape_lagrangian_at_knot(ipt, psi, dpsidxi, d2psidxi);
     
     // Premultiply the weights and the Jacobian
     double W = w * J;
     
     
     // Sum 'em
     double h_lubri=0.0;
     for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_dofs; k++)
        {
         h_lubri+=nodal_h_lubri(l, k) * psi(l, k);
        }
      }

      // Add contribution
     v_fluid+=h_lubri*W;
     
    } // End of loop over the integration points
  }
 

protected:

 /// hierher 
 void fill_in_contribution_to_residuals_lubri(Vector<double>& residuals)
  {

   // hierher
   // // Set up the initial conditions, if an IC pointer has been set
   // if (Solid_ic_pt != 0)
   // {
   //   fill_in_residuals_for_solid_ic(residuals);
   //   return;
   // }
   
   // Set the number of lagrangian coordinates
   const unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();
   
   // Find out how many nodes there are
   const unsigned n_node = nnode();
   
   // Find out how many positional dofs there are
   const unsigned n_position_type = nnodal_position_type();
   
   // Integer to store the local equation number
   int local_eqn = 0;
   
   
   // Set up memory for the shape functions:
   
   // # of nodes, # of positional dofs
   Shape psi(n_node, n_position_type);
   
   // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
   DShape dpsidxi(n_node, n_position_type, n_lagrangian);
   
   // # of nodes, # of positional dofs, # of derivs)
   DShape d2psidxi(n_node, n_position_type, n_lagrangian);
   
   // Set # of integration points
   const unsigned n_intpt = integral_pt()->nweight();
   
   // Get Physical Variables from Element

   // hierher
   double maybe_inverse_capillary=1.0; // 10.0;
   
   // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape functions w.r.t. Lagrangian coords
      double J = d2shape_lagrangian_at_knot(ipt, psi, dpsidxi, d2psidxi);

     // Initialise
     double h_lubri=0.0;
     double dh_lubri_dxi=0.0;
     double dh_dt=0.0;
     double curv=0.0;
     
     // Calculate values
     for (unsigned l = 0; l < n_node; l++)
      {
       // Loop over dof types
       for (unsigned k = 0; k < n_position_type; k++)
        {
         h_lubri+=nodal_h_lubri(l, k) * psi(l, k);

         
         // slope (take change into arclength into account
         dh_lubri_dxi+=nodal_h_lubri(l, k)*dpsidxi(l,k,0);

         // curvature of free surface, taking substrate (beam) curvature into
         // account hierher (i) linear (ii) ignore horizontal beam displacement!
         curv+=(nodal_h_lubri(l, k)+raw_nodal_position_gen(l,k,1))*
          d2psidxi(l,k,0);

         // Number of timsteps (past & present)
         const unsigned n_time =node_pt(l)->time_stepper_pt()->ntstorage();
         
         // Add the contributions to the time derivative
         for (unsigned t = 0; t < n_time; t++)
          {
           dh_dt+=node_pt(l)->time_stepper_pt()->weight(1,t)*
            nodal_h_lubri(t,l,k);
          }
        }
      }
     
     // Premultiply the weights and the Jacobian
     double W = w * J;
     
     
     // Loop over nodes
     for (unsigned j=0;j<n_node;j++)
      {
       // Loop over dof types
       for (unsigned k=0;k<n_position_type;k++)
        {
         local_eqn=nodal_local_eqn(j,k);
         // If it's not a boundary condition
         if (local_eqn >= 0)
          {
           residuals[local_eqn] +=
            (dh_dt*psi(j,k)+maybe_inverse_capillary*curv*
             (h_lubri*h_lubri*dh_lubri_dxi*dpsidxi(j,k,0)+
              h_lubri*h_lubri*h_lubri/3.0*d2psidxi(j,k,0)))*W;
           
            // (dh_dt*psi(j,k)+maybe_inverse_capillary*curv*d2psidxi(j,k,0))*W;
            
           // std::cout << "added: " << dh_dt << " "
           //           << psi(j,k) << " " 
           //           << maybe_inverse_capillary << " "
           //           << curv << " "
           //           << d2psidxi(j,k,0) << " " 
           //           << W << " "
           //           << " final: " << (dh_dt*psi(j,k)+maybe_inverse_capillary*curv*d2psidxi(j,k,0))*W
           //           << std::endl;
          }
        }
      }
     
    } // End of loop over the integration points
  }
 

private:

 /// Pointer to FSI parameter
 double* Q_pt;



 
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

 /// FSI parameter
 double Q_fsi=0.0;

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


 // Pin lubri (and reassign equation numbers)
 void pin_lubri()
  {
   unsigned nnod=Problem::mesh_pt()->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     mesh_pt()->node_pt(j)->pin(0);
     mesh_pt()->node_pt(j)->pin(1);
    }
   oomph_info << "Pinned lubri dofs: # of dofs = "
              << assign_eqn_numbers() << std::endl;
  }

 // Unpin lubri dofs (and reassign equation numbers)
 void unpin_lubri(const bool& no_flux_bc)
  {
   unsigned nnod=Problem::mesh_pt()->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     mesh_pt()->node_pt(j)->unpin(0);
     mesh_pt()->node_pt(j)->unpin(1);
    }

   // Now flux; height adjusts
   if (no_flux_bc)
    {
     mesh_pt()->boundary_node_pt(0,0)->pin(1);
     mesh_pt()->boundary_node_pt(1,0)->pin(1);
    }
   // Pinned ends; fluid drains out
   else
    {
     mesh_pt()->boundary_node_pt(0,0)->pin(0);
     mesh_pt()->boundary_node_pt(1,0)->pin(0);
    }
   
   oomph_info << "Unpinned lubri dofs: # of dofs = "
              << assign_eqn_numbers() << std::endl;
  }

 

 // hierher
 void set_initial_h_lubri(const double& h_mean, const double& h_amplitude)
  {
    // For now: Pin all lubri dofs and assign values
   // consistnt with constant film thickness
   unsigned nnod=Problem::mesh_pt()->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     double x=double(j)/double(nnod-1);
     double h=h_mean+h_amplitude*(1.0-cos(2.0*MathematicalConstants::Pi*x));
     // Value
     Problem::mesh_pt()->node_pt(j)->set_value(0,h);
     
     // Slope
     unsigned nel=mesh_pt()->nelement();
     double dh_ds=0.5*h_amplitude*2.0*MathematicalConstants::Pi/double(nel)*
      sin(2.0*MathematicalConstants::Pi*x);
    Problem::mesh_pt()->node_pt(j)->set_value(1,dh_ds);
  }
 
}


 /// Output stuff
 void doc_solution()
  {
   // Output file stream used for writing results
   ofstream file;
   
   // String used for the filename
   char filename[100]; 
   
   // Shape of beam and free surface
   sprintf(filename,"%s/beam%i.dat",
           Doc_info.directory().c_str(),
           Doc_info.number());
   file.open(filename);
   mesh_pt()->output(file,5);
   file.close();

   // Sum up energies and fluid volume
   unsigned nel=mesh_pt()->nelement();
   double kinetic_energy=0.0;
   double strain_energy=0.0;
   double v_lubri=0.0;
   for (unsigned e=0;e<nel;e++)
    {
     HermiteLubriBeamElement* el_pt=
      dynamic_cast<HermiteLubriBeamElement*>(
       mesh_pt()->element_pt(e));
     double el_strain_energy=0.0;
     double el_kin_energy=0.0;
     el_pt->get_energy(el_strain_energy,el_kin_energy);
     strain_energy+=el_strain_energy;
     kinetic_energy+=el_kin_energy;
     double el_v_lubri=0.0;
     el_pt->get_fluid_volume(el_v_lubri);
     v_lubri+=el_v_lubri;
    }

    // Write trace file
   Trace_file << Doc_info.number() << " "  // 1 
              << time_pt()->time() << " "  // 2
              << Global_Physical_Variables::P_ext  << " "  // 3
              << Global_Physical_Variables::Sigma0 << " "  // 4
              << Doc_node_pt->x(1) << " " // 5
              << Doc_node_pt->value(0) << " " // 6
              << kinetic_energy <<" " // 7
              << strain_energy <<" "  // 8
              << v_lubri << " "   // 9 
              << std::endl;
   
   // Bump
   Doc_info.number()++;
  }

 
 
private:

 /// Label for output
 DocInfo Doc_info;

 /// Trace file
 ofstream Trace_file;
 
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



   // hierher global parameter
   bool pin_lubri=false;
   if (pin_lubri)
    {
     mesh_pt()->boundary_node_pt(b,0)->pin(0);
    }
   // Symmetry/no flux:
   else
    {
     mesh_pt()->boundary_node_pt(b,0)->pin(1);
    }
  }


 // Assign uniform initial film thickness
 double h_initial=0.1; 
 unsigned nnod=Problem::mesh_pt()->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   // Value
   Problem::mesh_pt()->node_pt(j)->set_value(0,h_initial);

   // Slope
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
   elem_pt->q_pt() = &Global_Physical_Variables::Q_fsi;

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
 double pext_increment = -1.0e-7; // -0.00001;
 
 // Set initial values for control parameters 
 Global_Physical_Variables::P_ext = 0.0 - pext_increment;
 
 // Set the 2nd Piola Kirchhoff prestress
 Global_Physical_Variables::Sigma0=0.001;
 
 
 // Set output directory -- this function checks if the output
 // directory exists and issues a warning if it doesn't.
 Doc_info.set_directory("RESLT");
 
 // Open a trace file
 Trace_file.open("RESLT/trace_beam.dat");


 // hierher
 
 // // Write a header for the trace file
 // Trace_file <<  "VARIABLES=\"p_e_x_t\",\"sigma_0\""
 //            <<  ", \"d\"" << std::endl;


 // Pin lubri dofs during initial steady beam calculations (otherwise
 // free surface can fly away when using no flux bcs
 pin_lubri();

 // linear_solver_pt()=new FD_LU;

 // Doc initial condition/guess
 doc_solution();
 
    
 // STAGE 1: INFLATE, WITH GIVEN POSITIVE PRESTRESS
 //------------------------------------------------
 {

  // Switch off FSI
  Global_Physical_Variables::Q_fsi=0.0;
  
  // Set initial values for control parameters 
  Global_Physical_Variables::P_ext = 0.0 - pext_increment;
  
  // Loop over parameter increments
  unsigned nstep=10;
  for(unsigned i=1;i<=nstep;i++)
   {
    // Increment pressure
    Global_Physical_Variables::P_ext += pext_increment;
    
    oomph_info << "STAGE 1: Solving for p_ext sigma_0 = "
               << Global_Physical_Variables::P_ext << " " 
               << Global_Physical_Variables::Sigma0
               << std::endl;


    // DenseDoubleMatrix A;
    // DoubleVector b;
    // get_fd_jacobian(b,A);
    // A.sparse_indexed_output("FD_jac.dat");
    
    // get_jacobian(b,A);
    // A.sparse_indexed_output("jac.dat");

    // describe_dofs();
    //exit(0);
    
    // Solve the system
    steady_newton_solve();
    
    // Document the solution
    doc_solution();

    
   }
 }
 
 
 // STAGE 2: CHANGE SIGN OF PRESTRESS
 //----------------------------------
 {
  unsigned nstep=10;
  double d_sigma=
   2.0*Global_Physical_Variables::Sigma0/double(nstep-1);
  
  // Loop over parameter increments
  for(unsigned i=1;i<nstep;i++)
   {
    // Decrement pre-stress
    Global_Physical_Variables::Sigma0-=d_sigma;

    oomph_info << "STAGE 2: Solving for p_ext sigma_0 = "
               << Global_Physical_Variables::P_ext << " " 
               <<  Global_Physical_Variables::Sigma0
               << std::endl;
    
    // Solve the system
    steady_newton_solve();
    
    // Document the solution
    doc_solution();
    
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
    
    oomph_info << "STAGE 3: solving for p_ext sigma_0 = "
               << Global_Physical_Variables::P_ext << " " 
               <<  Global_Physical_Variables::Sigma0
               << std::endl;
    
    // Solve the system
    steady_newton_solve();
    
    // Document the solution
    doc_solution();

   }
 }



 //exit(0);
 
 // STAGE 4: TIMESTEP THE THING WITH FSI
 //-------------------------------------
 {
  // Set timestep
  double dt=0.01; // 1.0; 

  // Switch on FSI
  Global_Physical_Variables::Q_fsi=1.0e-4;
  
  // Set initial profile
  double h_mean=0.05; // 0.1;
  double h_amplitude=0.05;
  set_initial_h_lubri(h_mean,h_amplitude);

  // Unpin lubri dofs during initial steady beam calculationsc
  bool no_flux_bc=true;
  unpin_lubri(no_flux_bc);

  // Document the solution
  doc_solution();
  
  // Assign impulsive start
  assign_initial_values_impulsive(dt); // hierher try bypassing
  unsigned nstep=10000; // 1000;
  for (unsigned i=0;i<nstep;i++)
   {
    // Solve
    oomph_info << "STAGE 4: Doing unsteady solve for t = "
               << time_stepper_pt()->time_pt()->time()
               << std::endl;
    unsteady_newton_solve(dt);
    
    // Document the solution
    doc_solution();

        
   }
 }

 
} // end of parameter study

//========start_of_main================================================
/// Driver for beam (string under tension) test problem 
//=====================================================================
int main()
{

 // Set the non-dimensional thickness 
 Global_Physical_Variables::H=0.001; 
 
 // Set the length of domain // hierher get rid of this algotether
 double L = 1.0;

 // Number of elements (choose an even number if you want the control point 
 // to be located at the centre of the beam)
 unsigned n_element = 50; // 20;

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

