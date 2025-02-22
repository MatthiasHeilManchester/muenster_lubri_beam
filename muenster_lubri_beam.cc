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

//OOMPH-LIB includes
#include "generic.h"
#include "beam.h"
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;

using namespace oomph;

 
//========start_of_namespace========================
/// Namespace for parameters
//==================================================
namespace Global_Parameters
{

 // Beam parameters
 //-----------------

 /// Non-dimensional beam thickness
 double H_beam=0.05;

 /// 2nd Piola Kirchhoff pre-stress
 double Sigma0=0.0;

 /// Square of timescale ratio (i.e. non-dimensional density)  
 /// -- 1.0 for default value of scaling factor
 double Lambda_sq=1.0;

 /// Pressure load
 double P_ext=0.0;

 /// Load function: Apply a constant external pressure to the beam
 void load(const Vector<double>& xi, const Vector<double> &x,
           const Vector<double>& N, Vector<double>& load)
 {
  for(unsigned i=0;i<2;i++) {load[i] = -P_ext*N[i];}
 }

 /// Pressure during time-dependent simulation
 double P_ext_during_time_dependent_run=0.0;

 /// Time at which pressure is switched off during time-dependent
 /// simulation
 double T_switch_off_pressure_during_time_dependent_run=2.0;


 // FSI parameters
 //---------------
 
 /// FSI parameter (to be pointed to by element)
 double Q_fsi=0.0;

 /// Target FSI parameter for actual FSI run
 double Q_fsi_target=1.0e-4;



 // Lubrication theory parameters
 //------------------------------

 /// Scaled inverse capillary number
 double Scaled_inverse_capillary_number=1.0; 


 // Runtime parameters
 //-------------------

 /// Output directory
 std::string Dir_name="RESLT";

 /// Number of steady steps
 unsigned Nstep_steady=10;

 /// Max. pressure for steady solve
 double P_ext_max=-1.0e-6;

 /// Initial postive pre-stress
 double Sigma0_initial=0.005;
 
 /// Number of timesteps
 unsigned Ntstep=1000;

 /// Max. time for timestepping
 double T_max=1.0;

 

 // Parameters for validation of lubrication theory
 //-------------------------------------------------

 /// Validate solution of linearised equations?
 bool Use_linearised_flux=false;
 
 /// Mean film thickness for manufactured solution
 double H_lubri_mean_manufactured=0.1;

 /// Peak deviation of film thickness for manufactured solution
 double H_lubri_hat_manufactured=-0.05;

 /// Timescale for evolution of manufactured solution
 double T_lubri_manufactured=1.0;

 /// Manufactured solution for validation of lubrication theory
 double h_lubri_manufactured(const double& t, const double& x)
 {
  return H_lubri_mean_manufactured+
   H_lubri_hat_manufactured*
   cos(2.0*MathematicalConstants::Pi*x)*
   exp(-t/T_lubri_manufactured)*t*t;
 }

 /// Source term for manufactured solution
 double source_manufactured_solution(const double& t, const double& x)
 {
  double source=0.0;

  if (Use_linearised_flux)
   {
    source = 0.2e1 * H_lubri_hat_manufactured * cos(0.6283185308e1 *
    x) * t * exp(-t / T_lubri_manufactured) - H_lubri_hat_manufactured
    * cos(0.6283185308e1 * x) * t * t / T_lubri_manufactured * exp(-t
    / T_lubri_manufactured) + 0.5195151523e3 *
    Scaled_inverse_capillary_number * H_lubri_hat_manufactured *
    cos(0.6283185308e1 * x) * t * t * exp(-t / T_lubri_manufactured);
   }
  else
   {
    source= 0.2e1 * H_lubri_hat_manufactured * cos(0.6283185308e1 * x)
    * t * exp(-t / T_lubri_manufactured) - H_lubri_hat_manufactured *
    cos(0.6283185308e1 * x) * t * t / T_lubri_manufactured * exp(-t /
    T_lubri_manufactured) - 0.1558545457e4 *
    Scaled_inverse_capillary_number * pow(H_lubri_mean_manufactured +
    H_lubri_hat_manufactured * cos(0.6283185308e1 * x) * t * t *
    exp(-t / T_lubri_manufactured), 0.2e1) * H_lubri_hat_manufactured
    * H_lubri_hat_manufactured * pow(sin(0.6283185308e1 * x), 0.2e1) *
    pow(t, 0.4e1) * pow(exp(-t / T_lubri_manufactured), 0.2e1) +
    0.5195151523e3 * Scaled_inverse_capillary_number *
    pow(H_lubri_mean_manufactured + H_lubri_hat_manufactured *
    cos(0.6283185308e1 * x) * t * t * exp(-t / T_lubri_manufactured),
    0.3e1) * H_lubri_hat_manufactured * cos(0.6283185308e1 * x) * t *
    t * exp(-t / T_lubri_manufactured);
   }
  
  return source;
 }
 
} // end of namespace



/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


//=====================================================================
/// Upgraded Hermite Beam Element incorporates surface tension driven
/// thin film flow
//=====================================================================
class HermiteLubriBeamElement : public virtual HermiteBeamElement
{
 
public:
 
 /// Constructor: 
 HermiteLubriBeamElement() :
  Q_pt(0), Scaled_inverse_capillary_pt(0), Source_fct_pt(0),
  Manufactured_soln_fct_pt(0), Use_linearised_flux(false)
  {
  }


 ///  Required  # of `values' (pinned or actual dofs)
 /// at node n
 inline unsigned required_nvalue(const unsigned& n) const
  {
   return 2; 
  }

 
 /// Film thickness dof (local node j, type k)
 double nodal_h_lubri(const unsigned& j, const unsigned& k) const
  {
   return this->node_pt(j)->value(k);
  }
 
 /// Film thickness dof (local node j, type k) at previous time level t
 /// (t=0: present)
 double nodal_h_lubri(const unsigned& t, const unsigned& j, 
                      const unsigned& k) const
  {
   return this->node_pt(j)->value(t,k);
  }
 
 
 /// Film thickness at local coordinate s
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
         // account. All computed using (i) linearity (ii) ignoring
         // horizontal beam displacements (and stretching of beam!).
         curv+=(nodal_h_lubri(l, k)+raw_nodal_position_gen(l,k,1))*
          d2psidxi(l,k,0);

         // Number of timsteps (past & present)
         unsigned n_time =node_pt(l)->time_stepper_pt()->ntstorage();
         
         // Add the contributions to the time derivative
         for (unsigned t = 0; t < n_time; t++)
          {
           dh_dt+=node_pt(l)->time_stepper_pt()->weight(1, t)*
            nodal_h_lubri(t,l,k);
          }
        }
      }
      
     // Manufactured solution
     double t=node_pt(0)->time_stepper_pt()->time_pt()->time(); 
     double h_lubri_manufactured=0.0;
     if (Manufactured_soln_fct_pt!=0)
      {
       h_lubri_manufactured=Manufactured_soln_fct_pt(t,posn[0]);
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
             << h_lubri_manufactured << " " // 9
             << source(t,posn[0]) << " " // 10  
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


 /// Pointer to scaled inverse capillary number (const version)
 double* scaled_inverse_capillary_pt() const
  {
   return Scaled_inverse_capillary_pt;
  }
 

 /// Pointer to scaled inverse capillary number (read/write)
 double*& scaled_inverse_capillary_pt()
  {
   return Scaled_inverse_capillary_pt;
  }
 
 
 /// Scaled inverse capillary number 
 double scaled_inverse_capillary() const
  {
   if (Scaled_inverse_capillary_pt==0)
    {
     return 0.0;
    }
   return *Scaled_inverse_capillary_pt;
  }

 
 /// Function pointer to source function
 typedef double (*LubriSourceFctPt)(const double& t, const double& x);
 
 /// Pointer to source function (const version)
 LubriSourceFctPt source_fct_pt() const
  {
   return Source_fct_pt;
  }

 /// Pointer to source function (read/write)
 LubriSourceFctPt& source_fct_pt()
  {
   return Source_fct_pt;
  }
  
 /// Source function
 double source(const double& t, const double& x) const
  {
   if (Source_fct_pt==0)
    {
     return 0.0;
    }
   return Source_fct_pt(t,x);
  }


 /// Switch on linear flux (for validation)
 void enable_use_linearised_flux()
  {
   Use_linearised_flux=true;
  }
 
 /// Switch off linear flux
 void disable_use_linearised_flux()
  {
   Use_linearised_flux=false;
  }


 /// Function pointer to manufactured solution
 typedef double (*ManufacturedSolnFctPt)(const double& t, const double& x);
 
 /// Pointer to source function (const version)
 ManufacturedSolnFctPt manufactured_soln_fct_pt() const
  {
   return Manufactured_soln_fct_pt;
  }

 /// Pointer to source function (read/write)
 ManufacturedSolnFctPt& manufactured_soln_fct_pt()
  {
   return Manufactured_soln_fct_pt;
  }
  
 /// Manufactured solutionL h(t,x)
 double manufactured_soln_e(const double& t, const double& x) const
  {
   if (Manufactured_soln_fct_pt==0)
    {
     return 0.0;
    }
   return Manufactured_soln_fct_pt(t,x);
  }


 
 
 /// Fill in contribution to residuals: Combine the two single
 /// physics ones.
 void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
   fill_in_contribution_to_residuals_beam(residuals);
   fill_in_contribution_to_residuals_lubri(residuals);
  }

 /// Get the Jacobian and residuals. Overwrite the version in
 /// the beam class and reinstate brute-force finite differencing
 virtual void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian)
  {
   // Call FD version from base class
   SolidFiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
  }

 
 /// Overloaded version of the load vector function to incorporate lubrication
 /// pressure: Generic interface; pass number of integration point,
 /// Lagr. and Eulerian coordinate and normal vector and return the load
 /// vector.
 void load_vector(const unsigned& ipt,
                  const Vector<double>& xi,
                  const Vector<double>& x,
                  const Vector<double>& N,
                  Vector<double>& load)
  {
   
   // Call the version from the underlying beam element; this gets
   // the external pressure
   HermiteBeamElement::load_vector(ipt,xi,x,N,load);
   
   // Now add lubri traction
   //-----------------------
   
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
       // Curvature of free surface, taking substrate (beam) curvature into
       // account. All computed using (i) linearity (ii) ignoring
       // horizontal beam displacements (and stretching of beam!).
       curv+=(nodal_h_lubri(l, k)+raw_nodal_position_gen(l,k,1))*
        d2psidxi(l,k,0);
      }
    }

   // Add FSI load
   double q_fsi=q(); 
   load[0]+=q_fsi*curv*N[0];
   load[1]+=q_fsi*curv*N[1];
   
  }

 /// Volume (per unit depth) of fluid in the film
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

 
 /// Fill in contribution to residuals arising from lubrication theory
 void fill_in_contribution_to_residuals_lubri(Vector<double>& residuals)
  {

   // Version in base class provides this capability; haven't bothered
   // to reimplement it here.
   
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
   
   // Get parameters from Element

   // Cache inverse capillary number
   double inverse_capillary=scaled_inverse_capillary();

   // Get continuous time for evaluation of source function
   double t=node_pt(0)->time_stepper_pt()->time_pt()->time();
   
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
     double x=0.0;
     
     // Calculate values
     for (unsigned l = 0; l < n_node; l++)
      {
       // Loop over dof types
       for (unsigned k = 0; k < n_position_type; k++)
        {
         x+=raw_dnodal_position_gen_dt(0, l, k, 0) * psi(l, k);
         
         // filn thickness
         h_lubri+=nodal_h_lubri(l, k) * psi(l, k);
         
         // slope 
         dh_lubri_dxi+=nodal_h_lubri(l, k)*dpsidxi(l,k,0);
         
         // Curvature of free surface, taking substrate (beam) curvature into
         // account. Simplified: (i) linear (ii) ignore horizontal
         // beam displacement!
         curv+=(nodal_h_lubri(l, k)+raw_nodal_position_gen(l,k,1))*
          d2psidxi(l,k,0);
         
         // Number of timsteps (past & present)
         unsigned n_time =node_pt(l)->time_stepper_pt()->ntstorage();
         
         // Add the contributions to the time derivative
         for (unsigned t = 0; t < n_time; t++)
          {
           dh_dt+=node_pt(l)->time_stepper_pt()->weight(1,t)*
            nodal_h_lubri(t,l,k);
          }
        }
      }
          
     // Evaluate source function 
     double lubri_source=source(t,x);
     
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
           if (Use_linearised_flux)
            {
             residuals[local_eqn] +=
              ((dh_dt-lubri_source)*psi(j,k)+inverse_capillary*curv*
               1.0/3.0*d2psidxi(j,k,0))*W;
            }
           else
            {
             residuals[local_eqn] +=
              (
               (dh_dt-lubri_source)*psi(j,k)+
               inverse_capillary*curv*
               (h_lubri*h_lubri*dh_lubri_dxi * dpsidxi(j,k,0)+
                h_lubri*h_lubri*h_lubri/3.0  *d2psidxi(j,k,0) )
               )*W;
            }

          }
        }
      }
     
    } // End of loop over the integration points
  }
 

private:

 /// Pointer to FSI parameter
 double* Q_pt;

 /// Pointer to scaled inverse capillary number
 double* Scaled_inverse_capillary_pt;

 /// Pointer to source function
 LubriSourceFctPt Source_fct_pt;
 
 /// Pointer to manufactured solution (for validation)
 ManufacturedSolnFctPt Manufactured_soln_fct_pt;

 /// Use linearised flux (mainly for validatino)
 bool Use_linearised_flux;
      
};



/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////




//======start_of_problem_class==========================================
/// Beam problem object
//======================================================================
class MuensterLubriBeamProblem : public Problem
{
public:
 
 /// Constructor: The arguments are the number of elements
 MuensterLubriBeamProblem(const unsigned &n_elem);
 
 /// Conduct a parameter study
 void parameter_study();

 /// Validate lubrication theory
 void validate_lubri();

 /// Return pointer to the mesh (overloaded to return the
 /// actual Mesh type).
 OneDLagrangianMesh<HermiteLubriBeamElement>* mesh_pt() 
  {return dynamic_cast<OneDLagrangianMesh<HermiteLubriBeamElement>*>
    (Problem::mesh_pt());}

 /// No actions need to be performed after a solve
 void actions_after_newton_solve() {}

 /// No actions need to be performed before a solve
 void actions_before_newton_solve() {}

 /// Pin lubri dofs (and reassign equation numbers)
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

 /// Unpin lubri dofs (and reassign equation numbers)
 /// Boolean decides what type of BCs are to be applied
 void unpin_lubri(const bool& no_flux_bc)
  {
   unsigned nnod=Problem::mesh_pt()->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     mesh_pt()->node_pt(j)->unpin(0);
     mesh_pt()->node_pt(j)->unpin(1);
    }
   
   // No flux; height adjusts
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
 

 /// Set initial film profile: mean thickness h_mean
 /// perturbed by a cos (or sin) profile with amplitude
 /// h_amplitude
 void set_initial_h_lubri(const double& h_mean,
                          const double& h_amplitude,
                          const bool& use_cos)
  {
   // mapping for slope dofs
   double s_min=mesh_pt()->finite_element_pt(0)->s_min();
   double s_max=mesh_pt()->finite_element_pt(0)->s_max();
   double length=1.0;
   unsigned nel=mesh_pt()->nelement();
   double dx_ds=(length/double(nel))/(s_max-s_min);
   
   // visit nodes
   unsigned nnod=Problem::mesh_pt()->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     double x=double(j)/double(nnod-1);
     
     double h=0.0;
     if (use_cos)
      {
       h=h_mean+h_amplitude*cos(2.0*MathematicalConstants::Pi*x);
      }
     else
      {
       h=h_mean+h_amplitude*sin(MathematicalConstants::Pi*x);
      }
     
     // Value
     mesh_pt()->node_pt(j)->set_value(0,h);
     
     // Slope
     double dh_dx=0.0;
     if (use_cos)
      {
       dh_dx=-h_amplitude*2.0*MathematicalConstants::Pi*
        sin(2.0*MathematicalConstants::Pi*x);
      }
     else
      {
       dh_dx=h_amplitude*MathematicalConstants::Pi*
        cos(MathematicalConstants::Pi*x);
      }
      
     double dh_ds=dh_dx*dx_ds;
     mesh_pt()->node_pt(j)->set_value(1,dh_ds);
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
   double kinetic_energy=0.0;
   double strain_energy=0.0;
   double v_lubri=0.0;
   unsigned nel=mesh_pt()->nelement();
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
   Trace_file
    << Doc_info.number() << " "  // 1 
    << time_pt()->time() << " "  // 2
    << Global_Parameters::P_ext  << " "  // 3
    << Global_Parameters::Sigma0 << " "  // 4
    << Doc_node_pt->x(1) << " " // 5
    << Doc_node_pt->value(0) << " " // 6
    << kinetic_energy <<" " // 7
    << strain_energy <<" "  // 8
    << v_lubri << " "   // 9
    << Global_Parameters::h_lubri_manufactured(
     time_pt()->time(),
     Doc_node_pt->x(0)) << " " // 10
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

 /// Pointer to geometric object that represents the beam's undeformed shape
 GeomObject* Undef_beam_pt;

}; // end of problem class


//=============start_of_constructor=====================================
/// Constructor for elastic beam problem
//======================================================================
MuensterLubriBeamProblem::MuensterLubriBeamProblem(const unsigned &n_elem)
{

 // Taking away the pressure for negatively pre-stressed beam can be
 // a bit flakey...
 Problem::Max_newton_iterations=100;
 
 // Set the undeformed beam to be a straight line at y=0
 Undef_beam_pt=new StraightLine(0.0); 

 // Create the timestepper and add it to the Problem's collection of
 // timesteppers -- this creates the Problem's Time object.
 add_time_stepper_pt(new NewmarkBDF<2>); 

 
 // Create the (Lagrangian!) mesh, using the geometric object
 // Undef_beam_pt to specify the initial (Eulerian) position of the
 // nodes.
 double length=1.0;
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

   // For now apply symmetry/no flux conditions for lubrication theory
   // hierher add to talk
   mesh_pt()->boundary_node_pt(b,0)->pin(1);
  }


 // Assign uniform initial film thickness
 double h_mean=0.1;
 double h_amplitude=0.0;
 bool use_cos=true;
 set_initial_h_lubri(h_mean,h_amplitude,use_cos);
 
 //Find number of elements in the mesh
 unsigned n_element = mesh_pt()->nelement();
 
 //Loop over the elements to set physical parameters etc.
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   HermiteLubriBeamElement *elem_pt = 
    dynamic_cast<HermiteLubriBeamElement*>(mesh_pt()->element_pt(e));
   
   // Set physical parameters for each element:
   elem_pt->sigma0_pt() = &Global_Parameters::Sigma0;
   elem_pt->h_pt() = &Global_Parameters::H_beam;
   elem_pt->q_pt() = &Global_Parameters::Q_fsi;

   // Set the load Vector for each element
   elem_pt->load_vector_fct_pt() = &Global_Parameters::load;

   // Pass pointer to square of timescale ratio (non-dimensional density)
   elem_pt->lambda_sq_pt() = &Global_Parameters::Lambda_sq;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = Undef_beam_pt;

   // Set pointer to scaled inverse capillary number
   elem_pt->scaled_inverse_capillary_pt()=
    &Global_Parameters::Scaled_inverse_capillary_number;

   // Source function and manufactured solution for validation
   if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
    {
     elem_pt->source_fct_pt()=
      &Global_Parameters::source_manufactured_solution;
     elem_pt->manufactured_soln_fct_pt()=
      &Global_Parameters::h_lubri_manufactured;
     if (Global_Parameters::Use_linearised_flux)
      {
       elem_pt->enable_use_linearised_flux();
      }
     else
      {
       elem_pt->disable_use_linearised_flux();
      }
    }
   
  } // end of loop over elements

 // Choose node at which displacement is documented (halfway along -- provided
 // we have an odd number of nodes; complain if this is not the
 // case)
 unsigned n_nod=mesh_pt()->nnode();
 if (n_nod%2!=1)
  {
   cout << "Warning: Even number of nodes " << n_nod << std::endl;
   cout << "Comparison with exact solution will be misleading..."
        << std::endl;
  }
 Doc_node_pt=mesh_pt()->node_pt((n_nod+1)/2-1);
 
 // Assign the global and local equation numbers
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;

} // end of constructor









//=======start_of_parameter_study==========================================
/// Solver loop to perform parameter study
//=========================================================================
void MuensterLubriBeamProblem::parameter_study()
{
   
 // Set output directory -- this function checks if the output
 // directory exists and issues a warning if it doesn't.
 Doc_info.set_directory(Global_Parameters::Dir_name);
 
 // Open a trace file
 Trace_file.open(Global_Parameters::Dir_name+"/trace_beam.dat");

 // Switch linear solver in total desperation
 // linear_solver_pt()=new FD_LU;


 
    
 // STAGE 1: INFLATE, WITH GIVEN POSITIVE PRESTRESS
 //------------------------------------------------
 {

  // Pin lubri dofs during initial steady beam calculations (otherwise
  // free surface can fly away when using no flux bcs)
  pin_lubri();
  
  // Parameters for Inflation stage:
  double pext_increment = Global_Parameters::P_ext_max/
   double(Global_Parameters::Nstep_steady); // -1.0e-7; // hierher read in p_ext_max

  // Initial pressure
  Global_Parameters::P_ext = 0.0 - pext_increment;
  
  // Set the 2nd Piola Kirchhoff prestress (tensile)
  Global_Parameters::Sigma0=Global_Parameters::Sigma0_initial; // 0.005;  // hierher read in sigma0_max

  // Switch off FSI
  Global_Parameters::Q_fsi=0.0;
  
  // Set initial values for control parameters 
  Global_Parameters::P_ext = 0.0 - pext_increment;
  
  // Doc initial condition/guess
  doc_solution();
  
  // Loop over parameter increments
  for(unsigned i=1;i<=Global_Parameters::Nstep_steady;i++)
   {
    // Increment pressure
    Global_Parameters::P_ext += pext_increment;
    
    oomph_info << "STAGE 1: Solving for p_ext sigma_0 = "
               << Global_Parameters::P_ext << " " 
               << Global_Parameters::Sigma0
               << std::endl;


    
    // why is the Jacobian singular? Let's have a look...
    
    // DenseDoubleMatrix A;
    // DoubleVector b;
    // get_fd_jacobian(b,A);
    // A.sparse_indexed_output("FD_jac.dat");
    
    // get_jacobian(b,A);
    // A.sparse_indexed_output("jac.dat");

    // Ha! there's a zero row; which dofs is it associated with?
    // describe_dofs();
    
    // exit(0);
    
    // Solve the system
    steady_newton_solve();
    
    // Document the solution
    doc_solution();
    
   }
 }
 
 
 // STAGE 2: CHANGE SIGN OF PRESTRESS
 //----------------------------------
 {
  double d_sigma=
   2.0*Global_Parameters::Sigma0/double(Global_Parameters::Nstep_steady-1);
  
  // Loop over parameter increments
  for(unsigned i=1;i<Global_Parameters::Nstep_steady;i++)
   {
    // Decrement pre-stress
    Global_Parameters::Sigma0-=d_sigma;

    oomph_info << "STAGE 2: Solving for p_ext sigma_0 = "
               << Global_Parameters::P_ext << " " 
               <<  Global_Parameters::Sigma0
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
  double pext_increment=Global_Parameters::P_ext/
   double(Global_Parameters::Nstep_steady);
  for(unsigned i=1;i<=Global_Parameters::Nstep_steady;i++)
   {
    // Increment pressure
    Global_Parameters::P_ext -= pext_increment;
    
    oomph_info << "STAGE 3: solving for p_ext sigma_0 = "
               << Global_Parameters::P_ext << " " 
               <<  Global_Parameters::Sigma0
               << std::endl;
    
    // Solve the system
    steady_newton_solve();
    
    // Document the solution
    doc_solution();

   }
 }


 // Bail?
 if (CommandLineArgs::command_line_flag_has_been_set("--steady_only"))
  {
   oomph_info << "Bailing out after steady solves as requested" << std::endl;
   exit(0);
  }

 
 // STAGE 4: TIMESTEP THE THING WITH FSI
 //-------------------------------------
 {
  // Set timestep
  double dt=Global_Parameters::T_max/double(Global_Parameters::Ntstep);
 
  // Switch on FSI
  Global_Parameters::Q_fsi=
   Global_Parameters::Q_fsi_target;

  // Set pressure during time-dependent run (this is in addition to
  // any pressure from the fluid)
  Global_Parameters::P_ext=
   Global_Parameters::P_ext_during_time_dependent_run; 
    
  // Set initial profile 
  double h_mean=0.1;
  double h_amplitude=-0.05;
  bool use_cos=true;
  bool no_flux_bc=true;
 if (CommandLineArgs::command_line_flag_has_been_set("--pin_h_lubri_at_zero"))
  {
   use_cos=false;
   h_mean=0.0;
   h_amplitude=0.05;
   no_flux_bc=false;
  }
  set_initial_h_lubri(h_mean,h_amplitude,use_cos);

  // Unpin lubri dof
  unpin_lubri(no_flux_bc);

  // Document the solution
  doc_solution();
  
  // Assign impulsive start
  assign_initial_values_impulsive(dt); 
  for (unsigned i=0;i<Global_Parameters::Ntstep;i++)
   {
    // Solve
    oomph_info << "STAGE 4: Doing unsteady solve for t = "
               << time_stepper_pt()->time_pt()->time()
               << " for p_ext sigma_0 = "
               << Global_Parameters::P_ext << " " 
               <<  Global_Parameters::Sigma0
               << std::endl;
    unsteady_newton_solve(dt);
    
    // Document the solution
    doc_solution();
    
    // Switch off pressure
    if ((Global_Parameters::P_ext!=0.0)&&
        (time_stepper_pt()->time_pt()->time()>
         Global_Parameters::T_switch_off_pressure_during_time_dependent_run))
     {
      oomph_info << "Switching off pressure\n";
      Global_Parameters::P_ext=0.0;
     }
    
   }
 }

 
} // end of parameter study



//=======start_of_validate_lubri===========================================
/// Validate lubrication theory: Flat substrate, cos-shaped initial film
/// profile
//=========================================================================
void MuensterLubriBeamProblem::validate_lubri()
{
 
 // Set initial values for control parameters 
 Global_Parameters::P_ext = 0.0;
 
 // Set the 2nd Piola Kirchhoff prestress
 Global_Parameters::Sigma0=0.001;
 
 // Switch off FSI
 Global_Parameters::Q_fsi=0.0;
 
 // Set output directory -- this function checks if the output
 // directory exists and issues a warning if it doesn't.
 Doc_info.set_directory(Global_Parameters::Dir_name);
 
 // Open a trace file
 Trace_file.open(Global_Parameters::Dir_name+"/trace_beam.dat");
 
 // Set timestep
 double dt=Global_Parameters::T_max/double(Global_Parameters::Ntstep);
 
 // Set initial profile
 double h_mean=0.1;
 double h_amplitude=0.0;
 bool use_cos=true;
 set_initial_h_lubri(h_mean,h_amplitude,use_cos);
 
 // Unpin lubri dofs 
 bool no_flux_bc=true;
 unpin_lubri(no_flux_bc);
 
 // Document the solution
 doc_solution();
 
 // Assign impulsive start
 assign_initial_values_impulsive(dt); 
 for (unsigned i=0;i<Global_Parameters::Ntstep;i++)
  {
   // Solve
   oomph_info << "VALIDATION: Doing unsteady solve for t = "
              << time_stepper_pt()->time_pt()->time()
              << std::endl;
   unsteady_newton_solve(dt);
   
   // Document the solution
   doc_solution();   
  }

} // end of validation

//========start_of_main================================================
/// Driver 
//=====================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Do validation?
 CommandLineArgs::specify_command_line_flag("--validate");

 // Do validation with linearised solution?
 CommandLineArgs::specify_command_line_flag("--linearised_flux");

 // Name of output directory
 CommandLineArgs::specify_command_line_flag(
  "--dir_name",&Global_Parameters::Dir_name);
 
 // Number of elements (choose an even number if you want the control point 
 // to be located at the centre of the beam)
 unsigned n_element = 50;
 CommandLineArgs::specify_command_line_flag(
  "--nelement",&n_element);

 /// Initial postive pre-stress
 CommandLineArgs::specify_command_line_flag(
  "--sigma0_initial",&Global_Parameters::Sigma0_initial);
 
 // Number of timesteps
 CommandLineArgs::specify_command_line_flag(
  "--ntstep",&Global_Parameters::Ntstep);

 // Max. time for timestepping
 CommandLineArgs::specify_command_line_flag(
  "--t_max",&Global_Parameters::T_max);

 // Only steady computations in actual problem?
 CommandLineArgs::specify_command_line_flag("--steady_only");

 // Target FSI parameter
 CommandLineArgs::specify_command_line_flag(
  "--q_fsi_target",
  &Global_Parameters::Q_fsi_target);

  // Pressure during time-dependent run (this is in addition to
  // any pressure from the fluid)
 CommandLineArgs::specify_command_line_flag(
  "--p_ext_kick",
  &Global_Parameters::P_ext_during_time_dependent_run); 
 
 // Time at which additional pressure is switched off
 CommandLineArgs::specify_command_line_flag(
  "--t_switch_off_kick",
  &Global_Parameters::T_switch_off_pressure_during_time_dependent_run);
     
 // film pinned at zero (otherwise no flux/symmetry)
 CommandLineArgs::specify_command_line_flag("--pin_h_lubri_at_zero");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 // Linearised validation? 
 if (CommandLineArgs::command_line_flag_has_been_set
     ("--linearised_flux"))
  {
   Global_Parameters::Use_linearised_flux=true;
  }
 
   
 // Set the non-dimensional thickness 
 Global_Parameters::H_beam=0.001; 
 
 // Construst the problem
 MuensterLubriBeamProblem problem(n_element);
 
 // Validate lubrication theory?
 //-----------------------------
 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   problem.validate_lubri();
   exit(0);
  }

 // Conduct parameter study
 //------------------------
 problem.parameter_study();

} // end of main


