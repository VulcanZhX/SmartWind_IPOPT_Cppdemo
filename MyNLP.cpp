// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "MyNLP.hpp"
#include <cassert>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

using namespace Ipopt;

/* Constructor. */
MyNLP::MyNLP(const WindFarmOptimization& farmopt_cpy)
	: farmopt(farmopt_cpy)
{
	// No additional initialization required
}


MyNLP::~MyNLP()
{ }

bool MyNLP::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
   // The problem described in MyNLP.hpp has n variables
   n = MaxTurbines;

   // m equality constraint,
   m = 2;

   // nonzeros in the jacobian
   nnz_jac_g = m * n; // for simplicity, the jacobian matrix is assumed to be dense

   //// and 2 nonzeros in the hessian of the lagrangian
   //// (one in the hessian of the objective for x2,
   ////  and one in the hessian of the constraints for x1)
   //nnz_h_lag = 2;

   // We use the standard fortran index style for row/col entries
   index_style = C_STYLE;

   return true;
}

bool MyNLP::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
)
{
   // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
   // If desired, we could assert to make sure they are what we think they are.
   assert(n == MaxTurbines);
   assert(m == 2);

   // lower bound and upper bound on x
   for (Index i = 0; i < n; i++) {
	   x_l[i] = farmopt.yaw_lower; // x lower bound
	   x_u[i] = farmopt.yaw_upper; // x upper bound
   }

   // we have two inequality constraints. the lower bounds are zero
   g_l[0] = 0.0;
   g_l[1] = 0.0;

   return true;
}

bool MyNLP::get_starting_point(
   Index   n,
   bool    init_x,
   Number* x,
   bool    init_z,
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda
)
{
   // Here, we assume we only have starting values for x, if you code
   // your own NLP, you can provide starting values for the others if
   // you wish.
   assert(init_x == true);
   assert(init_z == false);
   assert(init_lambda == false);

   // we initialize x in bounds, in the upper right quadrant
   for (Index i = 0; i < n; i++) {
       x[i] = 0.0;
   }

   return true;
}

bool MyNLP::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
   // return the value of the objective function
	obj_value = objPower(const_cast<Number*>(x));

   return true;
}

bool MyNLP::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
   // return the gradient of the objective function grad_{x} f(x)

   // grad_{x1} f(x): x1 is not in the objective
   grad_f[0] = 0.0;

   // grad_{x2} f(x):
   Number x2 = x[1];
   grad_f[1] = -2.0 * (x2 - 2.0);

   return true;
}

bool MyNLP::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
   // return the value of the constraints: g(x)
   Number x1 = x[0];
   Number x2 = x[1];

   g[0] = -(x1 * x1 + x2 - 1.0);

   return true;
}

bool MyNLP::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
   if( values == NULL )
   {
      // return the structure of the jacobian of the constraints

      // element at 1,1: grad_{x1} g_{1}(x)
      iRow[0] = 1;
      jCol[0] = 1;

      // element at 1,2: grad_{x2} g_{1}(x)
      iRow[1] = 1;
      jCol[1] = 2;
   }
   else
   {
      // return the values of the jacobian of the constraints
      Number x1 = x[0];

      // element at 1,1: grad_{x1} g_{1}(x)
      values[0] = -2.0 * x1;

      // element at 1,2: grad_{x1} g_{1}(x)
      values[1] = -1.0;
   }

   return true;
}

//bool MyNLP::eval_h(
//   Index         n,
//   const Number* x,
//   bool          new_x,
//   Number        obj_factor,
//   Index         m,
//   const Number* lambda,
//   bool          new_lambda,
//   Index         nele_hess,
//   Index*        iRow,
//   Index*        jCol,
//   Number*       values
//)
//{
//   if( values == NULL )
//   {
//      // return the structure. This is a symmetric matrix, fill the lower left
//      // triangle only.
//
//      // element at 1,1: grad^2_{x1,x1} L(x,lambda)
//      iRow[0] = 1;
//      jCol[0] = 1;
//
//      // element at 2,2: grad^2_{x2,x2} L(x,lambda)
//      iRow[1] = 2;
//      jCol[1] = 2;
//
//      // Note: off-diagonal elements are zero for this problem
//   }
//   else
//   {
//      // return the values
//
//      // element at 1,1: grad^2_{x1,x1} L(x,lambda)
//      values[0] = -2.0 * lambda[0];
//
//      // element at 2,2: grad^2_{x2,x2} L(x,lambda)
//      values[1] = -2.0 * obj_factor;
//
//      // Note: off-diagonal elements are zero for this problem
//   }
//
//   return true;
//}

void MyNLP::finalize_solution(
   SolverReturn               status,
   Index                      n,
   const Number*              x,
   const Number*              z_L,
   const Number*              z_U,
   Index                      m,
   const Number*              g,
   const Number*              lambda,
   Number                     obj_value,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
)
{
   // here is where we would store the solution to variables, or write to a file, etc
   // so we could use the solution. Since the solution is displayed to the console,
   // we currently do nothing here.
}


// private methods imported

double MyNLP::objPower(Number* x) {
    // convert x to vector<double>
	std::vector<double> x_vec(MaxTurbines); // to be fixed
    for (Index i = 0; i < MaxTurbines; i++) {
        x_vec[i] = x[i];
    }
	// set yaw angles
	farmopt.setYawAngles(Eigen::Map<Eigen::VectorXd>(x_vec.data(), x_vec.size()));
    return -farmopt.getFarmPower();
}

// g(x) (constriants) implementation can be added here if needed

// grad_f and grad_g implementations can be added here if needed