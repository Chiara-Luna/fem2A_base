#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };
    
    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        if ( border ){
        for (int v =0; v<2; ++v) {
        vertices_.push_back(M.get_edge_vertex(i,v));
        }
        }
        else {
        for (int v =0; v<3; ++v) {
        vertices_.push_back(M.get_triangle_vertex(i,v));
        }
        }    
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        vertex r ;
        if (border_){ 
        r.x= (1- x_r.x) * vertices_[0].x + x_r.x* vertices_[1].x;
        r.y= (1- x_r.x) * vertices_[0].y + x_r.x* vertices_[1].y;
        }
        else {
        r.x= (1- x_r.x - x_r.y) * vertices_[0].x + x_r.x* vertices_[1].x +x_r.y* vertices_[2].x;
        r.y= (1- x_r.x - x_r.y) * vertices_[0].y + x_r.x* vertices_[1].y +x_r.y* vertices_[2].y;
        }
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        DenseMatrix J ;
        if (border_){ 
        J.set_size(2,1);
        J.set(0,0,-vertices_[0].x + vertices_[1].x);
        J.set(1,0,-vertices_[0].y + vertices_[1].y);
        }
        else {
        J.set_size(2,2);
        J.set(0,0,-vertices_[0].x + vertices_[1].x);
        J.set(1,0,-vertices_[0].y + vertices_[1].y);
        J.set(0,1,-vertices_[0].x + vertices_[2].x);
        J.set(1,1,-vertices_[0].y + vertices_[2].y);
        }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        double det;
        DenseMatrix J = ElementMapping::jacobian_matrix(x_r );
        if (border_){ 
        DenseMatrix T = J.transpose();
        det = sqrt(T.get(0,0)*J.get(0,0)+T.get(0,1)*J.get(1,0));
        }
        else {
        det = J.get(0,0)*J.get(1,1) - J.get(1,0)*J.get(0,1);
        }
        return det ;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order ){}

    int ShapeFunctions::nb_functions() const
    {
        int num = dim_+1;
        return num;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        float phi;
        if (i==0){
        phi = 1- x_r.x - x_r.y;
        }
        if (i==1){
        phi = x_r.x;
        }
        if (i==2){
        phi = x_r.y;
        }
        return phi;
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        Mesh mesh;
        mesh.load("data/square.mesh");
        ElementMapping element= ElementMapping(mesh, false, 4);
        vec2 g ;
        vec2 shape_coeff;
        DenseMatrix J = element.jacobian_matrix(x_r );
        
        // for triangle 
        if (dim_ ==2){
        DenseMatrix I = J.invert_2x2();
        DenseMatrix T = I.transpose();
        if (i==0){
        shape_coeff.x = -1;
        shape_coeff.y = -1;
        }
        if (i==1){
        shape_coeff.x = 1;
        shape_coeff.y = 0;
        }
        if (i==2){
        shape_coeff.x = 0;
        shape_coeff.y = 1;
        }
        g = T.mult_2x2_2( shape_coeff );
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
        Ke.set_size(3,3);
        for (int i = 0; i<3; ++i){
        	for (int j = 0; j<3; ++j){
        		double sum = 0;
        		for (int q = 0; q< quadrature.nb_points() ; ++q){
        			vertex pt = quadrature.point(q);
        			double scalar_prod = dot (reference_functions.evaluate_grad(i,pt), reference_functions.evaluate_grad(j,pt));
        			sum += quadrature.weight(q) * coefficient(pt) * scalar_prod * elt_mapping.jacobian(pt);
        }
        		Ke.set(i,j,sum);
        		}
        	}
    }

    void local_to_global_matrix(const Mesh& M,int t,const DenseMatrix& Ke,SparseMatrix& K )
    {
        for (int i = 0; i<3; ++i){
        	int I = M.get_triangle_vertex_index(t, i);
        	for (int j = 0; j<3; ++j){
        		double to_add = Ke.get(i,j);
        		int J = M.get_triangle_vertex_index(t, j );
        		K.add(I,J,to_add);
    		}
    	}
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (source term)" << '\n';
        for (int i = 0; i<2; ++i){
        	double sum = 0;
        	for (int q = 0; q< quadrature.nb_points() ; ++q){
        		vertex pt = quadrature.point(q);
        		vertex M_pt = elt_mapping.transform(pt);
        		sum += quadrature.weight(q) * reference_functions.evaluate(i, pt) * source(M_pt) * elt_mapping.jacobian(pt);
        		std::cout<< "sum "<< sum << std::endl;
        	}
        	Fe[i] = sum;
        	std::cout<< "Fe["<<i<<"] : "<< Fe[i]<< std::endl;
        }
    }

    /*void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        std::cout << "Fe -> F" << '\n';
        // TODO
    }*/

    void apply_dirichlet_boundary_conditions(
	const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::vector< bool > redondance (M.nb_vertices(),false);
        int penalty = 1000;
        for (int i=0; i< M.nb_edges(); ++i){
        	int j = M.get_edge_attribute(i);
        	for (int k = 0; k<2;++k){
        		int ind_vertice = M.get_edge_vertex_index(i,k);
        		if (attribute_is_dirichlet[j] && not redondance[ind_vertice]){
        			K.add(ind_vertice,ind_vertice,penalty);
        			F[ind_vertice]+= penalty* values[j];
        			redondance[ind_vertice] = true;
        		}
        	}
        }
    }

    /*void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }*/

}
