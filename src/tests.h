#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quadrature()
        {
        	/*Quadrature quad = Quadrature::get_quadrature(0, false);
         	std::cout<< quad.nb_points() << std::endl;
         	double sum =0;
         	for (int i= 0; i < quad.nb_points(); ++i){
         	std::cout<< quad.nb_points(i).x << std::endl;
         	std::cout<< quad.nb_points(i).y << std::endl;
         	std::cout<< quad.weight(i) << std::endl;
         	sum += quad.weight(i); 
         	}
         	std::cout<< sum << std::endl;*/
        	return true;
        }
        
        bool test_class_ElementMapping()
        {
        	Mesh mesh;
            	mesh.load("data/square.mesh");
        	ElementMapping element= ElementMapping(mesh, false, 4);
        	vertex x_r;
        	x_r.x = 0.2;
        	x_r.y = 0.4;
        	element.transform(x_r);
        	element.jacobian_matrix(x_r);
        	element.jacobian(x_r);
        	return true;
        }
        
        bool test_class_ShapeFunctions()
        {
        	ShapeFunctions func = ShapeFunctions(2,1);
        	func.nb_functions();
        	vertex x_r;
        	x_r.x = 0.2;
        	x_r.y = 0.4;
        	func.evaluate(1,x_r);
        	func.evaluate_grad( 1, x_r );
        	return true;
        }
        
                double unit_fct( vertex v )
        {
            return 1.;
        }
        
        bool test_assemble_triangle()
        {
        	Mesh mesh;
            	mesh.load("data/square.mesh");
        	ElementMapping elt_mapping = ElementMapping(mesh, false, 4);
        	ShapeFunctions reference_functions = ShapeFunctions(2,1);
        	
        	Quadrature quad = Quadrature::get_quadrature(4, false);
        	DenseMatrix Ke;
        	assemble_elementary_matrix(elt_mapping, reference_functions, quad, unit_fct, Ke );
        	
        	int ind_max = mesh.nb_triangles();
        	SparseMatrix K = SparseMatrix(ind_max*3);
        	local_to_global_matrix(mesh,4,Ke,K );
        	
        	std::vector< bool > attribute_is_dirichlet (mesh.nb_vertices(),true);
        	std::vector< double > values (mesh.nb_vertices(),1);
        	std::vector< double > F (mesh.nb_vertices(),0);
        	apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F);
        	return true;
        }
        	
    }
}
