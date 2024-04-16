#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh; 
            mesh.load(mesh_filename); //chargement maillage
            
            int ind_max = mesh.nb_triangles();
            SparseMatrix K = SparseMatrix(ind_max*3); //matrice globale
            std::vector< double > values (K.nb_rows(),0); //values g = x+y
            
            for (int element = 0; element< mesh.nb_triangles(); ++ element){
            	ElementMapping elt_mapping = ElementMapping(mesh, false, element); //element a calculer
            	ShapeFunctions reference_functions = ShapeFunctions(2,1); //calcul fonctions de base
            	Quadrature quad = Quadrature::get_quadrature(2, false); //quadrature de l'element
            	DenseMatrix Ke; 
            	assemble_elementary_matrix(elt_mapping, reference_functions, quad, unit_fct, Ke ); //matrice ke de l'element
            	local_to_global_matrix(mesh,4,Ke,K ); //ajout matrice element sur matrice globale
            	for (int point = 0; point<3;++point){ //ajoute valeur g=x+y dans vecteur values (condition de dirichlet simple)
            		vertex local_vertex = mesh.get_triangle_vertex(element,point);
            		vertex global_vertex= elt_mapping.transform(local_vertex);
            		int global_index = mesh.get_triangle_vertex_index(element,point);
            		values[global_index]= global_vertex.x +global_vertex.y;
            	}
            	} // fin travail par element
            	
            	std::vector< bool > attribute_is_dirichlet (K.nb_rows(),true); 
            	std::vector< double > F (K.nb_rows(),0); //dirichlet pur
            	apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F); //aplication cond de dirichlet sur K et F
            	std::vector< double > x; //creation vecteur solution
            	solve(K,F,x); // return true si simulation reussie / sauver solution?????
            }
        }

}
