#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

#define PI 3.14159265

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
        
        double sinus_bump_fct( vertex v )
        {
            return 2*pow(PI,2)*sin(PI *v.x)*sin(PI *v.y);
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
            mesh.load(mesh_filename+".mesh"); //chargement maillage
            
            //int ind_max = mesh.nb_triangles();
            SparseMatrix K = SparseMatrix(mesh.nb_vertices()); //matrice globale
            std::vector< double > values (K.nb_rows(),0); //values g = x+y
            
            for (int element = 0; element< mesh.nb_triangles(); ++ element){
            	ElementMapping elt_mapping = ElementMapping(mesh, false, element); //element a calculer
            	ShapeFunctions reference_functions = ShapeFunctions(2,1); //calcul fonctions de base
            	Quadrature quad = Quadrature::get_quadrature(2, false); //quadrature de l'element
            	DenseMatrix Ke; 
            	assemble_elementary_matrix(elt_mapping, reference_functions, quad, unit_fct, Ke ); //matrice ke de l'element
            	local_to_global_matrix(mesh,element,Ke,K); //ajout matrice element sur matrice globale	
            } // fin travail par element
            	
            for (int point = 0; point<mesh.nb_vertices(); ++point){
            	values[point]= mesh.get_vertex(point).x + mesh.get_vertex(point).y;
            }
            mesh.set_attribute(unit_fct,1,true);
            std::vector< bool > attribute_is_dirichlet (2,false);
            attribute_is_dirichlet[1] = true;
            std::vector< double > F (K.nb_rows(),0); //dirichlet pur
            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F); //aplication cond de dirichlet sur K et F
            std::vector< double > x; //creation vecteur solution
            solve(K,F,x); // return true si simulation reussie 
            save_solution(x,mesh_filename+".bb");
            }

        
        void source_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh; 
            mesh.load(mesh_filename+".mesh"); //chargement maillage
            
             //int ind_max = mesh.nb_triangles();
            SparseMatrix K = SparseMatrix(mesh.nb_vertices()); //matrice globale
            std::vector< double > values (K.nb_rows(),0); //values g = 0
            
            std::vector< double > F (K.nb_rows(),0); //terme source
            
            for (int element = 0; element< mesh.nb_triangles(); ++ element){
            	ElementMapping elt_mapping = ElementMapping(mesh, false, element); //element a calculer
            	ShapeFunctions reference_functions = ShapeFunctions(2,1); //calcul fonctions de base
            	Quadrature quad = Quadrature::get_quadrature(2, false); //quadrature de l'element
            	DenseMatrix Ke; 
            	assemble_elementary_matrix(elt_mapping, reference_functions, quad, unit_fct, Ke ); //matrice ke de l'element
            	local_to_global_matrix(mesh,element,Ke,K); //ajout matrice element sur matrice globale	
            	
            	double initialisation = 0;
            	std::vector< double > Fe (3, initialisation) ;
            	assemble_elementary_vector(elt_mapping, reference_functions, quad, unit_fct, Fe ); //vecteur fe de l'element
            	
            	local_to_global_vector(mesh, false, element, Fe, F); //ajout vecteur element a vecteur global
            	
            } // fin travail par element

            mesh.set_attribute(unit_fct,1,true);
            std::vector< bool > attribute_is_dirichlet (2,false);
            attribute_is_dirichlet[1] = true;

            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F); //aplication cond de dirichlet sur K et F
            std::vector< double > x; //creation vecteur solution
            solve(K,F,x); // return true si simulation reussie 
            save_solution(x,mesh_filename+".bb");
            
            
            }
            
        void sinus_bump_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh; 
            mesh.load(mesh_filename+".mesh"); //chargement maillage
            
             //int ind_max = mesh.nb_triangles();
            SparseMatrix K = SparseMatrix(mesh.nb_vertices()); //matrice globale
            std::vector< double > values (K.nb_rows(),0); //values g = 0
            
            std::vector< double > F (K.nb_rows(),0); //terme source
            
            for (int element = 0; element< mesh.nb_triangles(); ++ element){
            	ElementMapping elt_mapping = ElementMapping(mesh, false, element); //element a calculer
            	ShapeFunctions reference_functions = ShapeFunctions(2,1); //calcul fonctions de base
            	Quadrature quad = Quadrature::get_quadrature(2, false); //quadrature de l'element
            	DenseMatrix Ke; 
            	assemble_elementary_matrix(elt_mapping, reference_functions, quad, sinus_bump_fct, Ke ); //matrice ke de l'element
            	local_to_global_matrix(mesh,element,Ke,K); //ajout matrice element sur matrice globale	
            	
            	double initialisation = 0;
            	std::vector< double > Fe (3, initialisation) ;
            	assemble_elementary_vector(elt_mapping, reference_functions, quad, sinus_bump_fct, Fe ); //vecteur fe de l'element
            	
            	local_to_global_vector(mesh, false, element, Fe, F); //ajout vecteur element a vecteur global
            	
            } // fin travail par element
            
            
            mesh.set_attribute(unit_fct,1,true);
            std::vector< bool > attribute_is_dirichlet (2,false);
            attribute_is_dirichlet[1] = true;

            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F); //aplication cond de dirichlet sur K et F
            std::vector< double > x; //creation vecteur solution
            solve(K,F,x); // return true si simulation reussie
            save_solution(x,mesh_filename+".bb");
            
            }
            
        void sinus_bump_error_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh; 
            mesh.load(mesh_filename+".mesh"); //chargement maillage
            
             //int ind_max = mesh.nb_triangles();
            SparseMatrix K = SparseMatrix(mesh.nb_vertices()); //matrice globale
            std::vector< double > values (K.nb_rows(),0); //values g = 0
            
            std::vector< double > F (K.nb_rows(),0); //terme source
            
            std::vector< double > theorique (K.nb_rows(),0); //solution analytique
            
            for (int element = 0; element< mesh.nb_triangles(); ++ element){
            	ElementMapping elt_mapping = ElementMapping(mesh, false, element); //element a calculer
            	ShapeFunctions reference_functions = ShapeFunctions(2,1); //calcul fonctions de base
            	Quadrature quad = Quadrature::get_quadrature(2, false); //quadrature de l'element
            	DenseMatrix Ke; 
            	assemble_elementary_matrix(elt_mapping, reference_functions, quad, sinus_bump_fct, Ke ); //matrice ke de l'element
            	local_to_global_matrix(mesh,element,Ke,K); //ajout matrice element sur matrice globale	
            	
            	double initialisation = 0;
            	std::vector< double > Fe (3, initialisation) ;
            	assemble_elementary_vector(elt_mapping, reference_functions, quad, sinus_bump_fct, Fe ); //vecteur fe de l'element
            	
            	local_to_global_vector(mesh, false, element, Fe, F); //ajout vecteur element a vecteur global
            	
            } // fin travail par element
            
            for (int point = 0; point<mesh.nb_vertices(); ++point){ // calcul solution analytique 
            	theorique[point]= sin(PI*mesh.get_vertex(point).x) + sin(PI*mesh.get_vertex(point).y);
            }
            
            mesh.set_attribute(unit_fct,1,true);
            std::vector< bool > attribute_is_dirichlet (2,false);
            attribute_is_dirichlet[1] = true;

            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F); //aplication cond de dirichlet sur K et F
            std::vector< double > x; //creation vecteur solution
            solve(K,F,x); // return true si simulation reussie 
            
            std::vector< double > erreur (K.nb_rows(),0);
            for (int point = 0; point < mesh.nb_vertices(); ++point){
            	erreur[point] = abs(x[point]-theorique[point]); 
            }
            save_solution(erreur,mesh_filename+".bb");
            
            }
            
            
        /*void neumann_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh; 
            mesh.load(mesh_filename+".mesh"); //chargement maillage
            
             //int ind_max = mesh.nb_triangles();
            SparseMatrix K = SparseMatrix(mesh.nb_vertices()); //matrice globale
            std::vector< double > values (K.nb_rows(),0); //values g = 0
            
            
            for (int element = 0; element< mesh.nb_triangles(); ++ element){
            	ElementMapping elt_mapping = ElementMapping(mesh, false, element); //element a calculer
            	ShapeFunctions reference_functions = ShapeFunctions(2,1); //calcul fonctions de base
            	Quadrature quad = Quadrature::get_quadrature(2, false); //quadrature de l'element
            	DenseMatrix Ke; 
            	assemble_elementary_matrix(elt_mapping, reference_functions, quad, unit_fct, Ke ); //matrice ke de l'element
            	local_to_global_matrix(mesh,element,Ke,K); //ajout matrice element sur matrice globale	
            	
            	double initialisation = 0;
            	std::vector< double > Fe (3, initialisation) ;
            	assemble_elementary_vector(elt_mapping, reference_functions, quad, unit_fct, Fe ); //vecteur fe de l'element
            	
            	local_to_global_vector(mesh, false, element, Fe, F); //ajout vecteur element a vecteur global
            	
            } // fin travail par element
            
            
            mesh.set_attribute(unit_fct,1,true);
            std::vector< bool > attribute_is_dirichlet (2,false);
            attribute_is_dirichlet[1] = true;

            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F); //aplication cond de dirichlet sur K et F
            std::vector< double > x; //creation vecteur solution
            solve(K,F,x); // return true si simulation reussie 
            
            std::vector< double > erreur (K.nb_rows(),0);
            for (int point = 0; point < mesh.nb_vertices(); ++point){
            	erreur[point] = abs(x[point]-theorique[point]); 
            }
            save_solution(erreur,mesh_filename+".bb");
            
            }*/
            
        }
        
}

