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
        
        // for neumann problem
        
        double edge_dirichlet_fct( vertex v )
        {
            if (v.x == 1){
            	return 1.;
            }
            else {
            	return -1.;
            }
        }
        
        double edge_neumann_fct( vertex v )
        {
            if (v.x == 0){
            	return 1.;
            }
            else {
            	return -1.;
            }
        }
        
        double edge_null_neumann_fct( vertex v )
        {
            if (v.y == 1 || v.y == 0){
            	return 1.;
            }
            else {
            	return -1.;
            }
        }
        
        double cond_neumann_fct( vertex v )
        {
            return sin(PI*v.y);
        }
        
        // for mug problem
        
        double mug_dirichlet_fct( vertex v )
        {
            if ((v.x == 1 && v.y < 10 && v.y >1) || (v.x == 20 && v.y < 10 && v.y > 1) || (v.y == 1 && v.x < 20 && v.x > 1)){
            	return 1.;
            }
            else {
            	return -1.;
            }
        }
        
        double mug_neumann_fct( vertex v )
        {
            if ((v.x == 1 && v.y < 10 && v.y >1) || (v.x == 20 && v.y < 10 && v.y > 1) || (v.y == 1 && v.x < 20 && v.x > 1)){
            	return -1.;
            }
            else {
            	return 1.;
            }
        }
        
        double mug_cond_neumann_fct( vertex v )
        {
            return -0.1;
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
            
            
        void neumann_pb( const std::string& mesh_filename, bool verbose )
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
            
            mesh.set_attribute(edge_dirichlet_fct,1,true); // cond de dirichlet
            std::vector< bool > attribute_is_dirichlet (4,false);
            attribute_is_dirichlet[1] = true;
            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F); //aplication cond de dirichlet sur K et F
            
            mesh.set_attribute(edge_neumann_fct,2,true); // cond de neumann
            mesh.set_attribute(edge_null_neumann_fct,3,true); // cond de neumann nulle
            
            ShapeFunctions reference_functions = ShapeFunctions(1,1);
            Quadrature quad = Quadrature::get_quadrature(2, true);
            attribute_is_dirichlet[1] = false;
            for (int element = 0; element < mesh.nb_edges(); ++ element){
            	double initialisation = 0;
            	ElementMapping elt_mapping(mesh, true, element);
            	attribute_is_dirichlet[2] = true;
            	attribute_is_dirichlet[3] = false;
            	if ( attribute_is_dirichlet[mesh.get_edge_attribute(element)] ){
            		std::vector< double > Fe (2, initialisation) ;
            		assemble_elementary_neumann_vector(elt_mapping, reference_functions, quad, cond_neumann_fct, Fe );
        		local_to_global_vector(mesh, true, element, Fe, F);
            	}
            	attribute_is_dirichlet[2] = false;
            	attribute_is_dirichlet[3] = true;
            	if ( attribute_is_dirichlet[mesh.get_edge_attribute(element)] ){
            		std::vector< double > Fe (2, initialisation) ;
        		local_to_global_vector(mesh, true, element, Fe, F);
            	}
            }

            std::vector< double > x; //creation vecteur solution
            solve(K,F,x); // return true si simulation reussie
            save_solution(x,mesh_filename+".bb");
            
            }
            
            void mug_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh; 
            mesh.load(mesh_filename+".mesh"); //chargement maillage
            
             //int ind_max = mesh.nb_triangles();
            SparseMatrix K = SparseMatrix(mesh.nb_vertices()); //matrice globale
            std::vector< double > values (K.nb_rows(),100); //values g = 100
            
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
            	assemble_elementary_vector(elt_mapping, reference_functions, quad, zero_fct, Fe ); //vecteur fe de l'element
            	
            	local_to_global_vector(mesh, false, element, Fe, F); //ajout vecteur element a vecteur global
            	
            } // fin travail par element
            
            mesh.set_attribute(mug_dirichlet_fct,1,true); // cond de dirichlet
            std::vector< bool > attribute_is_dirichlet (3,false);
            attribute_is_dirichlet[1] = true;
            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet, values, K, F); //aplication cond de dirichlet sur K et F
            
            mesh.set_attribute(mug_neumann_fct,2,true); // cond de neumann
            
            ShapeFunctions reference_functions = ShapeFunctions(1,1);
            Quadrature quad = Quadrature::get_quadrature(2, true);
            attribute_is_dirichlet[1] = false;
            for (int element = 0; element < mesh.nb_edges(); ++ element){
            	double initialisation = 0;
            	ElementMapping elt_mapping(mesh, true, element);
            	attribute_is_dirichlet[2] = true;
            	if ( attribute_is_dirichlet[mesh.get_edge_attribute(element)] ){
            		std::vector< double > Fe (2, initialisation) ;
            		assemble_elementary_neumann_vector(elt_mapping, reference_functions, quad, mug_cond_neumann_fct, Fe );
        		local_to_global_vector(mesh, true, element, Fe, F);
            	}
            }

            std::vector< double > x; //creation vecteur solution
            solve(K,F,x); // return true si simulation reussie
            save_solution(x,mesh_filename+".bb");
            
            }
            
        }
        
}

