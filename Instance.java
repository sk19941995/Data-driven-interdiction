import java.io.File;
import java.math.BigDecimal;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import ilog.concert.IloConstraint;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

import java.math.RoundingMode;

public class Instance {
		
		int number_of_items;
		int number_of_lc;
		int number_of_fc;
		
		ArrayList<Double> mean = new ArrayList<Double>();	//mean of profits
		ArrayList<Double> std = new ArrayList<Double>(); //std of profits
		ArrayList<Double> min_cost = new ArrayList<Double>(); //lower bounds	
		ArrayList<Double> max_cost = new ArrayList<Double>(); //upper bounds
		
		ArrayList<ArrayList<Double>> F = new ArrayList<ArrayList<Double>>(); //follower's constraint matrix
		ArrayList<Double> f = new ArrayList<Double>(); //follower's right hand side vector
		
		ArrayList<ArrayList<Double>> H = new ArrayList<ArrayList<Double>>(); //leader's constraint matrix
		ArrayList<Double> h = new ArrayList<Double>(); //leader's right hand side vector
		
		double weight_l, weight_f;
		
		Instance(int n, int p, int d_f, double w_l, double w_f, int inverse, Random random) throws FileNotFoundException
		{	
			number_of_items = n;
			number_of_lc = p;
			number_of_fc = d_f;
			weight_l = w_l;
			weight_f = w_f;
			
			double lower_bound = 0.01;
			double upper_bound = 1.;
			
			for (int i = 0; i < p; i++)
			{
				ArrayList<Double> str = new ArrayList<Double>();
				double sum = 0;
				
				for (int j = 0; j < n; j++)
				{
					double val = random.nextDouble(lower_bound, upper_bound);
					
					str.add(val);
					sum += val;
				}
				
				H.add(str);
				h.add(w_l*sum);
		    }
				
			for (int i = 0; i < d_f; i++)
			{
				ArrayList<Double> str = new ArrayList<Double>();
				double sum = 0;
				
				for (int j = 0; j < n; j++)
				{
					double val = random.nextDouble(lower_bound, upper_bound);
					
					str.add(val);
					sum += val;
				}
				
				F.add(str);
				f.add(w_f*sum);
		    }
			
		   for (int i = 0; i < n; i++)
		   {
			     
			     double mean_01 = 0.5*(i + 1)*1./(n + 1);
			     double std_01 = 0.05 + 0.4*(i + 1)*1./(n + 1);
			     
			     if (inverse == 1) 
			     {
			    	 std_01 = 0.05 + (0.4 - 0.4*(i + 1)*1./(n + 1));
			     }
				 min_cost.add(lower_bound);
				 max_cost.add(upper_bound);
				 
				 mean.add(lower_bound + (upper_bound - lower_bound)*mean_01);
				 std.add((upper_bound - lower_bound)*std_01);
		   }
			  
		}
	
		Instance(int n, int p, double w_l, double w_f, Random random) throws FileNotFoundException
		{	
			number_of_items = n;
			number_of_lc = p;
			weight_l = w_l;
			weight_f = w_f;
			
			double lower_bound = 0.01;
			double upper_bound = 1.;
			
			for (int i = 0; i < p; i++)
			{
				ArrayList<Double> str = new ArrayList<Double>();
				double sum = 0;
				
				for (int j = 0; j < n; j++)
				{
					double val = random.nextDouble(lower_bound, upper_bound);
					
					str.add(val);
					sum += val;
				}
				
				H.add(str);
				h.add(w_l*sum);
		    }
			
			number_of_fc = 1;
			for (int i = 0; i < 1; i++)
			{
				ArrayList<Double> str = new ArrayList<Double>();
				double sum = 0;
				
				for (int j = 0; j < n; j++)
				{
					double val = 1.;
					
					str.add(val);
					sum += val;
				}
				
				F.add(str);
				f.add(Math.floor(0.2*sum));
		    }
		
			/*
			int num_nodes = (int)(d_f)/2;
			int num_edges = n;
			number_of_fc = num_edges;
			
			for (int i = 0; i < 2*num_nodes; i++)
			{
				ArrayList<Double> str = new ArrayList<Double>();
				for (int j = 0; j < num_edges; j++)
				{	
					str.add(0.);
				}
				
				F.add(str);
				f.add(1.);
		    }
			
			
			Map<ArrayList<Integer>, Integer> used = new HashMap<ArrayList<Integer>, Integer>(); 
			
			for (int j = 0; j < num_edges; j++)
			{
				int left = random.nextInt(0, num_nodes);
				int right = random.nextInt(num_nodes, 2*num_nodes);
				
				ArrayList<Integer> arc = new ArrayList<Integer>();
				arc.add(left);
				arc.add(right);
				
				if (used.containsKey(arc) == false)
				{
					F.get(left).set(j, 1.);
					F.get(right).set(j, 1.);
				}
				
				used.put(arc, 1);
			}*/
		   
			
		   for (int i = 0; i < n; i++)
		   {
			     double mean_01 = 0.5*(i + 1)*1./(n + 1);
			     double std_01 = 0.05 + 0.4*(i + 1)*1./(n + 1);
			     
				 min_cost.add(lower_bound);
				 max_cost.add(upper_bound);
				 
				 mean.add(lower_bound + (upper_bound - lower_bound)*mean_01);
				 std.add((upper_bound - lower_bound)*std_01);
		   }
			  
		}	
	public static ArrayList<ArrayList<Double>> SolveLeadersProblemFull(IloCplex cplex, Instance I, double epsilon_f, double alpha_l, double alpha_f, int k_l, int k_f, ArrayList<ArrayList<Double>> data_set_l, ArrayList<ArrayList<Double>> data_set_f,
					ArrayList<ArrayList<Double>> support_constraints, String type_of_risk_l, String type_of_risk_f,  ArrayList<ArrayList<Double>> delta_f, ArrayList<Double> average_l, ArrayList<Double> average_f, double M1)
		    {
			try 
		    {
		    	cplex.clearModel();
		    	cplex.setOut(null);
		    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
		    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
		    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
		            
		        double M = Double.MAX_VALUE;
		    	int number_of_sc = support_constraints.size() - 1; 
		    	
		    	ArrayList<ArrayList<IloNumVar>> nu_f = new ArrayList<ArrayList<IloNumVar>>();
		    	
		    	ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
		    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
		    	  
		    	ArrayList<IloNumVar> s_l = new ArrayList<IloNumVar>();
		    	ArrayList<IloNumVar> s_f = new ArrayList<IloNumVar>();
		    	
		    	//IloNumVar lambda_l = cplex.numVar(0, M, "lam_l");
		    	IloNumVar lambda_f = cplex.numVar(0, M, "lam_f");
		    	
		    	IloNumVar t_l = cplex.numVar(-M, M, "t_l");
		    	IloNumVar t_f = cplex.numVar(-M, M, "t_f");
		    	
		    	ArrayList<ArrayList<IloNumVar>> mu1_f = new ArrayList<ArrayList<IloNumVar>>();
		    	ArrayList<ArrayList<IloNumVar>> mu2_f = new ArrayList<ArrayList<IloNumVar>>();
		    	
		    	ArrayList<IloNumVar> beta1_f = new ArrayList<IloNumVar>();
		    	ArrayList<IloNumVar> beta2_f = new ArrayList<IloNumVar>();
		    	ArrayList<IloNumVar> gamma_f = new ArrayList<IloNumVar>();
		    	ArrayList<IloNumVar> prod_f = new ArrayList<IloNumVar>();
		    	
		    	
		    	for (int k = 0; k < k_l; k++)
		    	{
		    		s_l.add(cplex.numVar(0, M, "s_l_"+ k));
		    	}
		    	
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		s_f.add(cplex.numVar(0, M, "s_f_"+ k));
		    		nu_f.add(new ArrayList<IloNumVar>());
		    	
		    		for (int l = 0; l < number_of_sc; l++)
		    		{
		    			nu_f.get(k).add(cplex.numVar(0, M, "nu_f_"+ k + l));	
		    		}
		    	}
		    
		        for (int j = 0; j < I.number_of_items; j++)
		        {
		        	x.add(cplex.boolVar("x_" + j));
		            y.add(cplex.numVar(0, 1, "y_" + j));
		        	
		        	prod_f.add(cplex.numVar(0, M, "prod_" + j));
		        }
		    	
		        //Follower's dual variables
		        for (int j = 0; j < I.number_of_fc; j++)
		        {
		        	beta1_f.add(cplex.numVar(0, M, "beta1_f_"+ j));
		        }
		        
		        for (int i = 0; i < I.number_of_items; i++)
		        {
		        	beta2_f.add(cplex.numVar(0, M, "beta2_f_"+ i));
		        }
		        
		        for (int k = 0; k < k_f; k++)
		    	{
		    		gamma_f.add(cplex.numVar(0, M, "gamma_f_"+ k));
		    		
		    		mu1_f.add(new ArrayList<IloNumVar>());
		    		mu2_f.add(new ArrayList<IloNumVar>());
		    	
		    		 for (int j = 0; j < I.number_of_items; j++)
		    		{
		    			mu1_f.get(k).add(cplex.numVar(0, M, "mu1_f_"+ k + j));	
		    			mu2_f.get(k).add(cplex.numVar(0, M, "mu2_f_"+ k + j));
		    		}
		    	}
		        
		    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
		    	
		    	if ("Risk-averse".equals(type_of_risk_l))
		    	{   
		    		objective.addTerm(1., t_l);
		   
		    		for (int k = 0; k < k_l; k++)
		        	{
		        	     objective.addTerm(1./(k_l*(1. - alpha_l)), s_l.get(k));
		        	}
		    	}
		    	
		    	if ("Risk-neutral".equals(type_of_risk_l))
		    	{
		    		for (int i = 0; i < I.number_of_items; i++)
		    			objective.addTerm(average_l.get(i), y.get(i));
		    	  		
		    	}
		    	
		    	cplex.addMinimize(objective);
		    	
		    	//Constraints
		    	//Feasibility
		    	//X
		    	for (int j = 0; j < I.number_of_lc; j++)
		    	{
		    		IloLinearNumExpr str_1 = cplex.linearNumExpr(); 
		    		
		    		for (int i = 0; i < I.number_of_items; i++)
		        	{
		    			str_1.addTerm(I.H.get(j).get(i), x.get(i));
		        	}
		    		 
		    		cplex.addLe(str_1, I.h.get(j));
		         	
		        }
		    	
		    	 //Y
		    	for (int j = 0; j < I.number_of_fc; j++)
		    	{
		    		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
		    		
		    		for (int i = 0; i < I.number_of_items; i++)
		        	{
		    			str_2.addTerm(I.F.get(j).get(i), y.get(i));
		        	}
		    		 
		    		cplex.addLe(str_2, I.f.get(j));
		         	
		        }
		    	
		    	for (int j = 0; j < I.number_of_items; j++)
		    	{
		    		IloLinearNumExpr str = cplex.linearNumExpr(); 
		    		str.addTerm(1., y.get(j));
		    		str.addTerm(1., x.get(j));
		        	cplex.addLe(str, 1.);
		    	}
		        
		    	//
		     	/*for (int k = 0; k < k_l; k++)
		    		for (int j = 0; j < I.number_of_items; j++)
		    		{
		    			IloLinearNumExpr str1 = cplex.linearNumExpr();
		    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	    	
		    	    	str1.addTerm(-1., lambda_l);
		    	    	str2.addTerm(-1., lambda_l);
		    	    	
		    	    	str1.addTerm(1., y.get(j));
		    	    	str2.addTerm(-1., y.get(j));
		    	    	
		    	    	for (int l = 0; l < number_of_sc; l++)
		     	        {
		    	    		str1.addTerm(support_constraints.get(l).get(j), nu_l.get(k).get(l));
		    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_l.get(k).get(l));
		     	        }
		    	    	
		    	    	cplex.addLe(str1, 0.);
		    	    	cplex.addLe(str2, 0.);
		    		}*/
		    	
		     	//
		    	for (int k = 0; k < k_f; k++)
		    		for (int j = 0; j < I.number_of_items; j++)
		    		{
		    			IloLinearNumExpr str1 = cplex.linearNumExpr();
		    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	    	
		    	    	str1.addTerm(-1., lambda_f);
		    	    	str2.addTerm(-1., lambda_f);
		    	    	
		    	    	str1.addTerm(1., y.get(j));
		    	    	str2.addTerm(-1., y.get(j));
		    	    	
		    	    	for (int l = 0; l < number_of_sc; l++)
		     	        {
		    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(k).get(l));
		    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(k).get(l));
		     	        }
		    	    	
		    	    	cplex.addLe(str1, 0.);
		    	    	cplex.addLe(str2, 0.);
		    		}
		   	    
		    	//
		    	if ("Risk-averse".equals(type_of_risk_l))
		    	{   
		    		for (int k = 0; k < k_l; k++)
		    		{
		    			IloLinearNumExpr str1 = cplex.linearNumExpr();
		    			
		    			 for (int i = 0; i < I.number_of_items; i++)
		    			 {
		    				 str1.addTerm(data_set_l.get(i).get(k), y.get(i));
		    		     }
		    			
		    			str1.addTerm(-1., t_l);
		    			
		    			str1.addTerm(-1., s_l.get(k));
		    			cplex.addLe(str1, 0.);
		    			
		                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
		    			
		    			str2.addTerm(1., s_l.get(k));
		    			cplex.addGe(str2, 0);*/
		    		}
		    		
		    	}
		    	
		    	//
		    	if ("Risk-averse".equals(type_of_risk_f))
		    	{   
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			IloLinearNumExpr str1 = cplex.linearNumExpr();
		    			
		    			 for (int i = 0; i < I.number_of_items; i++)
		    			 {
		    				 str1.addTerm(-data_set_f.get(i).get(k), y.get(i));
		    		     }
		    			
		    			str1.addTerm(1., t_f);
		    			
		    
		        		for (int l = 0; l < number_of_sc; l++)
		        		{
		        			str1.addTerm(delta_f.get(k).get(l), nu_f.get(k).get(l));
		        		}
		    			
		    			str1.addTerm(-1., s_f.get(k));
		    			cplex.addLe(str1, 0.);
		    			
		                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
		    			
		    			str2.addTerm(1., s_f.get(k));
		    			cplex.addGe(str2, 0);*/
		    		}
		    		
		    	}
		    	
		    	//Dual constraints of the follower
		    	if ("Risk-averse".equals(type_of_risk_f))
		    	{
		    		//
			    	for (int i = 0; i < I.number_of_items; i++)
			    	 {	
			    		IloLinearNumExpr str1 = cplex.linearNumExpr();
						
			    		for (int j = 0; j < I.number_of_fc; j++)
			           	{
			       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
			           	}
			       		
			       		str1.addTerm(1., beta2_f.get(i));
			    		
			    		for (int k = 0; k < k_f; k++)
			    		{
			    			str1.addTerm(-data_set_f.get(i).get(k), gamma_f.get(k));
			    			str1.addTerm(-1., mu1_f.get(k).get(i));
			    			str1.addTerm(1., mu2_f.get(k).get(i));
			    		}
			        	
			    		cplex.addGe(str1, 0.);
			    	}
			    	
			    	//
			    	for (int k = 0; k < k_f; k++)
			    	{
			    		for (int l = 0; l < number_of_sc; l++)
			            {
			    			IloLinearNumExpr str1 = cplex.linearNumExpr();
				    		
			    			str1.addTerm(delta_f.get(k).get(l), gamma_f.get(k));
				    		
				    		for (int i = 0; i < I.number_of_items; i++)
					    	 {
				    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
				    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
					    	 }
				    		cplex.addGe(str1, 0.);
				    	} 
			    	}
			    	
			    	//
			    	IloLinearNumExpr str = cplex.linearNumExpr();
			    	for (int k = 0; k < k_f; k++)
			    	{
			    		IloLinearNumExpr str1 = cplex.linearNumExpr();
			    		
			    		str.addTerm(1., gamma_f.get(k));
			    		str1.addTerm(1., gamma_f.get(k));
			    		
			    		cplex.addLe(str1, 1./((1. - alpha_f) * k_f));
			    	}
			    	
			    	cplex.addEq(str, 1.);
		    	}
		    	
		    	if ("Risk-neutral".equals(type_of_risk_f))
		    	{   
		    		//
			    	for (int i = 0; i < I.number_of_items; i++)
			    	 {	
			    		IloLinearNumExpr str1 = cplex.linearNumExpr(-average_f.get(i));
						
			    		for (int j = 0; j < I.number_of_fc; j++)
			           	{
			       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
			           	}
			       		
			       		str1.addTerm(1., beta2_f.get(i));
			    		
			    		
			    		for (int k = 0; k < k_f; k++)
			    		{
			    			str1.addTerm(-1., mu1_f.get(k).get(i));
			    			str1.addTerm(1., mu2_f.get(k).get(i));
			    		}
			        	
			    		cplex.addGe(str1, 0.);
			    	}
			    	
			    	//
			    	for (int k = 0; k < k_f; k++)
			    	{
			    		for (int l = 0; l < number_of_sc; l++)
			            {
			    			IloLinearNumExpr str1 = cplex.linearNumExpr(delta_f.get(k).get(l)/(k_f + 0.));
				    		
				    		for (int i = 0; i < I.number_of_items; i++)
					    	 {
				    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
				    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
					    	 }
				    		cplex.addGe(str1, 0.);
				    	}
			    	}
		    	}
		    	
		    	//
		    	IloLinearNumExpr str3 = cplex.linearNumExpr();
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int i = 0; i < I.number_of_items; i++)
			    	 {
		    			str3.addTerm(1., mu1_f.get(k).get(i));
		    			str3.addTerm(1., mu2_f.get(k).get(i));
			    	 }
		    	}
		    	
		    	if ("Risk-averse".equals(type_of_risk_f))
		    		cplex.addLe(str3, epsilon_f/(1. - alpha_f));
		    	
		    	if ("Risk-neutral".equals(type_of_risk_f))
		    		cplex.addLe(str3, epsilon_f);
		    	
		    	
		    	//Strong duality
		    	if ("Risk-averse".equals(type_of_risk_f))
		    	{
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		IloLinearNumExpr str2 = cplex.linearNumExpr();
			    	
		    		
			    	for (int i = 0; i < I.number_of_items; i++)
			    	{
			    		str1.addTerm(1., prod_f.get(i));
			    		//str1.addTerm(-1., beta_f.get(i));
			    	}
			    	
			    	for (int j = 0; j < I.number_of_fc; j++)
			    	{
			    		str1.addTerm(I.f.get(j), beta1_f.get(j));
			    	}
			    	
			    	
			    	str2.addTerm(1., t_f);
			    	str2.addTerm(-epsilon_f/(1. - alpha_f), lambda_f);
			    	
			    	
			    	for (int k = 0; k < k_f; k++)
			    	{
			    		str2.addTerm(-1./(k_f * (1. - alpha_f)), s_f.get(k));
			    	}
			    	
			    	cplex.addEq(str1, str2);
		    	}
		    	
		    	if ("Risk-neutral".equals(type_of_risk_f))
		    	{
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		IloLinearNumExpr str2 = cplex.linearNumExpr();
			    	
			    	for (int i = 0; i < I.number_of_items; i++)
			    	{
			    		str1.addTerm(1., prod_f.get(i));
			    		//str1.addTerm(-1., beta_f.get(i));
			    	}
			        
			    	for (int j = 0; j < I.number_of_fc; j++)
			    	{
			    		str1.addTerm(I.f.get(j), beta1_f.get(j));
			    	}
		       		
			    	str2.addTerm(-epsilon_f, lambda_f);
			    	
			    	for (int i = 0; i < I.number_of_items; i++)
			    	{
			    		str2.addTerm(average_f.get(i), y.get(i));
			    	}
			    	
			    	for (int k = 0; k < k_f; k++)
		    			for (int l = 0; l < number_of_sc; l++)
		    		{
		    			str2.addTerm(-delta_f.get(k).get(l)/(k_f + 0.), nu_f.get(k).get(l));
		    		}
			    	
			    	cplex.addEq(str1, str2);
		    	}
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		IloLinearNumExpr str2 = cplex.linearNumExpr(-M1);
		    		IloLinearNumExpr str4 = cplex.linearNumExpr();
		    		
		    		str1.addTerm(1., prod_f.get(i));
		    		str1.addTerm(-1., beta2_f.get(i));
		        	
		    		str2.addTerm(1., prod_f.get(i));
		    		str2.addTerm(M1, x.get(i));
		        	
		    		str4.addTerm(1., prod_f.get(i));
		    		str4.addTerm(-1, beta2_f.get(i));
		    		str4.addTerm(M1, x.get(i));
		        	
		    		cplex.addLe(str1, 0.);
		    		cplex.addLe(str2, 0.);
		    		cplex.addGe(str4, 0.);
		    		
		    	}
		    	
		    	
		    	ArrayList<ArrayList<Double>> sol_opt = new ArrayList<ArrayList<Double>>();
		    	ArrayList<Double> y_opt = new ArrayList<Double>();
		    	ArrayList<Double> x_opt = new ArrayList<Double>();
		    
		    	if (cplex.solve())
		    	{   	  
		    	    for (int j = 0; j < I.number_of_items; j++)
		    	    {
		    	    	x_opt.add(cplex.getValue(x.get(j)));
		    			y_opt.add(cplex.getValue(y.get(j)));
		    	    }  
		    	    
		    	    x_opt.add(cplex.getObjValue());
		    	    
		    	    sol_opt.add(x_opt);
		    	    sol_opt.add(y_opt);
		    	    
		    	    /*double delta = 0;
		    	    for (int j = 0; j < I.number_of_items; j++)
		    	    {
		    	    	System.out.println(cplex.getValue(prod_f.get(j)));
		    	    	
		    	    	delta += cplex.getValue(prod_f.get(j)) - (1 - cplex.getValue(x.get(j))) * cplex.getValue(beta_f.get(j));
		    	    }
		    	    
		    	    System.out.println(delta);*/
		    	}
		    	else 
		        	System.out.println("No solution found");
		    	
		    	return sol_opt;
			}
		    catch (IloException Exc) { Exc.printStackTrace(); }
		    return new ArrayList<ArrayList<Double>>();
		}
	
	public static ArrayList<ArrayList<Double>> GenerateData(Instance I, int sample_size, Random rand, int fixed) 
    {
		ArrayList<ArrayList<Double>> data = new ArrayList<ArrayList<Double>>();
		
		Random rng = (fixed == 0) ? new Random() : rand;
    	
    	for (int j = 0; j < I.number_of_items; j++)
    	{  
	    	ArrayList<Double> samples = new ArrayList<Double>();
	    	
	    	double mean_01 = (I.mean.get(j) - I.min_cost.get(j))/(I.max_cost.get(j) - I.min_cost.get(j));
    		double std_01 = I.std.get(j)/(I.max_cost.get(j) - I.min_cost.get(j));
	    	
	    	for (int k = 0; k < sample_size; k++)
	    	{	    	
	    		double sample = -1;
	    		
	    		while ((sample <= 0) || (sample >= 1))
	    		{
	    			sample = mean_01 + std_01*rng.nextGaussian();
	    		}
	    		
	    		samples.add((I.max_cost.get(j) - I.min_cost.get(j))*sample + I.min_cost.get(j));
    		}
	    
	    	data.add(samples);
   
    	}
    	return data;
    }
	
	
	public static ArrayList<ArrayList<Double>> TakeSubsetColumns(ArrayList<ArrayList<Double>> data, int number_of_rows, int number_of_columns_reduced) 
    {
		ArrayList<ArrayList<Double>> data_2 = new ArrayList<ArrayList<Double>>();
    	
		for (int i = 0; i < number_of_rows; i++)
		{
			ArrayList<Double> samples = new ArrayList<Double>();
			
			for (int k = 0; k < number_of_columns_reduced; k++)
	    	{
				samples.add(data.get(i).get(k));
	    	}
			
			data_2.add(samples);
		}
    	return data_2;
    }
	
	public static ArrayList<ArrayList<Double>> TakeSubsetRows(ArrayList<ArrayList<Double>> data, int number_of_rows_reduced, int number_of_columns) 
    {
		ArrayList<ArrayList<Double>> data_2 = new ArrayList<ArrayList<Double>>();
    	
		for (int k = 0; k < number_of_rows_reduced; k++)
		{
			ArrayList<Double> samples = CopyDoubleVector1(data.get(k), number_of_columns);
			data_2.add(samples);
		}
    	return data_2;
    }
	
	public static ArrayList<Double> InitializeDoubleVector(int size, double element)
    {
    	ArrayList<Double> Vector = new ArrayList<Double>();
    	for (int i = 0; i < size; i++)
  		  Vector.add(element);
    	return Vector;
    }
	
	public static ArrayList<Double> CopyDoubleVector1(ArrayList<Double> vector, int size)
    {
    	ArrayList<Double> new_vector = new ArrayList<Double>();
    	for (int i = 0; i < size; i++)
  		  new_vector.add(vector.get(i));
    	return new_vector;
    }
	
	public static ArrayList<ArrayList<Double>> ConstructSupportConstraints(Instance I) 
    {
		ArrayList<ArrayList<Double>> support_constraints_left = new ArrayList<ArrayList<Double>>();
		ArrayList<Double> support_constraints_right = new ArrayList<Double>();
		
		for (int j = 0; j < I.number_of_items; j++)
		{

			support_constraints_left.add(InitializeDoubleVector(I.number_of_items, 0.));
			
			support_constraints_left.get(support_constraints_left.size() - 1).set(j, 1.);
			support_constraints_right.add(I.max_cost.get(j));
			
			support_constraints_left.add(InitializeDoubleVector(I.number_of_items, 0.));
			
			support_constraints_left.get(support_constraints_left.size() - 1).set(j, -1.);
			support_constraints_right.add(-I.min_cost.get(j));
		}
		
		support_constraints_left.add(support_constraints_right);

		return support_constraints_left;
    }
	
	public static double NominalCost(Instance I, ArrayList<Double> y)
	{	
		double cost = 0;
		for (int j = 0; j < I.number_of_items; j++) 
		{
 			cost += y.get(j)*I.mean.get(j);
		}
		
		return cost;
	}
	
	public static double AverageOfArray(ArrayList<Double> cost)
	{
		double average = 0;
		
		for (int i = 0; i < cost.size(); i++)
		{
			average += cost.get(i);
		}
		
		average = average/(cost.size() + 0.);
		return average;
	}
	
	public static double VarianceOfArray(ArrayList<Double> cost)
	{
		double average = AverageOfArray(cost);
		double mad = 0;
		
		for (int i = 0; i < cost.size(); i++)
		{
			mad += Math.abs(cost.get(i) - average);
		}
		
		mad = mad/(cost.size() + 0.);
		return mad;
	}
	
	public static double RightTailedCvaR(IloCplex cplex,  Instance I, ArrayList<ArrayList<Double>> data_set_all, int sample_size_all, double alpha, ArrayList<Double> y)
	{	
		double cost = -1;
		try 
	    {
		cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
        //cplex.setParam(IloCplex.Param.Preprocessing.RepeatPresolve, 0);
       
        double M = Double.MAX_VALUE;
        
        IloNumVar t = cplex.numVar(-M, M, "t");
        ArrayList<IloNumVar> z = new ArrayList<IloNumVar>();
        
        for (int k = 0; k < sample_size_all; k++)
    	{ 
        	z.add(cplex.numVar(-M, M, "z" + k));
    	}
        
        IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	objective.addTerm(1., t);
    		
    	for (int k = 0; k < sample_size_all; k++)
        {
        	objective.addTerm(1./(sample_size_all*(1 - alpha)), z.get(k));
        }
    	
        cplex.addMinimize(objective);
    	
    	//Constraints
        for (int k = 0; k < sample_size_all; k++)
    	{   
        	double sum = 0;
        	for (int i = 0; i < I.number_of_items; i++)
         	{
        		 sum += data_set_all.get(i).get(k)*y.get(i);
         	}
        	
        	IloLinearNumExpr str1 = cplex.linearNumExpr(); 
        	IloLinearNumExpr str2 = cplex.linearNumExpr(-sum); 
        
        	str1.addTerm(1., z.get(k));
        	str2.addTerm(1., z.get(k));
        	str2.addTerm(1., t);
        
        	cplex.addGe(str1, 0);
        	cplex.addGe(str2, 0);
    	}
        
        if (cplex.solve())
    	{   		
    	    cost = cplex.getObjValue();
    	}
    	else 
        	System.out.println("No solution found");
        
	    }
        catch (IloException Exc) { Exc.printStackTrace(); }
		return cost;
	}
	
		

	public static ArrayList<Double> SolveFollowersProblem(IloCplex cplex, Instance I, ArrayList<Double> x, double epsilon_f, double alpha_f, int sample_size, ArrayList<ArrayList<Double>> data_set,  ArrayList<ArrayList<Double>> support_constraints, String type_of_risk, ArrayList<ArrayList<Double>> delta_f, ArrayList<Double> average_f,  double binary_indicator)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
        
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	ArrayList<ArrayList<IloNumVar>> nu_f = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> s_f = new ArrayList<IloNumVar>();
    	IloNumVar lambda_f = cplex.numVar(0, M, "lam");
    	IloNumVar t_f = cplex.numVar(-M, M, "t");
    	
    	for (int k = 0; k < sample_size; k++)
    	{
    		s_f.add(cplex.numVar(0, M, "s_"+ k));
    		nu_f.add(new ArrayList<IloNumVar>());
    	
    		for (int l = 0; l < number_of_sc; l++)
    		{
    			nu_f.get(k).add(cplex.numVar(0, M, "nu_"+ k + l));	
    		}
    	}
    	
        for (int j = 0; j < I.number_of_items; j++)
        {
        	if (binary_indicator == 0) 
        		y.add(cplex.numVar(0, 1, "y_" + j));
        	else 
        		y.add(cplex.boolVar("y_" + j)); 
        }
        
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	if ("Risk-averse".equals(type_of_risk))
    	{   
    		objective.addTerm(1., t_f);
    		objective.addTerm(-epsilon_f/(1 - alpha_f), lambda_f);
    		
    		for (int k = 0; k < sample_size; k++)
        	{
        	     objective.addTerm(-1./(sample_size*(1 - alpha_f)), s_f.get(k));
        	}
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk))
    	{
    		objective.addTerm(-epsilon_f, lambda_f);
    		
    		for (int i = 0; i < I.number_of_items; i++)
    			objective.addTerm(average_f.get(i), y.get(i));
    		
    		for (int k = 0; k < sample_size; k++)
    			for (int l = 0; l < number_of_sc; l++)
    		{
    			objective.addTerm(-delta_f.get(k).get(l)/(sample_size + 0.), nu_f.get(k).get(l));
    		}
    		
    	}
    	
    	cplex.addMaximize(objective);
    	
    	//Constraints Fy + Lx <= f (b)
    	IloLinearNumExpr str = cplex.linearNumExpr(); 
    	
    	 //Y
    	for (int j = 0; j < I.number_of_fc; j++)
    	{
    		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_2.addTerm(I.F.get(j).get(i), y.get(i));
        	}
    		 
    		cplex.addLe(str_2, I.f.get(j));
         	
        }
    	
    	for (int j = 0; j < I.number_of_items; j++)
    	{
    		IloLinearNumExpr str1 = cplex.linearNumExpr(); 
    		str1.addTerm(1., y.get(j));
        	cplex.addLe(str1, 1 - x.get(j));
    	}
    
    	for (int k = 0; k < sample_size; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
    	    	
    	    	str1.addTerm(-1, lambda_f);
    	    	str2.addTerm(-1, lambda_f);
    	    	
    	    	str1.addTerm(1, y.get(j));
    	    	str2.addTerm(-1, y.get(j));
    	    	
    	    	for (int l = 0; l < number_of_sc; l++)
     	        {
    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(k).get(l));
    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(k).get(l));
     	        }
    	    	
    	    	cplex.addLe(str1, 0.);
    	    	cplex.addLe(str2, 0.);
    		}
   	    
    	if ("Risk-averse".equals(type_of_risk))
    	{   
    		for (int k = 0; k < sample_size; k++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			
    			 for (int i = 0; i < I.number_of_items; i++)
    			 {
    				 str1.addTerm(-data_set.get(i).get(k), y.get(i));
    		     }
    			
    			str1.addTerm(1., t_f);
    			
    
        		for (int l = 0; l < number_of_sc; l++)
        		{
        			str1.addTerm(delta_f.get(k).get(l), nu_f.get(k).get(l));
        		}
    			
    			str1.addTerm(-1., s_f.get(k));
    			cplex.addLe(str1, 0.);
    			
    			/*IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str2.addTerm(1., s_f.get(k));
    			cplex.addGe(str2, 0);*/
    		}
    		
    	}
    	
    	ArrayList<Double> y_opt = new ArrayList<Double>();
    	
    	if (cplex.solve())
    	{   		
    		/*for (int k = 0; k < sample_size; k++)
        	{
        		for (int l = 0; l < number_of_sc; l++)
        		{
        			System.out.println(cplex.getValue(nu_f.get(k).get(l)));	
        		}
        	}*/
    		
    		
    	    for (int j = 0; j < I.number_of_items; j++)
    	    {
    			y_opt.add(cplex.getValue(y.get(j)));
    	    }  
    	    
    	    y_opt.add(cplex.getObjValue());
    	}
    	else 
        	System.out.println("No solution found");
    	
    	return y_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}
	
	public static Double SolveFollowersProblemFixedY(IloCplex cplex, Instance I, ArrayList<Double> x, ArrayList<Double> y, double epsilon_f, double alpha_f, int sample_size, ArrayList<ArrayList<Double>> data_set,  ArrayList<ArrayList<Double>> support_constraints, String type_of_risk, ArrayList<ArrayList<Double>> delta_f, ArrayList<Double> average_f,  double binary_indicator)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
        
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	ArrayList<ArrayList<IloNumVar>> nu_f = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<IloNumVar> s_f = new ArrayList<IloNumVar>();
    	IloNumVar lambda_f = cplex.numVar(0, M, "lam");
    	IloNumVar t_f = cplex.numVar(-M, M, "t");
    	
    	for (int k = 0; k < sample_size; k++)
    	{
    		s_f.add(cplex.numVar(0, M, "s_"+ k));
    		nu_f.add(new ArrayList<IloNumVar>());
    	
    		for (int l = 0; l < number_of_sc; l++)
    		{
    			nu_f.get(k).add(cplex.numVar(0, M, "nu_"+ k + l));	
    		}
    	}
    	
        
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	if ("Risk-averse".equals(type_of_risk))
    	{   
    		objective.addTerm(1., t_f);
    		objective.addTerm(-epsilon_f/(1 - alpha_f), lambda_f);
    		
    		for (int k = 0; k < sample_size; k++)
        	{
        	     objective.addTerm(-1./(sample_size*(1 - alpha_f)), s_f.get(k));
        	}
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk))
    	{
    		
    		double sum = 0;
    		for (int i = 0; i < I.number_of_items; i++)
    			sum += average_f.get(i) * y.get(i);
    		
    		objective = cplex.linearNumExpr(sum);
    		
    		objective.addTerm(-epsilon_f, lambda_f);
    		
    		for (int k = 0; k < sample_size; k++)
    			for (int l = 0; l < number_of_sc; l++)
    		{
    			objective.addTerm(-delta_f.get(k).get(l)/(sample_size + 0.), nu_f.get(k).get(l));
    		}
    		
    	}
    	
    	cplex.addMaximize(objective);
    	
  
    	for (int k = 0; k < sample_size; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr(y.get(j));
    	    	IloLinearNumExpr str2 = cplex.linearNumExpr(-y.get(j));
    	    	
    	    	str1.addTerm(-1, lambda_f);
    	    	str2.addTerm(-1, lambda_f);
    	    	
    	    	for (int l = 0; l < number_of_sc; l++)
     	        {
    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(k).get(l));
    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(k).get(l));
     	        }
    	    	
    	    	cplex.addLe(str1, 0.);
    	    	cplex.addLe(str2, 0.);
    		}
   	    
    	if ("Risk-averse".equals(type_of_risk))
    	{   
    		for (int k = 0; k < sample_size; k++)
    		{
    			
    			double sum = 0;
    			 for (int i = 0; i < I.number_of_items; i++)
    			 {
    				 sum += -data_set.get(i).get(k) * y.get(i);
    		     }
    			
    			IloLinearNumExpr str1 = cplex.linearNumExpr(sum);
     			 
    			str1.addTerm(1., t_f);
    			
    
        		for (int l = 0; l < number_of_sc; l++)
        		{
        			str1.addTerm(delta_f.get(k).get(l), nu_f.get(k).get(l));
        		}
    			
    			str1.addTerm(-1., s_f.get(k));
    			cplex.addLe(str1, 0.);
    			
    			/*IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str2.addTerm(1., s_f.get(k));
    			cplex.addGe(str2, 0);*/
    		}
    		
    	}
    	
    	double z_opt = -1;
    	
    	if (cplex.solve())
    	{   		
    		
    	    z_opt = cplex.getObjValue();
    	}
    	else 
        	System.out.println("No solution found");
    	
    	return z_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return -1.;
}

	
	public static ArrayList<Double> SolveFollowersProblemFullInformation(IloCplex cplex, Instance I, ArrayList<Double> x, double alpha_f, int sample_size, ArrayList<ArrayList<Double>> data_set,  ArrayList<ArrayList<Double>> support_constraints, String type_of_risk, ArrayList<Double> average_f,  double binary_indicator) {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
       
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> s_f = new ArrayList<IloNumVar>();
    	IloNumVar t_f = cplex.numVar(-M, M, "t");
    	
    	for (int k = 0; k < sample_size; k++)
    	{
    		s_f.add(cplex.numVar(0, M, "s_"+ k));
    		
    	}
    	
        for (int j = 0; j < I.number_of_items; j++)
        {
        	if (binary_indicator == 0) 
        		y.add(cplex.numVar(0, 1, "y_" + j));
        	else 
        		y.add(cplex.boolVar("y_" + j)); 
        }
       
      
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	if ("Risk-averse".equals(type_of_risk))
    	{   
    		objective.addTerm(1., t_f);
    		
    		for (int k = 0; k < sample_size; k++)
        	{
        	     objective.addTerm(-1./(sample_size*(1 - alpha_f)), s_f.get(k));
        	}
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk))
    	{
    		
    		for (int i = 0; i < I.number_of_items; i++)
    			objective.addTerm(average_f.get(i), y.get(i));
    			
    	}
    	
    	cplex.addMaximize(objective);
    	
    	//Constraints Fy + Lx <= f (b)
    	IloLinearNumExpr str = cplex.linearNumExpr(); 
    	
    	 //Y
    	for (int j = 0; j < I.number_of_fc; j++)
    	{
    		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_2.addTerm(I.F.get(j).get(i), y.get(i));
        	}
    		 
    		cplex.addLe(str_2, I.f.get(j));
         	
        }
    	
    	for (int j = 0; j < I.number_of_items; j++)
    	{
    		IloLinearNumExpr str1 = cplex.linearNumExpr(); 
    		str1.addTerm(1., y.get(j));
        	cplex.addLe(str1, 1 - x.get(j));
    	}
    
    	
    	if ("Risk-averse".equals(type_of_risk))
    	{   
    		for (int k = 0; k < sample_size; k++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			
    			 for (int i = 0; i < I.number_of_items; i++)
    			 {
    				 str1.addTerm(-data_set.get(i).get(k), y.get(i));
    		     }
    			
    			str1.addTerm(1., t_f);
    			
    			str1.addTerm(-1., s_f.get(k));
    			cplex.addLe(str1, 0.);
    			
    			/*IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str2.addTerm(1., s_f.get(k));
    			cplex.addGe(str2, 0);*/
    		}
    		
    	}
    	
    	ArrayList<Double> y_opt = new ArrayList<Double>();
    	
    	if (cplex.solve())
    	{   		
    		/*for (int k = 0; k < sample_size; k++)
        	{
        		for (int l = 0; l < number_of_sc; l++)
        		{
        			System.out.println(cplex.getValue(nu_f.get(k).get(l)));	
        		}
        	}*/
    		
    		
    	    for (int j = 0; j < I.number_of_items; j++)
    	    {
    			y_opt.add(cplex.getValue(y.get(j)));
    	    }  
    	    
    	    y_opt.add(cplex.getObjValue());
    	}
    	else 
        	System.out.println("No solution found");
    	
    	return y_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}
	
	public static Double SolveFollowersProblemFullInformationFixedY(IloCplex cplex, Instance I, ArrayList<Double> x, ArrayList<Double> y, double alpha_f, int sample_size, ArrayList<ArrayList<Double>> data_set,  ArrayList<ArrayList<Double>> support_constraints, String type_of_risk, ArrayList<Double> average_f,  double binary_indicator) {
	try 
    {
		if ("Risk-neutral".equals(type_of_risk))
    	{
    		double sum = 0;
    		for (int i = 0; i < I.number_of_items; i++)
    			sum += average_f.get(i) * y.get(i);
    		
    		return sum; 	
    	}
		
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
        
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	ArrayList<IloNumVar> s_f = new ArrayList<IloNumVar>();
    	IloNumVar t_f = cplex.numVar(-M, M, "t");
    	
    	for (int k = 0; k < sample_size; k++)
    	{
    		s_f.add(cplex.numVar(0, M, "s_"+ k));
    		
    	}
    	
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	if ("Risk-averse".equals(type_of_risk))
    	{   
    		objective.addTerm(1., t_f);
    		
    		for (int k = 0; k < sample_size; k++)
        	{
        	     objective.addTerm(-1./(sample_size*(1 - alpha_f)), s_f.get(k));
        	}
    	}
    	
    	cplex.addMaximize(objective);
    	
    	if ("Risk-averse".equals(type_of_risk))
    	{   
    		for (int k = 0; k < sample_size; k++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			
    			double sum = 0;
    			for (int i = 0; i < I.number_of_items; i++)
    			{
    				 sum += -data_set.get(i).get(k) * y.get(i);
    		    }
    			
    			str1.addTerm(1., t_f);
    			str1.addTerm(-1., s_f.get(k));
    			
    			cplex.addLe(str1, -sum);
    			
    			/*IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str2.addTerm(1., s_f.get(k));
    			cplex.addGe(str2, 0);*/
    		}
    		
    	}
    	
    	double z_opt = -1; 
    	
    	if (cplex.solve())
    	{   		
    	    z_opt = cplex.getObjValue();
    	}
    	else 
        	System.out.println("No solution found");
    	
    	return z_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return -1.;
}

	
		public static ArrayList<ArrayList<Double>> SolveLeadersProblemFullFixedX(IloCplex cplex, Instance I, ArrayList<Double> x, double epsilon_f, double alpha_l, double alpha_f, int k_l, int k_f, ArrayList<ArrayList<Double>> data_set_l, ArrayList<ArrayList<Double>> data_set_f,
				ArrayList<ArrayList<Double>> support_constraints, String type_of_risk_l, String type_of_risk_f,  ArrayList<ArrayList<Double>> delta_f, ArrayList<Double> average_l, ArrayList<Double> average_f)
	    {
		try 
	    {
	    	cplex.clearModel();
	    	cplex.setOut(null);
	    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
	    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
	    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
	         	
	        double M = Double.MAX_VALUE;
	    	int number_of_sc = support_constraints.size() - 1; 
	    	
	    	ArrayList<ArrayList<IloNumVar>> nu_f = new ArrayList<ArrayList<IloNumVar>>();
	    	
	    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
	    	
	    	ArrayList<IloNumVar> s_l = new ArrayList<IloNumVar>();
	    	ArrayList<IloNumVar> s_f = new ArrayList<IloNumVar>();
	    	
	    	//IloNumVar lambda_l = cplex.numVar(0, M, "lam_l");
	    	IloNumVar lambda_f = cplex.numVar(0, M, "lam_f");
	    	
	    	IloNumVar t_l = cplex.numVar(-M, M, "t_l");
	    	IloNumVar t_f = cplex.numVar(-M, M, "t_f");
	    	
	    	ArrayList<ArrayList<IloNumVar>> mu1_f = new ArrayList<ArrayList<IloNumVar>>();
	    	ArrayList<ArrayList<IloNumVar>> mu2_f = new ArrayList<ArrayList<IloNumVar>>();
	    	
	    	ArrayList<IloNumVar> beta1_f = new ArrayList<IloNumVar>();
	    	ArrayList<IloNumVar> beta2_f = new ArrayList<IloNumVar>();
	    	ArrayList<IloNumVar> gamma_f = new ArrayList<IloNumVar>();
	    	//ArrayList<IloNumVar> prod_f = new ArrayList<IloNumVar>();
	    	
	    	
	    	for (int k = 0; k < k_l; k++)
	    	{
	    		s_l.add(cplex.numVar(0, M, "s_l_"+ k));
	    	}
	    	
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		s_f.add(cplex.numVar(0, M, "s_f_"+ k));
	    		nu_f.add(new ArrayList<IloNumVar>());
	    	
	    		for (int l = 0; l < number_of_sc; l++)
	    		{
	    			nu_f.get(k).add(cplex.numVar(0, M, "nu_f_"+ k + l));	
	    		}
	    	}
	    
	        for (int j = 0; j < I.number_of_items; j++)
	        {
	        	y.add(cplex.numVar(0, 1, "y_" + j));  	
	        	//prod_f.add(cplex.numVar(0, M, "prod_" + j));
	        }
	    	
	        //Follower's dual variables
	        for (int j = 0; j < I.number_of_fc; j++)
	        {
	        	beta1_f.add(cplex.numVar(0, M, "beta1_f_"+ j));
	        }
	        
	        for (int i = 0; i < I.number_of_items; i++)
	        {
	        	beta2_f.add(cplex.numVar(0, M, "beta2_f_"+ i));
	        }
	        
	        for (int k = 0; k < k_f; k++)
	    	{
	    		gamma_f.add(cplex.numVar(0, M, "gamma_f_"+ k));
	    		
	    		mu1_f.add(new ArrayList<IloNumVar>());
	    		mu2_f.add(new ArrayList<IloNumVar>());
	    	
	    		 for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			mu1_f.get(k).add(cplex.numVar(0, M, "mu1_f_"+ k + j));	
	    			mu2_f.get(k).add(cplex.numVar(0, M, "mu2_f_"+ k + j));
	    		}
	    	}
	        
	    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
	    	
	    	if ("Risk-averse".equals(type_of_risk_l))
	    	{   
	    		objective.addTerm(1., t_l);
	   
	    		for (int k = 0; k < k_l; k++)
	        	{
	        	     objective.addTerm(1./(k_l*(1. - alpha_l)), s_l.get(k));
	        	}
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_l))
	    	{
	    		for (int i = 0; i < I.number_of_items; i++)
	    			objective.addTerm(average_l.get(i), y.get(i));
	    	  		
	    	}
	    	
	    	cplex.addMinimize(objective);
	    	
	    	//Constraints
	    	//Feasibility
	    	for (int j = 0; j < I.number_of_fc; j++)
	    	{
	    		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
	    		
	    		for (int i = 0; i < I.number_of_items; i++)
	        	{
	    			str_2.addTerm(I.F.get(j).get(i), y.get(i));
	        	}
	    		 
	    		cplex.addLe(str_2, I.f.get(j));
	         	
	        }
	    	
	    	for (int j = 0; j < I.number_of_items; j++)
	    	{
	    		IloLinearNumExpr str = cplex.linearNumExpr(); 
	    		str.addTerm(1., y.get(j));
	    	    cplex.addLe(str, 1. - x.get(j));
	    	}
	        
	    	//
	     	/*for (int k = 0; k < k_l; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	    	
	    	    	str1.addTerm(-1., lambda_l);
	    	    	str2.addTerm(-1., lambda_l);
	    	    	
	    	    	str1.addTerm(1., y.get(j));
	    	    	str2.addTerm(-1., y.get(j));
	    	    	
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(support_constraints.get(l).get(j), nu_l.get(k).get(l));
	    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_l.get(k).get(l));
	     	        }
	    	    	
	    	    	cplex.addLe(str1, 0.);
	    	    	cplex.addLe(str2, 0.);
	    		}*/
	    	
	     	//
	    	for (int k = 0; k < k_f; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	    	
	    	    	str1.addTerm(-1., lambda_f);
	    	    	str2.addTerm(-1., lambda_f);
	    	    	
	    	    	str1.addTerm(1., y.get(j));
	    	    	str2.addTerm(-1., y.get(j));
	    	    	
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(k).get(l));
	    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(k).get(l));
	     	        }
	    	    	
	    	    	cplex.addLe(str1, 0.);
	    	    	cplex.addLe(str2, 0.);
	    		}
	   	    
	    	//
	    	if ("Risk-averse".equals(type_of_risk_l))
	    	{   
	    		for (int k = 0; k < k_l; k++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    			
	    			 for (int i = 0; i < I.number_of_items; i++)
	    			 {
	    				 str1.addTerm(data_set_l.get(i).get(k), y.get(i));
	    		     }
	    			
	    			str1.addTerm(-1., t_l);
	    			
	    			str1.addTerm(-1., s_l.get(k));
	    			cplex.addLe(str1, 0.);
	    			
	                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
	    			
	    			str2.addTerm(1., s_l.get(k));
	    			cplex.addGe(str2, 0);*/
	    		}
	    		
	    	}
	    	
	    	//
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{   
	    		for (int k = 0; k < k_f; k++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    			
	    			 for (int i = 0; i < I.number_of_items; i++)
	    			 {
	    				 str1.addTerm(-data_set_f.get(i).get(k), y.get(i));
	    		     }
	    			
	    			str1.addTerm(1., t_f);
	    			
	    
	        		for (int l = 0; l < number_of_sc; l++)
	        		{
	        			str1.addTerm(delta_f.get(k).get(l), nu_f.get(k).get(l));
	        		}
	    			
	    			str1.addTerm(-1., s_f.get(k));
	    			cplex.addLe(str1, 0.);
	    			
	                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
	    			
	    			str2.addTerm(1., s_f.get(k));
	    			cplex.addGe(str2, 0);*/
	    		}
	    		
	    	}
	    	
	    	//Dual constraints of the follower
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{
	    		//
		    	for (int i = 0; i < I.number_of_items; i++)
		    	 {	
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
					
		    		for (int j = 0; j < I.number_of_fc; j++)
		           	{
		       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
		           	}
		       		
		       		str1.addTerm(1., beta2_f.get(i));
		    		
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			str1.addTerm(-data_set_f.get(i).get(k), gamma_f.get(k));
		    			str1.addTerm(-1., mu1_f.get(k).get(i));
		    			str1.addTerm(1., mu2_f.get(k).get(i));
		    		}
		        	
		    		cplex.addGe(str1, 0.);
		    	}
		    	
		    	//
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int l = 0; l < number_of_sc; l++)
		            {
		    			IloLinearNumExpr str1 = cplex.linearNumExpr();
			    		
		    			str1.addTerm(delta_f.get(k).get(l), gamma_f.get(k));
			    		
			    		for (int i = 0; i < I.number_of_items; i++)
				    	 {
			    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
			    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
				    	 }
			    		cplex.addGe(str1, 0.);
			    	} 
		    	}
		    	
		    	//
		    	IloLinearNumExpr str = cplex.linearNumExpr();
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		
		    		str.addTerm(1., gamma_f.get(k));
		    		str1.addTerm(1., gamma_f.get(k));
		    		
		    		cplex.addLe(str1, 1./((1. - alpha_f) * k_f));
		    	}
		    	
		    	cplex.addEq(str, 1.);
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    	{   
	    		//
		    	for (int i = 0; i < I.number_of_items; i++)
		    	 {	
		    		IloLinearNumExpr str1 = cplex.linearNumExpr(-average_f.get(i));
					
		    		for (int j = 0; j < I.number_of_fc; j++)
		           	{
		       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
		           	}
		       		
		       		str1.addTerm(1., beta2_f.get(i));
		    		
		    		
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			str1.addTerm(-1., mu1_f.get(k).get(i));
		    			str1.addTerm(1., mu2_f.get(k).get(i));
		    		}
		        	
		    		cplex.addGe(str1, 0.);
		    	}
		    	
		    	//
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int l = 0; l < number_of_sc; l++)
		            {
		    			IloLinearNumExpr str1 = cplex.linearNumExpr(delta_f.get(k).get(l)/(k_f + 0.));
			    		
			    		for (int i = 0; i < I.number_of_items; i++)
				    	 {
			    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
			    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
				    	 }
			    		cplex.addGe(str1, 0.);
			    	}
		    	}
	    	}
	    	
	    	//
	    	IloLinearNumExpr str3 = cplex.linearNumExpr();
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		for (int i = 0; i < I.number_of_items; i++)
		    	 {
	    			str3.addTerm(1., mu1_f.get(k).get(i));
	    			str3.addTerm(1., mu2_f.get(k).get(i));
		    	 }
	    	}
	    	
	    	if ("Risk-averse".equals(type_of_risk_f))
	    		cplex.addLe(str3, epsilon_f/(1. - alpha_f));
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    		cplex.addLe(str3, epsilon_f);
	    	
	    	
	    	//Strong duality 
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	
	    		
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str1.addTerm(1 - x.get(i), beta2_f.get(i));
		    	}
		    	
		    	for (int j = 0; j < I.number_of_fc; j++)
		    	{
		    		str1.addTerm(I.f.get(j), beta1_f.get(j));
		    	}
		    	
		    	str2.addTerm(1., t_f);
		    	str2.addTerm(-epsilon_f/(1. - alpha_f), lambda_f);
		    	
		    	
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		str2.addTerm(-1./(k_f * (1. - alpha_f)), s_f.get(k));
		    	}
		    	
		    	cplex.addEq(str1, str2);
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str1.addTerm(1 - x.get(i), beta2_f.get(i));
		        }
		        
		    	for (int j = 0; j < I.number_of_fc; j++)
		    	{
		    		str1.addTerm(I.f.get(j), beta1_f.get(j));
		    	}
		    	
		    	str2.addTerm(-epsilon_f, lambda_f);
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str2.addTerm(average_f.get(i), y.get(i));
		    	}
		    	
		    	for (int k = 0; k < k_f; k++)
	    			for (int l = 0; l < number_of_sc; l++)
	    		{
	    			str2.addTerm(-delta_f.get(k).get(l)/(k_f + 0.), nu_f.get(k).get(l));
	    		}
		    	
		    	cplex.addEq(str1, str2);
	    	}
	    	
	    	/*for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr(-M1);
	    		IloLinearNumExpr str4 = cplex.linearNumExpr();
	    		
	    		str1.addTerm(1., prod_f.get(i));
	    		str1.addTerm(-1., beta_f.get(i));
	        	
	    		str2.addTerm(1., prod_f.get(i));
	    		str2.addTerm(M1, x.get(i));
	        	
	    		str4.addTerm(1., prod_f.get(i));
	    		str4.addTerm(-1, beta_f.get(i));
	    		str4.addTerm(M1, x.get(i));
	        	
	    		cplex.addLe(str1, 0.);
	    		cplex.addLe(str2, 0.);
	    		cplex.addGe(str4, 0.);
	    		
	    	}*/
	    	
	    	
	    	ArrayList<ArrayList<Double>> sol_opt = new ArrayList<ArrayList<Double>>();
	    	ArrayList<Double> y_opt = new ArrayList<Double>();
	    	ArrayList<Double> x_opt = new ArrayList<Double>();
	    	
	    	
	    	if (cplex.solve())
	    	{   	
	
	    	    /*for (int j = 0; j < I.number_of_items; j++)
	    	    {
	    	    	System.out.println((1 - x.get(j)) * cplex.getValue(beta_f.get(j)));
	    	    }*/
	    	      		
	    	    for (int j = 0; j < I.number_of_items; j++)
	    	    {
	    	    	//x_opt.add(cplex.getValue(x.get(j)));
	    			y_opt.add(cplex.getValue(y.get(j)));
	    	    }  
	    	    
	    	    x_opt.add(cplex.getObjValue());
	    	    
	    	    sol_opt.add(x_opt);
	    	    sol_opt.add(y_opt);
	    	}
	    	else 
	        	System.out.println("No solution found");
	    	
	    	return sol_opt;
		}
	    catch (IloException Exc) { Exc.printStackTrace(); }
	    return new ArrayList<ArrayList<Double>>();
	}
	
		
	public static ArrayList<ArrayList<Double>> SolveLeadersProblem(IloCplex cplex, Instance I, double epsilon_l, double epsilon_f, double alpha_l, double alpha_f, int k_l, int k_f, ArrayList<ArrayList<Double>> data_set_l, ArrayList<ArrayList<Double>> data_set_f,
			ArrayList<ArrayList<Double>> support_constraints, String type_of_risk_l, String type_of_risk_f,  ArrayList<ArrayList<Double>> delta_l,  ArrayList<ArrayList<Double>> delta_f, ArrayList<Double> average_l, ArrayList<Double> average_f, double M1)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
    	cplex.setParam(IloCplex.Param.TimeLimit, 3600.0);
       
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	ArrayList<ArrayList<IloNumVar>> nu_l = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> nu_f = new ArrayList<ArrayList<IloNumVar>>();
    	
    	ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
    	
    	ArrayList<IloNumVar> s_l = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> s_f = new ArrayList<IloNumVar>();
    	
    	IloNumVar lambda_l = cplex.numVar(0, M, "lam_l");
    	IloNumVar lambda_f = cplex.numVar(0, M, "lam_f");
    	
    	IloNumVar t_l = cplex.numVar(-M, M, "t_l");
    	IloNumVar t_f = cplex.numVar(-M, M, "t_f");
    	
    	
    	ArrayList<ArrayList<IloNumVar>> mu1_f = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> mu2_f = new ArrayList<ArrayList<IloNumVar>>();
    	
    	ArrayList<IloNumVar> beta1_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> beta2_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> gamma_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> prod_f = new ArrayList<IloNumVar>();
    	
    	for (int k = 0; k < k_l; k++)
    	{
    		s_l.add(cplex.numVar(0, M, "s_l_"+ k));
    		nu_l.add(new ArrayList<IloNumVar>());
    	
    		for (int l = 0; l < number_of_sc; l++)
    		{
    			nu_l.get(k).add(cplex.numVar(0, M, "nu_l_"+ k + l));	
    		}
    	}
    	
    	for (int k = 0; k < k_f; k++)
    	{
    		s_f.add(cplex.numVar(0, M, "s_f_"+ k));
    		nu_f.add(new ArrayList<IloNumVar>());
    	
    		for (int l = 0; l < number_of_sc; l++)
    		{
    			nu_f.get(k).add(cplex.numVar(0, M, "nu_f_"+ k + l));	
    		}
    	}
    
        for (int j = 0; j < I.number_of_items; j++)
        {
        	x.add(cplex.boolVar("x_" + j));
        	y.add(cplex.numVar(0, 1, "y_" + j));
        	
        	prod_f.add(cplex.numVar(0, M, "prod_" + j));
        }
    	
        //Follower's dual variables
        for (int j = 0; j < I.number_of_fc; j++)
        {
        	beta1_f.add(cplex.numVar(0, M, "beta1_f_"+ j));
        }
        
        for (int i = 0; i < I.number_of_items; i++)
        {
        	beta2_f.add(cplex.numVar(0, M, "beta2_f_"+ i));
        }
        
        for (int k = 0; k < k_f; k++)
    	{
    		gamma_f.add(cplex.numVar(0, M, "gamma_f_"+ k));
    		
    		mu1_f.add(new ArrayList<IloNumVar>());
    		mu2_f.add(new ArrayList<IloNumVar>());
    	
    		 for (int j = 0; j < I.number_of_items; j++)
    		{
    			mu1_f.get(k).add(cplex.numVar(0, M, "mu1_f_"+ k + j));	
    			mu2_f.get(k).add(cplex.numVar(0, M, "mu2_f_"+ k + j));
    		}
    	}
        
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	if ("Risk-averse".equals(type_of_risk_l))
    	{   
    		objective.addTerm(1., t_l);
    		objective.addTerm(epsilon_l/(1. - alpha_l), lambda_l);
    		
    		for (int k = 0; k < k_l; k++)
        	{
        	     objective.addTerm(1./(k_l*(1. - alpha_l)), s_l.get(k));
        	}
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk_l))
    	{
    		objective.addTerm(epsilon_l, lambda_l);
    		
    		for (int i = 0; i < I.number_of_items; i++)
    			objective.addTerm(average_l.get(i), y.get(i));
    		
    		
    		for (int k = 0; k < k_l; k++)
    		{
    			for (int l = 0; l < number_of_sc; l++)
	    		{
	    			objective.addTerm(delta_l.get(k).get(l)/(k_l + 0.), nu_l.get(k).get(l));
	    		}
    		
    			/*for (int i = 0; i < I.number_of_items; i++)
       			 {
       				 objective.addTerm(data_set_l.get(i).get(k)/(k_l + 0.), y.get(i));
       		     }*/
    		}
    			
    		
    	}
    	
    	cplex.addMinimize(objective);
    	
    	//Constraints
    	//Feasibility
    	//X
    	for (int j = 0; j < I.number_of_lc; j++)
    	{
    		IloLinearNumExpr str_1 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_1.addTerm(I.H.get(j).get(i), x.get(i));
        	}
    		 
    		cplex.addLe(str_1, I.h.get(j));
         	
        }
    	
    	 //Y
    	for (int j = 0; j < I.number_of_fc; j++)
    	{
    		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_2.addTerm(I.F.get(j).get(i), y.get(i));
        	}
    		 
    		cplex.addLe(str_2, I.f.get(j));
         	
        }
   
    	// y <= 1 - x 
    	for (int i = 0; i < I.number_of_items; i++)
    	{
    		IloLinearNumExpr str = cplex.linearNumExpr(); 
    		str.addTerm(1., y.get(i));
    		str.addTerm(1., x.get(i));
        	cplex.addLe(str, 1.);
    	}
        
    	//
     	for (int k = 0; k < k_l; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
    	    	
    	    	str1.addTerm(-1., lambda_l);
    	    	str2.addTerm(-1., lambda_l);
    	    	
    	    	str1.addTerm(-1., y.get(j));
    	    	str2.addTerm(1., y.get(j));
    	    	
    	    	for (int l = 0; l < number_of_sc; l++)
     	        {
    	    		str1.addTerm(support_constraints.get(l).get(j), nu_l.get(k).get(l));
    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_l.get(k).get(l));
     	        }
    	    	
    	    	cplex.addLe(str1, 0.);
    	    	cplex.addLe(str2, 0.);
    		}
    	
     	//
    	for (int k = 0; k < k_f; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
    	    	
    	    	str1.addTerm(-1., lambda_f);
    	    	str2.addTerm(-1., lambda_f);
    	    	
    	    	str1.addTerm(1., y.get(j));
    	    	str2.addTerm(-1., y.get(j));
    	    	
    	    	for (int l = 0; l < number_of_sc; l++)
     	        {
    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(k).get(l));
    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(k).get(l));
     	        }
    	    	
    	    	cplex.addLe(str1, 0.);
    	    	cplex.addLe(str2, 0.);
    		}
   	    
    	//
    	if ("Risk-averse".equals(type_of_risk_l))
    	{   
    		for (int k = 0; k < k_l; k++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			
    			 for (int i = 0; i < I.number_of_items; i++)
    			 {
    				 str1.addTerm(data_set_l.get(i).get(k), y.get(i));
    		     }
    			
    			str1.addTerm(-1., t_l);
    			
    
        		for (int l = 0; l < number_of_sc; l++)
        		{
        			str1.addTerm(delta_l.get(k).get(l), nu_l.get(k).get(l));
        		}
    			
    			str1.addTerm(-1., s_l.get(k));
    			cplex.addLe(str1, 0.);
    			
                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str2.addTerm(1., s_l.get(k));
    			cplex.addGe(str2, 0);*/
    		}
    		
    	}
    	
    	//
    	if ("Risk-averse".equals(type_of_risk_f))
    	{   
    		for (int k = 0; k < k_f; k++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			
  
    			 for (int i = 0; i < I.number_of_items; i++)
    			 {
    				 str1.addTerm(-data_set_f.get(i).get(k), y.get(i));
    		     }
    			
    			str1.addTerm(1., t_f);
    			
    
        		for (int l = 0; l < number_of_sc; l++)
        		{
        			str1.addTerm(delta_f.get(k).get(l), nu_f.get(k).get(l));
        		}
    			
    			str1.addTerm(-1., s_f.get(k));
    			cplex.addLe(str1, 0.);
    			
                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str2.addTerm(1., s_f.get(k));
    			cplex.addGe(str2, 0);*/
    		}
    		
    	}
    	
    	//Dual constraints of the follower
    	if ("Risk-averse".equals(type_of_risk_f))
    	{
    		//
	    	for (int i = 0; i < I.number_of_items; i++)
	    	 {	
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
				
	    		for (int j = 0; j < I.number_of_fc; j++)
	           	{
	       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
	           	}
	       		
	       		str1.addTerm(1., beta2_f.get(i));
	    		
	    		for (int k = 0; k < k_f; k++)
	    		{
	    			str1.addTerm(-data_set_f.get(i).get(k), gamma_f.get(k));
	    			str1.addTerm(-1., mu1_f.get(k).get(i));
	    			str1.addTerm(1., mu2_f.get(k).get(i));
	    		}
	        	
	    		cplex.addGe(str1, 0.);
	    	}
	    	
	    	//
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		for (int l = 0; l < number_of_sc; l++)
	            {
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		
	    			str1.addTerm(delta_f.get(k).get(l), gamma_f.get(k));
		    		
		    		for (int i = 0; i < I.number_of_items; i++)
			    	 {
		    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
		    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
			    	 }
		    		cplex.addGe(str1, 0.);
		    	} 
	    	}
	    	
	    	//
	    	IloLinearNumExpr str = cplex.linearNumExpr();
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		
	    		str.addTerm(1., gamma_f.get(k));
	    		str1.addTerm(1., gamma_f.get(k));
	    		
	    		cplex.addLe(str1, 1./((1. - alpha_f) * k_f));
	    	}
	    	
	    	cplex.addEq(str, 1.);
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk_f))
    	{   
    		
    		
    		//
	    	for (int i = 0; i < I.number_of_items; i++)
	    	 {	
	    		IloLinearNumExpr str1 = cplex.linearNumExpr(-average_f.get(i));
				
	    		for (int j = 0; j < I.number_of_fc; j++)
	           	{
	       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
	           	}
	       		
	       		str1.addTerm(1., beta2_f.get(i));
	    		
	    		for (int k = 0; k < k_f; k++)
	    		{
	    			str1.addTerm(-1., mu1_f.get(k).get(i));
	    			str1.addTerm(1., mu2_f.get(k).get(i));
	    		}
	        	
	    		cplex.addGe(str1, 0.);
	    	}
	    	
	    	//
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		for (int l = 0; l < number_of_sc; l++)
	            {
	    			IloLinearNumExpr str1 = cplex.linearNumExpr(delta_f.get(k).get(l)/(k_f + 0.));
		    		
		    		for (int i = 0; i < I.number_of_items; i++)
			    	 {
		    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
		    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
			    	 }
		    		cplex.addGe(str1, 0.);
		    	}
	    	}
    	}
    	
    	//
    	IloLinearNumExpr str3 = cplex.linearNumExpr();
    	for (int k = 0; k < k_f; k++)
    	{
    		for (int i = 0; i < I.number_of_items; i++)
	    	 {
    			str3.addTerm(1., mu1_f.get(k).get(i));
    			str3.addTerm(1., mu2_f.get(k).get(i));
	    	 }
    	}
    	
    	if ("Risk-averse".equals(type_of_risk_f))
    		cplex.addLe(str3, epsilon_f/(1. - alpha_f));
    	
    	if ("Risk-neutral".equals(type_of_risk_f))
    		cplex.addLe(str3, epsilon_f);
    	
    	
    	//Strong duality 
    	if ("Risk-averse".equals(type_of_risk_f))
    	{
    		IloLinearNumExpr str1 = cplex.linearNumExpr();
    		IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	
    		
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		str1.addTerm(1., prod_f.get(i));
	    	}
	    	
	    	for (int j = 0; j < I.number_of_fc; j++)
	    	{
	    		str1.addTerm(I.f.get(j), beta1_f.get(j));
	    	}
	    	
	    	str2.addTerm(1., t_f);
	    	str2.addTerm(-epsilon_f/(1. - alpha_f), lambda_f);
	    	
	    	
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		str2.addTerm(-1./(k_f * (1. - alpha_f)), s_f.get(k));
	    	}
	    	
	    	cplex.addEq(str1, str2);
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk_f))
    	{
    		IloLinearNumExpr str1 = cplex.linearNumExpr();
    		IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		str1.addTerm(1., prod_f.get(i));
	    		
	    	}
	        
	    	for (int j = 0; j < I.number_of_fc; j++)
	    	{
	    		str1.addTerm(I.f.get(j), beta1_f.get(j));
	    	}
	    	
	    	str2.addTerm(-epsilon_f, lambda_f);
	    	
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		str2.addTerm(average_f.get(i), y.get(i));
	    	}
	    	
	    	for (int k = 0; k < k_f; k++)
    			for (int l = 0; l < number_of_sc; l++)
    		{
    			str2.addTerm(-delta_f.get(k).get(l)/(k_f + 0.), nu_f.get(k).get(l));
    		}
	    	
	    	cplex.addEq(str1, str2);
    	}
    	
    	for (int i = 0; i < I.number_of_items; i++)
    	{
    		IloLinearNumExpr str1 = cplex.linearNumExpr();
    		IloLinearNumExpr str2 = cplex.linearNumExpr(-M1);
    		IloLinearNumExpr str4 = cplex.linearNumExpr();
    		
    		str1.addTerm(1., prod_f.get(i));
    		str1.addTerm(-1., beta2_f.get(i));
        	
    		str2.addTerm(1., prod_f.get(i));
    		str2.addTerm(M1, x.get(i));
        	
    		str4.addTerm(1., prod_f.get(i));
    		str4.addTerm(-1, beta2_f.get(i));
    		str4.addTerm(M1, x.get(i));
        	
    		cplex.addLe(str1, 0.);
    		cplex.addLe(str2, 0.);
    		cplex.addGe(str4, 0.);
    	}
    	
    	
    	ArrayList<ArrayList<Double>> sol_opt = new ArrayList<ArrayList<Double>>();
    	ArrayList<Double> y_opt = new ArrayList<Double>();
    	ArrayList<Double> x_opt = new ArrayList<Double>();
    	
    	
    	if (cplex.solve())
    	{   	    		 		
    	    for (int j = 0; j < I.number_of_items; j++)
    	    {
    	    	x_opt.add(cplex.getValue(x.get(j)));
    			y_opt.add(cplex.getValue(y.get(j)));
    	    }  
    	    
    	    x_opt.add(cplex.getObjValue());
    	    
    	    sol_opt.add(x_opt);
    	    sol_opt.add(y_opt);
    	}
    	else 
        	System.out.println("No solution found");
    	
    	if (cplex.getStatus() != IloCplex.Status.Optimal)
    	{
    		 return new ArrayList<ArrayList<Double>>();
    	}
    	
    	return sol_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<ArrayList<Double>>();
}
	
	public static ArrayList<ArrayList<Double>> SolveLeadersProblemNominal(IloCplex cplex, Instance I,   
		 double M1)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
        
        double M = Double.MAX_VALUE;
    	
    	
    	ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
    	
    	ArrayList<IloNumVar> beta1_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> beta2_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> prod_f = new ArrayList<IloNumVar>();
    	
        for (int i = 0; i < I.number_of_items; i++)
        {
        	x.add(cplex.boolVar("x_" + i));
        	y.add(cplex.numVar(0, 1, "y_" + i));
        	prod_f.add(cplex.numVar(0, M, "prod_" + i));
        }
    	
        //Follower's dual variables
        for (int j = 0; j < I.number_of_fc; j++)
        {
        	beta1_f.add(cplex.numVar(0, M, "beta1_f_"+ j));
        }
        
        for (int i = 0; i < I.number_of_items; i++)
        {
        	beta2_f.add(cplex.numVar(0, M, "beta2_f_"+ i));
        }
        
      
        IloLinearNumExpr objective = cplex.linearNumExpr();
      
        for (int i = 0; i < I.number_of_items; i++)
        {
        	objective.addTerm(I.mean.get(i), y.get(i));
        }
    		
    	cplex.addMinimize(objective);
    	
    	//Constraints
    	//X
    	for (int j = 0; j < I.number_of_lc; j++)
    	{
    		IloLinearNumExpr str_1 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_1.addTerm(I.H.get(j).get(i), x.get(i));
        	}
    		 
    		cplex.addLe(str_1, I.h.get(j));
         	
        }
    	
    	 //Y
    	for (int j = 0; j < I.number_of_fc; j++)
    	{
    		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_2.addTerm(I.F.get(j).get(i), y.get(i));
        	}
    		 
    		cplex.addLe(str_2, I.f.get(j));
         	
        }
   
    	// y <= 1 - x 
    	for (int i = 0; i < I.number_of_items; i++)
    	{
    		IloLinearNumExpr str = cplex.linearNumExpr(); 
    		str.addTerm(1., y.get(i));
    		str.addTerm(1., x.get(i));
        	cplex.addLe(str, 1.);
    	}
        
    
    	//Dual constraints of the follower
    	for (int i = 0; i < I.number_of_items; i++)
    	 {	
    		IloLinearNumExpr str1 = cplex.linearNumExpr(-I.mean.get(i));
			
    		for (int j = 0; j < I.number_of_fc; j++)
        	{
    			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
        	}
    		
    		str1.addTerm(1., beta2_f.get(i));
 
    		cplex.addGe(str1, 0.);
    	}
	    
    	//Strong duality 
    
    	IloLinearNumExpr str1 = cplex.linearNumExpr();
    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	
    	for (int j = 0; j < I.number_of_fc; j++)
    	{
    		str1.addTerm(I.f.get(j), beta1_f.get(j));
    	}
    	
    	for (int i = 0; i < I.number_of_items; i++)
    	{
    		str1.addTerm(1., prod_f.get(i));
    
    	}
    
    	for (int i = 0; i < I.number_of_items; i++)
    	{
    		str2.addTerm(I.mean.get(i), y.get(i));
    	}
    
    	cplex.addEq(str1, str2);
   
    	
    	for (int i = 0; i < I.number_of_items; i++)
    	{
    		IloLinearNumExpr str3 = cplex.linearNumExpr();
    		IloLinearNumExpr str4 = cplex.linearNumExpr(-M1);
    		IloLinearNumExpr str5 = cplex.linearNumExpr();
    		
    		str3.addTerm(1., prod_f.get(i));
    		str3.addTerm(-1., beta2_f.get(i));
        	
    		str4.addTerm(1., prod_f.get(i));
    		str4.addTerm(M1, x.get(i));
        	
    		str5.addTerm(1., prod_f.get(i));
    		str5.addTerm(-1, beta2_f.get(i));
    		str5.addTerm(M1, x.get(i));
        	
    		cplex.addLe(str3, 0.);
    		cplex.addLe(str4, 0.);
    		cplex.addGe(str5, 0.);	
    	}
    	
    	
    	ArrayList<ArrayList<Double>> sol_opt = new ArrayList<ArrayList<Double>>();
    	ArrayList<Double> y_opt = new ArrayList<Double>();
    	ArrayList<Double> x_opt = new ArrayList<Double>();
    	
    	
    	if (cplex.solve())
    	{   	    		 		
    	    for (int j = 0; j < I.number_of_items; j++)
    	    {
    	    	x_opt.add(cplex.getValue(x.get(j)));
    			y_opt.add(cplex.getValue(y.get(j)));
    	    }  
    	    
    	    x_opt.add(cplex.getObjValue());
    	    
    	    sol_opt.add(x_opt);
    	    sol_opt.add(y_opt);
    	}
    	else 
        	System.out.println("No solution found");
    	
    	return sol_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<ArrayList<Double>>();
}
	
	
/*	public static ArrayList<ArrayList<Double>> SolveLeadersProblemModified(IloCplex cplex, Instance I, double epsilon_l, double epsilon_f, int k_l, int k_f, ArrayList<ArrayList<Double>> data_set_l, ArrayList<ArrayList<Double>> data_set_f,
			ArrayList<ArrayList<Double>> support_constraints, String type_of_risk_l, String type_of_risk_f,  ArrayList<ArrayList<Double>> delta_l,  ArrayList<ArrayList<Double>> delta_f, ArrayList<Double> average_l, ArrayList<Double> average_f, double M1, double binary_indicator)
	{
		try 
	    {
	    	cplex.clearModel();
	    	cplex.setOut(null);
	    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
	    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
	    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
	        //cplex.setParam(IloCplex.Param.Preprocessing.RepeatPresolve, 0);
	         
	        
	        double M = Double.MAX_VALUE;
	    	int number_of_sc = support_constraints.size() - 1; 
	    	
	    	ArrayList<ArrayList<IloNumVar>> nu_l = new ArrayList<ArrayList<IloNumVar>>();
	    	ArrayList<ArrayList<IloNumVar>> nu_f = new ArrayList<ArrayList<IloNumVar>>();
	    	
	    	ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
	    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
	    	
	    	IloNumVar lambda_l = cplex.numVar(0, M, "lam_l");
	    	IloNumVar lambda_f = cplex.numVar(0, M, "lam_f");
	    	
	    	IloNumVar phi_f = cplex.numVar(0, M, "phi");
	    	
	    	ArrayList<ArrayList<IloNumVar>> mu1_f = new ArrayList<ArrayList<IloNumVar>>();
	    	ArrayList<ArrayList<IloNumVar>> mu2_f = new ArrayList<ArrayList<IloNumVar>>();
	    	
	    	ArrayList<IloNumVar> beta_f = new ArrayList<IloNumVar>();
	    	ArrayList<IloNumVar> prod_f = new ArrayList<IloNumVar>();
	    	
	    	
	    	
	    	for (int k = 0; k < k_l; k++)
	    	{
	    		nu_l.add(new ArrayList<IloNumVar>());
	    	
	    		for (int l = 0; l < number_of_sc; l++)
	    		{
	    			nu_l.get(k).add(cplex.numVar(0, M, "nu_l_"+ k + l));	
	    		}
	    	}
	    	
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		nu_f.add(new ArrayList<IloNumVar>());
	    	
	    		for (int l = 0; l < number_of_sc; l++)
	    		{
	    			nu_f.get(k).add(cplex.numVar(0, M, "nu_f_"+ k + l));	
	    		}
	    	}
	    
	        for (int j = 0; j < I.number_of_items; j++)
	        {
	        	x.add(cplex.boolVar("x_" + j));
	        	
	        	if (binary_indicator == 0) 
	        		y.add(cplex.numVar(0, 1, "y_" + j));
	        	else y.add(cplex.boolVar("y_" + j));
	        	
	        	prod_f.add(cplex.numVar(0, M, "prod_" + j));
	        }
	    	
	        //Follower's dual variables
	        for (int j = 0; j < I.number_of_items; j++)
	        {
	        	beta_f.add(cplex.numVar(0, M, "beta_f_"+ j));
	        }
	        
	        for (int k = 0; k < k_f; k++)
	    	{
	    		
	    		mu1_f.add(new ArrayList<IloNumVar>());
	    		mu2_f.add(new ArrayList<IloNumVar>());
	    	
	    		 for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			mu1_f.get(k).add(cplex.numVar(0, M, "mu1_f_"+ k + j));	
	    			mu2_f.get(k).add(cplex.numVar(0, M, "mu2_f_"+ k + j));
	    		}
	    	}
	        
	    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
	    	
	    	if ("Risk-averse".equals(type_of_risk_l))
	    	{   
	    		objective.addTerm(-epsilon_l, lambda_l);
	    		
	    		for (int k = 0; k < k_l; k++)
	        	{
	    			for (int i = 0; i < I.number_of_items; i++)
	    			 {
	    				 objective.addTerm(data_set_l.get(i).get(k)/(k_l + 0.), y.get(i));
	    		     }
	    			
	        		for (int l = 0; l < number_of_sc; l++)
	        		{
	        			objective.addTerm(-delta_l.get(k).get(l)/(k_l + 0.), nu_l.get(k).get(l));
	        		}
	    
	        	}
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_l))
	    	{
	    		objective.addTerm(-epsilon_l, lambda_l);
	    		
	    		for (int i = 0; i < I.number_of_items; i++)
	    			objective.addTerm(average_l.get(i), y.get(i));
	    		
	    		
	    		for (int k = 0; k < k_l; k++)
	    		{
	    			for (int l = 0; l < number_of_sc; l++)
		    		{
		    			objective.addTerm(-delta_l.get(k).get(l)/(k_l + 0.), nu_l.get(k).get(l));
		    		}
	    		
	    			//for (int i = 0; i < I.number_of_items; i++)
	       			// {
	       			//	 objective.addTerm(data_set_l.get(i).get(k)/(k_l + 0.), y.get(i));
	       		    // }
	    		}
	    			
	    		
	    	}
	    	
	    	cplex.addMaximize(objective);
	    	
	    	//Constraints
	    	//Feasibility
	    	IloLinearNumExpr str_1 = cplex.linearNumExpr(); 
	    	IloLinearNumExpr str_2 = cplex.linearNumExpr();
	    	
	    	for (int j = 0; j < I.number_of_items; j++)
	    	{
	    		str_1.addTerm(1., x.get(j));
	    		str_2.addTerm(1., y.get(j));
	    	}
	    	
	    	cplex.addLe(str_1, I.weight_l);
	    	cplex.addGe(str_2, I.weight_f);
	    	
	    	for (int j = 0; j < I.number_of_items; j++)
	    	{
	    		IloLinearNumExpr str = cplex.linearNumExpr(); 
	    		str.addTerm(1., y.get(j));
	    		str.addTerm(1., x.get(j));
	        	cplex.addLe(str, 1.);
	    	}
	        
	    	//
	     	for (int k = 0; k < k_l; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	    	
	    	    	str1.addTerm(-1., lambda_l);
	    	    	str2.addTerm(-1., lambda_l);
	    	    	
	    	    	str1.addTerm(1., y.get(j));
	    	    	str2.addTerm(-1., y.get(j));
	    	    	
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(support_constraints.get(l).get(j), nu_l.get(k).get(l));
	    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_l.get(k).get(l));
	     	        }
	    	    	
	    	    	cplex.addLe(str1, 0.);
	    	    	cplex.addLe(str2, 0.);
	    		}
	    	
	     	//
	    	for (int k = 0; k < k_f; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	    	
	    	    	str1.addTerm(-1., lambda_f);
	    	    	str2.addTerm(-1., lambda_f);
	    	    	
	    	    	str1.addTerm(-1., y.get(j));
	    	    	str2.addTerm(1., y.get(j));
	    	    	
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(k).get(l));
	    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(k).get(l));
	     	        }
	    	    	
	    	    	cplex.addLe(str1, 0.);
	    	    	cplex.addLe(str2, 0.);
	    		}
	   	    	    		
	    	//Dual constraints of the follower
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{
	    		//
		    	for (int i = 0; i < I.number_of_items; i++)
		    	 {	
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
					
		    		str1.addTerm(1., beta_f.get(i));
		    		str1.addTerm(-1., phi_f);
		    		
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			str1.addTerm(1., mu1_f.get(k).get(i));
		    			str1.addTerm(-1., mu2_f.get(k).get(i));
		    		}
		    		
		    		cplex.addGe(str1, -average_f.get(i));
		    	}
		    	
		    	//
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int l = 0; l < number_of_sc; l++)
		            {
		    			IloLinearNumExpr str1 = cplex.linearNumExpr();
			    	
			    		for (int i = 0; i < I.number_of_items; i++)
				    	 {
			    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
			    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
				    	 }
			    		cplex.addGe(str1, -delta_f.get(k).get(l)/(k_f + 0.));
			    	} 
		    	}

	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    	{   
	    		//
		    	for (int i = 0; i < I.number_of_items; i++)
		    	 {	
		    		IloLinearNumExpr str1 = cplex.linearNumExpr(average_f.get(i));
					
		    		str1.addTerm(1., beta_f.get(i));
		    		str1.addTerm(-1., phi_f);
		    		
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			str1.addTerm(1., mu1_f.get(k).get(i));
		    			str1.addTerm(-1., mu2_f.get(k).get(i));
		    		}
		        	
		    		cplex.addGe(str1, 0.);
		    	}
		    	
		    	//
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int l = 0; l < number_of_sc; l++)
		            {
		    			IloLinearNumExpr str1 = cplex.linearNumExpr(delta_f.get(k).get(l)/(k_f + 0.));
			    		
			    		for (int i = 0; i < I.number_of_items; i++)
				    	 {
			    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
			    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
				    	 }
			    		cplex.addGe(str1, 0.);
			    	}
		    	}
	    	}
	    	
	    	//
	    	IloLinearNumExpr str3 = cplex.linearNumExpr();
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		for (int i = 0; i < I.number_of_items; i++)
		    	 {
	    			str3.addTerm(1., mu1_f.get(k).get(i));
	    			str3.addTerm(1., mu2_f.get(k).get(i));
		    	 }
	    	}
	    	
	    	if ("Risk-averse".equals(type_of_risk_f))
	    		cplex.addLe(str3, epsilon_f);
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    		cplex.addLe(str3, epsilon_f);
	    	
	    	
	    	//Strong duality 
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	
	    		
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str1.addTerm(-1., prod_f.get(i));
		    		//str1.addTerm(-1., beta_f.get(i));
		    	}
		    	str1.addTerm(I.weight_f, phi_f);
		    	
		    	str2.addTerm(epsilon_f, lambda_f);
		    	
		    	
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int i = 0; i < I.number_of_items; i++)
	    			 {
	    				 str2.addTerm(data_set_f.get(i).get(k)/(k_f + 0.), y.get(i));
	    		     }
	    			
	    			
	    
	        		for (int l = 0; l < number_of_sc; l++)
	        		{
	        			str2.addTerm(delta_f.get(k).get(l)/(k_f + 0.), nu_f.get(k).get(l));
	        		}
		    	}
		    	
		    	cplex.addEq(str1, str2);
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str1.addTerm(-1., prod_f.get(i));
		    		//str1.addTerm(-1., beta_f.get(i));
		    	}
		        
		    	str1.addTerm(I.weight_f, phi_f);
		    	
		    	str2.addTerm(epsilon_f, lambda_f);
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str2.addTerm(average_f.get(i), y.get(i));
		    	}
		    	
		    	for (int k = 0; k < k_f; k++)
	    			for (int l = 0; l < number_of_sc; l++)
	    		{
	    			str2.addTerm(delta_f.get(k).get(l)/(k_f + 0.), nu_f.get(k).get(l));
	    		}
		    	
		    	cplex.addEq(str1, str2);
	    	}
	    	
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr(-M1);
	    		IloLinearNumExpr str4 = cplex.linearNumExpr();
	    		
	    		str1.addTerm(1., prod_f.get(i));
	    		str1.addTerm(-1., beta_f.get(i));
	        	
	    		str2.addTerm(1., prod_f.get(i));
	    		str2.addTerm(M1, x.get(i));
	        	
	    		str4.addTerm(1., prod_f.get(i));
	    		str4.addTerm(-1, beta_f.get(i));
	    		str4.addTerm(M1, x.get(i));
	        	
	    		cplex.addLe(str1, 0.);
	    		cplex.addLe(str2, 0.);
	    		cplex.addGe(str4, 0.);
	    		
	    	}
	    	
	    	
	    	ArrayList<ArrayList<Double>> sol_opt = new ArrayList<ArrayList<Double>>();
	    	ArrayList<Double> y_opt = new ArrayList<Double>();
	    	ArrayList<Double> x_opt = new ArrayList<Double>();
	    	
	    	
	    	if (cplex.solve())
	    	{   		
	    	    for (int j = 0; j < I.number_of_items; j++)
	    	    {
	    	    	x_opt.add(cplex.getValue(x.get(j)));
	    			y_opt.add(cplex.getValue(y.get(j)));
	    	    }  
	    	    
	    	    x_opt.add(cplex.getObjValue());
	    	    
	    	    sol_opt.add(x_opt);
	    	    sol_opt.add(y_opt);
	    	}
	    	else 
	        	System.out.println("No solution found");
	    	
	    	return sol_opt;
		}
	    catch (IloException Exc) { Exc.printStackTrace(); }
	    return new ArrayList<ArrayList<Double>>();
}*/
	

	public static ArrayList<Double> SolveLeadersProblemScenarios(IloCplex cplex, Instance I, int number_of_scenarios, double epsilon_l, double epsilon_f, double alpha_l, double alpha_f, int k_l, int k_f, ArrayList<ArrayList<Double>> data_set_l, Map<Integer, ArrayList<ArrayList<Double>>> data_sets_f,
			ArrayList<ArrayList<Double>> support_constraints, String type_of_risk_l, String type_of_risk_f,  ArrayList<ArrayList<Double>> delta_l,  Map<Integer, ArrayList<ArrayList<Double>>> deltas_f, ArrayList<Double> average_l, Map<Integer, ArrayList<Double>> averages_f, double M1)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
    	cplex.setParam(IloCplex.Param.TimeLimit, 3600.0);
        
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	IloNumVar z = cplex.numVar(-M, M, "z");
    	
    	Map<Integer, ArrayList<ArrayList<IloNumVar>>> nu_l = new HashMap<Integer, ArrayList<ArrayList<IloNumVar>>>();
    	Map<Integer, ArrayList<ArrayList<IloNumVar>>> nu_f = new HashMap<Integer, ArrayList<ArrayList<IloNumVar>>>();
    	
    	ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
    	Map<Integer, ArrayList<IloNumVar>> y = new HashMap<Integer, ArrayList<IloNumVar>>();
    	
    	Map<Integer, ArrayList<IloNumVar>> s_l = new HashMap<Integer, ArrayList<IloNumVar>>();
    	Map<Integer, ArrayList<IloNumVar>> s_f = new HashMap<Integer, ArrayList<IloNumVar>>();
        
    	Map<Integer, IloNumVar> lambda_l = new HashMap<Integer, IloNumVar>(); //cplex.numVar(0, M, "lam_l");
    	Map<Integer, IloNumVar> lambda_f = new HashMap<Integer, IloNumVar>();
    	
    	for (int r = 0; r < number_of_scenarios; r++)
    	{
    		lambda_l.put(r, cplex.numVar(0, M, "lam_l" + r));
    		lambda_f.put(r, cplex.numVar(0, M, "lam_f" + r));
    	}
    		
    	Map<Integer, IloNumVar> t_l = new HashMap<Integer, IloNumVar>(); //cplex.numVar(0, M, "lam_l");
    	Map<Integer, IloNumVar> t_f = new HashMap<Integer, IloNumVar>();
    	
    	for (int r = 0; r < number_of_scenarios; r++)
    	{
    		t_l.put(r, cplex.numVar(-M, M, "t_l" + r));
    		t_f.put(r, cplex.numVar(-M, M, "t_f" + r));
    	}
    		
    	Map<Integer, ArrayList<ArrayList<IloNumVar>>> mu1_f = new HashMap<Integer, ArrayList<ArrayList<IloNumVar>>>();
    	Map<Integer, ArrayList<ArrayList<IloNumVar>>> mu2_f = new HashMap<Integer, ArrayList<ArrayList<IloNumVar>>>();
    	
    	Map<Integer, ArrayList<IloNumVar>> beta1_f = new HashMap<Integer, ArrayList<IloNumVar>>();
    	Map<Integer, ArrayList<IloNumVar>> beta2_f = new HashMap<Integer, ArrayList<IloNumVar>>();
    	Map<Integer, ArrayList<IloNumVar>> gamma_f = new HashMap<Integer, ArrayList<IloNumVar>>();
    	Map<Integer, ArrayList<IloNumVar>> prod_f = new HashMap<Integer, ArrayList<IloNumVar>>();
    	
    	for (int j = 0; j < I.number_of_items; j++)
        	x.add(cplex.boolVar("x_" + j));
    	
    	for (int r = 0; r < number_of_scenarios; r++)
    	{
    		ArrayList<IloNumVar> s_lr = new ArrayList<IloNumVar>();
    		ArrayList<ArrayList<IloNumVar>> nu_lr = new ArrayList<ArrayList<IloNumVar>>();
    	
	    	for (int k = 0; k < k_l; k++)
	    	{
	    		
	    		s_lr.add(cplex.numVar(0, M, "s_l_"+ k + r));
	    		nu_lr.add(new ArrayList<IloNumVar>());
	    	
	    		for (int l = 0; l < number_of_sc; l++)
	    		{
	    			nu_lr.get(k).add(cplex.numVar(0, M, "nu_l_"+ k + l + r));	
	    		}
	    	}
	    	
	    	s_l.put(r, s_lr);
    		nu_l.put(r, nu_lr);
	    	
    		ArrayList<IloNumVar> s_fr = new ArrayList<IloNumVar>();
    		ArrayList<ArrayList<IloNumVar>> nu_fr = new ArrayList<ArrayList<IloNumVar>>();
    		
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		s_fr.add(cplex.numVar(0, M, "s_f_"+ k + r));
	    		nu_fr.add(new ArrayList<IloNumVar>());
	    	
	    		for (int l = 0; l < number_of_sc; l++)
	    		{
	    			nu_fr.get(k).add(cplex.numVar(0, M, "nu_f_"+ k + l + r));	
	    		}
	    	}
	    	
	    	s_f.put(r, s_fr);
    		nu_f.put(r, nu_fr);
	    
    		ArrayList<IloNumVar> y_r = new ArrayList<IloNumVar>();
        	ArrayList<IloNumVar> prod_fr = new ArrayList<IloNumVar>();
    		
	        for (int j = 0; j < I.number_of_items; j++)
	        {	
	        	y_r.add(cplex.numVar(0, 1, "y_" + j + r));
	        	prod_fr.add(cplex.numVar(0, M, "prod_" + j + r));
	        }
	        
	        y.put(r, y_r);
    		prod_f.put(r, prod_fr);
	    	
    		ArrayList<IloNumVar> beta1_fr = new ArrayList<IloNumVar>();
    		ArrayList<IloNumVar> beta2_fr = new ArrayList<IloNumVar>();
	        
    		//Follower's dual variables
    		  for (int j = 0; j < I.number_of_fc; j++)
    	      {
    			    beta1_fr.add(cplex.numVar(0, M, "beta1_f_"+ j + r));
    	      }
    	        
    	      for (int i = 0; i < I.number_of_items; i++)
    	      {
    	    	  beta2_fr.add(cplex.numVar(0, M, "beta2_f_"+ i + r));
    	      }
    		
	        
	        beta1_f.put(r, beta1_fr);
	        beta2_f.put(r, beta2_fr);
	        
	        ArrayList<IloNumVar> gamma_fr = new ArrayList<IloNumVar>();
	        ArrayList<ArrayList<IloNumVar>> mu1_fr = new ArrayList<ArrayList<IloNumVar>>();
	        ArrayList<ArrayList<IloNumVar>> mu2_fr = new ArrayList<ArrayList<IloNumVar>>();
    		
	        for (int k = 0; k < k_f; k++)
	    	{
	    		gamma_fr.add(cplex.numVar(0, M, "gamma_f_"+ k + r));
	    		
	    		mu1_fr.add(new ArrayList<IloNumVar>());
	    		mu2_fr.add(new ArrayList<IloNumVar>());
	    	
	    		 for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			mu1_fr.get(k).add(cplex.numVar(0, M, "mu1_f_"+ k + j + r));	
	    			mu2_fr.get(k).add(cplex.numVar(0, M, "mu2_f_"+ k + j + r));
	    		}
	    	}
	        
	        gamma_f.put(r, gamma_fr);
    		mu1_f.put(r, mu1_fr);
    		mu2_f.put(r, mu2_fr);
    	}
    	
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	objective.addTerm(1., z);
    	
    	cplex.addMinimize(objective);
    	
    	//Constraints
    	//Scenario-based constraints for z
    	for (int r = 0; r < number_of_scenarios; r++)
    	{ 
    		
    		IloLinearNumExpr str = cplex.linearNumExpr(); 
        	str.addTerm(-1., z);
    		
	    	if ("Risk-averse".equals(type_of_risk_l))
	    	{   
	    		str.addTerm(1., t_l.get(r));
	    		str.addTerm(epsilon_l/(1. - alpha_l), lambda_l.get(r));
	    		
	    		for (int k = 0; k < k_l; k++)
	        	{
	        	     str.addTerm(1./(k_l*(1. - alpha_l)), s_l.get(r).get(k));
	        	}
	    	}
    	
	    	if ("Risk-neutral".equals(type_of_risk_l))
	    	{
	    		str.addTerm(epsilon_l, lambda_l.get(r));
	    		
	    		for (int i = 0; i < I.number_of_items; i++)
	    			str.addTerm(average_l.get(i), y.get(r).get(i));
	    		
	    		for (int k = 0; k < k_l; k++)
	    			for (int l = 0; l < number_of_sc; l++)
	    		{
	    			str.addTerm(delta_l.get(k).get(l)/(k_l + 0.), nu_l.get(r).get(k).get(l));
	    		}
	    		
	    	}
	    	
	    	cplex.addLe(str, 0.);
    	}
    	
   
    	//Feasibility
    	//X
    	for (int j = 0; j < I.number_of_lc; j++)
    	{
    		IloLinearNumExpr str_1 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_1.addTerm(I.H.get(j).get(i), x.get(i));
        	}
    		 
    		cplex.addLe(str_1, I.h.get(j));
         	
        }
    	
    	//Y
    	for (int r = 0; r < number_of_scenarios; r++)
    	{
    		for (int j = 0; j < I.number_of_fc; j++)
        	{
        		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
        		
        		for (int i = 0; i < I.number_of_items; i++)
            	{
        			str_2.addTerm(I.F.get(j).get(i), y.get(r).get(i));
            	}
        		 
        		cplex.addLe(str_2, I.f.get(j));
             	
            }
	    	
	    	for (int j = 0; j < I.number_of_items; j++)
	    	{
	    		IloLinearNumExpr str = cplex.linearNumExpr(); 
	    		str.addTerm(1., y.get(r).get(j));
	    		str.addTerm(1., x.get(j));
	        	cplex.addLe(str, 1.);
	    	}
    	}
    
    	//
    	for (int r = 0; r < number_of_scenarios; r++)
    	{	
	     	for (int k = 0; k < k_l; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	    	
	    	    	str1.addTerm(-1., lambda_l.get(r));
	    	    	str2.addTerm(-1., lambda_l.get(r));
	    	    	
	    	    	str1.addTerm(-1., y.get(r).get(j));
	    	    	str2.addTerm(1., y.get(r).get(j));
	    	    	
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(support_constraints.get(l).get(j), nu_l.get(r).get(k).get(l));
	    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_l.get(r).get(k).get(l));
	     	        }
	    	    	
	    	    	cplex.addLe(str1, 0.);
	    	    	cplex.addLe(str2, 0.);
	    		}
    	
	     	//
	    	for (int k = 0; k < k_f; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	    	
	    	    	str1.addTerm(-1., lambda_f.get(r));
	    	    	str2.addTerm(-1., lambda_f.get(r));
	    	    	
	    	    	str1.addTerm(1., y.get(r).get(j));
	    	    	str2.addTerm(-1., y.get(r).get(j));
	    	    	
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(r).get(k).get(l));
	    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(r).get(k).get(l));
	     	        }
	    	    	
	    	    	cplex.addLe(str1, 0.);
	    	    	cplex.addLe(str2, 0.);
	    		}
    	
	    	//
	    	if ("Risk-averse".equals(type_of_risk_l))
	    	{   
	    		for (int k = 0; k < k_l; k++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    			
	    			 for (int i = 0; i < I.number_of_items; i++)
	    			 {
	    				 str1.addTerm(data_set_l.get(i).get(k), y.get(r).get(i));
	    		     }
	    			
	    			str1.addTerm(-1., t_l.get(r));
	    			
	    
	        		for (int l = 0; l < number_of_sc; l++)
	        		{
	        			str1.addTerm(delta_l.get(k).get(l), nu_l.get(r).get(k).get(l));
	        		}
	    			
	    			str1.addTerm(-1., s_l.get(r).get(k));
	    			cplex.addLe(str1, 0.);
	    			
	                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
	    			
	    			str2.addTerm(1., s_l.get(k));
	    			cplex.addGe(str2, 0);*/
	    		}
	    		
	    	}
	    	
	    	//
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{   
	    		for (int k = 0; k < k_f; k++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    			
	    			 for (int i = 0; i < I.number_of_items; i++)
	    			 {
	    				 str1.addTerm(-data_sets_f.get(r).get(i).get(k), y.get(r).get(i));
	    		     }
	    			
	    			str1.addTerm(1., t_f.get(r));
	    			
	    
	        		for (int l = 0; l < number_of_sc; l++)
	        		{
	        			str1.addTerm(deltas_f.get(r).get(k).get(l), nu_f.get(r).get(k).get(l));
	        		}
	    			
	    			str1.addTerm(-1., s_f.get(r).get(k));
	    			cplex.addLe(str1, 0.);
	    			
	                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
	    			
	    			str2.addTerm(1., s_f.get(k));
	    			cplex.addGe(str2, 0);*/
	    		}
	    		
	    	}
	    	
	    	//Dual constraints of the follower
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{
	    		//
		    	for (int i = 0; i < I.number_of_items; i++)
		    	 {	
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
					
		    		for (int j = 0; j < I.number_of_fc; j++)
		           	{
		       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(r).get(j));
		           	}
		       		
		       		str1.addTerm(1., beta2_f.get(r).get(i));
		    		
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			str1.addTerm(-data_sets_f.get(r).get(i).get(k), gamma_f.get(r).get(k));
		    			str1.addTerm(-1., mu1_f.get(r).get(k).get(i));
		    			str1.addTerm(1., mu2_f.get(r).get(k).get(i));
		    		}
		        	
		    		cplex.addGe(str1, 0.);
		    	}
		    	
		    	//
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int l = 0; l < number_of_sc; l++)
		            {
		    			IloLinearNumExpr str1 = cplex.linearNumExpr();
			    		
		    			str1.addTerm(deltas_f.get(r).get(k).get(l), gamma_f.get(r).get(k));
			    		
			    		for (int i = 0; i < I.number_of_items; i++)
				    	 {
			    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(r).get(k).get(i));
			    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(r).get(k).get(i));
				    	 }
			    		cplex.addGe(str1, 0.);
			    	} 
		    	}
		    	
		    	//
		    	IloLinearNumExpr str = cplex.linearNumExpr();
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		
		    		str.addTerm(1., gamma_f.get(r).get(k));
		    		str1.addTerm(1., gamma_f.get(r).get(k));
		    		
		    		cplex.addLe(str1, 1./((1. - alpha_f) * k_f));
		    	}
		    	
		    	cplex.addEq(str, 1.);
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    	{   
	    		//
		    	for (int i = 0; i < I.number_of_items; i++)
		    	 {	
		    		IloLinearNumExpr str1 = cplex.linearNumExpr(-averages_f.get(r).get(i));
					
		    		for (int j = 0; j < I.number_of_fc; j++)
		           	{
		       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(r).get(j));
		           	}
		       		
		       		str1.addTerm(1., beta2_f.get(r).get(i));
		    		
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			str1.addTerm(-1., mu1_f.get(r).get(k).get(i));
		    			str1.addTerm(1., mu2_f.get(r).get(k).get(i));
		    		}
		        	
		    		cplex.addGe(str1, 0.);
		    	}
		    	
		    	//
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int l = 0; l < number_of_sc; l++)
		            {
		    			IloLinearNumExpr str1 = cplex.linearNumExpr(deltas_f.get(r).get(k).get(l)/(k_f + 0.));
			    		
			    		for (int i = 0; i < I.number_of_items; i++)
				    	 {
			    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(r).get(k).get(i));
			    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(r).get(k).get(i));
				    	 }
			    		cplex.addGe(str1, 0.);
			    	}
		    	}
	    	}
	    	
	    	//
	    	IloLinearNumExpr str3 = cplex.linearNumExpr();
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		for (int i = 0; i < I.number_of_items; i++)
		    	 {
	    			str3.addTerm(1., mu1_f.get(r).get(k).get(i));
	    			str3.addTerm(1., mu2_f.get(r).get(k).get(i));
		    	 }
	    	}
	    	
	    	if ("Risk-averse".equals(type_of_risk_f))
	    		cplex.addLe(str3, epsilon_f/(1. - alpha_f));
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    		cplex.addLe(str3, epsilon_f);
	    	
	    	
	    	//Strong duality 
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	
	    		
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str1.addTerm(1., prod_f.get(r).get(i));
		    	}

		    	for (int j = 0; j < I.number_of_fc; j++)
		    	{
		    		str1.addTerm(I.f.get(j), beta1_f.get(r).get(j));
		    	}
		    	
		    	str2.addTerm(1., t_f.get(r));
		    	str2.addTerm(-epsilon_f/(1. - alpha_f), lambda_f.get(r));
		    	
		    	
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		str2.addTerm(-1./(k_f * (1. - alpha_f)), s_f.get(r).get(k));
		    	}
		    	
		    	cplex.addEq(str1, str2);
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str1.addTerm(1., prod_f.get(r).get(i));
		    		//str1.addTerm(-1., beta_f.get(i));
		    	}
		        
		    	for (int j = 0; j < I.number_of_fc; j++)
		    	{
		    		str1.addTerm(I.f.get(j), beta1_f.get(r).get(j));
		    	}
		    	
		    	str2.addTerm(-epsilon_f, lambda_f.get(r));
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str2.addTerm(averages_f.get(r).get(i), y.get(r).get(i));
		    	}
		    	
		    	for (int k = 0; k < k_f; k++)
	    			for (int l = 0; l < number_of_sc; l++)
	    		{
	    			str2.addTerm(-deltas_f.get(r).get(k).get(l)/(k_f + 0.), nu_f.get(r).get(k).get(l));
	    		}
		    	
		    	cplex.addEq(str1, str2);
	    	}
	    	
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr(-M1);
	    		IloLinearNumExpr str4 = cplex.linearNumExpr();
	    		
	    		str1.addTerm(1., prod_f.get(r).get(i));
	    		str1.addTerm(-1., beta2_f.get(r).get(i));
	        	
	    		str2.addTerm(1., prod_f.get(r).get(i));
	    		str2.addTerm(M1, x.get(i));
	        	
	    		str4.addTerm(1., prod_f.get(r).get(i));
	    		str4.addTerm(-1, beta2_f.get(r).get(i));
	    		str4.addTerm(M1, x.get(i));
	        	
	    		cplex.addLe(str1, 0.);
	    		cplex.addLe(str2, 0.);
	    		cplex.addGe(str4, 0.);
	    		
	    	}
    	}
    	
    	//ArrayList<ArrayList<Double>> sol_opt = new ArrayList<ArrayList<Double>>();
    	Map<Integer, ArrayList<Double>> y_opt = new HashMap<Integer, ArrayList<Double>>();
    	ArrayList<Double> x_opt = new ArrayList<Double>();
    	
    	
    	if (cplex.solve())
    	{   	
    	    for (int j = 0; j < I.number_of_items; j++)
    	    {
    	    	x_opt.add(cplex.getValue(x.get(j)));
    	    }  
    	    
    	    //System.out.println(x_opt);
    		
    	    x_opt.add(cplex.getObjValue());
    	    
    	    
    	    for (int r = 0; r < number_of_scenarios; r++)
            {
            	ArrayList<Double> y_r = new ArrayList<Double>();
            	
            	for (int j = 0; j < I.number_of_items; j++)
            		y_r.add(cplex.getValue(y.get(r).get(j)));
            	
            	y_opt.put(r, y_r);
            }
            
            //System.out.println(y_opt);
    	}
    	else 
        	System.out.println("No solution found");
    	
    	if (cplex.getStatus() != IloCplex.Status.Optimal)
    	{
    		 return new ArrayList<Double>();
    	}
    	
    	return x_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}

	public static ArrayList<Double> SolveLeadersProblemScenarios2(IloCplex cplex, Instance I, int number_of_scenarios, double epsilon_l, Map<Integer, Double> epsilons_f, double alpha_l, Map<Integer, Double> alphas_f, int k_l, int k_f, ArrayList<ArrayList<Double>> data_set_l, Map<Integer, ArrayList<ArrayList<Double>>> data_sets_f,
			ArrayList<ArrayList<Double>> support_constraints, String type_of_risk_l, String type_of_risk_f,  ArrayList<ArrayList<Double>> delta_l,  Map<Integer, ArrayList<ArrayList<Double>>> deltas_f, ArrayList<Double> average_l, Map<Integer, ArrayList<Double>> averages_f, double M1)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
    	cplex.setParam(IloCplex.Param.TimeLimit, 3600.0);
        
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	IloNumVar z = cplex.numVar(-M, M, "z");
    	
    	Map<Integer, ArrayList<ArrayList<IloNumVar>>> nu_l = new HashMap<Integer, ArrayList<ArrayList<IloNumVar>>>();
    	Map<Integer, ArrayList<ArrayList<IloNumVar>>> nu_f = new HashMap<Integer, ArrayList<ArrayList<IloNumVar>>>();
    	
    	ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
    	Map<Integer, ArrayList<IloNumVar>> y = new HashMap<Integer, ArrayList<IloNumVar>>();
    	
    	Map<Integer, ArrayList<IloNumVar>> s_l = new HashMap<Integer, ArrayList<IloNumVar>>();
    	Map<Integer, ArrayList<IloNumVar>> s_f = new HashMap<Integer, ArrayList<IloNumVar>>();
        
    	Map<Integer, IloNumVar> lambda_l = new HashMap<Integer, IloNumVar>(); //cplex.numVar(0, M, "lam_l");
    	Map<Integer, IloNumVar> lambda_f = new HashMap<Integer, IloNumVar>();
    	
    	for (int r = 0; r < number_of_scenarios; r++)
    	{
    		lambda_l.put(r, cplex.numVar(0, M, "lam_l" + r));
    		lambda_f.put(r, cplex.numVar(0, M, "lam_f" + r));
    	}
    		
    	Map<Integer, IloNumVar> t_l = new HashMap<Integer, IloNumVar>(); //cplex.numVar(0, M, "lam_l");
    	Map<Integer, IloNumVar> t_f = new HashMap<Integer, IloNumVar>();
    	
    	for (int r = 0; r < number_of_scenarios; r++)
    	{
    		t_l.put(r, cplex.numVar(-M, M, "t_l" + r));
    		t_f.put(r, cplex.numVar(-M, M, "t_f" + r));
    	}
    		
    	Map<Integer, ArrayList<ArrayList<IloNumVar>>> mu1_f = new HashMap<Integer, ArrayList<ArrayList<IloNumVar>>>();
    	Map<Integer, ArrayList<ArrayList<IloNumVar>>> mu2_f = new HashMap<Integer, ArrayList<ArrayList<IloNumVar>>>();
    	
    	Map<Integer, ArrayList<IloNumVar>> beta1_f = new HashMap<Integer, ArrayList<IloNumVar>>();
    	Map<Integer, ArrayList<IloNumVar>> beta2_f = new HashMap<Integer, ArrayList<IloNumVar>>();
    	Map<Integer, ArrayList<IloNumVar>> gamma_f = new HashMap<Integer, ArrayList<IloNumVar>>();
    	Map<Integer, ArrayList<IloNumVar>> prod_f = new HashMap<Integer, ArrayList<IloNumVar>>();
    	
    	for (int j = 0; j < I.number_of_items; j++)
        	x.add(cplex.boolVar("x_" + j));
    	
    	for (int r = 0; r < number_of_scenarios; r++)
    	{
    		ArrayList<IloNumVar> s_lr = new ArrayList<IloNumVar>();
    		ArrayList<ArrayList<IloNumVar>> nu_lr = new ArrayList<ArrayList<IloNumVar>>();
    	
	    	for (int k = 0; k < k_l; k++)
	    	{
	    		
	    		s_lr.add(cplex.numVar(0, M, "s_l_"+ k + r));
	    		nu_lr.add(new ArrayList<IloNumVar>());
	    	
	    		for (int l = 0; l < number_of_sc; l++)
	    		{
	    			nu_lr.get(k).add(cplex.numVar(0, M, "nu_l_"+ k + l + r));	
	    		}
	    	}
	    	
	    	s_l.put(r, s_lr);
    		nu_l.put(r, nu_lr);
	    	
    		ArrayList<IloNumVar> s_fr = new ArrayList<IloNumVar>();
    		ArrayList<ArrayList<IloNumVar>> nu_fr = new ArrayList<ArrayList<IloNumVar>>();
    		
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		s_fr.add(cplex.numVar(0, M, "s_f_"+ k + r));
	    		nu_fr.add(new ArrayList<IloNumVar>());
	    	
	    		for (int l = 0; l < number_of_sc; l++)
	    		{
	    			nu_fr.get(k).add(cplex.numVar(0, M, "nu_f_"+ k + l + r));	
	    		}
	    	}
	    	
	    	s_f.put(r, s_fr);
    		nu_f.put(r, nu_fr);
	    
    		ArrayList<IloNumVar> y_r = new ArrayList<IloNumVar>();
        	ArrayList<IloNumVar> prod_fr = new ArrayList<IloNumVar>();
    		
	        for (int j = 0; j < I.number_of_items; j++)
	        {	
	        	y_r.add(cplex.numVar(0, 1, "y_" + j + r));
	        	prod_fr.add(cplex.numVar(0, M, "prod_" + j + r));
	        }
	        
	        y.put(r, y_r);
    		prod_f.put(r, prod_fr);
	    	
    		ArrayList<IloNumVar> beta1_fr = new ArrayList<IloNumVar>();
    		ArrayList<IloNumVar> beta2_fr = new ArrayList<IloNumVar>();
	        
    		//Follower's dual variables
    		  for (int j = 0; j < I.number_of_fc; j++)
    	      {
    			    beta1_fr.add(cplex.numVar(0, M, "beta1_f_"+ j + r));
    	      }
    	        
    	      for (int i = 0; i < I.number_of_items; i++)
    	      {
    	    	  beta2_fr.add(cplex.numVar(0, M, "beta2_f_"+ i + r));
    	      }
    		
	        
	        beta1_f.put(r, beta1_fr);
	        beta2_f.put(r, beta2_fr);
	        
	        ArrayList<IloNumVar> gamma_fr = new ArrayList<IloNumVar>();
	        ArrayList<ArrayList<IloNumVar>> mu1_fr = new ArrayList<ArrayList<IloNumVar>>();
	        ArrayList<ArrayList<IloNumVar>> mu2_fr = new ArrayList<ArrayList<IloNumVar>>();
    		
	        for (int k = 0; k < k_f; k++)
	    	{
	    		gamma_fr.add(cplex.numVar(0, M, "gamma_f_"+ k + r));
	    		
	    		mu1_fr.add(new ArrayList<IloNumVar>());
	    		mu2_fr.add(new ArrayList<IloNumVar>());
	    	
	    		 for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			mu1_fr.get(k).add(cplex.numVar(0, M, "mu1_f_"+ k + j + r));	
	    			mu2_fr.get(k).add(cplex.numVar(0, M, "mu2_f_"+ k + j + r));
	    		}
	    	}
	        
	        gamma_f.put(r, gamma_fr);
    		mu1_f.put(r, mu1_fr);
    		mu2_f.put(r, mu2_fr);
    	}
    	
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	objective.addTerm(1., z);
    	
    	cplex.addMinimize(objective);
    	
    	//Constraints
    	//Scenario-based constraints for z
    	for (int r = 0; r < number_of_scenarios; r++)
    	{ 
    		
    		IloLinearNumExpr str = cplex.linearNumExpr(); 
        	str.addTerm(-1., z);
    		
	    	if ("Risk-averse".equals(type_of_risk_l))
	    	{   
	    		str.addTerm(1., t_l.get(r));
	    		str.addTerm(epsilon_l/(1. - alpha_l), lambda_l.get(r));
	    		
	    		for (int k = 0; k < k_l; k++)
	        	{
	        	     str.addTerm(1./(k_l*(1. - alpha_l)), s_l.get(r).get(k));
	        	}
	    	}
    	
	    	if ("Risk-neutral".equals(type_of_risk_l))
	    	{
	    		str.addTerm(epsilon_l, lambda_l.get(r));
	    		
	    		for (int i = 0; i < I.number_of_items; i++)
	    			str.addTerm(average_l.get(i), y.get(r).get(i));
	    		
	    		for (int k = 0; k < k_l; k++)
	    			for (int l = 0; l < number_of_sc; l++)
	    		{
	    			str.addTerm(delta_l.get(k).get(l)/(k_l + 0.), nu_l.get(r).get(k).get(l));
	    		}
	    		
	    	}
	    	
	    	cplex.addLe(str, 0.);
    	}
    	
   
    	//Feasibility
    	//X
    	for (int j = 0; j < I.number_of_lc; j++)
    	{
    		IloLinearNumExpr str_1 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_1.addTerm(I.H.get(j).get(i), x.get(i));
        	}
    		 
    		cplex.addLe(str_1, I.h.get(j));
         	
        }
    	
    	//Y
    	for (int r = 0; r < number_of_scenarios; r++)
    	{
    		for (int j = 0; j < I.number_of_fc; j++)
        	{
        		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
        		
        		for (int i = 0; i < I.number_of_items; i++)
            	{
        			str_2.addTerm(I.F.get(j).get(i), y.get(r).get(i));
            	}
        		 
        		cplex.addLe(str_2, I.f.get(j));
             	
            }
	    	
	    	for (int j = 0; j < I.number_of_items; j++)
	    	{
	    		IloLinearNumExpr str = cplex.linearNumExpr(); 
	    		str.addTerm(1., y.get(r).get(j));
	    		str.addTerm(1., x.get(j));
	        	cplex.addLe(str, 1.);
	    	}
    	}
    
    	//
    	for (int r = 0; r < number_of_scenarios; r++)
    	{	
	     	for (int k = 0; k < k_l; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	    	
	    	    	str1.addTerm(-1., lambda_l.get(r));
	    	    	str2.addTerm(-1., lambda_l.get(r));
	    	    	
	    	    	str1.addTerm(-1., y.get(r).get(j));
	    	    	str2.addTerm(1., y.get(r).get(j));
	    	    	
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(support_constraints.get(l).get(j), nu_l.get(r).get(k).get(l));
	    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_l.get(r).get(k).get(l));
	     	        }
	    	    	
	    	    	cplex.addLe(str1, 0.);
	    	    	cplex.addLe(str2, 0.);
	    		}
    	
	     	//
	    	for (int k = 0; k < k_f; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	    	
	    	    	str1.addTerm(-1., lambda_f.get(r));
	    	    	str2.addTerm(-1., lambda_f.get(r));
	    	    	
	    	    	str1.addTerm(1., y.get(r).get(j));
	    	    	str2.addTerm(-1., y.get(r).get(j));
	    	    	
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(r).get(k).get(l));
	    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(r).get(k).get(l));
	     	        }
	    	    	
	    	    	cplex.addLe(str1, 0.);
	    	    	cplex.addLe(str2, 0.);
	    		}
    	
	    	//
	    	if ("Risk-averse".equals(type_of_risk_l))
	    	{   
	    		for (int k = 0; k < k_l; k++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    			
	    			 for (int i = 0; i < I.number_of_items; i++)
	    			 {
	    				 str1.addTerm(data_set_l.get(i).get(k), y.get(r).get(i));
	    		     }
	    			
	    			str1.addTerm(-1., t_l.get(r));
	    			
	    
	        		for (int l = 0; l < number_of_sc; l++)
	        		{
	        			str1.addTerm(delta_l.get(k).get(l), nu_l.get(r).get(k).get(l));
	        		}
	    			
	    			str1.addTerm(-1., s_l.get(r).get(k));
	    			cplex.addLe(str1, 0.);
	    			
	                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
	    			
	    			str2.addTerm(1., s_l.get(k));
	    			cplex.addGe(str2, 0);*/
	    		}
	    		
	    	}
	    	
	    	//
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{   
	    		for (int k = 0; k < k_f; k++)
	    		{
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
	    			
	    			 for (int i = 0; i < I.number_of_items; i++)
	    			 {
	    				 str1.addTerm(-data_sets_f.get(r).get(i).get(k), y.get(r).get(i));
	    		     }
	    			
	    			str1.addTerm(1., t_f.get(r));
	    			
	    
	        		for (int l = 0; l < number_of_sc; l++)
	        		{
	        			str1.addTerm(deltas_f.get(r).get(k).get(l), nu_f.get(r).get(k).get(l));
	        		}
	    			
	    			str1.addTerm(-1., s_f.get(r).get(k));
	    			cplex.addLe(str1, 0.);
	    			
	                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
	    			
	    			str2.addTerm(1., s_f.get(k));
	    			cplex.addGe(str2, 0);*/
	    		}
	    		
	    	}
	    	
	    	//Dual constraints of the follower
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{
	    		//
		    	for (int i = 0; i < I.number_of_items; i++)
		    	 {	
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
					
		    		for (int j = 0; j < I.number_of_fc; j++)
		           	{
		       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(r).get(j));
		           	}
		       		
		       		str1.addTerm(1., beta2_f.get(r).get(i));
		    		
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			str1.addTerm(-data_sets_f.get(r).get(i).get(k), gamma_f.get(r).get(k));
		    			str1.addTerm(-1., mu1_f.get(r).get(k).get(i));
		    			str1.addTerm(1., mu2_f.get(r).get(k).get(i));
		    		}
		        	
		    		cplex.addGe(str1, 0.);
		    	}
		    	
		    	//
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int l = 0; l < number_of_sc; l++)
		            {
		    			IloLinearNumExpr str1 = cplex.linearNumExpr();
			    		
		    			str1.addTerm(deltas_f.get(r).get(k).get(l), gamma_f.get(r).get(k));
			    		
			    		for (int i = 0; i < I.number_of_items; i++)
				    	 {
			    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(r).get(k).get(i));
			    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(r).get(k).get(i));
				    	 }
			    		cplex.addGe(str1, 0.);
			    	} 
		    	}
		    	
		    	//
		    	IloLinearNumExpr str = cplex.linearNumExpr();
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		
		    		str.addTerm(1., gamma_f.get(r).get(k));
		    		str1.addTerm(1., gamma_f.get(r).get(k));
		    		
		    		cplex.addLe(str1, 1./((1. - alphas_f.get(r)) * k_f));
		    	}
		    	
		    	cplex.addEq(str, 1.);
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    	{   
	    		//
		    	for (int i = 0; i < I.number_of_items; i++)
		    	 {	
		    		IloLinearNumExpr str1 = cplex.linearNumExpr(-averages_f.get(r).get(i));
					
		    		for (int j = 0; j < I.number_of_fc; j++)
		           	{
		       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(r).get(j));
		           	}
		       		
		       		str1.addTerm(1., beta2_f.get(r).get(i));
		    		
		    		for (int k = 0; k < k_f; k++)
		    		{
		    			str1.addTerm(-1., mu1_f.get(r).get(k).get(i));
		    			str1.addTerm(1., mu2_f.get(r).get(k).get(i));
		    		}
		        	
		    		cplex.addGe(str1, 0.);
		    	}
		    	
		    	//
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		for (int l = 0; l < number_of_sc; l++)
		            {
		    			IloLinearNumExpr str1 = cplex.linearNumExpr(deltas_f.get(r).get(k).get(l)/(k_f + 0.));
			    		
			    		for (int i = 0; i < I.number_of_items; i++)
				    	 {
			    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(r).get(k).get(i));
			    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(r).get(k).get(i));
				    	 }
			    		cplex.addGe(str1, 0.);
			    	}
		    	}
	    	}
	    	
	    	//
	    	IloLinearNumExpr str3 = cplex.linearNumExpr();
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		for (int i = 0; i < I.number_of_items; i++)
		    	 {
	    			str3.addTerm(1., mu1_f.get(r).get(k).get(i));
	    			str3.addTerm(1., mu2_f.get(r).get(k).get(i));
		    	 }
	    	}
	    	
	    	if ("Risk-averse".equals(type_of_risk_f))
	    		cplex.addLe(str3, epsilons_f.get(r)/(1. - alphas_f.get(r)));
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    		cplex.addLe(str3, epsilons_f.get(r));
	    	
	    	
	    	//Strong duality 
	    	if ("Risk-averse".equals(type_of_risk_f))
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	
	    		
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str1.addTerm(1., prod_f.get(r).get(i));
		    	}

		    	for (int j = 0; j < I.number_of_fc; j++)
		    	{
		    		str1.addTerm(I.f.get(j), beta1_f.get(r).get(j));
		    	}
		    	
		    	str2.addTerm(1., t_f.get(r));
		    	str2.addTerm(-epsilons_f.get(r)/(1. - alphas_f.get(r)), lambda_f.get(r));
		    	
		    	
		    	for (int k = 0; k < k_f; k++)
		    	{
		    		str2.addTerm(-1./(k_f * (1. - alphas_f.get(r))), s_f.get(r).get(k));
		    	}
		    	
		    	cplex.addEq(str1, str2);
	    	}
	    	
	    	if ("Risk-neutral".equals(type_of_risk_f))
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr();
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str1.addTerm(1., prod_f.get(r).get(i));
		    		//str1.addTerm(-1., beta_f.get(i));
		    	}
		        
		    	for (int j = 0; j < I.number_of_fc; j++)
		    	{
		    		str1.addTerm(I.f.get(j), beta1_f.get(r).get(j));
		    	}
		    	
		    	str2.addTerm(-epsilons_f.get(r), lambda_f.get(r));
		    	
		    	for (int i = 0; i < I.number_of_items; i++)
		    	{
		    		str2.addTerm(averages_f.get(r).get(i), y.get(r).get(i));
		    	}
		    	
		    	for (int k = 0; k < k_f; k++)
	    			for (int l = 0; l < number_of_sc; l++)
	    		{
	    			str2.addTerm(-deltas_f.get(r).get(k).get(l)/(k_f + 0.), nu_f.get(r).get(k).get(l));
	    		}
		    	
		    	cplex.addEq(str1, str2);
	    	}
	    	
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		IloLinearNumExpr str2 = cplex.linearNumExpr(-M1);
	    		IloLinearNumExpr str4 = cplex.linearNumExpr();
	    		
	    		str1.addTerm(1., prod_f.get(r).get(i));
	    		str1.addTerm(-1., beta2_f.get(r).get(i));
	        	
	    		str2.addTerm(1., prod_f.get(r).get(i));
	    		str2.addTerm(M1, x.get(i));
	        	
	    		str4.addTerm(1., prod_f.get(r).get(i));
	    		str4.addTerm(-1, beta2_f.get(r).get(i));
	    		str4.addTerm(M1, x.get(i));
	        	
	    		cplex.addLe(str1, 0.);
	    		cplex.addLe(str2, 0.);
	    		cplex.addGe(str4, 0.);
	    		
	    	}
    	}
    	
    	//ArrayList<ArrayList<Double>> sol_opt = new ArrayList<ArrayList<Double>>();
    	Map<Integer, ArrayList<Double>> y_opt = new HashMap<Integer, ArrayList<Double>>();
    	ArrayList<Double> x_opt = new ArrayList<Double>();
    	
    	
    	if (cplex.solve())
    	{   	
    	    for (int j = 0; j < I.number_of_items; j++)
    	    {
    	    	x_opt.add(cplex.getValue(x.get(j)));
    	    }  
    	    
    	    //System.out.println(x_opt);
    		
    	    x_opt.add(cplex.getObjValue());
    	    
    	    
    	    for (int r = 0; r < number_of_scenarios; r++)
            {
            	ArrayList<Double> y_r = new ArrayList<Double>();
            	
            	for (int j = 0; j < I.number_of_items; j++)
            		y_r.add(cplex.getValue(y.get(r).get(j)));
            	
            	y_opt.put(r, y_r);
            }
            
            //System.out.println(y_opt);
    	}
    	else 
        	System.out.println("No solution found");
    	
    	if (cplex.getStatus() != IloCplex.Status.Optimal)
    	{
    		 return new ArrayList<Double>();
    	}
    	
    	return x_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}
	
	
	public static ArrayList<Double> SolveMasterProblemAmbiguityFree(IloCplex cplex, IloCplex cplex2, Instance I, double lower_bound_initial, double upper_bound_initial, double alpha_l, int k_l, ArrayList<ArrayList<Double>> data_set_l, double M1, double time_0)
    {
	try 
    {
		cplex2.clearModel();
    	cplex2.setOut(null);
    	cplex2.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex2.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex2.setParam(IloCplex.Param.WorkMem, 2048);
    	cplex2.setParam(IloCplex.Param.TimeLimit, 3600.0);
        
        double lower_bound = lower_bound_initial, upper_bound = upper_bound_initial;
        
        double M = Double.MAX_VALUE;
    	
    	ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
    	//ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
    	IloNumVar z = cplex2.numVar(-M, M, "z");
    
    	ArrayList<ArrayList<Double>> c_gamma = new ArrayList<ArrayList<Double>>();
    	ArrayList<ArrayList<IloNumVar>> beta1_gamma = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> beta2_gamma = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> prod_gamma = new ArrayList<ArrayList<IloNumVar>>();
    	
    	for (int j = 0; j < I.number_of_items; j++)
        {
        	x.add(cplex2.boolVar("x_" + j));
        }
     	
        IloLinearNumExpr objective = cplex2.linearNumExpr(); 
    	
    	objective.addTerm(1., z);
        	
    	cplex2.addMinimize(objective);
    	
    	//X
    	for (int j = 0; j < I.number_of_lc; j++)
    	{
    		IloLinearNumExpr str_1 = cplex2.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_1.addTerm(I.H.get(j).get(i), x.get(i));
        	}
    		 
    		cplex2.addLe(str_1, I.h.get(j));
         	
        }
    	
    	ArrayList<Double> x_opt = new ArrayList<Double>();    	
    	
    	//Compute initial lower bound at x = 0, initialize Gamma
    	ArrayList<Double> x_curr = InitializeDoubleVector(I.number_of_items, 0.);
    	ArrayList<ArrayList<Double>> Gamma = new ArrayList<ArrayList<Double>>();
    	
    	int curr = 0;
    	while ((Math.abs(upper_bound - lower_bound) > Math.pow(10, -3)))
    	{
    		ArrayList<Double> sol_f = Instance.SolveFollowersProblemPessimisticAmbiguityFree(cplex, I, x_curr, alpha_l, k_l, data_set_l);
        	
    		if (sol_f.size() == 0)
    		{
    			return new ArrayList<Double>();
       	    }
    		
    		if (sol_f.get(k_l) < upper_bound)
    		{
    			upper_bound = sol_f.get(k_l);
    	        x_opt = CopyDoubleVector1(x_curr, I.number_of_items);		
    		}
    		
    		ArrayList<Double> new_gamma = CopyDoubleVector1(sol_f, k_l);
        	Gamma.add(new_gamma);
    		
        	//Compute c_gamma
        	c_gamma.add(new ArrayList<Double>());
            
        	for (int i = 0; i < I.number_of_items; i++)
        	{
        		double c_current = 0;
    	    	
        		for (int k = 0; k < k_l; k++)
    	    	{
    	    		c_current += 1./((int)Math.round(k_l * (1. - alpha_l)))*data_set_l.get(i).get(k)*Gamma.get(curr).get(k);
    	    	}
        		
        		c_gamma.get(curr).add(c_current);
        	}
        	
        	beta1_gamma.add(new ArrayList<IloNumVar>());
        	beta2_gamma.add(new ArrayList<IloNumVar>());
        	prod_gamma.add(new ArrayList<IloNumVar>());
         	
        	for (int j = 0; j < I.number_of_fc; j++)
            {
                beta1_gamma.get(curr).add(cplex2.numVar(0, M, "beta1_"+ curr + j));
            }
             
        	for (int j = 0; j < I.number_of_items; j++)
            {
     			beta2_gamma.get(curr).add(cplex2.numVar(0, M, "beta2_"+ curr + j));	
     			prod_gamma.get(curr).add(cplex2.numVar(0, M, "prod_"+ curr + j));
     		}
        	 
        	//Dual constraints
        	for (int i = 0; i < I.number_of_items; i++)
    	   	 {	
    	   		IloLinearNumExpr str1 = cplex2.linearNumExpr(-c_gamma.get(curr).get(i));
    					
	    		for (int j = 0; j < I.number_of_fc; j++)
	           	{
	       			str1.addTerm(I.F.get(j).get(i), beta1_gamma.get(curr).get(j));
	           	}
	       		
	       		str1.addTerm(1., beta2_gamma.get(curr).get(i));

    	   		cplex2.addGe(str1, 0.);
    	   	}
        	
        	IloLinearNumExpr str1 = cplex2.linearNumExpr();
    		str1.addTerm(1., z);
    		
        	for (int i = 0; i < I.number_of_items; i++)
        	{
        		str1.addTerm(-1., prod_gamma.get(curr).get(i));
        	}
        	
        	for (int j = 0; j < I.number_of_fc; j++)
	    	{
	    		str1.addTerm(-I.f.get(j), beta1_gamma.get(curr).get(j));
	    	}
        	
        	cplex2.addGe(str1, 0.);
    	    
        	for (int i = 0; i < I.number_of_items; i++)
        	{
        		IloLinearNumExpr str2 = cplex2.linearNumExpr();
        		IloLinearNumExpr str3 = cplex2.linearNumExpr(-M1);
        		IloLinearNumExpr str4 = cplex2.linearNumExpr();
        		
        		str2.addTerm(1., prod_gamma.get(curr).get(i));
        		str2.addTerm(-1., beta2_gamma.get(curr).get(i));
            	
        		str3.addTerm(1., prod_gamma.get(curr).get(i));
        		str3.addTerm(M1, x.get(i));
            	
        		str4.addTerm(1., prod_gamma.get(curr).get(i));
        		str4.addTerm(-1, beta2_gamma.get(curr).get(i));
        		str4.addTerm(M1, x.get(i));
            	
        		cplex2.addLe(str2, 0.);
        		cplex2.addLe(str3, 0.);
        		cplex2.addGe(str4, 0.);
        	}

        	x_curr = new ArrayList<Double>();
        	
        	if (cplex2.solve())
        	{   	
        	    for (int j = 0; j < I.number_of_items; j++)
        	    {
        	    	x_curr.add(cplex2.getValue(x.get(j)));
        			
        	    }  
        	    lower_bound = cplex2.getObjValue();
        	}
        	else 
            	System.out.println("No solution found");
      
        	if ((cplex2.getStatus() != IloCplex.Status.Optimal) || ((System.currentTimeMillis() - time_0)*Math.pow(10, -3) > 3600.))
        	{
       		 return new ArrayList<Double>();
       	    }
        	
        	curr++;
    	}
    	
        x_opt.add(lower_bound);
    	x_opt.add(upper_bound);
    	return x_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}
	
	public static ArrayList<Double> SolveMasterProblemRiskNeutral(IloCplex cplex, IloCplex cplex2, Instance I, double lower_bound_initial, double upper_bound_initial, ArrayList<ArrayList<Double>> support_constraints, double epsilon_l, int k_l, ArrayList<ArrayList<Double>> data_set_l, ArrayList<ArrayList<Double>> delta_l, ArrayList<Double> average_l, double M_objective, double binary_indicator, double time_0)
    {
	try 
    {
		cplex2.clearModel();
    	cplex2.setOut(null);
    	cplex2.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex2.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex2.setParam(IloCplex.Param.WorkMem, 2048);
    	cplex2.setParam(IloCplex.Param.TimeLimit, 3600.0);
        
        double lower_bound = lower_bound_initial, upper_bound = upper_bound_initial;
        
        double M = Double.MAX_VALUE;
        int number_of_sc = support_constraints.size() - 1; 
        
    	ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
    	IloNumVar z = cplex2.numVar(-M, M, "z");
    	
    	ArrayList<ArrayList<ArrayList<IloNumVar>>> nu_y = new ArrayList<ArrayList<ArrayList<IloNumVar>>>();
    	ArrayList<IloNumVar> lambda_y = new ArrayList<IloNumVar>(); 
    	
    	for (int j = 0; j < I.number_of_items; j++)
        {
        	x.add(cplex2.boolVar("x_" + j));
        }
     	
        IloLinearNumExpr objective = cplex2.linearNumExpr(); 
    	
    	objective.addTerm(1., z);
        	
    	cplex2.addMinimize(objective);
    	
    	//X
    	for (int j = 0; j < I.number_of_lc; j++)
    	{
    		IloLinearNumExpr str_1 = cplex2.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_1.addTerm(I.H.get(j).get(i), x.get(i));
        	}
    		 
    		cplex2.addLe(str_1, I.h.get(j));
         	
        }
    	
        ArrayList<Double> x_opt = new ArrayList<Double>();    	
    	
    	//Compute initial lower bound at x = 0, initialize Gamma
    	ArrayList<Double> x_curr = InitializeDoubleVector(I.number_of_items, 0.);
    	ArrayList<ArrayList<Double>> Y = new ArrayList<ArrayList<Double>>();
         
    	int curr = 0;
    	
    	while ((Math.abs(upper_bound - lower_bound) > Math.pow(10, -3)))
    	{
    		ArrayList<Double> sol_f = Instance.SolveFollowersProblemPessimisticRiskNeutral(cplex, I, x_curr, support_constraints, epsilon_l, k_l, data_set_l, average_l, binary_indicator);
        	
    		if (sol_f.size() == 0)
    		{
    			return new ArrayList<Double>();
       	    }
    		
    		if (sol_f.get(I.number_of_items) < upper_bound)
    		{
    			upper_bound = sol_f.get(I.number_of_items);
    	        x_opt = CopyDoubleVector1(x_curr, I.number_of_items);		
    		}
    		
    		ArrayList<Double> new_y = CopyDoubleVector1(sol_f, I.number_of_items);
    		
    		//ArrayList<Double> sol_2 = Instance.SolveFollowersProblemPessimisticRiskNeutralFixedY(cplex, I, new_y, support_constraints, epsilon_l, k_l, data_set_l, average_l, binary_indicator);
    		
        	Y.add(new_y);
    		
	        nu_y.add(new ArrayList<ArrayList<IloNumVar>>());
	        lambda_y.add(cplex2.numVar(0, M, "lam_" + curr));
	    	
	    	for (int k = 0; k < k_l; k++)
	    	{   
	    		nu_y.get(curr).add(new ArrayList<IloNumVar>());
	    		for (int l = 0; l < number_of_sc; l++)
	    		{
	    			nu_y.get(curr).get(k).add(cplex2.numVar(0, M, "nu_"+ curr + k + l));
	    		}
	    	}
    		
	    	//
	    	double sum = 0;
	    	for (int i = 0; i < I.number_of_items; i++)
		   	{	
		   		sum += average_l.get(i) * Y.get(curr).get(i);
		   	}
	    	
	    	IloLinearNumExpr str = cplex2.linearNumExpr(sum);
			
	    	str.addTerm(-1., z);
			str.addTerm(epsilon_l, lambda_y.get(curr));
			
			for (int i = 0; i < I.number_of_items; i++)
		   	{	
		   		str.addTerm(-M_objective * Y.get(curr).get(i), x.get(i));
		   	}
			
			for (int k = 0; k < k_l; k++)
				for (int l = 0; l < number_of_sc; l++)
			{
				str.addTerm(delta_l.get(k).get(l)/(k_l + 0.), nu_y.get(curr).get(k).get(l));
			}
			
	    	cplex2.addLe(str, 0.);
	    	
	    	//
	    	for (int k = 0; k < k_l; k++)
	    		for (int j = 0; j < I.number_of_items; j++)
	    		{
	    			IloLinearNumExpr str1 = cplex2.linearNumExpr(Y.get(curr).get(j));
	    	    	IloLinearNumExpr str2 = cplex2.linearNumExpr(-Y.get(curr).get(j));
	    	    	
	    	    	str1.addTerm(-1., lambda_y.get(curr));
	    	    	str2.addTerm(-1., lambda_y.get(curr));
	    	    		
	    	    	for (int l = 0; l < number_of_sc; l++)
	     	        {
	    	    		str1.addTerm(-support_constraints.get(l).get(j), nu_y.get(curr).get(k).get(l));
	    	    		str2.addTerm(support_constraints.get(l).get(j), nu_y.get(curr).get(k).get(l));
	     	        }
	    	    	
	    	    	cplex2.addLe(str1, 0.);
	    	    	cplex2.addLe(str2, 0.);
	    		}
	    	
            x_curr = new ArrayList<Double>();
        	ArrayList<Double> rhs_opt = new ArrayList<Double>();
        	if (cplex2.solve())
        	{   	
        	    for (int j = 0; j < I.number_of_items; j++)
        	    {
        	    	x_curr.add(cplex2.getValue(x.get(j)));
        	    }  
        	    lower_bound = cplex2.getObjValue();
        	    /*
        	    for (int l = 0; l < Y.size(); l++)
        	    {
        	    	double sum1 = 0;
        	    	sum1 += epsilon_l*cplex2.getValue(lambda_y.get(l));
        	    	
        	    	for (int i = 0; i < I.number_of_items; i++)
	    		   	{	
        	    		sum1 += average_l.get(i)*Y.get(l).get(i);
    	    			
	    		   		sum1 += -M_objective * Y.get(l).get(i) * cplex2.getValue(x.get(i));
	    		   	}
	    			
	    			for (int k = 0; k < k_l; k++)
	    				for (int r = 0; r < number_of_sc; r++)
	    			{
	    				sum1 += delta_l.get(k).get(r)/(k_l + 0.) * cplex2.getValue(nu_y.get(l).get(k).get(r));
	    			}		
        	    	
        	    	rhs_opt.add(sum1);
        	    }
        	    */
        	}
        	else 
            	System.out.println("No solution found");
             
        	if ((cplex2.getStatus() != IloCplex.Status.Optimal) || ((System.currentTimeMillis() - time_0)*Math.pow(10, -3) > 3600.0))
        	{
        		return new ArrayList<Double>();
       	    }
  
        	curr++;
    	}  
    	x_opt.add(lower_bound);
    	x_opt.add(upper_bound);
    		
    	return x_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}
	
	
	public static ArrayList<Double> SolveFollowersProblemPessimisticAmbiguityFree(IloCplex cplex, Instance I, ArrayList<Double> x, double alpha_l, int k_l, ArrayList<ArrayList<Double>> data_set_l)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
    	cplex.setParam(IloCplex.Param.TimeLimit, 3600.0);
        
        double M = Double.MAX_VALUE;
    	
    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
    	
    	ArrayList<IloNumVar> gamma_l = new ArrayList<IloNumVar>();
    	ArrayList<ArrayList<IloNumVar>> prod_gamma_y = new ArrayList<ArrayList<IloNumVar>>();
    	
    	
    	for (int k = 0; k < k_l; k++)
    	{
    		gamma_l.add(cplex.boolVar("gamma_l_"+ k));
    		
    		prod_gamma_y.add(new ArrayList<IloNumVar>());
        
    		for (int j = 0; j < I.number_of_items; j++)
            {
    			prod_gamma_y.get(k).add(cplex.numVar(0, M, "prod_"+ k + j));
            }
    	}
    	
        for (int j = 0; j < I.number_of_items; j++)
        {
        	y.add(cplex.numVar(0, 1, "y_" + j));
        }
    	
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
		for (int k = 0; k < k_l; k++)
			for (int j = 0; j < I.number_of_items; j++)
            {
				objective.addTerm(1./((int)Math.round(k_l * (1. - alpha_l)))*data_set_l.get(j).get(k), prod_gamma_y.get(k).get(j));
    	    }
	
		cplex.addMaximize(objective);
    	
    	//Constraints
    	//Feasibility
    	 //Y
    	for (int j = 0; j < I.number_of_fc; j++)
    	{
    		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_2.addTerm(I.F.get(j).get(i), y.get(i));
        	}
    		 
    		cplex.addLe(str_2, I.f.get(j));
         	
        }
    	
    	for (int j = 0; j < I.number_of_items; j++)
    	{
    		IloLinearNumExpr str = cplex.linearNumExpr(); 
    		str.addTerm(1., y.get(j));
    		cplex.addLe(str, 1. - x.get(j));
    	}
        
        //
		IloLinearNumExpr str = cplex.linearNumExpr();
		
		for (int k = 0; k < k_l; k++)
			str.addTerm(1., gamma_l.get(k));
		
		cplex.addEq(str, (int)Math.round(k_l * (1. - alpha_l)));
	
        //
		for (int k = 0; k < k_l; k++)
    		for (int i = 0; i < I.number_of_items; i++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
        		IloLinearNumExpr str2 = cplex.linearNumExpr(1.);
        		IloLinearNumExpr str3 = cplex.linearNumExpr();
        		
        		str1.addTerm(1., prod_gamma_y.get(k).get(i));
        		str1.addTerm(-1., y.get(i));
        		
        		str2.addTerm(1., prod_gamma_y.get(k).get(i));
        		str2.addTerm(-1., y.get(i));
        		str2.addTerm(-1., gamma_l.get(k));
        		
        		str3.addTerm(1., prod_gamma_y.get(k).get(i));
        		str3.addTerm(-1., gamma_l.get(k));
        		
        		cplex.addLe(str1, 0.);
        		cplex.addGe(str2, 0.);
        		cplex.addLe(str3, 0.);
    		}
        
    	ArrayList<Double> gamma_opt = new ArrayList<Double>();
    	ArrayList<Double> y_opt = new ArrayList<Double>();
    	
    	//System.out.println(cplex.getNcols());
    	//System.out.println(cplex.getNrows());
    	
    	if (cplex.solve())
    	{   	
    	
    	    for (int j = 0; j < k_l; j++)
    	    {
    	    	gamma_opt.add(cplex.getValue(gamma_l.get(j)));
    	    }  
    	    
    	    gamma_opt.add(cplex.getObjValue());
    	
    	}
    	else 
        	System.out.println("No solution found");
    	
    	if (cplex.getStatus() != IloCplex.Status.Optimal)
    	{
    		 return new ArrayList<Double>();
    	}
    	
    	return gamma_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}
    
	public static ArrayList<Double> SolveFollowersProblemPessimisticRiskNeutralFixedY(IloCplex cplex, Instance I, ArrayList<Double> y, ArrayList<ArrayList<Double>> support_constraints, double epsilon_l, int k_l, ArrayList<ArrayList<Double>> data_set_l, ArrayList<Double> average_l, double binary_indicator)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
        
        double M = Double.MAX_VALUE;
        int number_of_sc = support_constraints.size() - 1; 
        
    	ArrayList<ArrayList<IloNumVar>> xi_l = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> eta_l = new ArrayList<ArrayList<IloNumVar>>();
    	
    
    	for (int k = 0; k < k_l; k++)
    	{
    		xi_l.add(new ArrayList<IloNumVar>());
    		eta_l.add(new ArrayList<IloNumVar>());
    	
    		for (int j = 0; j < I.number_of_items; j++)
            {
    			xi_l.get(k).add(cplex.numVar(-M, M, "xi_"+ k + j));
    			eta_l.get(k).add(cplex.numVar(0, M, "eta_"+ k + j));
            }
    	}
    		
    	double sum = 0;
    	for (int j = 0; j < I.number_of_items; j++)	
		{
			sum += average_l.get(j) * y.get(j);
		}
    	
    	IloLinearNumExpr objective = cplex.linearNumExpr(sum); 
    	
    	
		for (int j = 0; j < I.number_of_items; j++)	
		{
			for (int k = 0; k < k_l; k++)
            {
				objective.addTerm(-1./k_l*y.get(j), xi_l.get(k).get(j));
    	    }
		}
			
    	cplex.addMaximize(objective);
    	
    	//Constraints
    	//Feasibility
    	for (int k = 0; k < k_l; k++)
    		for (int l = 0; l < number_of_sc; l++)
	    	{
	    		IloLinearNumExpr str = cplex.linearNumExpr();
	    		
	    		double sum1 = 0;
	    		
	    		for (int j = 0; j < I.number_of_items; j++)
	        	{
	    			str.addTerm(-support_constraints.get(l).get(j), xi_l.get(k).get(j));
	    			sum1 += support_constraints.get(l).get(j) * data_set_l.get(j).get(k);
	        	}
	    		
	    		cplex.addLe(str, support_constraints.get(number_of_sc).get(l) - sum1); 
	    	}
    	
 
    	//
    	for (int k = 0; k < k_l; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str1.addTerm(-1., eta_l.get(k).get(j));
    			str1.addTerm(-1., xi_l.get(k).get(j));
    			
    			str2.addTerm(1., xi_l.get(k).get(j));
    			str2.addTerm(-1., eta_l.get(k).get(j));
    			
    			cplex.addLe(str1, 0.);
    		    cplex.addLe(str2, 0.);
    		}
    	
    	
    	//
    	IloLinearNumExpr str = cplex.linearNumExpr(); 
    	
    	for (int k = 0; k < k_l; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    			str.addTerm(1./k_l, eta_l.get(k).get(j));
		
		cplex.addLe(str, epsilon_l);
	
    	ArrayList<Double> opt_val = new ArrayList<Double>();
    	ArrayList<ArrayList<Double>> xi_opt = new ArrayList<ArrayList<Double>>();
    	ArrayList<ArrayList<Double>> eta_opt = new ArrayList<ArrayList<Double>>();
    	
    	if (cplex.solve())
    	{   	
    		/*
    		for (int k = 0; k < k_l; k++)
    		{
    			xi_opt.add(new ArrayList<Double>());
    			eta_opt.add(new ArrayList<Double>());
    			
        		for (int i = 0; i < I.number_of_items; i++)
        		{
        			xi_opt.get(k).add(cplex.getValue(xi_l.get(k).get(i)));
        			eta_opt.get(k).add(cplex.getValue(eta_l.get(k).get(i)));
        		}
    		}
            
    		System.out.println(xi_opt);
    		System.out.println(eta_opt);*/
    	    opt_val.add(cplex.getObjValue());
    	}
    	else 
        	System.out.println("No solution found");
    	
    	return opt_val;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}
	
	public static ArrayList<Double> SolveFollowersProblemPessimisticRiskNeutral(IloCplex cplex, Instance I, ArrayList<Double> x, ArrayList<ArrayList<Double>> support_constraints, double epsilon_l, int k_l, ArrayList<ArrayList<Double>> data_set_l, ArrayList<Double> average_l, double binary_indicator)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
    	cplex.setParam(IloCplex.Param.TimeLimit, 3600.0);
    	
        double M = Double.MAX_VALUE;
        int number_of_sc = support_constraints.size() - 1; 
        
    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
    	
    	ArrayList<ArrayList<IloNumVar>> xi_l = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> eta_l = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> prod_xi_y = new ArrayList<ArrayList<IloNumVar>>();
    	
    
    	for (int k = 0; k < k_l; k++)
    	{
    		xi_l.add(new ArrayList<IloNumVar>());
    		eta_l.add(new ArrayList<IloNumVar>());
    		prod_xi_y.add(new ArrayList<IloNumVar>());
            
    		for (int j = 0; j < I.number_of_items; j++)
            {
    			double xi_min = Math.max(-k_l*epsilon_l, -I.max_cost.get(j) + I.min_cost.get(j));
    			double xi_max = Math.min(k_l*epsilon_l, I.max_cost.get(j) - I.min_cost.get(j));
    	    	
    			xi_l.get(k).add(cplex.numVar(-M, M, "xi_"+ k + j));
    			eta_l.get(k).add(cplex.numVar(0, M, "eta_"+ k + j));
    			prod_xi_y.get(k).add(cplex.numVar(xi_min, xi_max, "prod_"+ k + j));
            }
    	}
    	
        for (int j = 0; j < I.number_of_items; j++)
        {
        	y.add(cplex.boolVar("y_" + j));
        }
    	
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
		for (int j = 0; j < I.number_of_items; j++)	
		{
			objective.addTerm(average_l.get(j), y.get(j));

			for (int k = 0; k < k_l; k++)
            {
				objective.addTerm(-1./k_l, prod_xi_y.get(k).get(j));
    	    }
		}
			
    	cplex.addMaximize(objective);
    	
    	//Constraints
    	//Feasibility
	   	 //Y
	   	for (int j = 0; j < I.number_of_fc; j++)
	   	{
	   		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
	   		
	   		for (int i = 0; i < I.number_of_items; i++)
	       	{
	   			str_2.addTerm(I.F.get(j).get(i), y.get(i));
	       	}
	   		 
	   		cplex.addLe(str_2, I.f.get(j));	
	       }
    	
    	for (int j = 0; j < I.number_of_items; j++)
    	{
    		IloLinearNumExpr str = cplex.linearNumExpr(); 
    		str.addTerm(1., y.get(j));
    		cplex.addLe(str, 1. - x.get(j));
    	}
        
    	//
    	for (int k = 0; k < k_l; k++)
    		for (int l = 0; l < number_of_sc; l++)
	    	{
	    		IloLinearNumExpr str = cplex.linearNumExpr();
	    		
	    		double sum = 0;
	    		
	    		for (int j = 0; j < I.number_of_items; j++)
	        	{
	    			str.addTerm(-support_constraints.get(l).get(j), xi_l.get(k).get(j));
	    			sum += support_constraints.get(l).get(j) * data_set_l.get(j).get(k);
	        	}
	    		
	    		cplex.addLe(str, support_constraints.get(number_of_sc).get(l) - sum); 
	    	}
    	
    	for (int k = 0; k < k_l; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str1.addTerm(-1., eta_l.get(k).get(j));
    			str1.addTerm(-1., xi_l.get(k).get(j));
    			
    			str2.addTerm(1., xi_l.get(k).get(j));
    			str2.addTerm(-1., eta_l.get(k).get(j));
    			
    			cplex.addLe(str1, 0.);
    		    cplex.addLe(str2, 0.);
    		}
    	
    	
    	//
    	IloLinearNumExpr str = cplex.linearNumExpr(); 
    	
    	for (int k = 0; k < k_l; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    			str.addTerm(1./k_l, eta_l.get(k).get(j));
		
		cplex.addLe(str, epsilon_l);
	    
		//products xi * y
		for (int k = 0; k < k_l; k++)
    		for (int i = 0; i < I.number_of_items; i++)
    		{	
    			double xi_min = Math.max(-k_l*epsilon_l, -I.max_cost.get(i) + I.min_cost.get(i));
    			double xi_max = Math.min(k_l*epsilon_l, I.max_cost.get(i) - I.min_cost.get(i));
    	    	
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
        		IloLinearNumExpr str2 = cplex.linearNumExpr();
        		IloLinearNumExpr str3 = cplex.linearNumExpr();
        		IloLinearNumExpr str4 = cplex.linearNumExpr();
        		
        		
        		str1.addTerm(1., prod_xi_y.get(k).get(i));
        		str1.addTerm(-xi_max, y.get(i));
        		
        		str2.addTerm(1., prod_xi_y.get(k).get(i));
        		str2.addTerm(-xi_min, y.get(i));
        			
        		str3.addTerm(1., prod_xi_y.get(k).get(i));
        		str3.addTerm(-xi_min, y.get(i));
        		str3.addTerm(-1., xi_l.get(k).get(i));
        		
        		str4.addTerm(1., prod_xi_y.get(k).get(i));
        		str4.addTerm(-xi_max, y.get(i));
        		str4.addTerm(-1., xi_l.get(k).get(i));
        		
        		cplex.addLe(str1, 0.);
        		cplex.addGe(str2, 0.);
        		cplex.addLe(str3, -xi_min);
        		cplex.addGe(str4, -xi_max);
    		}
 	
    	ArrayList<Double> y_opt = new ArrayList<Double>();
        ArrayList<ArrayList<Double>> prod_opt = new ArrayList<ArrayList<Double>>();
        ArrayList<ArrayList<Double>> xi_opt = new ArrayList<ArrayList<Double>>();
        ArrayList<ArrayList<Double>> eta_opt = new ArrayList<ArrayList<Double>>();
        
    	if (cplex.solve())
    	{     		
    		for (int j = 0; j < I.number_of_items; j++)
    	    {
    	    	y_opt.add(cplex.getValue(y.get(j)));
    	    }  
    	    
    		/*
    		for (int k = 0; k < k_l; k++)
    		{
    			prod_opt.add(new ArrayList<Double>());
    			xi_opt.add(new ArrayList<Double>());
    			eta_opt.add(new ArrayList<Double>());
    			
        		for (int i = 0; i < I.number_of_items; i++)
        		{
        			prod_opt.get(k).add(cplex.getValue(prod_xi_y.get(k).get(i)));
        			xi_opt.get(k).add(cplex.getValue(xi_l.get(k).get(i)));
        			eta_opt.get(k).add(cplex.getValue(eta_l.get(k).get(i)));
        		}
    		}
    
    		System.out.println(prod_opt);
    		System.out.println(xi_opt);
    		System.out.println(eta_opt);*/
    	    y_opt.add(cplex.getObjValue());   	    
    	}
    	else 
        	System.out.println("No solution found");
    	
    	if (cplex.getStatus() != IloCplex.Status.Optimal)
    	{
    		 return new ArrayList<Double>();
    	}
    	
    	//System.out.println(cplex.getModel());
    	return y_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}
	
	public static ArrayList<ArrayList<Double>> SolveLeadersProblemFixedX(IloCplex cplex, Instance I, ArrayList<Double> x, double epsilon_l, double epsilon_f, double alpha_l, double alpha_f, int k_l, int k_f, ArrayList<ArrayList<Double>> data_set_l, ArrayList<ArrayList<Double>> data_set_f,
			ArrayList<ArrayList<Double>> support_constraints, String type_of_risk_l, String type_of_risk_f,  ArrayList<ArrayList<Double>> delta_l,  ArrayList<ArrayList<Double>> delta_f, ArrayList<Double> average_l, ArrayList<Double> average_f)
    {
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
        
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	ArrayList<ArrayList<IloNumVar>> nu_l = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> nu_f = new ArrayList<ArrayList<IloNumVar>>();
    	
    	//ArrayList<IloNumVar> x = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> y = new ArrayList<IloNumVar>();
    	
    	ArrayList<IloNumVar> s_l = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> s_f = new ArrayList<IloNumVar>();
    	
    	IloNumVar lambda_l = cplex.numVar(0, M, "lam_l");
    	IloNumVar lambda_f = cplex.numVar(0, M, "lam_f");
    	
    	IloNumVar t_l = cplex.numVar(-M, M, "t_l");
    	IloNumVar t_f = cplex.numVar(-M, M, "t_f");
    	
    	ArrayList<ArrayList<IloNumVar>> mu1_f = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> mu2_f = new ArrayList<ArrayList<IloNumVar>>();
    	
    	ArrayList<IloNumVar> beta1_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> beta2_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> gamma_f = new ArrayList<IloNumVar>();
    	//ArrayList<IloNumVar> prod_f = new ArrayList<IloNumVar>();
    	
    	
    	for (int k = 0; k < k_l; k++)
    	{
    		s_l.add(cplex.numVar(0, M, "s_l_"+ k));
    		nu_l.add(new ArrayList<IloNumVar>());
    	
    		for (int l = 0; l < number_of_sc; l++)
    		{
    			nu_l.get(k).add(cplex.numVar(0, M, "nu_l_"+ k + l));	
    		}
    	}
    	
    	for (int k = 0; k < k_f; k++)
    	{
    		s_f.add(cplex.numVar(0, M, "s_f_"+ k));
    		nu_f.add(new ArrayList<IloNumVar>());
    	
    		for (int l = 0; l < number_of_sc; l++)
    		{
    			nu_f.get(k).add(cplex.numVar(0, M, "nu_f_"+ k + l));	
    		}
    	}
    
        for (int j = 0; j < I.number_of_items; j++)
        {
        		y.add(cplex.numVar(0, 1, "y_" + j));
        }
    	
        //Follower's dual variables
        for (int j = 0; j < I.number_of_fc; j++)
        {
        	beta1_f.add(cplex.numVar(0, M, "beta1_f_"+ j));
        }
        
        for (int i = 0; i < I.number_of_items; i++)
        {
        	beta2_f.add(cplex.numVar(0, M, "beta2_f_"+ i));
        }
        
        for (int k = 0; k < k_f; k++)
    	{
    		gamma_f.add(cplex.numVar(0, M, "gamma_f_"+ k));
    		
    		mu1_f.add(new ArrayList<IloNumVar>());
    		mu2_f.add(new ArrayList<IloNumVar>());
    	
    		 for (int j = 0; j < I.number_of_items; j++)
    		{
    			mu1_f.get(k).add(cplex.numVar(0, M, "mu1_f_"+ k + j));	
    			mu2_f.get(k).add(cplex.numVar(0, M, "mu2_f_"+ k + j));
    		}
    	}
        
    	IloLinearNumExpr objective = cplex.linearNumExpr(); 
    	
    	if ("Risk-averse".equals(type_of_risk_l))
    	{   
    		objective.addTerm(1., t_l);
    		objective.addTerm(epsilon_l/(1. - alpha_l), lambda_l);
    		
    		for (int k = 0; k < k_l; k++)
        	{
        	     objective.addTerm(1./(k_l*(1. - alpha_l)), s_l.get(k));
        	}
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk_l))
    	{
    		objective.addTerm(epsilon_l, lambda_l);
    		
    		for (int i = 0; i < I.number_of_items; i++)
    			objective.addTerm(average_l.get(i), y.get(i));
    		
    		for (int k = 0; k < k_l; k++)
    			for (int l = 0; l < number_of_sc; l++)
    		{
    			objective.addTerm(delta_l.get(k).get(l)/(k_l + 0.), nu_l.get(k).get(l));
    		}
    		
    	}
    	
    	cplex.addMinimize(objective);
    	
    	//Constraints
    	//Feasibility
    	//IloLinearNumExpr str_1 = cplex.linearNumExpr(); 
    	for (int j = 0; j < I.number_of_fc; j++)
    	{
    		IloLinearNumExpr str_2 = cplex.linearNumExpr(); 
    		
    		for (int i = 0; i < I.number_of_items; i++)
        	{
    			str_2.addTerm(I.F.get(j).get(i), y.get(i));
        	}
    		 
    		cplex.addLe(str_2, I.f.get(j));
         	
        }
    	
    	for (int j = 0; j < I.number_of_items; j++)
    	{
    		IloLinearNumExpr str = cplex.linearNumExpr(); 
    		str.addTerm(1., y.get(j));
        	cplex.addLe(str, 1. - x.get(j));
    	}
        
    	//
     	for (int k = 0; k < k_l; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
    	    	
    	    	str1.addTerm(-1., lambda_l);
    	    	str2.addTerm(-1., lambda_l);
    	    	
    	    	str1.addTerm(-1., y.get(j));
    	    	str2.addTerm(1., y.get(j));
    	    	
    	    	for (int l = 0; l < number_of_sc; l++)
     	        {
    	    		str1.addTerm(support_constraints.get(l).get(j), nu_l.get(k).get(l));
    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_l.get(k).get(l));
     	        }
    	    	
    	    	cplex.addLe(str1, 0.);
    	    	cplex.addLe(str2, 0.);
    		}
    	
     	//
    	for (int k = 0; k < k_f; k++)
    		for (int j = 0; j < I.number_of_items; j++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    	    	IloLinearNumExpr str2 = cplex.linearNumExpr();
    	    	
    	    	str1.addTerm(-1., lambda_f);
    	    	str2.addTerm(-1., lambda_f);
    	    	
    	    	str1.addTerm(1., y.get(j));
    	    	str2.addTerm(-1., y.get(j));
    	    	
    	    	for (int l = 0; l < number_of_sc; l++)
     	        {
    	    		str1.addTerm(support_constraints.get(l).get(j), nu_f.get(k).get(l));
    	    		str2.addTerm(-support_constraints.get(l).get(j), nu_f.get(k).get(l));
     	        }
    	    	
    	    	cplex.addLe(str1, 0.);
    	    	cplex.addLe(str2, 0.);
    		}
   	    
    	//
    	if ("Risk-averse".equals(type_of_risk_l))
    	{   
    		for (int k = 0; k < k_l; k++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			
    			 for (int i = 0; i < I.number_of_items; i++)
    			 {
    				 str1.addTerm(data_set_l.get(i).get(k), y.get(i));
    		     }
    			
    			str1.addTerm(-1., t_l);
    			
    
        		for (int l = 0; l < number_of_sc; l++)
        		{
        			str1.addTerm(delta_l.get(k).get(l), nu_l.get(k).get(l));
        		}
    			
    			str1.addTerm(-1., s_l.get(k));
    			cplex.addLe(str1, 0.);
    			
                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str2.addTerm(1., s_l.get(k));
    			cplex.addGe(str2, 0);*/
    		}
    		
    	}
    	
    	//
    	if ("Risk-averse".equals(type_of_risk_f))
    	{   
    		for (int k = 0; k < k_f; k++)
    		{
    			IloLinearNumExpr str1 = cplex.linearNumExpr();
    			
    			 for (int i = 0; i < I.number_of_items; i++)
    			 {
    				 str1.addTerm(-data_set_f.get(i).get(k), y.get(i));
    		     }
    			
    			str1.addTerm(1., t_f);
    			
    
        		for (int l = 0; l < number_of_sc; l++)
        		{
        			str1.addTerm(delta_f.get(k).get(l), nu_f.get(k).get(l));
        		}
    			
    			str1.addTerm(-1., s_f.get(k));
    			cplex.addLe(str1, 0.);
    			
                /*IloLinearNumExpr str2 = cplex.linearNumExpr();
    			
    			str2.addTerm(1., s_f.get(k));
    			cplex.addGe(str2, 0);*/
    		}
    		
    	}
    	
    	//Dual constraints of the follower
    	if ("Risk-averse".equals(type_of_risk_f))
    	{
    		//
	    	for (int i = 0; i < I.number_of_items; i++)
	    	 {	
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
				
	    		for (int j = 0; j < I.number_of_fc; j++)
	           	{
	       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
	           	}
	       		
	       		str1.addTerm(1., beta2_f.get(i));

	    		
	    		for (int k = 0; k < k_f; k++)
	    		{
	    			str1.addTerm(-data_set_f.get(i).get(k), gamma_f.get(k));
	    			str1.addTerm(-1., mu1_f.get(k).get(i));
	    			str1.addTerm(1., mu2_f.get(k).get(i));
	    		}
	        	
	    		cplex.addGe(str1, 0.);
	    	}
	    	
	    	//
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		for (int l = 0; l < number_of_sc; l++)
	            {
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		
	    			str1.addTerm(delta_f.get(k).get(l), gamma_f.get(k));
		    		
		    		for (int i = 0; i < I.number_of_items; i++)
			    	 {
		    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
		    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
			    	 }
		    		cplex.addGe(str1, 0.);
		    	} 
	    	}
	    	
	    	//
	    	IloLinearNumExpr str = cplex.linearNumExpr();
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		
	    		str.addTerm(1., gamma_f.get(k));
	    		str1.addTerm(1., gamma_f.get(k));
	    		
	    		cplex.addLe(str1, 1./((1. - alpha_f) * k_f));
	    	}
	    	
	    	cplex.addEq(str, 1.);
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk_f))
    	{   
    		//
	    	for (int i = 0; i < I.number_of_items; i++)
	    	 {	
	    		IloLinearNumExpr str1 = cplex.linearNumExpr(-average_f.get(i));
				
	    		for (int j = 0; j < I.number_of_fc; j++)
	           	{
	       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
	           	}
	       		
	       		str1.addTerm(1., beta2_f.get(i));

	    		
	    		for (int k = 0; k < k_f; k++)
	    		{
	    			str1.addTerm(-1., mu1_f.get(k).get(i));
	    			str1.addTerm(1., mu2_f.get(k).get(i));
	    		}
	        	
	    		cplex.addGe(str1, 0.);
	    	}
	    	
	    	//
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		for (int l = 0; l < number_of_sc; l++)
	            {
	    			IloLinearNumExpr str1 = cplex.linearNumExpr(delta_f.get(k).get(l)/(k_f + 0.));
		    		
		    		for (int i = 0; i < I.number_of_items; i++)
			    	 {
		    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
		    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
			    	 }
		    		cplex.addGe(str1, 0.);
		    	}
	    	}
    	}
    	
    	//
    	IloLinearNumExpr str3 = cplex.linearNumExpr();
    	for (int k = 0; k < k_f; k++)
    	{
    		for (int i = 0; i < I.number_of_items; i++)
	    	 {
    			str3.addTerm(1., mu1_f.get(k).get(i));
    			str3.addTerm(1., mu2_f.get(k).get(i));
	    	 }
    	}
    	
    	if ("Risk-averse".equals(type_of_risk_f))
    		cplex.addLe(str3, epsilon_f/(1. - alpha_f));
    	
    	if ("Risk-neutral".equals(type_of_risk_f))
    		cplex.addLe(str3, epsilon_f);
    	
    	
    	////Feasibility
    	if ("Risk-averse".equals(type_of_risk_f))
    	{
    		IloLinearNumExpr str1 = cplex.linearNumExpr();
    		IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	
    		
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		str1.addTerm(1 - x.get(i), beta2_f.get(i));
	    	}
	    	
	    	for (int j = 0; j < I.number_of_fc; j++)
	    	{
	    		str1.addTerm(I.f.get(j), beta1_f.get(j));
	    	}
	    	
	    	str2.addTerm(1., t_f);
	    	str2.addTerm(-epsilon_f/(1. - alpha_f), lambda_f);
	    	
	    	
	    	for (int k = 0; k < k_f; k++)
	    	{
	    		str2.addTerm(-1./(k_f * (1. - alpha_f)), s_f.get(k));
	    	}
	    	
	    	cplex.addEq(str1, str2);
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk_f))
    	{
    		IloLinearNumExpr str1 = cplex.linearNumExpr();
    		IloLinearNumExpr str2 = cplex.linearNumExpr();
	    	
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		str1.addTerm(-x.get(i), beta2_f.get(i));
	    		str1.addTerm(1., beta2_f.get(i));
	    	}
	        
	    	for (int j = 0; j < I.number_of_fc; j++)
	    	{
	    		str1.addTerm(I.f.get(j), beta1_f.get(j));
	    	}
	    	
	    	str2.addTerm(-epsilon_f, lambda_f);
	    	
	    	for (int i = 0; i < I.number_of_items; i++)
	    	{
	    		str2.addTerm(average_f.get(i), y.get(i));
	    	}
	    	
	    	for (int k = 0; k < k_f; k++)
    			for (int l = 0; l < number_of_sc; l++)
    		{
    			str2.addTerm(-delta_f.get(k).get(l)/(k_f + 0.), nu_f.get(k).get(l));
    		}
	    	
	    	cplex.addEq(str1, str2);
    	}
    	
    	ArrayList<ArrayList<Double>> sol_opt = new ArrayList<ArrayList<Double>>();
    	ArrayList<Double> y_opt = new ArrayList<Double>();
    	ArrayList<Double> x_opt = new ArrayList<Double>();
    	
    	
    	if (cplex.solve())
    	{   		
    		
    	    for (int j = 0; j < I.number_of_items; j++)
    	    {
    	    	//x_opt.add(cplex.getValue(x.get(j)));
    			y_opt.add(cplex.getValue(y.get(j)));
    	    }  
    	    
    	    x_opt.add(cplex.getObjValue());
    	    
    	    sol_opt.add(x_opt);
    	    sol_opt.add(y_opt);
    	}
    	else 
        	System.out.println("No solution found");
    	
    	return sol_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<ArrayList<Double>>();
}

	public static ArrayList<ArrayList<Double>> ComputeDelta(Instance I, ArrayList<ArrayList<Double>> data_set, ArrayList<ArrayList<Double>> support_constraints)
	{    
	   ArrayList<ArrayList<Double>> delta_f = new  ArrayList<ArrayList<Double>>();
	        
	   int sample_size = data_set.get(0).size();
	   int number_of_sc = support_constraints.size() - 1; 
	    	
       for (int k = 0; k < sample_size; k++)
       {
           delta_f.add(new ArrayList<Double>());  
	        
           for (int l = 0; l < number_of_sc; l++)
	        {
	        	double sum = 0;
	        	
	        	for (int i = 0; i < I.number_of_items; i++)
	        	{
	        		sum += support_constraints.get(l).get(i) * data_set.get(i).get(k);
	        	}
	        		
	        	delta_f.get(k).add(support_constraints.get(number_of_sc).get(l) - sum);
	        }
   	}
	return delta_f;
	}
	
	public static ArrayList<ArrayList<Double>> JoinDataSets(ArrayList<ArrayList<Double>> data_set_1, ArrayList<ArrayList<Double>> data_set_2, int number_of_rows) 
    {
		ArrayList<ArrayList<Double>> data_set = new ArrayList<ArrayList<Double>>();
    	
		int sample_size_1 = data_set_1.get(0).size();
		int sample_size_2 = data_set_2.get(0).size();
	
		
		for (int i = 0; i < number_of_rows; i++)
		{
			data_set.add(new ArrayList<Double>());
			
			for (int j = 0; j < sample_size_1; j++)
				data_set.get(data_set.size() - 1).add(data_set_1.get(i).get(j));
			
			for (int j = 0; j < sample_size_2; j++)
				data_set.get(data_set.size() - 1).add(data_set_2.get(i).get(j));
		}
		
    	return data_set;
    }

	public static ArrayList<ArrayList<Double>> GenerateNewDataSet(Instance I, ArrayList<ArrayList<Double>> data_set_true, int number_of_rows, int k_lf, double noise_level, ArrayList<ArrayList<Double>> shifts, ArrayList<ArrayList<Double>> missed) 
    {
		ArrayList<ArrayList<Double>> data_set = new ArrayList<ArrayList<Double>>();
    	
		int sample_size = data_set_true.get(0).size();
	    
		for (int i = 0; i < number_of_rows; i++)
		{
			data_set.add(new ArrayList<Double>());
		
			//the common part
			for (int j = 0; j < k_lf; j++)
				data_set.get(data_set.size() - 1).add(data_set_true.get(i).get(j));
			
			//the new part
			for (int j = k_lf; j < sample_size; j++)
			{
				double rand = missed.get(i).get(j);
			    if (rand < 0.5)
				{   
					double elem = data_set_true.get(i).get(j);
				    
					/*double l0 = I.min_cost.get(i);
					double u0 = I.max_cost.get(i);
				
					double lb = -1;
					double ub = -1;
					
					int num = (int)(1./(2*noise_level));
					
					for (int l = 0; l < num; l++)
					{
						double left = l0 + (l + 0.)/num*(u0 - l0);
						double right = l0 + (l + 1.)/num*(u0 - l0);
						
						if (elem >= left && elem <= right)
						{
							lb = left;
							ub = right;
							break;
						}
					}*/
					
					double a = shifts.get(i).get(j);
					double lb = Math.max(I.min_cost.get(i), elem - a*noise_level);
					double ub = Math.min(I.max_cost.get(i), elem + (1 - a)*noise_level);
					
					data_set.get(data_set.size() - 1).add(ThreadLocalRandom.current().nextDouble(lb, ub));
				}
			    else
				{
			    	data_set.get(data_set.size() - 1).add(data_set_true.get(i).get(j));	
				}
			}
		}
		
    	return data_set;
    }

	public static ArrayList<Double> SolveFollowersProblem2(IloCplex cplex, Instance I, ArrayList<Double> x, double epsilon_f, double alpha_f, int sample_size, ArrayList<ArrayList<Double>> data_set,  ArrayList<ArrayList<Double>> support_constraints, String type_of_risk, ArrayList<ArrayList<Double>> delta_f, ArrayList<Double> average_f,  double binary_indicator){
	try 
    {
    	cplex.clearModel();
    	cplex.setOut(null);
    	cplex.setParam(IloCplex.IntParam.VarSel, 4);
    	cplex.setParam(IloCplex.BooleanParam.MemoryEmphasis, true);
    	cplex.setParam(IloCplex.Param.WorkMem, 2048);
       
        double M = Double.MAX_VALUE;
    	int number_of_sc = support_constraints.size() - 1; 
    	
    	ArrayList<ArrayList<IloNumVar>> mu1_f = new ArrayList<ArrayList<IloNumVar>>();
    	ArrayList<ArrayList<IloNumVar>> mu2_f = new ArrayList<ArrayList<IloNumVar>>();
    	
    	ArrayList<IloNumVar> beta1_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> beta2_f = new ArrayList<IloNumVar>();
    	ArrayList<IloNumVar> gamma_f = new ArrayList<IloNumVar>();
    	//ArrayList<IloNumVar> prod_f = new ArrayList<IloNumVar>();
    	
    	
    	 for (int j = 0; j < I.number_of_fc; j++)
         {
         	beta1_f.add(cplex.numVar(0, M, "beta1_f_"+ j));
         }
         
         for (int i = 0; i < I.number_of_items; i++)
         {
         	beta2_f.add(cplex.numVar(0, M, "beta2_f_"+ i));
         }
        
        for (int k = 0; k < sample_size; k++)
    	{
    		gamma_f.add(cplex.numVar(0, M, "gamma_f_"+ k));
    		
    		mu1_f.add(new ArrayList<IloNumVar>());
    		mu2_f.add(new ArrayList<IloNumVar>());
    	
    	    for (int j = 0; j < I.number_of_items; j++)
    		{
    			mu1_f.get(k).add(cplex.numVar(0, M, "mu1_f_"+ k + j));	
    			mu2_f.get(k).add(cplex.numVar(0, M, "mu2_f_"+ k + j));
    		}
    	}
        
    	IloLinearNumExpr objective = cplex.linearNumExpr(0.); 
    	
    	
    	for (int i = 0; i < I.number_of_items; i++)
    	{
    		//str1.addTerm(1., prod_f.get(i));
    		objective.addTerm(-x.get(i) + 1., beta2_f.get(i));
    	}
    	
    	for (int j = 0; j < I.number_of_fc; j++)
    	{
    		objective.addTerm(I.f.get(j), beta1_f.get(j));
    	}
    	
    	cplex.addMinimize(objective);
    	
    	//Constraints
    		
    	//Dual constraints of the follower
    	if ("Risk-averse".equals(type_of_risk))
    	{
    		//
	    	for (int i = 0; i < I.number_of_items; i++)
	    	 {	
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
				
	    		for (int j = 0; j < I.number_of_fc; j++)
	           	{
	       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
	           	}
	       		
	       		str1.addTerm(1., beta2_f.get(i));
	    		
	    		for (int k = 0; k < sample_size; k++)
	    		{
	    			str1.addTerm(-data_set.get(i).get(k), gamma_f.get(k));
	    			str1.addTerm(-1., mu1_f.get(k).get(i));
	    			str1.addTerm(1., mu2_f.get(k).get(i));
	    		}
	        	
	    		cplex.addGe(str1, 0.);
	    	}
	    	
	    	//
	    	for (int k = 0; k < sample_size; k++)
	    	{
	    		for (int l = 0; l < number_of_sc; l++)
	            {
	    			IloLinearNumExpr str1 = cplex.linearNumExpr();
		    		
	    			str1.addTerm(delta_f.get(k).get(l), gamma_f.get(k));
		    		
		    		for (int i = 0; i < I.number_of_items; i++)
			    	 {
		    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
		    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
			    	 }
		    		cplex.addGe(str1, 0.);
		    	} 
	    	}
	    	
	    	//
	    	IloLinearNumExpr str0 = cplex.linearNumExpr();
	    	for (int k = 0; k < sample_size; k++)
	    	{
	    		IloLinearNumExpr str1 = cplex.linearNumExpr();
	    		
	    		str0.addTerm(1., gamma_f.get(k));
	    		str1.addTerm(1., gamma_f.get(k));
	    		
	    		cplex.addLe(str1, 1./((1. - alpha_f) * sample_size));
	    	}
	    	
	    	cplex.addEq(str0, 1.);
    	}
    	
    	if ("Risk-neutral".equals(type_of_risk))
    	{   
    		//
	    	for (int i = 0; i < I.number_of_items; i++)
	    	 {
	    		IloLinearNumExpr str1 = cplex.linearNumExpr(-average_f.get(i));
				
	    		for (int j = 0; j < I.number_of_fc; j++)
	           	{
	       			str1.addTerm(I.F.get(j).get(i), beta1_f.get(j));
	           	}
	       		
	       		str1.addTerm(1., beta2_f.get(i));
	    		
	    		
	    		for (int k = 0; k < sample_size; k++)
	    		{
	    			str1.addTerm(-1., mu1_f.get(k).get(i));
	    			str1.addTerm(1., mu2_f.get(k).get(i));
	    		}
	        	
	    		cplex.addGe(str1, 0.);
	    	}
	    	
	    	//
	    	for (int k = 0; k < sample_size; k++)
	    	{
	    		for (int l = 0; l < number_of_sc; l++)
	            {
	    			IloLinearNumExpr str1 = cplex.linearNumExpr(delta_f.get(k).get(l)/(sample_size + 0.));
		    		
		    		for (int i = 0; i < I.number_of_items; i++)
			    	 {
		    			str1.addTerm(-support_constraints.get(l).get(i), mu1_f.get(k).get(i));
		    			str1.addTerm(support_constraints.get(l).get(i), mu2_f.get(k).get(i));
			    	 }
		    		cplex.addGe(str1, 0.);
		    	}
	    	}
    	}
    	
    	//
    	IloLinearNumExpr str3 = cplex.linearNumExpr();
    	for (int k = 0; k < sample_size; k++)
    	{
    		for (int i = 0; i < I.number_of_items; i++)
	    	 {
    			str3.addTerm(1., mu1_f.get(k).get(i));
    			str3.addTerm(1., mu2_f.get(k).get(i));
	    	 }
    	}
    	
    	if ("Risk-averse".equals(type_of_risk))
    		cplex.addLe(str3, epsilon_f/(1. - alpha_f));
    	
    	if ("Risk-neutral".equals(type_of_risk))
    		cplex.addLe(str3, epsilon_f);
    	
    	ArrayList<Double> y_opt = new ArrayList<Double>();
    	
    	if (cplex.solve())
    	{   		
    	    y_opt.add(cplex.getObjValue());
    	}
    	else 
        	System.out.println("No solution found");
    	
    	return y_opt;
	}
    catch (IloException Exc) { Exc.printStackTrace(); }
    return new ArrayList<Double>();
}

	public static ArrayList<Double> Round(ArrayList<Double> x)
	{
		ArrayList<Double> x_round = InitializeDoubleVector(x.size(), 0.);
		
		for (int i = 0; i < x.size(); i++)
			x_round.set(i, (double)(Math.round(x.get(i))));
		return x_round;
	}
	
}
