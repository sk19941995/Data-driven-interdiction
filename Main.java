import java.io.IOException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


import ilog.concert.IloException;
import ilog.cplex.IloCplex;
public class Main {

	  public static void main(String[] args) throws IloException, IOException {
		  
			IloCplex cplex = new IloCplex();
			IloCplex cplex2 = new IloCplex();
		            
	        int sample_size_all = 1000; //sample size for computing exact CVaR
	        
	    	int k_l = -1, k_f = -1; //sample sizes
	    	double const_l = -1., const_f = -1.;
	    	double epsilon_l = -1., epsilon_f = -1.; //Wasserstein radii
	        double M1 = -1;
	        
	        int num_of_experiment = -1;
	        int binary_indicator = 0;
	        
//***********************************************************************************
	        //Generate test instances:
	        Map<Integer, Instance> random_instances = new HashMap<Integer, Instance>();
	        Map<Integer, Instance> random_instances_u = new HashMap<Integer, Instance>();
	        Map<Integer, Instance> random_instances_inv = new HashMap<Integer, Instance>();
	        
	        
	        long globalSeed = 12345L; 
	        Random globalRandom = new Random(globalSeed);
	        
	        double alpha_l = 0.95, alpha_f = 0.95;
	        
	        int n = 10;
	        int p = 1;
	        int d_f = 10;
	        
	        double w_l = 0.4, w_f = 0.4;
	        
	        int number_of_instances_1 = 10;
	        int number_of_instances_2 = 10;
	           
	  
	        for (int i = 0; i < number_of_instances_1; i++)
	        {
	        	long instanceSeed = globalRandom.nextLong();
	        	Random random = new Random(instanceSeed);
	        
	        	Instance I = new Instance(n, p, d_f, w_l, w_f, 0, random); 
	        	Instance I_inv = new Instance(n, p, d_f, w_l, w_f, 1, random);
	        	Instance I_u = new Instance(n, p, w_l, w_f, random);
	        	random_instances.put(i, I);
	        	random_instances_inv.put(i, I_inv);
	        	random_instances_u.put(i, I_u);
	        }
	        
	        
//***********************************************************************************
	        //Experiment 0 (solve test instances)
	        //num_of_experiment = 0;
	        
	        if (num_of_experiment == 0)
	        {
	        	System.out.println("Experiment 0: \n");
	        	
	        	Map<ArrayList<Integer>, Double> times = new HashMap<ArrayList<Integer>, Double>();
			       
	        	for (int t1 = 0; t1 < number_of_instances_1; t1++)
			        {
				        Instance I = random_instances.get(t1);
				        
				        ArrayList<Integer> instance = new ArrayList<Integer>();
			        	instance.add(t1);
				        
				        double start = System.currentTimeMillis();
						
			            ArrayList<ArrayList<Double>> sol = Instance.SolveLeadersProblemNominal(cplex, I, I.number_of_items + 1);
					     
			            double finish = System.currentTimeMillis();
							
			            System.out.println(sol);
			            times.put(instance, (finish-start)*Math.pow(10, -3));
			        }
	        	System.out.println(times);
	        }
	        	 
		 
//***********************************************************************************	         
        //Experiment 1 (convergence follower)
        //num_of_experiment = 1;
        if (num_of_experiment == 1)
        {
        	System.out.println("Experiment 1: \n");
        	int error_objective_num = 0;
        	
        	int[] sample_sizes = {5, 10, 15, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200, 250, 300};
        	k_l = 300;
        	
	        const_l = 0.1; const_f = 0.5;
	        epsilon_l = const_l/Math.sqrt(k_l);
	        
	        Map<ArrayList<Integer>, Double> relative_losses1 = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> relative_losses2 = new HashMap<ArrayList<Integer>, Double>();
	        
	        Map<ArrayList<Integer>, Double> times1 = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> times2 = new HashMap<ArrayList<Integer>, Double>();
	        
	        
	        for (int t1 = 0; t1 < number_of_instances_1; t1++)
	        {
	        
		        Instance I = random_instances.get(t1);
		    
		        //Generate support constraints
		        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
		        int number_of_sc = support_constraints.size() - 1; 
		    	     
		        for (int t2 = 0; t2 < number_of_instances_2; t2++)
		        {    
			        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
		        	System.out.println(number_of_instance);	
			        
		        	//Generate data
			        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
			        ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
			         
			        //Compute delta_l and average_l
			        ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);   
			        ArrayList<Double> average_l = new ArrayList<Double>();
			        
			        for (int i = 0; i < I.number_of_items; i++)
			    	{
			        	average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
			    	}
			        
			        //Main parameter loop
			        int num_param = 1;
			        for (int l = 0; l < sample_sizes.length; l++)
			        {
			        	 k_f = sample_sizes[l];
			        	 ArrayList<Integer> instance = new ArrayList<Integer>();
			        	 instance.add(number_of_instance);
			             instance.add(num_param);
			             
			             System.out.println(num_param);
			             num_param ++;
			             
			             epsilon_f = const_f/Math.sqrt(k_f);
			             
			             ArrayList<ArrayList<Double>> data_set_f = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_f); 
			             ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l, k_f, number_of_sc);      
			             
			             ArrayList<Double> average_f = new ArrayList<Double>();      
			             ArrayList<Double> average_all = new ArrayList<Double>();
			             
			             for (int i = 0; i < I.number_of_items; i++)
					    	{
					        	average_f.add(Instance.AverageOfArray(data_set_f.get(i)));	
					        	average_all.add(Instance.AverageOfArray(data_set_all.get(i)));	
					    	} 
			             
			             /*ArrayList<Double> x_check = Instance.InitializeDoubleVector(I.number_of_items, 0.);
			             
			             ArrayList<Double> sol_1 = Instance.SolveFollowersProblem(cplex, I, x_check, epsilon_f, alpha_f, k_f, data_set_f, support_constraints, "Risk-averse", delta_f, average_f);
			             ArrayList<Double> sol_2 = Instance.SolveFollowersProblem2(cplex, I, x_check, epsilon_f, alpha_f, k_f, data_set_f, support_constraints, "Risk-averse", delta_f, average_f);
				           
			             System.out.println(sol_1.get(I.number_of_items));
			             System.out.println(sol_2.get(0));*/
			             
			             M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
			             
			             double start = System.currentTimeMillis();
							
			             ArrayList<ArrayList<Double>> sol1 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
					    			support_constraints, "Risk-averse", "Risk-averse",  delta_l, delta_f, average_l, average_f, M1);
					     
			             double finish = System.currentTimeMillis();
							
			             times1.put(instance, (finish-start)*Math.pow(10, -3));
			             
			             start = System.currentTimeMillis();
							
			             ArrayList<ArrayList<Double>> sol2 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
					    			support_constraints, "Risk-averse", "Risk-neutral",  delta_l, delta_f, average_l, average_f, M1);
					     
			             finish = System.currentTimeMillis();
							
			             times2.put(instance, (finish-start)*Math.pow(10, -3));
			             
			             ArrayList<Double> x_opt1 = Instance.Round(Instance.CopyDoubleVector1(sol1.get(0), I.number_of_items)), y_opt1 = sol1.get(1);
			             ArrayList<Double> x_opt2 = Instance.Round(Instance.CopyDoubleVector1(sol2.get(0), I.number_of_items)), y_opt2 = sol2.get(1);				             
			             		     
			             double numerator1 = Instance.SolveFollowersProblemFullInformationFixedY(cplex, I, x_opt1, y_opt1, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-averse", average_all, binary_indicator);
			             double numerator2 = Instance.SolveFollowersProblemFullInformationFixedY(cplex, I, x_opt2, y_opt2, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-neutral", average_all, binary_indicator);
			             
			             ArrayList<Double> full_information_solution1 = Instance.SolveFollowersProblemFullInformation(cplex, I, x_opt1, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-averse", average_all, binary_indicator);
			             ArrayList<Double> full_information_solution2 = Instance.SolveFollowersProblemFullInformation(cplex, I, x_opt2, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-neutral", average_all, binary_indicator);
			             
			             
			             double denominator1 = full_information_solution1.get(I.number_of_items);
			             double denominator2 = full_information_solution2.get(I.number_of_items);
			            
			             if ((numerator1 > denominator1 + Math.pow(10, -3)) || (numerator2 > denominator2 + Math.pow(10, -3)))
	            			{
	            				error_objective_num ++;
	            			}
			             
			             relative_losses1.put(instance, numerator1/denominator1);    	
			             relative_losses2.put(instance, numerator2/denominator2);    
			        }
		        }
		    }          
	      
	        //Compute mean and variance
	        ArrayList<Double> mean_loss1 = new ArrayList<Double>();
	        ArrayList<Double> var_loss1 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_loss2 = new ArrayList<Double>();
	        ArrayList<Double> var_loss2 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_time1 = new ArrayList<Double>();
	        ArrayList<Double> var_time1 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_time2 = new ArrayList<Double>();
	        ArrayList<Double> var_time2 = new ArrayList<Double>();
	        
	        int num_param = 1;
	        for (int l = 0; l < sample_sizes.length; l++)
	        {
	        	k_f = sample_sizes[l];
	        	ArrayList<Double> relative_losses_param1 = new ArrayList<Double>();
	        	ArrayList<Double> times_param1 = new ArrayList<Double>();
	        	
	        	ArrayList<Double> relative_losses_param2 = new ArrayList<Double>();
	        	ArrayList<Double> times_param2 = new ArrayList<Double>();
	        	
	        	
	        	for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
	        	{
	        		ArrayList<Integer> instance = new ArrayList<Integer>();
		            instance.add(t);
		            instance.add(num_param);
		            
		            relative_losses_param1.add(relative_losses1.get(instance));    
		            times_param1.add(times1.get(instance));  
		            
		            relative_losses_param2.add(relative_losses2.get(instance));    
		            times_param2.add(times2.get(instance)); 
	        	}
	        	
	        	double mean1 = Instance.AverageOfArray(relative_losses_param1);
	        	double var1 = Instance.VarianceOfArray(relative_losses_param1);
	        	
	        	double mean2 = Instance.AverageOfArray(relative_losses_param2);
	        	double var2 = Instance.VarianceOfArray(relative_losses_param2);
	        	
	        	double meant1 = Instance.AverageOfArray(times_param1);
	        	double vart1 = Instance.VarianceOfArray(times_param1);
	        	
	        	double meant2 = Instance.AverageOfArray(times_param2);
	        	double vart2 = Instance.VarianceOfArray(times_param2);
	        	
	        	mean_loss1.add(mean1);
	        	var_loss1.add(var1); 
	        	
	        	mean_loss2.add(mean2);
	        	var_loss2.add(var2); 
	        	
	        	mean_time1.add(meant1);
	        	var_time1.add(vart1);
	        	
	        	mean_time2.add(meant2);
	        	var_time2.add(vart2);
	        	
	        	num_param ++;
	        }
	        
	        //Output
	        FileWriter fileWriter = new FileWriter("Experiment_1.txt");
	       
	        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
	        {
	        	printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
	        	printWriter.println();
	        	
	        	printWriter.printf("Parameters: n = %d, k_l = %d, const_l = %.2f, const_f = %.2f \n", n, k_l, const_l, const_f); 
	        	printWriter.printf("Number of objective errors: %d \n", error_objective_num); 
	        	
	        	printWriter.println();
	        	
	        	num_param = 1;
				
	        	 for (int l = 0; l < sample_sizes.length; l++)
			    {
			        k_f = sample_sizes[l];
					printWriter.printf("k_f = %d : %.3f (%.2f) %.2f (%.2f) | %.3f (%.2f) %.2f (%.2f) \n", k_f, mean_loss1.get(num_param - 1), var_loss1.get(num_param - 1),
						 mean_time1.get(num_param - 1), var_time1.get(num_param - 1), mean_loss2.get(num_param - 1), var_loss2.get(num_param - 1),
						 mean_time2.get(num_param - 1), var_time2.get(num_param - 1));			        
					num_param++; 
				}
				  
				printWriter.flush();
				printWriter.close();
			}
	         
        }
      
//***********************************************************************************	         	        
        //Experiment 2 (convergence leader)        
        //num_of_experiment = 2;
        if (num_of_experiment == 2)
        { 
        	System.out.println("Experiment 2: \n");
        	int error_objective_num = 0;
		    
	      
	        //Leader's convergence depends on epsilon_f (the faster the larger epsilon_f is)
        	int[] sample_sizes = {10, 15, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200, 250, 300};
        	k_f = 10;
        	int k_l_max = 300;
        	
	        const_l = 0.5; const_f = 0.1;
	        epsilon_f = const_f/Math.sqrt(k_f);
	       
	        
	        Map<ArrayList<Integer>, Double> relative_losses1 = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> relative_losses2 = new HashMap<ArrayList<Integer>, Double>();
	        
	        Map<ArrayList<Integer>, Double> times1 = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> times2 = new HashMap<ArrayList<Integer>, Double>();
	        
	        for (int t1 = 0; t1 < number_of_instances_1; t1++)
	        {
	        
		        Instance I = random_instances.get(t1);
		        
		        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
		        int number_of_sc = support_constraints.size() - 1; 
		        
		        for (int t2 = 0; t2 < number_of_instances_2; t2++)
		        {    
			        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
		        	System.out.println(number_of_instance);	
			        
			        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
			        
			        ArrayList<ArrayList<Double>> data_set_l_full = Instance.GenerateData(I, k_l_max); 
			        ArrayList<ArrayList<Double>> delta_l_full = Instance.ComputeDelta(I, data_set_l_full, support_constraints);
			           
			        ArrayList<ArrayList<Double>> data_set_f = Instance.TakeSubsetColumns(data_set_l_full, I.number_of_items, k_f); 
			        ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l_full, k_f, number_of_sc);  
		              
			        ArrayList<Double> average_all = new ArrayList<Double>(); 
			        ArrayList<Double> average_f = new ArrayList<Double>();      
				     
			        for (int i = 0; i < I.number_of_items; i++)
			    	{
			    	   average_all.add(Instance.AverageOfArray(data_set_all.get(i)));	
			    	   average_f.add(Instance.AverageOfArray(data_set_f.get(i)));	
			    	} 
			        
			        M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
                    
			        ArrayList<ArrayList<Double>> full_information_solution1 = Instance.SolveLeadersProblemFull(cplex, I, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f,
    			    		support_constraints, "Risk-averse", "Risk-averse", delta_f, average_all, average_f, M1);
    			    
                    ArrayList<ArrayList<Double>> full_information_solution2 = Instance.SolveLeadersProblemFull(cplex, I, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f,
    			    		support_constraints, "Risk-neutral", "Risk-averse", delta_f, average_all, average_f, M1);
    			       
                    
			        //Main parameter loop
                    int num_param = 1;
                    for (int l = 0; l < sample_sizes.length; l++)
    		        {
    		        	 k_l = sample_sizes[l];
			        	 ArrayList<Integer> instance = new ArrayList<Integer>();
			             instance.add(number_of_instance);
			             instance.add(num_param);
			             
			             System.out.println(num_param);
			             num_param ++; 
			             
			             epsilon_l = const_l/Math.sqrt(k_l);
			             
			             ArrayList<ArrayList<Double>> data_set_l = Instance.TakeSubsetColumns(data_set_l_full, I.number_of_items, k_l); 
			             ArrayList<ArrayList<Double>> delta_l = Instance.TakeSubsetRows(delta_l_full, k_l, number_of_sc);      
			             ArrayList<Double> average_l = new ArrayList<Double>();      
			             
			             for (int i = 0; i < I.number_of_items; i++)
					    	{
					        	average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
					    	} 
			             
			             //M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
			             
			             double start = System.currentTimeMillis();
							   
			             ArrayList<ArrayList<Double>> sol1 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
					    			support_constraints, "Risk-averse", "Risk-averse", delta_l, delta_f, average_l, average_f, M1);
					     
			             double finish = System.currentTimeMillis();
						
		             
                         times1.put(instance, (finish-start)*Math.pow(10, -3));
			             
			             start = System.currentTimeMillis();
							
			             ArrayList<ArrayList<Double>> sol2 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
					    			support_constraints, "Risk-neutral", "Risk-averse",  delta_l, delta_f, average_l, average_f, M1);
					     
			             finish = System.currentTimeMillis();
							
			             times2.put(instance, (finish-start)*Math.pow(10, -3));
			             
			  
			             ArrayList<Double> x_opt1 = Instance.Round(Instance.CopyDoubleVector1(sol1.get(0), I.number_of_items));
			             ArrayList<Double> x_opt2 = Instance.Round(Instance.CopyDoubleVector1(sol2.get(0), I.number_of_items));
			             
			             
                         /*ArrayList<Double> solf = Instance.SolveFollowersProblem(cplex, I, x_opt1, epsilon_f, alpha_f, k_f, data_set_f,  support_constraints, "Risk-averse", delta_f, average_f, binary_indicator);
                         ArrayList<Double> solf2 = Instance.SolveFollowersProblem(cplex, I, x_opt1, 100., alpha_f, k_f, data_set_f,  support_constraints, "Risk-averse", delta_f, average_f, binary_indicator);
                         
                         
			             System.out.println(solf.get(I.number_of_items));
			             System.out.println(solf2.get(I.number_of_items));*/
			             
			             ArrayList<ArrayList<Double>> sol_x_opt1 = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_opt1, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f,
  				  				support_constraints, "Risk-averse", "Risk-averse", delta_f, average_all, average_f);
  				         
			             ArrayList<ArrayList<Double>> sol_x_opt2 = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_opt2, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f,
	  				  				support_constraints, "Risk-neutral", "Risk-averse", delta_f, average_all, average_f);
	  				           	
			             
			             double numerator1 = sol_x_opt1.get(0).get(0);
			             double numerator2 = sol_x_opt2.get(0).get(0);
			            
			             double denominator1 = full_information_solution1.get(0).get(I.number_of_items);
			             double denominator2 = full_information_solution2.get(0).get(I.number_of_items);
			             
			             if ((numerator1 < denominator1 - Math.pow(10, -3)) || (numerator2 < denominator2 - Math.pow(10, -3)))
	            			{
	            				error_objective_num ++;
	            			}
			             
			             relative_losses1.put(instance, numerator1/denominator1);    	
			             relative_losses2.put(instance, numerator2/denominator2);    
			        }
		        }
		    }          
	      
	        //Compute mean and variance
	        ArrayList<Double> mean_loss1 = new ArrayList<Double>();
	        ArrayList<Double> var_loss1 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_loss2 = new ArrayList<Double>();
	        ArrayList<Double> var_loss2 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_time1 = new ArrayList<Double>();
	        ArrayList<Double> var_time1 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_time2 = new ArrayList<Double>();
	        ArrayList<Double> var_time2 = new ArrayList<Double>();
	        
	        int num_param = 1;
	        for (int l = 0; l < sample_sizes.length; l++)
	        {
	        	k_l = sample_sizes[l];
	        	ArrayList<Double> relative_losses_param1 = new ArrayList<Double>();
	        	ArrayList<Double> times_param1 = new ArrayList<Double>();
	        	
	        	ArrayList<Double> relative_losses_param2 = new ArrayList<Double>();
	        	ArrayList<Double> times_param2 = new ArrayList<Double>();
	        	
	        	
	        	for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
	        	{
	        		ArrayList<Integer> instance = new ArrayList<Integer>();
		            instance.add(t);
		            instance.add(num_param);
		            
		            relative_losses_param1.add(relative_losses1.get(instance));    
		            times_param1.add(times1.get(instance));  
		            
		            relative_losses_param2.add(relative_losses2.get(instance));    
		            times_param2.add(times2.get(instance)); 
	        	}
	        	
	        	double mean1 = Instance.AverageOfArray(relative_losses_param1);
	        	double var1 = Instance.VarianceOfArray(relative_losses_param1);
	        	
	        	double mean2 = Instance.AverageOfArray(relative_losses_param2);
	        	double var2 = Instance.VarianceOfArray(relative_losses_param2);
	        	
	        	double meant1 = Instance.AverageOfArray(times_param1);
	        	double vart1 = Instance.VarianceOfArray(times_param1);
	        	
	        	double meant2 = Instance.AverageOfArray(times_param2);
	        	double vart2 = Instance.VarianceOfArray(times_param2);
	        	
	        	mean_loss1.add(mean1);
	        	var_loss1.add(var1); 
	        	
	        	mean_loss2.add(mean2);
	        	var_loss2.add(var2); 
	        	
	        	mean_time1.add(meant1);
	        	var_time1.add(vart1);
	        	
	        	mean_time2.add(meant2);
	        	var_time2.add(vart2);
	        	
	        	num_param ++;
	        }
	        
	        FileWriter fileWriter = new FileWriter("Experiment_2.txt");
		       
	        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
	        {
	        	printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
	        	printWriter.println();
	        	
	        	printWriter.printf("Parameters: n = %d, k_f = %d, const_l = %.2f, const_f = %.2f \n", n, k_f, const_l, const_f); 
	        	printWriter.printf("Number of objective errors: %d \n", error_objective_num); 
	        	
	        	printWriter.println();
	        	
	        	num_param = 1;
				
	        	 for (int l = 0; l < sample_sizes.length; l++)
			     {
			        k_l = sample_sizes[l];
			        printWriter.printf("k_l = %d : %.3f (%.2f) %.2f (%.2f) | %.3f (%.2f) %.2f (%.2f) \n", k_l, mean_loss1.get(num_param - 1), var_loss1.get(num_param - 1),
							 mean_time1.get(num_param - 1), var_time1.get(num_param - 1), mean_loss2.get(num_param - 1), var_loss2.get(num_param - 1),
							 mean_time2.get(num_param - 1), var_time2.get(num_param - 1));		
			        num_param++;
				 }
				  
				printWriter.flush();
				printWriter.close();
			}
        }
	        
//***********************************************************************************	         	   
        //Experiment 3 (epsilon follower)
        //num_of_experiment = 3;
        if (num_of_experiment == 3)
        {
        	System.out.println("Experiment 3 (epsilon follower): \n");
        	int error_objective_num = 0;
        	
	        k_l = 30;
	        k_f = 30;
	        const_l = 0.1;
	        epsilon_l = const_l/Math.sqrt(k_l);
	        
	        double[] gammas = {0.0001, 0.0002, 0.0003, 0.0004, 0.0007, 0.0011, 0.0018, 0.0029, 0.0046, 0.0075,
	        		0.0121, 0.0196, 0.0316, 0.0511, 0.0825, 0.1, 0.1334, 0.2154, 0.3481, 0.5623,
	        		0.9085, 1.4678, 2.3714, 3.8312, 6.1897, 10.0};
	        
	        Map<ArrayList<Integer>, Double> relative_losses1 = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> relative_losses2 = new HashMap<ArrayList<Integer>, Double>();
	        
	        Map<ArrayList<Integer>, Double> times1 = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> times2 = new HashMap<ArrayList<Integer>, Double>();
	        
	        for (int t1 = 0; t1 < number_of_instances_1; t1++)
	        {
	        
		        Instance I = random_instances.get(t1);
		        
		        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
		        int number_of_sc = support_constraints.size() - 1; 
		       
		        for (int t2 = 0; t2 < number_of_instances_2; t2++)
		        {    
			        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
		        	System.out.println(number_of_instance);	
			        
			        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
			      
			        ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
			        ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
			        ArrayList<Double> average_l = new ArrayList<Double>();
			        
			        for (int i = 0; i < I.number_of_items; i++)
			    	{
			        	average_l.add(Instance.AverageOfArray(data_set_l.get(i)));		
			    	}
			        
			        ArrayList<ArrayList<Double>> data_set_f = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_f); 
			        ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l, k_f, number_of_sc);      
			        ArrayList<Double> average_f = new ArrayList<Double>();      
			        ArrayList<Double> average_all = new ArrayList<Double>();     
			        
			        for (int i = 0; i < I.number_of_items; i++)
			    	{
			        	average_f.add(Instance.AverageOfArray(data_set_f.get(i)));	
			        	average_all.add(Instance.AverageOfArray(data_set_all.get(i)));
			    	} 
			        
			        //Main parameter loop
			        int num_param = 1;
			        for (int l = 0; l < gammas.length; l++)
				     {
				         double gamma = gammas[l];
			        	 ArrayList<Integer> instance = new ArrayList<Integer>();
			             instance.add(number_of_instance);
			             instance.add(num_param);
			             
			             System.out.println(num_param);
			             num_param ++;
			             
			             epsilon_f = gamma/Math.sqrt(k_f);
				        	
                         M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
			             System.out.println(M1);
			             
			             double start = System.currentTimeMillis();
							
			             ArrayList<ArrayList<Double>> sol1 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
					    			support_constraints, "Risk-averse", "Risk-averse",  delta_l, delta_f, average_l, average_f, M1);
					     
			             double finish = System.currentTimeMillis();
							
			             times1.put(instance, (finish-start)*Math.pow(10, -3));
			             
			             start = System.currentTimeMillis();
							
			             ArrayList<ArrayList<Double>> sol2 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
					    			support_constraints, "Risk-averse", "Risk-neutral",  delta_l, delta_f, average_l, average_f, M1);
					     
			             finish = System.currentTimeMillis();
							
			             times2.put(instance, (finish-start)*Math.pow(10, -3));
			             
			             ArrayList<Double> x_opt1 = Instance.Round(Instance.CopyDoubleVector1(sol1.get(0), I.number_of_items)), y_opt1 = sol1.get(1);
			             ArrayList<Double> x_opt2 = Instance.Round(Instance.CopyDoubleVector1(sol2.get(0), I.number_of_items)), y_opt2 = sol2.get(1);				             
			             		     
			             double numerator1 = Instance.SolveFollowersProblemFullInformationFixedY(cplex, I, x_opt1, y_opt1, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-averse", average_all, binary_indicator);
			             double numerator2 = Instance.SolveFollowersProblemFullInformationFixedY(cplex, I, x_opt2, y_opt2, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-neutral", average_all, binary_indicator);
			             
			             ArrayList<Double> full_information_solution1 = Instance.SolveFollowersProblemFullInformation(cplex, I, x_opt1, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-averse", average_all, binary_indicator);
			             ArrayList<Double> full_information_solution2 = Instance.SolveFollowersProblemFullInformation(cplex, I, x_opt2, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-neutral", average_all, binary_indicator);
			             
			             
			             double denominator1 = full_information_solution1.get(I.number_of_items);
			             double denominator2 = full_information_solution2.get(I.number_of_items);
			            
			             if ((numerator1 > denominator1 + Math.pow(10, -3)) || (numerator2 > denominator2 + Math.pow(10, -3)))
	            			{
	            				error_objective_num ++;
	            			}
			             
			             relative_losses1.put(instance, numerator1/denominator1);    	
			             relative_losses2.put(instance, numerator2/denominator2);   
			        }
		        }
		    }          
	      
	      //Compute mean and variance
	        ArrayList<Double> mean_loss1 = new ArrayList<Double>();
	        ArrayList<Double> var_loss1 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_loss2 = new ArrayList<Double>();
	        ArrayList<Double> var_loss2 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_time1 = new ArrayList<Double>();
	        ArrayList<Double> var_time1 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_time2 = new ArrayList<Double>();
	        ArrayList<Double> var_time2 = new ArrayList<Double>();
	        
	        int num_param = 1;
	        for (int l = 0; l < gammas.length; l++)
	        {
	        	ArrayList<Double> relative_losses_param1 = new ArrayList<Double>();
	        	ArrayList<Double> times_param1 = new ArrayList<Double>();
	        	
	        	ArrayList<Double> relative_losses_param2 = new ArrayList<Double>();
	        	ArrayList<Double> times_param2 = new ArrayList<Double>();
	        	
	        	
	        	for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
	        	{
	        		ArrayList<Integer> instance = new ArrayList<Integer>();
		            instance.add(t);
		            instance.add(num_param);
		            
		            relative_losses_param1.add(relative_losses1.get(instance));    
		            times_param1.add(times1.get(instance));  
		            
		            relative_losses_param2.add(relative_losses2.get(instance));    
		            times_param2.add(times2.get(instance)); 
	        	}
	        	
	        	double mean1 = Instance.AverageOfArray(relative_losses_param1);
	        	double var1 = Instance.VarianceOfArray(relative_losses_param1);
	        	
	        	double mean2 = Instance.AverageOfArray(relative_losses_param2);
	        	double var2 = Instance.VarianceOfArray(relative_losses_param2);
	        	
	        	double meant1 = Instance.AverageOfArray(times_param1);
	        	double vart1 = Instance.VarianceOfArray(times_param1);
	        	
	        	double meant2 = Instance.AverageOfArray(times_param2);
	        	double vart2 = Instance.VarianceOfArray(times_param2);
	        	
	        	mean_loss1.add(mean1);
	        	var_loss1.add(var1); 
	        	
	        	mean_loss2.add(mean2);
	        	var_loss2.add(var2); 
	        	
	        	mean_time1.add(meant1);
	        	var_time1.add(vart1);
	        	
	        	mean_time2.add(meant2);
	        	var_time2.add(vart2);
	        	
	        	num_param ++;
	        }
	        
	        FileWriter fileWriter = new FileWriter("Experiment_3.txt");
		       
	        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
	        {
	        	printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
	        	printWriter.println();
	        	
	        	printWriter.printf("Parameters: n = %d, k_l = %d, k_f = %d, const_l = %.2f \n", n, k_l, k_f, const_l); 
	        	printWriter.printf("Number of objective errors: %d \n", error_objective_num); 
	        	
	        	printWriter.println();
	        	
	        	num_param = 0;
				
	        	for (int l = 0; l < gammas.length; l++)
				    {
	        		    double gamma = gammas[l];
			        	num_param ++;
			        	printWriter.printf("gamma = %.4f : %.2f (%.2f) %.2f (%.2f) | %.2f (%.2f) %.2f (%.2f) \n", gamma, mean_loss1.get(num_param - 1), var_loss1.get(num_param - 1),
								 mean_time1.get(num_param - 1), var_time1.get(num_param - 1), mean_loss2.get(num_param - 1), var_loss2.get(num_param - 1),
								 mean_time2.get(num_param - 1), var_time2.get(num_param - 1));		
			        }
				  
				printWriter.flush();
				printWriter.close();
			}
        }
    
	        
//***********************************************************************************	         	    
      //Experiment 4 (epsilon leader)
      //num_of_experiment = 4;   
      if (num_of_experiment == 4)
  	        { 
  	        	System.out.println("Experiment 4: \n");
  			    int error_objective_num = 0;
  			
  		        k_l = 30;
		        k_f = 30;
		        const_f = 0.1;
		        epsilon_f = const_f/Math.sqrt(k_f);
		           
		        double[] gammas = {0.0001, 0.0002, 0.0003, 0.0004, 0.0007, 0.0011, 0.0018, 0.0029, 0.0046, 0.0075,
		        		0.0121, 0.0196, 0.0316, 0.0511, 0.0825, 0.1, 0.1334, 0.2154, 0.3481, 0.5623,
		        		0.9085, 1.4678, 2.3714, 3.8312, 6.1897, 10.0};
		        //double[] gammas = {0.0};
		        
		        Map<ArrayList<Integer>, Double> relative_losses1 = new HashMap<ArrayList<Integer>, Double>();
		        Map<ArrayList<Integer>, Double> relative_losses2 = new HashMap<ArrayList<Integer>, Double>();
		        
		        Map<ArrayList<Integer>, Double> times1 = new HashMap<ArrayList<Integer>, Double>();
		        Map<ArrayList<Integer>, Double> times2 = new HashMap<ArrayList<Integer>, Double>();
		        
  		        for (int t1 = 0; t1 < number_of_instances_1; t1++)
  		        {
  		        
  			        Instance I = random_instances.get(t1);
  			        
  			        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
			        int number_of_sc = support_constraints.size() - 1; 
			    	   
  			        for (int t2 = 0; t2 < number_of_instances_2; t2++)
  			        {    
  				        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
  			        	System.out.println(number_of_instance);	
  				        
  				        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
  				       
  				        ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
  				        ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
  				       
  				        //ArrayList<ArrayList<Double>> data_set_f = Instance.GenerateData(I, k_f);
  				        ArrayList<ArrayList<Double>> data_set_f = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_f); 
  				        ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l, k_f, number_of_sc);  
  			           
  				       //Compute average for leader and follower
  				       ArrayList<Double> average_l = new ArrayList<Double>(); 
  				       ArrayList<Double> average_f = new ArrayList<Double>();   
  				       ArrayList<Double> average_all = new ArrayList<Double>();      
  				       ArrayList<Double> variance_all = new ArrayList<Double>(); 
  					     
  				       for (int i = 0; i < I.number_of_items; i++)
  				    	{
  				    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
  				    	    average_f.add(Instance.AverageOfArray(data_set_f.get(i)));	
  				    	    average_all.add(Instance.AverageOfArray(data_set_all.get(i)));
  				    	    variance_all.add(Instance.VarianceOfArray(data_set_all.get(i)));	
  				    	} 
  				        	
  				       //System.out.print(average_all);
  				       //System.out.print(variance_all);
  				       M1 = Math.max(I.number_of_items, 1. + epsilon_f/(1 - alpha_f));
  				       
  				       ArrayList<ArrayList<Double>> full_information_solution1 = Instance.SolveLeadersProblemFull(cplex, I, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f,
  	    			    		support_constraints, "Risk-averse", "Risk-averse",  delta_f, average_all, average_f, M1);
  	    			   
  				       ArrayList<ArrayList<Double>> full_information_solution2 = Instance.SolveLeadersProblemFull(cplex, I, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f,
  	    			    		support_constraints, "Risk-neutral", "Risk-averse",  delta_f, average_all, average_f, M1);
  	    			       
  				       //Main parameter loop
  	                   int num_param = 1;
  	                   for (int l = 0; l < gammas.length; l++)
					    {
		        		     double gamma = gammas[l];
  				        	 ArrayList<Integer> instance = new ArrayList<Integer>();
  				             instance.add(number_of_instance);
  				             instance.add(num_param);
  				             
  				             System.out.println(num_param);
  				             num_param ++;
  				             
  				             epsilon_l = gamma/Math.sqrt(k_l);
				             
				             double start = System.currentTimeMillis();
								   
				             ArrayList<ArrayList<Double>> sol1 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
						    			support_constraints, "Risk-averse", "Risk-averse", delta_l, delta_f, average_l, average_f, M1);
						     
				             double finish = System.currentTimeMillis();
							
                             times1.put(instance, (finish-start)*Math.pow(10, -3));
				             
				             start = System.currentTimeMillis();
								
				             ArrayList<ArrayList<Double>> sol2 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
						    			support_constraints, "Risk-neutral", "Risk-averse",  delta_l, delta_f, average_l, average_f, M1);
						     
				             finish = System.currentTimeMillis();
								
				             times2.put(instance, (finish-start)*Math.pow(10, -3));
				             
				  
				             ArrayList<Double> x_opt1 = Instance.Round(Instance.CopyDoubleVector1(sol1.get(0), I.number_of_items));
				             ArrayList<Double> x_opt2 = Instance.Round(Instance.CopyDoubleVector1(sol2.get(0), I.number_of_items));
				             
				             ArrayList<ArrayList<Double>> sol_x_opt1 = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_opt1, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f,
	  				  				support_constraints, "Risk-averse", "Risk-averse", delta_f, average_all, average_f);
	  				         
				             ArrayList<ArrayList<Double>> sol_x_opt2 = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_opt2, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f,
		  				  				support_constraints, "Risk-neutral", "Risk-averse", delta_f, average_all, average_f);
		  				           	
				             
				             double numerator1 = sol_x_opt1.get(0).get(0);
				             double numerator2 = sol_x_opt2.get(0).get(0);
				            
				             double denominator1 = full_information_solution1.get(0).get(I.number_of_items);
				             double denominator2 = full_information_solution2.get(0).get(I.number_of_items);
				             
				             if ((numerator1 < denominator1 - Math.pow(10, -3)) || (numerator2 < denominator2 - Math.pow(10, -3)))
		            			{
		            				error_objective_num ++;
		            			}
				             
				             relative_losses1.put(instance, numerator1/denominator1);    	
				             relative_losses2.put(instance, numerator2/denominator2); 	        
  				        }
  			        }
  			    }          
  		       
	  		    //Compute mean and variance
		        ArrayList<Double> mean_loss1 = new ArrayList<Double>();
		        ArrayList<Double> var_loss1 = new ArrayList<Double>();
		        
		        ArrayList<Double> mean_loss2 = new ArrayList<Double>();
		        ArrayList<Double> var_loss2 = new ArrayList<Double>();
		        
		        ArrayList<Double> mean_time1 = new ArrayList<Double>();
		        ArrayList<Double> var_time1 = new ArrayList<Double>();
		        
		        ArrayList<Double> mean_time2 = new ArrayList<Double>();
		        ArrayList<Double> var_time2 = new ArrayList<Double>();
		        
		        int num_param = 1;
		        for (int l = 0; l < gammas.length; l++)
		        {
		        	ArrayList<Double> relative_losses_param1 = new ArrayList<Double>();
		        	ArrayList<Double> times_param1 = new ArrayList<Double>();
		        	
		        	ArrayList<Double> relative_losses_param2 = new ArrayList<Double>();
		        	ArrayList<Double> times_param2 = new ArrayList<Double>();
		        	
		        	
		        	for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
		        	{
		        		ArrayList<Integer> instance = new ArrayList<Integer>();
			            instance.add(t);
			            instance.add(num_param);
			            
			            relative_losses_param1.add(relative_losses1.get(instance));    
			            times_param1.add(times1.get(instance));  
			            
			            relative_losses_param2.add(relative_losses2.get(instance));    
			            times_param2.add(times2.get(instance)); 
		        	}
		        	
		        	double mean1 = Instance.AverageOfArray(relative_losses_param1);
		        	double var1 = Instance.VarianceOfArray(relative_losses_param1);
		        	
		        	double mean2 = Instance.AverageOfArray(relative_losses_param2);
		        	double var2 = Instance.VarianceOfArray(relative_losses_param2);
		        	
		        	double meant1 = Instance.AverageOfArray(times_param1);
		        	double vart1 = Instance.VarianceOfArray(times_param1);
		        	
		        	double meant2 = Instance.AverageOfArray(times_param2);
		        	double vart2 = Instance.VarianceOfArray(times_param2);
		        	
		        	mean_loss1.add(mean1);
		        	var_loss1.add(var1); 
		        	
		        	mean_loss2.add(mean2);
		        	var_loss2.add(var2); 
		        	
		        	mean_time1.add(meant1);
		        	var_time1.add(vart1);
		        	
		        	mean_time2.add(meant2);
		        	var_time2.add(vart2);
		        	
		        	num_param ++;
		        }
	  		       
	  		    FileWriter fileWriter = new FileWriter("Experiment_4.txt");
			       
		        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
		        {
		        	printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
		        	printWriter.println();
		        	
		        	printWriter.printf("Parameters: n = %d, k_l = %d, k_f = %d, const_f = %.2f \n", n, k_l, k_f, const_f); 
		        	printWriter.printf("Number of objective errors: %d \n", error_objective_num); 
		        	
		        	printWriter.println();
		        	
		        	num_param = 0;
					
		        	for (int l = 0; l < gammas.length; l++)
				    {
	        		    double gamma = gammas[l];
			        	num_param ++;
			        	printWriter.printf("gamma = %.4f : %.3f (%.3f) %.2f (%.2f) | %.3f (%.3f) %.2f (%.2f) \n", gamma, mean_loss1.get(num_param - 1), var_loss1.get(num_param - 1),
								 mean_time1.get(num_param - 1), var_time1.get(num_param - 1), mean_loss2.get(num_param - 1), var_loss2.get(num_param - 1),
								 mean_time2.get(num_param - 1), var_time2.get(num_param - 1));		
			        }
					  
					printWriter.flush();
					printWriter.close();
				}    
	  	 }
	       
      
//***********************************************************************************	         	  	        
  	    //Experiment 5 (risk preferences)
        //num_of_experiment = 5;
        if (num_of_experiment == 5)
  	        { 
  	        	System.out.println("Experiment 5: \n");
  	       
  		        const_l = 0.1; const_f = 0.1;
  		        
  		        int[] sample_sizes = {30, 300};
  		        
  		        Map<ArrayList<Integer>, Double> relative_losses = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> times = new HashMap<ArrayList<Integer>, Double>();
  		        
  		        for (int t1 = 0; t1 < number_of_instances_1; t1++)
  		        {
  		     
  			        Instance I = random_instances.get(t1);
  			        
  			        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
			        int number_of_sc = support_constraints.size() - 1; 
			    	 
  			        for (int t2 = 0; t2 < number_of_instances_2; t2++)
  			        {    
  				        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
  			        	System.out.println(number_of_instance);	
  				        
  				        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
  				       			        	                      
			  		    //Main parameter loop
  				        int num_param = 1;
  				  	for (int l = 0; l < sample_sizes.length; l++)
					    	{
					    		k_l = sample_sizes[l];
					    	    k_f = sample_sizes[l];
			  			        
					    	    epsilon_l = const_l/Math.sqrt(k_l);
					    	    epsilon_f = const_f/Math.sqrt(k_f);
					    	     
					    	    M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
					    	    
			  			        ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
		  				        ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
			           
		  				        ArrayList<ArrayList<Double>> data_set_f = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_f); 
		  				        ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l, k_f, number_of_sc);  
		  			           
		  				        ArrayList<Double> average_l = new ArrayList<Double>(); 
		  				        ArrayList<Double> average_f = new ArrayList<Double>();  
		  				        ArrayList<Double> average_all = new ArrayList<Double>(); 
		  				        ArrayList<Double> variance_all = new ArrayList<Double>(); 
		  					     
		  				        for (int i = 0; i < I.number_of_items; i++)
		  				    	{
		  				    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
		  				    	    average_f.add(Instance.AverageOfArray(data_set_f.get(i)));	
		  				    	    average_all.add(Instance.AverageOfArray(data_set_all.get(i)));	
		  				    	    variance_all.add(Instance.VarianceOfArray(data_set_all.get(i)));	
		  				    	}
		  			            
		  				        //System.out.println(average_all);
		  				        //System.out.println(variance_all);
		  				        
			  				    for (int risk = 0; risk < 4; risk++)
			  				    {	
			  				    	System.out.printf("(%d, %d) \n", sample_sizes[l], risk);
			  				    	
			  				    	String type_of_risk_l = "Risk-neutral", type_of_risk_f = "Risk-neutral";
		  				             
		  				             if ((risk == 1) || (risk == 3))
		  				             {
		  				            	type_of_risk_f = "Risk-averse";
		  				             }
		  				             
		  				             if ((risk == 2) || (risk == 3))
		  				             {
		  				            	type_of_risk_l = "Risk-averse";
		  				             }
			  				    	
		  				             double start = System.currentTimeMillis();
  								   
		  				             ArrayList<ArrayList<Double>> sol = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
		  						    			support_constraints, type_of_risk_l, type_of_risk_f, delta_l, delta_f, average_l, average_f, M1);
		  						     
		  				             double finish = System.currentTimeMillis();
		  							  
		  				             ArrayList<Double> x_opt = Instance.Round(Instance.CopyDoubleVector1(sol.get(0), I.number_of_items)), y_opt = sol.get(1);
		  				            
		  				             //System.out.printf("%d \n", risk);
		  				             //System.out.println(y_opt);
		  				             
			  				    	 for (int model = 0; model < 3; model++)   
							         {
			  				    		ArrayList<Integer> instance = new ArrayList<Integer>();
			  				            instance.add(number_of_instance);
			  				    		instance.add(num_param);
			  				            num_param ++;
			  				            
			  				            double numerator = -1;
			  				            
			  				            if (model == 0)
			  				            {
			  				            	 numerator = Instance.SolveFollowersProblemFullInformationFixedY(cplex, I, x_opt, y_opt, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-averse", average_all, binary_indicator);
			  				            }
			  				             
			  				             if (model == 2)
			  				             {
			  				            	 numerator = Instance.RightTailedCvaR(cplex, I, data_set_all, sample_size_all, alpha_l, y_opt);
			  				             }
			  				            
			  				             if (model == 1)
			  				             {
			  				            	numerator = Instance.SolveFollowersProblemFullInformationFixedY(cplex, I, x_opt, y_opt, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-neutral", average_all, binary_indicator);
			  				             }
			  				            
			  				             relative_losses.put(instance, numerator);
			  				             times.put(instance, (finish-start)*Math.pow(10, -3));
			  				        }
			  			     }
  			        }
  		        }
  		      }		      
	          ArrayList<Double> mean_loss = new ArrayList<Double>();
	          ArrayList<Double> var_loss = new ArrayList<Double>();
	        
	          ArrayList<Double> mean_time = new ArrayList<Double>();
	          ArrayList<Double> var_time = new ArrayList<Double>();
	        
  		      int num_param = 1;
  		    
  		      for (int l = 0; l < sample_sizes.length; l++)
			        for (int risk = 0; risk < 4; risk++)
		  		    	for (int model = 0; model < 3; model++)   
				        {
		  		        	ArrayList<Double> relative_losses_param = new ArrayList<Double>();
		  		        	ArrayList<Double> times_param = new ArrayList<Double>();
		  		        	
		  		        	
		  		        	for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
		  		        	{
		  		        		ArrayList<Integer> instance = new ArrayList<Integer>();
		  			            instance.add(t);
		  			            instance.add(num_param);
		  			            
		  			            relative_losses_param.add(relative_losses.get(instance));    
		  			            times_param.add(times.get(instance));  
		  		        	}
		  		        	
		  		        	double mean = Instance.AverageOfArray(relative_losses_param);
		  		        	double var = Instance.VarianceOfArray(relative_losses_param);
		  		        	
		  		        	double mean2 = Instance.AverageOfArray(times_param);
		  		        	double var2 = Instance.VarianceOfArray(times_param);
		  		        	
		  		        	mean_loss.add(mean);
		  		        	var_loss.add(var); 
		  		        	
		  		        	mean_time.add(mean2);
		  		        	var_time.add(var2);
		  		        	
		  		        	num_param ++;
		  		        }
		  		   
	  		  
	  		    FileWriter fileWriter = new FileWriter("Experiment_5.txt");
		       
		        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
		        {
		        	printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
		        	printWriter.println();
		        	
		        	printWriter.printf("Parameters: n = %d, const_l = %.2f, const_f = %.2f \n", n, const_l, const_f); 
		        	printWriter.println();
		        	
		        	num_param = 1;
		        	
		        	for (int l = 0; l < sample_sizes.length; l++)
					    {	
		        		  printWriter.printf("sample_size = %d \n", sample_sizes[l]);	 
			  		      for (int risk = 0; risk < 4; risk++)
			  		      {
			  		    	   if (risk == 0) printWriter.printf("NN: ");
		  		        	   if (risk == 1) printWriter.printf("NA: ");
		  		        	   if (risk == 2) printWriter.printf("AN: ");
		  		        	   if (risk == 3) printWriter.printf("AA: ");
		  		        	   
			  		           for (int model = 0; model < 3; model++)   
						       {    
			  		        	    printWriter.printf("%.2f (%.2f) ", mean_loss.get(num_param - 1), var_loss.get(num_param - 1));
			  		        	    num_param ++;
						       }
			  		         printWriter.println();
			  		      }
			  		    printWriter.println();
			  	    }
					
					printWriter.flush();
					printWriter.close();
				}     
  	        }

//***********************************************************************************
	    //Experiment 5 (risk preferences)
	    //num_of_experiment = 51;
        if (num_of_experiment == 51)
  	        { 
  	        	System.out.println("Experiment 5a: \n");
  	       
  		        const_l = 0.1; const_f = 0.1;
  		        
  		        int[] sample_sizes = {30, 300};
  		        
  		        Map<ArrayList<Integer>, Double> relative_losses = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> times = new HashMap<ArrayList<Integer>, Double>();
  		        
  		        for (int t1 = 0; t1 < number_of_instances_1; t1++)
  		        {
  		     
  			        Instance I = random_instances_inv.get(t1);
  			        
  			        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
			        int number_of_sc = support_constraints.size() - 1; 
			    	 
  			        for (int t2 = 0; t2 < number_of_instances_2; t2++)
  			        {    
  				        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
  			        	System.out.println(number_of_instance);	
  				        
  				        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
  				       			        	                      
			  		    //Main parameter loop
  				        int num_param = 1;
  				  	for (int l = 0; l < sample_sizes.length; l++)
					    	{
					    		k_l = sample_sizes[l];
					    	    k_f = sample_sizes[l];
			  			        
					    	    epsilon_l = const_l/Math.sqrt(k_l);
					    	    epsilon_f = const_f/Math.sqrt(k_f);
					    	     
					    	    M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
					    	    
			  			        ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
		  				        ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
			           
		  				        ArrayList<ArrayList<Double>> data_set_f = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_f); 
		  				        ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l, k_f, number_of_sc);  
		  			           
		  				        ArrayList<Double> average_l = new ArrayList<Double>(); 
		  				        ArrayList<Double> average_f = new ArrayList<Double>();  
		  				        ArrayList<Double> average_all = new ArrayList<Double>(); 
		  				        ArrayList<Double> variance_all = new ArrayList<Double>(); 
		  					     
		  				        for (int i = 0; i < I.number_of_items; i++)
		  				    	{
		  				    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
		  				    	    average_f.add(Instance.AverageOfArray(data_set_f.get(i)));	
		  				    	    average_all.add(Instance.AverageOfArray(data_set_all.get(i)));	
		  				    	    variance_all.add(Instance.VarianceOfArray(data_set_all.get(i)));	
		  				    	}
		  			            
		  				        //System.out.println(average_all);
		  				        //System.out.println(variance_all);
		  				        
			  				    for (int risk = 0; risk < 4; risk++)
			  				    {	
			  				    	System.out.printf("(%d, %d) \n", sample_sizes[l], risk);
			  				    	
			  				    	String type_of_risk_l = "Risk-neutral", type_of_risk_f = "Risk-neutral";
		  				             
		  				             if ((risk == 1) || (risk == 3))
		  				             {
		  				            	type_of_risk_f = "Risk-averse";
		  				             }
		  				             
		  				             if ((risk == 2) || (risk == 3))
		  				             {
		  				            	type_of_risk_l = "Risk-averse";
		  				             }
			  				    	
		  				             double start = System.currentTimeMillis();
  								   
		  				             ArrayList<ArrayList<Double>> sol = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
		  						    			support_constraints, type_of_risk_l, type_of_risk_f, delta_l, delta_f, average_l, average_f, M1);
		  						     
		  				             double finish = System.currentTimeMillis();
		  							  
		  				             ArrayList<Double> x_opt = Instance.Round(Instance.CopyDoubleVector1(sol.get(0), I.number_of_items)), y_opt = sol.get(1);
		  				            
		  				             //System.out.printf("%d \n", risk);
		  				             //System.out.println(y_opt);
		  				             
			  				    	 for (int model = 0; model < 3; model++)   
							         {
			  				    		ArrayList<Integer> instance = new ArrayList<Integer>();
			  				            instance.add(number_of_instance);
			  				    		instance.add(num_param);
			  				            num_param ++;
			  				            
			  				            double numerator = -1;
			  				            
			  				            if (model == 0)
			  				            {
			  				            	 numerator = Instance.SolveFollowersProblemFullInformationFixedY(cplex, I, x_opt, y_opt, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-averse", average_all, binary_indicator);
			  				            }
			  				             
			  				             if (model == 2)
			  				             {
			  				            	 numerator = Instance.RightTailedCvaR(cplex, I, data_set_all, sample_size_all, alpha_l, y_opt);
			  				             }
			  				            
			  				             if (model == 1)
			  				             {
			  				            	numerator = Instance.SolveFollowersProblemFullInformationFixedY(cplex, I, x_opt, y_opt, alpha_f, sample_size_all, data_set_all, support_constraints, "Risk-neutral", average_all, binary_indicator);
			  				             }
			  				            
			  				             relative_losses.put(instance, numerator);
			  				             times.put(instance, (finish-start)*Math.pow(10, -3));
			  				        }
			  			     }
  			        }
  		        }
  		      }		      
	          ArrayList<Double> mean_loss = new ArrayList<Double>();
	          ArrayList<Double> var_loss = new ArrayList<Double>();
	        
	          ArrayList<Double> mean_time = new ArrayList<Double>();
	          ArrayList<Double> var_time = new ArrayList<Double>();
	        
  		      int num_param = 1;
  		    
  		      for (int l = 0; l < sample_sizes.length; l++)
			        for (int risk = 0; risk < 4; risk++)
		  		    	for (int model = 0; model < 3; model++)   
				        {
		  		        	ArrayList<Double> relative_losses_param = new ArrayList<Double>();
		  		        	ArrayList<Double> times_param = new ArrayList<Double>();
		  		        	
		  		        	
		  		        	for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
		  		        	{
		  		        		ArrayList<Integer> instance = new ArrayList<Integer>();
		  			            instance.add(t);
		  			            instance.add(num_param);
		  			            
		  			            relative_losses_param.add(relative_losses.get(instance));    
		  			            times_param.add(times.get(instance));  
		  		        	}
		  		        	
		  		        	double mean = Instance.AverageOfArray(relative_losses_param);
		  		        	double var = Instance.VarianceOfArray(relative_losses_param);
		  		        	
		  		        	double mean2 = Instance.AverageOfArray(times_param);
		  		        	double var2 = Instance.VarianceOfArray(times_param);
		  		        	
		  		        	mean_loss.add(mean);
		  		        	var_loss.add(var); 
		  		        	
		  		        	mean_time.add(mean2);
		  		        	var_time.add(var2);
		  		        	
		  		        	num_param ++;
		  		        }
		  		   
	  		  
	  		    FileWriter fileWriter = new FileWriter("Experiment_5a.txt");
		       
		        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
		        {
		        	printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
		        	printWriter.println();
		        	
		        	printWriter.printf("Parameters: n = %d, const_l = %.2f, const_f = %.2f \n", n, const_l, const_f); 
		        	printWriter.println();
		        	
		        	num_param = 1;
		        	
		        	for (int l = 0; l < sample_sizes.length; l++)
					    {	
		        		  printWriter.printf("sample_size = %d \n", sample_sizes[l]);	 
			  		      for (int risk = 0; risk < 4; risk++)
			  		      {
			  		    	   if (risk == 0) printWriter.printf("NN: ");
		  		        	   if (risk == 1) printWriter.printf("NA: ");
		  		        	   if (risk == 2) printWriter.printf("AN: ");
		  		        	   if (risk == 3) printWriter.printf("AA: ");
		  		        	   
			  		           for (int model = 0; model < 3; model++)   
						       {    
			  		        	    printWriter.printf("%.2f (%.2f) ", mean_loss.get(num_param - 1), var_loss.get(num_param - 1));
			  		        	    num_param ++;
						       }
			  		         printWriter.println();
			  		      }
			  		    printWriter.println();
			  	    }
					
					printWriter.flush();
					printWriter.close();
				}     
  	        }
	        
//***********************************************************************************
        //Experiment 6 (in-sample and out-sample analysis of the basic model and the semi-pessimistic approximations, number of scenarios selection)  
  	    //num_of_experiment = 6;
  	    if (num_of_experiment == 6)
        { 
        	System.out.println("Experiment 6: \n");
            int error_objective_num = 0;
            
	        k_l = 30;
	        k_f = 30;
	        
	        int k_lf = 20;
	        
	        const_l = 0.1; 
	        const_f = 0.1;
	        
	        epsilon_l = const_l/Math.sqrt(k_f);
	        epsilon_f = const_f/Math.sqrt(k_f);
	        
	      
	        int number_of_scenarios = 10;
	        double max_noise = 0.5;
	        
	        Map<ArrayList<Integer>, Double> relative_losses_sp_in = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> relative_losses_sp_out = new HashMap<ArrayList<Integer>, Double>();
	        
	        ArrayList<Double> relative_losses_b_out = new ArrayList<Double>();
	        ArrayList<Double> relative_losses_b_in = new ArrayList<Double>();
	        ArrayList<Double> relative_losses_true_out = new ArrayList<Double>();
	       
	        Map<ArrayList<Integer>, Double> times_sp = new HashMap<ArrayList<Integer>, Double>();
	        ArrayList<Double> times_b = new ArrayList<Double>();
	        
	        for (int t1 = 0; t1 < number_of_instances_1; t1++)
	        {
		        Instance I = random_instances.get(t1);
		        
		        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
		        int number_of_sc = support_constraints.size() - 1; 
		        
		        
		        ArrayList<ArrayList<Double>> shifts = new ArrayList<ArrayList<Double>>();
  		        ArrayList<ArrayList<Double>> missed = new ArrayList<ArrayList<Double>>();
  		        
  		        for (int i = 0; i < I.number_of_items; i++)
  		        {
  		        	shifts.add(new ArrayList<Double>());
  		        	missed.add(new ArrayList<Double>());
  		        	
  		            for (int j = 0; j < k_f; j++)
  		            {
  		            	shifts.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
  		            	missed.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
  		            }
  		        }
		        
		        for (int t2 = 0; t2 < number_of_instances_2; t2++)
		        {    
			        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
		        	System.out.println(number_of_instance);	
			        
			        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
			        	                      
	    			ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
	    			ArrayList<ArrayList<Double>> data_set_part_f = Instance.GenerateData(I, k_f - k_lf); 
	    			ArrayList<ArrayList<Double>> data_set_common = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_lf); 
	  
	    			ArrayList<ArrayList<Double>> data_set_f_true = Instance.JoinDataSets(data_set_common, data_set_part_f, I.number_of_items);
	    			
	    			ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
	    			//ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l, k_lf, number_of_sc);  
			        ArrayList<ArrayList<Double>> delta_f_true = Instance.ComputeDelta(I, data_set_f_true, support_constraints);
			           
	  				//Compute average for leader
	  				ArrayList<Double> average_l = new ArrayList<Double>(); 
	  				ArrayList<Double> average_f = new ArrayList<Double>(); 
	  				ArrayList<Double> average_f_true = new ArrayList<Double>(); 
	  				ArrayList<Double> average_all = new ArrayList<Double>(); 
	  				
	  				for (int i = 0; i < I.number_of_items; i++)
  				    {
  				    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
  				    	    average_f.add(Instance.AverageOfArray(data_set_common.get(i)));	
  				    	    average_f_true.add(Instance.AverageOfArray(data_set_f_true.get(i)));	
  				    	    average_all.add(Instance.AverageOfArray(data_set_all.get(i)));	
  				    }
	  				
	  				 //Solve true formulation
	  				M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
	  				
	  				ArrayList<ArrayList<Double>> sol_true = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f_true,
				    			support_constraints, "Risk-averse", "Risk-averse",  delta_l, delta_f_true, average_l, average_f_true, M1);
				    
	  				//Solve basic formulation
	  				double start = System.currentTimeMillis();
	  				
	  				ArrayList<ArrayList<Double>> sol_b = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_l,
				    			support_constraints, "Risk-averse", "Risk-averse",  delta_l, delta_l, average_l, average_l, M1);
				    
	  				double finish = System.currentTimeMillis();
	  				
	  				times_b.add((finish-start)*Math.pow(10, -3));
			  		
	  			    //Compute in-sample performance for the basic formulation
	  			    double z_b = sol_b.get(0).get(I.number_of_items);	
	  			    double z_true = sol_true.get(0).get(I.number_of_items);	
	  			    
	  			    relative_losses_b_in.add(z_b/z_true);
	  			    
	  			    //Compute out-of-sample performance for the basic formulations
	  			    ArrayList<Double> x_true = Instance.CopyDoubleVector1(sol_true.get(0), I.number_of_items);  
	  			    ArrayList<Double> x_b = Instance.CopyDoubleVector1(sol_b.get(0), I.number_of_items); 
	  			    
	  			    ArrayList<ArrayList<Double>> sol_x_true = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_true, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
			  				support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true);
		        
	  			    
	  			    ArrayList<ArrayList<Double>> sol_x_b = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_b, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
			  				support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true);
		        
	  			    double numerator_true = sol_x_true.get(0).get(0);
		            double numerator_b = sol_x_b.get(0).get(0);
		            
		            ArrayList<ArrayList<Double>> full_information_solution = Instance.SolveLeadersProblemFull(cplex, I, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
	    			    		support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true, M1);
	  				
	  			    double denominator = full_information_solution.get(0).get(I.number_of_items);	
	  			   
	  			    if ((numerator_b < denominator - Math.pow(10, -3)) || ((numerator_true < denominator - Math.pow(10, -3))))
	          			{
	          				error_objective_num ++;
	          			}
	  			    
	  			    relative_losses_true_out.add(numerator_true/denominator);
	  			    relative_losses_b_out.add(numerator_b/denominator);
	  			  
		  			//Main parameter loop
		  			int num_param = 1;	
		  			
		  			 for (double noise_level = 0.2; noise_level <= max_noise; noise_level += 0.3) 
		  			{
		  				Map<Integer, ArrayList<ArrayList<Double>>> data_sets_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
		  				Map<Integer, ArrayList<ArrayList<Double>>> deltas_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
		  				Map<Integer, ArrayList<Double>> averages_f = new HashMap<Integer, ArrayList<Double>>();
		  				
		  				for (int r = 0; r < number_of_scenarios; r += 1)
		  				{
		  					ArrayList<ArrayList<Double>> data_set_f_r = Instance.GenerateNewDataSet(I, data_set_f_true, I.number_of_items, k_lf, noise_level, shifts, missed);
		  					ArrayList<ArrayList<Double>> delta_f_r = Instance.ComputeDelta(I, data_set_f_r, support_constraints);  
		  					
		  					ArrayList<Double> average_f_r = new ArrayList<Double>();
		  					
		  					for (int i = 0; i < I.number_of_items; i++)
		  				    {
		  				    	    average_f_r.add(Instance.AverageOfArray(data_set_f_r.get(i)));	
		  				    }
		  					
		  					averages_f.put(r, average_f_r);
		  					data_sets_f.put(r, data_set_f_r);
		  					deltas_f.put(r, delta_f_r);
		  				}
		  				  				
		  				for (int r = 0; r < number_of_scenarios; r += 1)
		  				{    
		  					ArrayList<Integer> instance = new ArrayList<Integer>();
 				            instance.add(number_of_instance);
 				            instance.add(num_param);
 				             
 				            System.out.println(num_param);
 				             
 				            num_param ++;
 				            
 				            M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
 			  				
				  		    start = System.currentTimeMillis();
									        
				  		    ArrayList<Double> sol_sp =  Instance.SolveLeadersProblemScenarios(cplex, I, r + 1, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_sets_f,
			            			support_constraints, "Risk-averse",  "Risk-averse", delta_l,  deltas_f, average_l, averages_f, M1);
					        
				  		    finish = System.currentTimeMillis();
				  		    
				  		    times_sp.put(instance, (finish-start)*Math.pow(10, -3));
				  		
				  		    //Compute in-sample performance
				  		    double z_sp = sol_sp.get(I.number_of_items);	
				  		   
				  		    relative_losses_sp_in.put(instance, z_sp/z_true);
				  		    
				  		    //Compute out-of-sample performance 
				  		    ArrayList<Double> x_sp = Instance.CopyDoubleVector1(sol_sp, I.number_of_items); 
				  		    
				  		    ArrayList<ArrayList<Double>> sol_x_sp = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_sp, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
	  				  				support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true);
				             
				            double numerator_sp = sol_x_sp.get(0).get(0);
				           
				            if (numerator_sp < denominator - Math.pow(10, -3))
		          			{
		          				error_objective_num ++;
		          			}
				            
				  		    relative_losses_sp_out.put(instance, numerator_sp/denominator);
				  		    
		  				} 
		  			}
	  			}
	  			
	      }
		  		      
	      ArrayList<Double> min_loss_in = new ArrayList<Double>();
          
	      ArrayList<Double> mean_loss_in = new ArrayList<Double>();
          ArrayList<Double> var_loss_in = new ArrayList<Double>();
          
          ArrayList<Double> mean_loss_out = new ArrayList<Double>();
          ArrayList<Double> var_loss_out = new ArrayList<Double>();
          
          ArrayList<Double> mean_time = new ArrayList<Double>();
          ArrayList<Double> var_time = new ArrayList<Double>();
          	
          int num_param = 0;
          
          for (double noise_level = 0.2; noise_level <= max_noise; noise_level += 0.3) 
        	  for (int r = 0; r < number_of_scenarios; r += 1)
        	  { 
        		  num_param++;
        		  ArrayList<Double> relative_losses_param_in = new ArrayList<Double>();
        		  ArrayList<Double> relative_losses_param_out = new ArrayList<Double>();
  		          
  		          ArrayList<Double> times_param = new ArrayList<Double>();
  		          
  		          for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
  		        	{
  		        		ArrayList<Integer> instance = new ArrayList<Integer>();
  			            instance.add(t);
  			            instance.add(num_param);
  			            
  			            relative_losses_param_in.add(relative_losses_sp_in.get(instance));  
  			            times_param.add(times_sp.get(instance));
  			             
  			            relative_losses_param_out.add(relative_losses_sp_out.get(instance));        
  		        	}
        		  
  		         double mean_in = Instance.AverageOfArray(relative_losses_param_in);
  		         double var_in = Instance.VarianceOfArray(relative_losses_param_in);
  	  		     
  		         Collections.sort(relative_losses_param_in); 
	  		    
  		         double index = 0.05*(number_of_instances_1*number_of_instances_2 - 1);
  		         int lowerIndex = (int) Math.floor(index);
  		         int upperIndex = (int) Math.ceil(index);
  		         
  		         double min_in = relative_losses_param_in.get(lowerIndex) + (index - lowerIndex) * (relative_losses_param_in.get(upperIndex) - relative_losses_param_in.get(lowerIndex));
  		         
  	  		     double mean_out = Instance.AverageOfArray(relative_losses_param_out);
	             double var_out = Instance.VarianceOfArray(relative_losses_param_out);
  		     
  	  		     double mean_t = Instance.AverageOfArray(times_param);
	  		     double var_t = Instance.VarianceOfArray(times_param);  
	  		     
	  		     min_loss_in.add(min_in);
  		  		 
	  		     mean_loss_in.add(mean_in);
  	  		     var_loss_in.add(var_in); 
  	  		     
  	  		     mean_loss_out.add(mean_out);
  		         var_loss_out.add(var_out); 
  		     
  	  		     mean_time.add(mean_t);
	  		     var_time.add(var_t);
	  	  	}
	  			   	
  		   
          FileWriter fileWriter = new FileWriter("Experiment_6.txt");
          PrintWriter printWriter = new PrintWriter(fileWriter);

          printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
          printWriter.println();
          
          printWriter.printf("Parameters: n = %d, number_of_scenarios = %d, k_l = %d, k_f = %d, const_l = %.2f, const_f = %.2f \n", n, number_of_scenarios, k_l, k_f, const_l, const_f); 
          printWriter.printf("Number of objective errors: %d \n", error_objective_num); 
        	
          printWriter.println();
          
          printWriter.printf("True formulation: \n");
          printWriter.printf("Out-of-sample performance: %.3f (%.3f) \n", Instance.AverageOfArray(relative_losses_true_out), Instance.VarianceOfArray(relative_losses_true_out));
          
          printWriter.println();
          
          printWriter.printf("Basic formulation: \n");
          printWriter.printf("Out-of-sample performance: %.3f (%.3f) \n", Instance.AverageOfArray(relative_losses_b_out), Instance.VarianceOfArray(relative_losses_b_out));
          printWriter.printf("In-sample performance: %.3f (%.3f) \n", Instance.AverageOfArray(relative_losses_b_in), Instance.VarianceOfArray(relative_losses_b_in));
          printWriter.printf("Time: %.2f (%.2f) \n", Instance.AverageOfArray(times_b), Instance.VarianceOfArray(times_b));
          
          printWriter.println();
          
          printWriter.printf("Average in-sample performance of the semi-pessimistic approximation: \n");
          
          num_param = 1;
          for (double noise_level = 0.2; noise_level <= max_noise; noise_level += 0.3) 
          {
        	  printWriter.printf("noise level = %.2f: ", noise_level);
              for (int r = 0; r < number_of_scenarios; r += 1) 
              {
                 printWriter.printf("%.3f (%.3f) ", mean_loss_in.get(num_param - 1), var_loss_in.get(num_param - 1));
                 num_param++;
              }
              printWriter.println();  
          }
          
          printWriter.println();
          
          printWriter.printf("Worst-case in-sample performance of the semi-pessimistic approximation: \n");
          
          num_param = 1;
          for (double noise_level = 0.2; noise_level <= max_noise; noise_level += 0.3) 
          {
        	  printWriter.printf("noise level = %.2f: ", noise_level);
              for (int r = 0; r < number_of_scenarios; r += 1) 
              {
                  printWriter.printf("%.3f ", min_loss_in.get(num_param - 1));
                  num_param++;
              }
              printWriter.println();  
          }
          
          printWriter.println();
         
          printWriter.printf("Average out-of-sample performance of the semi-pessimistic approximation: \n");
          
          num_param = 1;
          for (double noise_level = 0.2; noise_level <= max_noise; noise_level += 0.3) 
          {
        	  printWriter.printf("noise level = %.2f: ", noise_level);
              for (int r = 0; r < number_of_scenarios; r += 1) 
              {
                 printWriter.printf("%.3f (%.3f) ", mean_loss_out.get(num_param - 1), var_loss_out.get(num_param - 1));
                 num_param++;
              }
              printWriter.println();  
          }
          
          printWriter.println();
         
          
          printWriter.printf("Solution times: \n");
          
          num_param = 1;
          for (double noise_level = 0.2; noise_level <= max_noise; noise_level += 0.3) 
          {
        	  printWriter.printf("noise level = %.2f: ", noise_level);
              for (int r = 0; r < number_of_scenarios; r += 1) 
              {
                  printWriter.printf("%.2f (%.2f) ", mean_time.get(num_param - 1), var_time.get(num_param - 1));
                  num_param++;
              }
              printWriter.println();  
          }
          printWriter.flush();
          printWriter.close();
      }   
        
	        
//***********************************************************************************	         	  	        
	        //Experiment 7 (in-sample and out-sample analysis, sample-average leader)  
	  	    num_of_experiment = 7;
	  	    if (num_of_experiment == 7)
  	        { 
  	        	System.out.println("Experiment 7: \n");
  	            int error_objective_num = 0;
  	            
  	            double alpha_l_mod = 0.9;
  		        k_l = 30;
  		        k_f = 30;
  		        
		        const_l = 0.; 
		        double const_f_true = 0.1;
		        
		        epsilon_l = const_l/Math.sqrt(k_f);
		        double epsilon_f_true = const_f_true/Math.sqrt(k_f);
		        
		        double alpha_f_true = 0.95;
		        
		        int[] sample_sizes = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
		        
		        int number_of_scenarios = 5;
		        double max_noise = 0.5;
		        
		        Map<ArrayList<Integer>, Double> relative_losses_b_in = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> relative_losses_sp_in = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> relative_losses_p_in = new HashMap<ArrayList<Integer>, Double>();
		        
  		        Map<ArrayList<Integer>, Double> relative_losses_b_out = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> relative_losses_sp_out = new HashMap<ArrayList<Integer>, Double>(); 
  		        Map<ArrayList<Integer>, Double> relative_losses_p_out = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> relative_losses_true_out = new HashMap<ArrayList<Integer>, Double>();
  		        
  		        Map<ArrayList<Integer>, Double> times_b = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> times_sp = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> times_p = new HashMap<ArrayList<Integer>, Double>();
 		       
  		        
  		        for (int t1 = 0; t1 < number_of_instances_1; t1++)
  		        {
  			        Instance I = random_instances.get(t1);
  			        
  			        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
  			        int number_of_sc = support_constraints.size() - 1; 
  			        
  			        ArrayList<ArrayList<Double>> shifts = new ArrayList<ArrayList<Double>>();
	  		        ArrayList<ArrayList<Double>> missed = new ArrayList<ArrayList<Double>>();
	  		        
	  		        for (int i = 0; i < I.number_of_items; i++)
	  		        {
	  		        	shifts.add(new ArrayList<Double>());
	  		        	missed.add(new ArrayList<Double>());
	  		        	
	  		            for (int j = 0; j < k_f; j++)
	  		            {
	  		            	shifts.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
	  		            	missed.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
	  		            }
	  		        }
  			        
  			        for (int t2 = 0; t2 < number_of_instances_2; t2++)
  			        {    
  				        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
  			        	System.out.println(number_of_instance);	
  				        
  				        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
  				        	                      
		    			ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
		    			ArrayList<ArrayList<Double>> data_set_to_add = Instance.GenerateData(I, k_f); 
		    			
		    			ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
  				           
		  				//Compute average for leader
		  				ArrayList<Double> average_l = new ArrayList<Double>(); 
		  				ArrayList<Double> average_all = new ArrayList<Double>(); 
		  				
		  				for (int i = 0; i < I.number_of_items; i++)
	  				    {
	  				    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));		
	  				    	    average_all.add(Instance.AverageOfArray(data_set_all.get(i)));	
	  				    }
		  				
		  			    //Solve pessimistic formulation
		  				double lower_bound = 0.;
		  				double upper_bound = I.number_of_items; 
		  				M1 = Math.max(I.number_of_items, 1.);
		  			 
		  			    double start_p = System.currentTimeMillis();
						
		  			    ArrayList<Double> sol_p = Instance.SolveMasterProblemAmbiguityFree(cplex, cplex2, I, lower_bound, upper_bound, alpha_l_mod, k_l, data_set_l, M1, start_p);
		  				
		  			    double finish_p = System.currentTimeMillis();
						
			  			//Main parameter loop
			  			int num_param = 1;	
			  			
			  			for (int l = 0; l < sample_sizes.length; l++)
			  			{
			  				ArrayList<Integer> instance = new ArrayList<Integer>();
 				            instance.add(number_of_instance);
 				            instance.add(num_param);
 				             
 				            System.out.println(num_param);
 				             
 				            num_param ++;
			  				
			  				int k_lf = sample_sizes[l];
			  			
			  				ArrayList<ArrayList<Double>> data_set_part_f = Instance.TakeSubsetColumns(data_set_to_add, I.number_of_items, k_f - k_lf);
			    			ArrayList<ArrayList<Double>> data_set_common = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_lf); 
			    			ArrayList<ArrayList<Double>> data_set_f_true = Instance.JoinDataSets(data_set_common, data_set_part_f, I.number_of_items);
			    			ArrayList<ArrayList<Double>> delta_f_true = Instance.ComputeDelta(I, data_set_f_true, support_constraints);
			    			
			    			ArrayList<Double> average_f_true = new ArrayList<Double>(); 
			    			ArrayList<Double> average_f = new ArrayList<Double>(); 
					  		 
			    			for (int i = 0; i < I.number_of_items; i++)
		  				    {
		  				       average_f_true.add(Instance.AverageOfArray(data_set_f_true.get(i)));
		  				       average_f.add(Instance.AverageOfArray(data_set_common.get(i)));	
		  				    }
			  				
			    			//Solve true formulation
			  				M1 = Math.max(I.number_of_items, 1 + epsilon_f_true/(1 - alpha_f_true));
			  				
			  				ArrayList<ArrayList<Double>> sol_true = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f_true, alpha_l_mod, alpha_f_true, k_l, k_f, data_set_l, data_set_f_true,
						    			support_constraints, "Risk-averse", "Risk-averse",  delta_l, delta_f_true, average_l, average_f_true, M1);
						    
			    			//Solve basic formulation
			  				double start = System.currentTimeMillis();
			  				
			  				ArrayList<ArrayList<Double>> sol_b = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f_true, alpha_l_mod, alpha_f_true, k_l, k_f, data_set_l, data_set_l,
						    			support_constraints, "Risk-averse", "Risk-averse",  delta_l, delta_l, average_l, average_l, M1);
						    
			  				double finish = System.currentTimeMillis();
			  				
			  	
			  				times_b.put(instance, (finish-start)*Math.pow(10, -3));
			  				times_p.put(instance, (finish_p-start_p)*Math.pow(10, -3));
						  		
			  			    //Compute in-sample performance of basic and pessimistic formulations
			  				double z_true = sol_true.get(0).get(I.number_of_items);	
			  				double z_b = sol_b.get(0).get(I.number_of_items);	
			  			    double z_p = sol_p.get(I.number_of_items);
			  			    
			  			    relative_losses_b_in.put(instance, z_b/z_true);
			  			    relative_losses_p_in.put(instance, z_p/z_true);
			  			    
			  			    
			  			    //Compute out-of-sample performance for the basic and pessimistic formulations
			  			    ArrayList<Double> x_true = Instance.CopyDoubleVector1(sol_true.get(0), I.number_of_items);  
			  			    ArrayList<Double> x_b = Instance.CopyDoubleVector1(sol_b.get(0), I.number_of_items); 
			  			    ArrayList<Double> x_p = Instance.CopyDoubleVector1(sol_p, I.number_of_items); 
			  			   
			  			    ArrayList<ArrayList<Double>> sol_x_true = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_true, epsilon_f_true, alpha_l_mod, alpha_f_true, sample_size_all, k_f, data_set_all, data_set_f_true,
					  				support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true);
			  			    
			  			    ArrayList<ArrayList<Double>> sol_x_b = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_b, epsilon_f_true, alpha_l_mod, alpha_f_true, sample_size_all, k_f, data_set_all, data_set_f_true,
					  				support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true);
				        
				            ArrayList<ArrayList<Double>> sol_x_p = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_p, epsilon_f_true, alpha_l_mod, alpha_f_true, sample_size_all, k_f, data_set_all, data_set_f_true,
		  				  				support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true);
					           
				            double numerator_true = sol_x_true.get(0).get(0);
				            double numerator_b = sol_x_b.get(0).get(0);
				            double numerator_p = sol_x_p.get(0).get(0);
				            
				            ArrayList<ArrayList<Double>> full_information_solution = Instance.SolveLeadersProblemFull(cplex, I, epsilon_f_true, alpha_l_mod, alpha_f_true, sample_size_all, k_f, data_set_all, data_set_f_true,
			    			    		support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true, M1);
			  				
			  			    double denominator = full_information_solution.get(0).get(I.number_of_items);	
			  			   
			  			    if ((numerator_b < denominator - Math.pow(10, -3)) || (numerator_p < denominator - Math.pow(10, -3)) || (numerator_true < denominator - Math.pow(10, -3)) || (z_p < z_true - Math.pow(10, -3)))
			          			{
			          				error_objective_num ++;
			          			}
			  			    
			  			    relative_losses_true_out.put(instance, numerator_true/denominator);
			  			    relative_losses_b_out.put(instance, numerator_b/denominator);
			  			    relative_losses_p_out.put(instance, numerator_p/denominator);
			  			    
			  			    //Solve semi-pessimistic formulation 
			  			    Map<Integer, ArrayList<ArrayList<Double>>> data_sets_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
			  				Map<Integer, ArrayList<ArrayList<Double>>> deltas_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
			  				Map<Integer, ArrayList<Double>> averages_f = new HashMap<Integer, ArrayList<Double>>();
			  				Map<Integer, Double> epsilons_f = new HashMap<Integer, Double>();
			  				Map<Integer, Double> alphas_f = new HashMap<Integer, Double>();
			  				
			  				for (int r = 0; r < number_of_scenarios; r += 1)
			  				{
			  					ArrayList<ArrayList<Double>> data_set_f_r = Instance.GenerateNewDataSet(I, data_set_f_true, I.number_of_items, k_lf, max_noise, shifts, missed);
			  					ArrayList<ArrayList<Double>> delta_f_r = Instance.ComputeDelta(I, data_set_f_r, support_constraints);  
			  					
			  					ArrayList<Double> average_f_r = new ArrayList<Double>();
			  					
			  					for (int i = 0; i < I.number_of_items; i++)
			  				    {
			  				    	    average_f_r.add(Instance.AverageOfArray(data_set_f_r.get(i)));	
			  				    }
			  					
			  					averages_f.put(r, average_f_r);
			  					data_sets_f.put(r, data_set_f_r);
			  					deltas_f.put(r, delta_f_r);
			  					
			  				    epsilons_f.put(r, const_f_true/Math.sqrt(k_f));
			  					alphas_f.put(r, alpha_f_true);
			  					
			  					//epsilons_f.put(r, ThreadLocalRandom.current().nextDouble(0.05, 0.15)/Math.sqrt(k_f));
			  					//alphas_f.put(r, ThreadLocalRandom.current().nextDouble(0.9, 0.99));
			  				}
			  				
			  				M1 = Math.max(I.number_of_items, 1 + epsilon_f_true/(1 - alpha_f_true));
			  				
			  			    //M1 = Math.max(I.number_of_items, 1 + 0.15/(1 - 0.99));
			  				
					  		start = System.currentTimeMillis();
										        
					  		ArrayList<Double> sol_sp =  Instance.SolveLeadersProblemScenarios2(cplex, I, number_of_scenarios, epsilon_l, epsilons_f, alpha_l_mod, alphas_f, k_l, k_f, data_set_l, data_sets_f,
				            			support_constraints, "Risk-averse",  "Risk-averse", delta_l,  deltas_f, average_l, averages_f, M1);
						        
					  		finish = System.currentTimeMillis();
					  		    
					  		times_sp.put(instance, (finish-start)*Math.pow(10, -3));
					  		
					  		//Compute in-sample performance
					  		double z_sp = sol_sp.get(I.number_of_items);	
					  		   
					  		relative_losses_sp_in.put(instance, z_sp/z_true);
					  		    
					  		//Compute out-of-sample performance 
					  	    ArrayList<Double> x_sp = Instance.CopyDoubleVector1(sol_sp, I.number_of_items); 
					  		    
					  	    ArrayList<ArrayList<Double>> sol_x_sp = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_sp, epsilon_f_true, alpha_l_mod, alpha_f_true, sample_size_all, k_f, data_set_all, data_set_f_true,
		  				  				support_constraints, "Risk-averse", "Risk-averse", delta_f_true, average_all, average_f_true);
					             
					        double numerator_sp = sol_x_sp.get(0).get(0);
					           
					        if (numerator_sp < denominator - Math.pow(10, -3))
			          			{
			          				error_objective_num ++;
			          			}
					            
					  		relative_losses_sp_out.put(instance, numerator_sp/denominator);
					  		    
			  			} 
			  		}
		  		}
		  	
  		      ArrayList<Double> mean_loss_true_out = new ArrayList<Double>();
	          ArrayList<Double> var_loss_true_out = new ArrayList<Double>();  
  		        
  		      ArrayList<Double> min_loss_sp_in = new ArrayList<Double>();
	         
  		      ArrayList<Double> mean_loss_b_in = new ArrayList<Double>();
	          ArrayList<Double> var_loss_b_in = new ArrayList<Double>();
  		      
  		      ArrayList<Double> mean_loss_sp_in = new ArrayList<Double>();
	          ArrayList<Double> var_loss_sp_in = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_loss_p_in = new ArrayList<Double>();
	          ArrayList<Double> var_loss_p_in = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_loss_b_out = new ArrayList<Double>();
	          ArrayList<Double> var_loss_b_out = new ArrayList<Double>();
	         
	          ArrayList<Double> mean_loss_sp_out = new ArrayList<Double>();
	          ArrayList<Double> var_loss_sp_out = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_loss_p_out = new ArrayList<Double>();
	          ArrayList<Double> var_loss_p_out = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_b_time = new ArrayList<Double>();
	          ArrayList<Double> var_b_time = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_sp_time = new ArrayList<Double>();
	          ArrayList<Double> var_sp_time = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_p_time = new ArrayList<Double>();
	          ArrayList<Double> var_p_time = new ArrayList<Double>();
	          	
	          int num_param = 1;
	          
		      for (int l = 0; l < sample_sizes.length; l++)
		        {
		    	    ArrayList<Double> relative_losses_param_true_out = new ArrayList<Double>();
		      
		    	    ArrayList<Double> relative_losses_param_b_in = new ArrayList<Double>();
		    	    ArrayList<Double> relative_losses_param_b_out = new ArrayList<Double>();
  		      
		    	    ArrayList<Double> relative_losses_param_sp_in = new ArrayList<Double>();
        		    ArrayList<Double> relative_losses_param_sp_out = new ArrayList<Double>();
  		          
        		    ArrayList<Double> relative_losses_param_p_in = new ArrayList<Double>();
        		    ArrayList<Double> relative_losses_param_p_out = new ArrayList<Double>();
  		          
        		    ArrayList<Double> times_param_b = new ArrayList<Double>();
        		    ArrayList<Double> times_param_sp = new ArrayList<Double>();
        		    ArrayList<Double> times_param_p = new ArrayList<Double>();
    		            
		        	for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
		        	{
		        		ArrayList<Integer> instance = new ArrayList<Integer>();
			            instance.add(t);
			            instance.add(num_param);
			            
			            relative_losses_param_true_out.add(relative_losses_true_out.get(instance));
			            
			            relative_losses_param_b_in.add(relative_losses_b_in.get(instance));
			            relative_losses_param_b_out.add(relative_losses_b_out.get(instance));
			            relative_losses_param_sp_in.add(relative_losses_sp_in.get(instance));
			            relative_losses_param_sp_out.add(relative_losses_sp_out.get(instance));
			            relative_losses_param_p_in.add(relative_losses_p_in.get(instance));
			            relative_losses_param_p_out.add(relative_losses_p_out.get(instance));
			            
			            times_param_b.add(times_b.get(instance));
			            times_param_sp.add(times_sp.get(instance));
			            times_param_p.add(times_p.get(instance));   
			            
		        	}
		        	
		        	double mean_true_out = Instance.AverageOfArray(relative_losses_param_true_out);
		        	double var_true_out = Instance.VarianceOfArray(relative_losses_param_true_out);
		        	
		        	double mean_b_in = Instance.AverageOfArray(relative_losses_param_b_in);
		        	double var_b_in = Instance.VarianceOfArray(relative_losses_param_b_in);
		        	
		        	double mean_b_out = Instance.AverageOfArray(relative_losses_param_b_out);
		        	double var_b_out = Instance.VarianceOfArray(relative_losses_param_b_out);
		        	
		        	double mean_sp_in = Instance.AverageOfArray(relative_losses_param_sp_in);
		        	double var_sp_in = Instance.VarianceOfArray(relative_losses_param_sp_in);
		        	
		        	double mean_sp_out = Instance.AverageOfArray(relative_losses_param_sp_out);
		        	double var_sp_out = Instance.VarianceOfArray(relative_losses_param_sp_out);
		        	
		        	double mean_p_in = Instance.AverageOfArray(relative_losses_param_p_in);
		        	double var_p_in = Instance.VarianceOfArray(relative_losses_param_p_in);
		        	
		        	double mean_p_out = Instance.AverageOfArray(relative_losses_param_p_out);
		        	double var_p_out = Instance.VarianceOfArray(relative_losses_param_p_out);
		        	
		        	double mean_t_b = Instance.AverageOfArray(times_param_b);
		        	double var_t_b = Instance.VarianceOfArray(times_param_b);
		        	
		        	double mean_t_sp = Instance.AverageOfArray(times_param_sp);
		        	double var_t_sp = Instance.VarianceOfArray(times_param_sp);
		        	
		        	double mean_t_p = Instance.AverageOfArray(times_param_p);
		        	double var_t_p = Instance.VarianceOfArray(times_param_p);
		        	
		        	Collections.sort(relative_losses_param_sp_in); 
		  		    
	  		        double index = 0.05*(number_of_instances_1*number_of_instances_2 - 1);
	  		        int lowerIndex = (int) Math.floor(index);
	  		        int upperIndex = (int) Math.ceil(index);
	  		         
	  		        double min_sp = relative_losses_param_sp_in.get(lowerIndex) + (index - lowerIndex) * (relative_losses_param_sp_in.get(upperIndex) - relative_losses_param_sp_in.get(lowerIndex));
	  		         
	  		        mean_loss_true_out.add(mean_true_out);
		        	var_loss_true_out.add(var_true_out);
	  		        
	  		        mean_loss_b_in.add(mean_b_in);
		        	var_loss_b_in.add(var_b_in);
		        	
	  		        mean_loss_b_out.add(mean_b_out);
		        	var_loss_b_out.add(var_b_out);
		        	
		        	mean_loss_sp_in.add(mean_sp_in);
		        	var_loss_sp_in.add(var_sp_in);
		        	
		        	mean_loss_p_in.add(mean_p_in);
		        	var_loss_p_in.add(var_p_in);
		        	
		        	mean_loss_sp_out.add(mean_sp_out);
		        	var_loss_sp_out.add(var_sp_out);
		        	
		        	mean_loss_p_out.add(mean_p_out);
		        	var_loss_p_out.add(var_p_out);
		        	
		        	min_loss_sp_in.add(min_sp);
		        	
		        	mean_b_time.add(mean_t_b);
		        	mean_sp_time.add(mean_t_sp);
		        	mean_p_time.add(mean_t_p);
		        	var_b_time.add(var_t_b);
		        	var_sp_time.add(var_t_sp);
		        	var_p_time.add(var_t_p);
		        	
		        	num_param ++;
		        }

	          FileWriter fileWriter = new FileWriter("Experiment_7.txt");
	          PrintWriter printWriter = new PrintWriter(fileWriter);

	          printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
	          printWriter.println();
	          
	          printWriter.printf("Parameters: n = %d, number_of_scenarios = %d, k_l = %d, k_f = %d, const_l = %.2f, const_f = %.2f \n", n, number_of_scenarios, k_l, k_f, const_l, const_f); 
	          printWriter.printf("Number of objective errors: %d \n", error_objective_num); 
	        	
	          printWriter.println();
	          
	          //True basic formulation
	          printWriter.printf("True basic formulation: \n");
	          printWriter.println();
              printWriter.printf("Out-of-sample performance: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.3f (%.3f) ", mean_loss_true_out.get(num_param - 1), var_loss_true_out.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();
	          printWriter.println();
	          
	          //Basic formulation
	          printWriter.printf("Basic formulation: \n");
	          printWriter.println();
	          
              printWriter.printf("In-sample performance: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.3f (%.3f) ", mean_loss_b_in.get(num_param - 1), var_loss_b_in.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();  
	          
	          printWriter.printf("Out-of-sample performance: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.3f (%.3f) ", mean_loss_b_out.get(num_param - 1), var_loss_b_out.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();  
	          
		  	  printWriter.printf("Time: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.2f (%.2f) ", mean_b_time.get(num_param - 1), var_b_time.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();  
	          printWriter.println();
	          
	          printWriter.printf("Semi-pessimistic formulation: \n");
	          printWriter.println();
	          
	          printWriter.printf("In-sample performance: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.3f (%.3f) ", mean_loss_sp_in.get(num_param - 1), var_loss_sp_in.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();
	          
              printWriter.printf("Worst-case performance: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.3f ", min_loss_sp_in.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();
	          
              printWriter.printf("Out-of-sample performance: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.3f (%.3f) ", mean_loss_sp_out.get(num_param - 1), var_loss_sp_out.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();
	          
	          printWriter.printf("Time: ");
	         
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.2f (%.2f) ", mean_sp_time.get(num_param - 1), var_sp_time.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();  
	          printWriter.println();  
	          
	          printWriter.printf("Pessimistic formulation: \n");
	          printWriter.println();
	          
	          printWriter.printf("In-sample performance: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.3f (%.3f) ", mean_loss_p_in.get(num_param - 1), var_loss_p_in.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();
	              
              printWriter.printf("Out-of-sample performance: ");
	          
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.3f (%.3f) ", mean_loss_p_out.get(num_param - 1), var_loss_p_out.get(num_param - 1));
	                 num_param++;
	          }
	          
	          printWriter.println();
	          
	          printWriter.printf("Time: ");
	         
	          num_param = 1;
	          for (int l = 0; l < sample_sizes.length; l++)
	 		  {
	                 printWriter.printf("%.2f (%.2f) ", mean_p_time.get(num_param - 1), var_p_time.get(num_param - 1));
	                 num_param++;
	          }
	          
	         
	          printWriter.flush();
	          printWriter.close();
  	      }  

//**************************************************************************************  
	  	  //Experiment 8 (in-sample and out-sample analysis, risk-neutral leader, small epsilon_f)  //don't forget to put binary indicator = 0
		  // num_of_experiment = 8;
	  	    if (num_of_experiment == 8)
  	        { 
  	        	System.out.println("Experiment 8: \n");
  	            int error_objective_num = 0;
  	            
  	            k_l = 30;
  		        k_f = 30;
  		      
		        const_l = 0.1; 
		        const_f = 0.1;
		     
		        epsilon_l = const_l/Math.sqrt(k_f);
		        epsilon_f = const_f/Math.sqrt(k_f);
		        
		        int[] sample_sizes = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
	        	
		        int number_of_scenarios = 5;
		        double max_noise = 0.5;
		        
		        Map<ArrayList<Integer>, Double> relative_losses_true_out = new HashMap<ArrayList<Integer>, Double>();
  		        
		        Map<ArrayList<Integer>, Double> relative_losses_b_in = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> relative_losses_sp_in = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> relative_losses_p_in = new HashMap<ArrayList<Integer>, Double>();
		        
  		        Map<ArrayList<Integer>, Double> relative_losses_b_out = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> relative_losses_sp_out = new HashMap<ArrayList<Integer>, Double>(); 
  		        Map<ArrayList<Integer>, Double> relative_losses_p_out = new HashMap<ArrayList<Integer>, Double>();
		           
  		        Map<ArrayList<Integer>, Double> times_b = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> times_sp = new HashMap<ArrayList<Integer>, Double>();
  		        Map<ArrayList<Integer>, Double> times_p = new HashMap<ArrayList<Integer>, Double>();
 		       
  		        
	  		    for (int t1 = 0; t1 < number_of_instances_1; t1++)
	  		        {
	  			        Instance I = random_instances_u.get(t1);
	  			      
	  			        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
	  			        int number_of_sc = support_constraints.size() - 1; 
	  			        
	  			        ArrayList<ArrayList<Double>> shifts = new ArrayList<ArrayList<Double>>();
		  		        ArrayList<ArrayList<Double>> missed = new ArrayList<ArrayList<Double>>();
		  		        
		  		        for (int i = 0; i < I.number_of_items; i++)
		  		        {
		  		        	shifts.add(new ArrayList<Double>());
		  		        	missed.add(new ArrayList<Double>());
		  		        	
		  		            for (int j = 0; j < k_f; j++)
		  		            {
		  		            	shifts.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
		  		            	missed.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
		  		            }
		  		        }
	  			        
	  			        for (int t2 = 0; t2 < number_of_instances_2; t2++)
	  			        {    
	  				        int number_of_instance = t1*number_of_instances_2 + t2 + 1;
	  			        	System.out.println(number_of_instance);	
	  				         
	  				        ArrayList<ArrayList<Double>> data_set_all = Instance.GenerateData(I, sample_size_all); 
	  				        	                      
			    			ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
			    			ArrayList<ArrayList<Double>> data_set_to_add = Instance.GenerateData(I, k_f); 
			    			ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
	  				           
			  				//Compute average for leader
			  				ArrayList<Double> average_l = new ArrayList<Double>(); 
			  				ArrayList<Double> average_all = new ArrayList<Double>(); 
			  				
			  				for (int i = 0; i < I.number_of_items; i++)
		  				    {
		  				    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));		
		  				    	    average_all.add(Instance.AverageOfArray(data_set_all.get(i)));	
		  				    }
			  				
			  			    //Solve pessimistic formulation
			  				double lower_bound = 0.;
			  				double upper_bound = I.number_of_items; 
			  				M1 = Math.max(I.number_of_items, epsilon_f/(1 - alpha_f));
			  			    double M_objective = I.number_of_items;
			  			    
			  			    double start_p = System.currentTimeMillis();
							
			  			    ArrayList<Double> sol_p = Instance.SolveMasterProblemRiskNeutral(cplex, cplex2, I, lower_bound, upper_bound, support_constraints, epsilon_l, k_l, data_set_l, delta_l, average_l, M_objective, binary_indicator, start_p);
			  				
			  			    double finish_p = System.currentTimeMillis();
			  			     	
				  			//Main parameter loop
				  			int num_param = 1;	
				  			
				  			for (int ind = 0; ind < 2; ind++)
				  				for (int l = 0; l < sample_sizes.length; l++)
					  			{
					  		        String type_of_risk_f = "Risk-averse";
					  		        
					  		        if (ind == 0) type_of_risk_f = "Risk-neutral"; 
				  				
					  				ArrayList<Integer> instance = new ArrayList<Integer>();
		 				            instance.add(number_of_instance);
		 				            instance.add(num_param);
					  				
		 				            System.out.println(num_param);
						             
						            num_param ++;
		 				            
					  				int k_lf = sample_sizes[l];
					  
					  				ArrayList<ArrayList<Double>> data_set_part_f = Instance.TakeSubsetColumns(data_set_to_add, I.number_of_items, k_f - k_lf);
					    			ArrayList<ArrayList<Double>> data_set_common = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_lf); 
					    			ArrayList<ArrayList<Double>> data_set_f_true = Instance.JoinDataSets(data_set_common, data_set_part_f, I.number_of_items);
					    			ArrayList<ArrayList<Double>> delta_f_true = Instance.ComputeDelta(I, data_set_f_true, support_constraints);
					    			
					    			ArrayList<Double> average_f_true = new ArrayList<Double>(); 
					    			ArrayList<Double> average_f = new ArrayList<Double>(); 
							  		 
					    			for (int i = 0; i < I.number_of_items; i++)
				  				    {
				  				       average_f_true.add(Instance.AverageOfArray(data_set_f_true.get(i)));
				  				       average_f.add(Instance.AverageOfArray(data_set_common.get(i)));	
				  				    }
					  				
					  			    //Solve true formulation
					  				M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
					  				
					  				ArrayList<ArrayList<Double>> sol_true = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f_true,
								    			support_constraints, "Risk-neutral", type_of_risk_f,  delta_l, delta_f_true, average_l, average_f_true, M1);
								    
					    			
					  			    //Solve basic formulation
					  				double start_b = System.currentTimeMillis();
					  				
					  				ArrayList<ArrayList<Double>> sol_b = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_l,
								    			support_constraints, "Risk-neutral", type_of_risk_f,  delta_l, delta_l, average_l, average_l, M1);
								    
					  				double finish_b = System.currentTimeMillis();
					  				
					  				times_b.put(instance, (finish_b-start_b)*Math.pow(10, -3));
					  				times_p.put(instance, (finish_p-start_p)*Math.pow(10, -3));
					  			   
					  				//Compute in-sample performance of the pessimistic formulation
					  				double z_true = sol_true.get(0).get(I.number_of_items);	
					  				double z_b = sol_b.get(0).get(I.number_of_items);	
					  			    double z_p = sol_p.get(I.number_of_items);
					  			    
					  			    relative_losses_b_in.put(instance, z_b/z_true);
					  			    relative_losses_p_in.put(instance, z_p/z_true);
					  			   
					  			   
					  			    //Compute out-of-sample performance for the basic and pessimistic formulations
					  			    ArrayList<Double> x_b = Instance.CopyDoubleVector1(sol_b.get(0), I.number_of_items); 
					  			    ArrayList<Double> x_p = Instance.CopyDoubleVector1(sol_p, I.number_of_items); 
					  			    ArrayList<Double> x_true = Instance.CopyDoubleVector1(sol_true.get(0), I.number_of_items); 

					  			    ArrayList<ArrayList<Double>> sol_x_true = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_true, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
							  				support_constraints, "Risk-neutral", type_of_risk_f, delta_f_true, average_all, average_f_true);
					  			    
					  			    ArrayList<ArrayList<Double>> sol_x_b = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_b, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
							  				support_constraints, "Risk-neutral", type_of_risk_f, delta_f_true, average_all, average_f_true);
						        
						            ArrayList<ArrayList<Double>> sol_x_p = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_p, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
				  				  				support_constraints, "Risk-neutral", type_of_risk_f, delta_f_true, average_all, average_f_true);
							           
						            double numerator_true = sol_x_true.get(0).get(0);
						            double numerator_b = sol_x_b.get(0).get(0);
						            double numerator_p = sol_x_p.get(0).get(0);
					  			    
						            ArrayList<ArrayList<Double>> full_information_solution = Instance.SolveLeadersProblemFull(cplex, I, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
					    			    		support_constraints, "Risk-neutral", type_of_risk_f, delta_f_true, average_all, average_f_true, M1);
					  				
					  			    double denominator = full_information_solution.get(0).get(I.number_of_items);	
					  			   
					  			     if ((numerator_b < denominator - Math.pow(10, -3)) || (numerator_p < denominator - Math.pow(10, -3)) || (numerator_true < denominator - Math.pow(10, -3)) || (z_p < z_true - Math.pow(10, -3)))
				          			{
				          				error_objective_num ++;
				          			}
				  			    
				  			        relative_losses_true_out.put(instance, numerator_true/denominator);
					  				relative_losses_b_out.put(instance, numerator_b/denominator);
						  			relative_losses_p_out.put(instance, numerator_p/denominator);
					  			    
						  			//Solve semi-pessimistic formulation
					  			    Map<Integer, ArrayList<ArrayList<Double>>> data_sets_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
					  				Map<Integer, ArrayList<ArrayList<Double>>> deltas_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
					  				Map<Integer, ArrayList<Double>> averages_f = new HashMap<Integer, ArrayList<Double>>();
					  				
					  				for (int r = 0; r < number_of_scenarios; r += 1)
					  				{
					  					ArrayList<ArrayList<Double>> data_set_f_r = Instance.GenerateNewDataSet(I, data_set_f_true, I.number_of_items, k_lf, max_noise, shifts, missed);
					  					ArrayList<ArrayList<Double>> delta_f_r = Instance.ComputeDelta(I, data_set_f_r, support_constraints);  
					  					
					  					ArrayList<Double> average_f_r = new ArrayList<Double>();
					  					
					  					for (int i = 0; i < I.number_of_items; i++)
					  				    {
					  				    	    average_f_r.add(Instance.AverageOfArray(data_set_f_r.get(i)));	
					  				    }
					  					
					  					averages_f.put(r, average_f_r);
					  					data_sets_f.put(r, data_set_f_r);
					  					deltas_f.put(r, delta_f_r);
					  				}
					  				  
					  				
					  					
					  				double start = System.currentTimeMillis();
												        
							  		ArrayList<Double> sol_sp =  Instance.SolveLeadersProblemScenarios(cplex, I, number_of_scenarios, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_sets_f,
						            			support_constraints, "Risk-neutral", type_of_risk_f, delta_l,  deltas_f, average_l, averages_f, M1);
								        
							  		double finish = System.currentTimeMillis();
							  		    
							  		times_sp.put(instance, (finish-start)*Math.pow(10, -3));
							  		
							  		
							  		
							  		//Compute in-sample performance
							  		double z_sp = sol_sp.get(I.number_of_items);	
							  		   
							  		relative_losses_sp_in.put(instance, z_sp/z_true);
							  		    
							  		//Compute out-of-sample performance 
							  	    ArrayList<Double> x_sp = Instance.CopyDoubleVector1(sol_sp, I.number_of_items); 
							  		    
							  	    ArrayList<ArrayList<Double>> sol_x_sp = Instance.SolveLeadersProblemFullFixedX(cplex, I, x_sp, epsilon_f, alpha_l, alpha_f, sample_size_all, k_f, data_set_all, data_set_f_true,
				  				  				support_constraints, "Risk-neutral", type_of_risk_f, delta_f_true, average_all, average_f_true);
							             
							        double numerator_sp = sol_x_sp.get(0).get(0);
							           
							        if (numerator_sp < denominator - Math.pow(10, -3))
					          			{
					          				error_objective_num ++;
					          			}
							            
							  		relative_losses_sp_out.put(instance, numerator_sp/denominator);
							  		     
					  	    } 
				  		}
			  		}
  		        	 	   
	  		  ArrayList<Double> mean_loss_true_out = new ArrayList<Double>();
	          ArrayList<Double> var_loss_true_out = new ArrayList<Double>();  
  		        
  		      ArrayList<Double> min_loss_sp_in = new ArrayList<Double>();
	         
  		      ArrayList<Double> mean_loss_b_in = new ArrayList<Double>();
	          ArrayList<Double> var_loss_b_in = new ArrayList<Double>();
  		      
  		      ArrayList<Double> mean_loss_sp_in = new ArrayList<Double>();
	          ArrayList<Double> var_loss_sp_in = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_loss_p_in = new ArrayList<Double>();
	          ArrayList<Double> var_loss_p_in = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_loss_b_out = new ArrayList<Double>();
	          ArrayList<Double> var_loss_b_out = new ArrayList<Double>();
	         
	          ArrayList<Double> mean_loss_sp_out = new ArrayList<Double>();
	          ArrayList<Double> var_loss_sp_out = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_loss_p_out = new ArrayList<Double>();
	          ArrayList<Double> var_loss_p_out = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_b_time = new ArrayList<Double>();
	          ArrayList<Double> var_b_time = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_sp_time = new ArrayList<Double>();
	          ArrayList<Double> var_sp_time = new ArrayList<Double>();
	          
	          ArrayList<Double> mean_p_time = new ArrayList<Double>();
	          ArrayList<Double> var_p_time = new ArrayList<Double>();
	          	
	          int num_param = 1;
	          
	          for (int ind = 0; ind < 2; ind++)
	        	  for (int l = 0; l < sample_sizes.length; l++)		 
	  		  {
	        		  ArrayList<Double> relative_losses_param_true_out = new ArrayList<Double>();
	    		      
			    	    ArrayList<Double> relative_losses_param_b_in = new ArrayList<Double>();
			    	    ArrayList<Double> relative_losses_param_b_out = new ArrayList<Double>();
	  		      
			    	    ArrayList<Double> relative_losses_param_sp_in = new ArrayList<Double>();
	        		    ArrayList<Double> relative_losses_param_sp_out = new ArrayList<Double>();
	  		          
	        		    ArrayList<Double> relative_losses_param_p_in = new ArrayList<Double>();
	        		    ArrayList<Double> relative_losses_param_p_out = new ArrayList<Double>();
	  		          
	        		    ArrayList<Double> times_param_b = new ArrayList<Double>();
	        		    ArrayList<Double> times_param_sp = new ArrayList<Double>();
	        		    ArrayList<Double> times_param_p = new ArrayList<Double>();
	    		            
			        	for (int t = 1; t <= number_of_instances_1*number_of_instances_2; t++)
			        	{
			        		ArrayList<Integer> instance = new ArrayList<Integer>();
				            instance.add(t);
				            instance.add(num_param);
				            
				            relative_losses_param_true_out.add(relative_losses_true_out.get(instance));
				            
				            relative_losses_param_b_in.add(relative_losses_b_in.get(instance));
				            relative_losses_param_b_out.add(relative_losses_b_out.get(instance));
				            relative_losses_param_sp_in.add(relative_losses_sp_in.get(instance));
				            relative_losses_param_sp_out.add(relative_losses_sp_out.get(instance));
				            relative_losses_param_p_in.add(relative_losses_p_in.get(instance));
				            relative_losses_param_p_out.add(relative_losses_p_out.get(instance));
				            
				            times_param_b.add(times_b.get(instance));
				            times_param_sp.add(times_sp.get(instance));
				            times_param_p.add(times_p.get(instance));   
				            
			        	}
			        	
			        	double mean_true_out = Instance.AverageOfArray(relative_losses_param_true_out);
			        	double var_true_out = Instance.VarianceOfArray(relative_losses_param_true_out);
			        	
			        	double mean_b_in = Instance.AverageOfArray(relative_losses_param_b_in);
			        	double var_b_in = Instance.VarianceOfArray(relative_losses_param_b_in);
			        	
			        	double mean_b_out = Instance.AverageOfArray(relative_losses_param_b_out);
			        	double var_b_out = Instance.VarianceOfArray(relative_losses_param_b_out);
			        	
			        	double mean_sp_in = Instance.AverageOfArray(relative_losses_param_sp_in);
			        	double var_sp_in = Instance.VarianceOfArray(relative_losses_param_sp_in);
			        	
			        	double mean_sp_out = Instance.AverageOfArray(relative_losses_param_sp_out);
			        	double var_sp_out = Instance.VarianceOfArray(relative_losses_param_sp_out);
			        	
			        	double mean_p_in = Instance.AverageOfArray(relative_losses_param_p_in);
			        	double var_p_in = Instance.VarianceOfArray(relative_losses_param_p_in);
			        	
			        	double mean_p_out = Instance.AverageOfArray(relative_losses_param_p_out);
			        	double var_p_out = Instance.VarianceOfArray(relative_losses_param_p_out);
			        	
			        	double mean_t_b = Instance.AverageOfArray(times_param_b);
			        	double var_t_b = Instance.VarianceOfArray(times_param_b);
			        	
			        	double mean_t_sp = Instance.AverageOfArray(times_param_sp);
			        	double var_t_sp = Instance.VarianceOfArray(times_param_sp);
			        	
			        	double mean_t_p = Instance.AverageOfArray(times_param_p);
			        	double var_t_p = Instance.VarianceOfArray(times_param_p);
			        	
			        	Collections.sort(relative_losses_param_sp_in); 
			  		    
		  		        double index = 0.05*(number_of_instances_1*number_of_instances_2 - 1);
		  		        int lowerIndex = (int) Math.floor(index);
		  		        int upperIndex = (int) Math.ceil(index);
		  		         
		  		        double min_sp = relative_losses_param_sp_in.get(lowerIndex) + (index - lowerIndex) * (relative_losses_param_sp_in.get(upperIndex) - relative_losses_param_sp_in.get(lowerIndex));
		  		         
		  		        mean_loss_true_out.add(mean_true_out);
			        	var_loss_true_out.add(var_true_out);
		  		        
		  		        mean_loss_b_in.add(mean_b_in);
			        	var_loss_b_in.add(var_b_in);
			        	
		  		        mean_loss_b_out.add(mean_b_out);
			        	var_loss_b_out.add(var_b_out);
			        	
			        	mean_loss_sp_in.add(mean_sp_in);
			        	var_loss_sp_in.add(var_sp_in);
			        	
			        	mean_loss_p_in.add(mean_p_in);
			        	var_loss_p_in.add(var_p_in);
			        	
			        	mean_loss_sp_out.add(mean_sp_out);
			        	var_loss_sp_out.add(var_sp_out);
			        	
			        	mean_loss_p_out.add(mean_p_out);
			        	var_loss_p_out.add(var_p_out);
			        	
			        	min_loss_sp_in.add(min_sp);
			        	
			        	mean_b_time.add(mean_t_b);
			        	mean_sp_time.add(mean_t_sp);
			        	mean_p_time.add(mean_t_p);
			        	var_b_time.add(var_t_b);
			        	var_sp_time.add(var_t_sp);
			        	var_p_time.add(var_t_p);
			        	
			        	num_param ++;
		        }

	          FileWriter fileWriter = new FileWriter("Experiment_8.txt");
	          PrintWriter printWriter = new PrintWriter(fileWriter);

	          printWriter.printf("%d distributions, %d data-sets \n", number_of_instances_1, number_of_instances_2); 
	          printWriter.println();
	          
	          printWriter.printf("Parameters: n = %d, number_of_scenarios = %d, k_l = %d, k_f = %d, const_l = %.2f, const_f = %.2f \n", n, number_of_scenarios, k_l, k_f, const_l, const_f); 
	          printWriter.printf("Number of objective errors: %d \n", error_objective_num); 
	        	
	          printWriter.println();
	          
	          //True basic formulation
	          printWriter.printf("True basic formulation: \n");
	          printWriter.println();
              printWriter.printf("Out-of-sample performance: ");
	          
              num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	        	  
	        	  for (int l = 0; l < sample_sizes.length; l++)
		 		  {
	                 printWriter.printf("%.3f (%.3f) ", mean_loss_true_out.get(num_param - 1), var_loss_true_out.get(num_param - 1));
	                 num_param++;
		 		  }
	        	  
	        	  printWriter.println();  
	          }
	          
	          printWriter.println();
	          
	          //Basic formulation
	          printWriter.printf("Basic formulation: \n");
	          printWriter.println();
	          
	          printWriter.printf("In-sample performance: ");
	          printWriter.println();
	          
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	        		 
		          for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_loss_b_in.get(num_param - 1), var_loss_b_in.get(num_param - 1));
		                 num_param++;
		          }
		          
		          printWriter.println();  
	          }
	          
	          printWriter.println();
				
	          
	          printWriter.printf("Out-of-sample performance: ");
	          printWriter.println();
	          
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	        		 
		          for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_loss_b_out.get(num_param - 1), var_loss_b_out.get(num_param - 1));
		                 num_param++;
		          }
		          
		          printWriter.println();  
	          }
	          
	          printWriter.println();
				
		  	  printWriter.printf("Time: ");
		  	  
		  	  printWriter.println();
		  	
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	       
		          for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_b_time.get(num_param - 1), var_b_time.get(num_param - 1));
		                 num_param++;
		          }
	              printWriter.println();  
	          }
	          
	          printWriter.println();
	          
	          printWriter.printf("Semi-pessimistic formulation: \n");
	          printWriter.println();
	          
	          printWriter.printf("In-sample performance: ");
	          printWriter.println();
			  	
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	       
	              for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_loss_sp_in.get(num_param - 1), var_loss_sp_in.get(num_param - 1));
		                 num_param++;
		          }
	              printWriter.println();
	          }
	          printWriter.println();
				 
              printWriter.printf("Worst-case performance: ");
              printWriter.println();
			  	
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	       
	        	  for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f ", min_loss_sp_in.get(num_param - 1));
		                 num_param++;
		          }
	              printWriter.println();
	          }
	          printWriter.println();
				 
              printWriter.printf("Out-of-sample performance: ");
	          
              printWriter.println();
			  	
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	       
	        	  for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_loss_sp_out.get(num_param - 1), var_loss_sp_out.get(num_param - 1));
		                 num_param++;
		          }
	          
	              printWriter.println();
	          }
	          printWriter.println();
				 
	          printWriter.printf("Time: ");
	         
	          printWriter.println();
			  	
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	       
	        	  for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_sp_time.get(num_param - 1), var_sp_time.get(num_param - 1));
		                 num_param++;
		          }
	          
	              printWriter.println(); 
	          }
	          printWriter.println();  
	          
	          printWriter.printf("Pessimistic formulation: \n");
	          printWriter.println();
	          
	          printWriter.printf("In-sample performance: ");
	          
	          printWriter.println();
			  
	          num_param = 1;
	          
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	       
	        	  for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_loss_p_in.get(num_param - 1), var_loss_p_in.get(num_param - 1));
		                 num_param++;
		          }
	          
	              printWriter.println();
	          }
	          printWriter.println();
			  	 
              printWriter.printf("Out-of-sample performance: ");
	          
              printWriter.println();
			  	
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	       
	        	  for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_loss_p_out.get(num_param - 1), var_loss_p_out.get(num_param - 1));
		                 num_param++;
		          }
	          
	              printWriter.println();
	          }
	          printWriter.println();
				  
	          printWriter.printf("Time: ");
	         
	          printWriter.println();
			  	
	          num_param = 1;
	          for (int ind = 0; ind < 2; ind++)	
	          {
	        	  if (ind == 0) printWriter.printf("Risk-neutral follower: ");
	        	  else printWriter.printf("Risk-averse follower: ");
	       
	        	  for (int l = 0; l < sample_sizes.length; l++)
		 		  {
		                 printWriter.printf("%.3f (%.3f) ", mean_p_time.get(num_param - 1), var_p_time.get(num_param - 1));
		                 num_param++;
		          }
	        	  printWriter.println();
				  
	          }
	         
			  
	          printWriter.flush();
	          printWriter.close();
  	      }   
	  	    
//***************************************************************************************
	  	//Experiment 9 (solution times k_l, k_f)        
	    //num_of_experiment = 9;     
        if (num_of_experiment == 9)
        { 
        	System.out.println("Experiment 9: \n");
         
	        int[] sample_sizes = {5, 10, 15, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200, 250, 300};
	        int max_index = sample_sizes.length;

	        int number_of_scenarios = 10;
        	double max_noise = 0.2;
        	
	        Map<ArrayList<Integer>, Double> times_b = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> times_sp_5 = new HashMap<ArrayList<Integer>, Double>();
	        Map<ArrayList<Integer>, Double> times_sp_10 = new HashMap<ArrayList<Integer>, Double>();
	        
	        for (int t1 = 0; t1 < number_of_instances_1; t1++)
	        {  
	        	int number_of_instance = t1 + 1;
		        System.out.println(number_of_instance);	
		        
		        Instance I = random_instances.get(t1);
		        
		        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
		        int number_of_sc = support_constraints.size() - 1; 
		        
		        ArrayList<ArrayList<Double>> shifts = new ArrayList<ArrayList<Double>>();
  		        ArrayList<ArrayList<Double>> missed = new ArrayList<ArrayList<Double>>();
  		        
  		        for (int i = 0; i < I.number_of_items; i++)
  		        {
  		        	shifts.add(new ArrayList<Double>());
  		        	missed.add(new ArrayList<Double>());
  		        	
  		            for (int j = 0; j < k_f; j++)
  		            {
  		            	shifts.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
  		            	missed.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
  		            }
  		        }
		        
		      
		        //Main parameter loop
                int num_param = 1;
                for (int l = 0; l < max_index; l++)
		        {
		        	
		        	 ArrayList<Integer> instance = new ArrayList<Integer>();
		             instance.add(number_of_instance);
		             instance.add(num_param);
		             
		             System.out.println(num_param);
		             num_param ++; 
		            
		             k_l = sample_sizes[l];
		        	 k_f = sample_sizes[l];
		             int k_lf = (int)(2./3.*k_f);
		        	 
		             const_l = 0.1;
		        	 const_f = 0.1; 
		        
		        	 
		        	 epsilon_l = const_l/Math.sqrt(k_l);
		        	 epsilon_f = const_l/Math.sqrt(k_l);
		        	 
		        	 ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
			         ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
			        
		             ArrayList<ArrayList<Double>> data_set_f = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_f); 
				     ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l, k_f, number_of_sc);  
			              
				     ArrayList<Double> average_l = new ArrayList<Double>();    
				     ArrayList<Double> average_f = new ArrayList<Double>();
				     
				       for (int i = 0; i < I.number_of_items; i++)
				    	{
				    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
				    	    average_f.add(Instance.AverageOfArray(data_set_f.get(i)));	
				    	} 
		             
		             M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
		             
		             double start = System.currentTimeMillis();
						   
		             ArrayList<ArrayList<Double>> sol1 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
				    			support_constraints, "Risk-averse", "Risk-averse", delta_l, delta_f, average_l, average_f, M1);
				     
		             double finish = System.currentTimeMillis();
					
		             if (sol1.size() == 0)
		  			 {
		  				 max_index = l;
		  				 break;
		  			 }
		             
                     times_b.put(instance, (finish-start)*Math.pow(10, -3));
		             
                     //Solve semi-pessimistic formulation
	  			     Map<Integer, ArrayList<ArrayList<Double>>> data_sets_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
	  				 Map<Integer, ArrayList<ArrayList<Double>>> deltas_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
	  				 Map<Integer, ArrayList<Double>> averages_f = new HashMap<Integer, ArrayList<Double>>();
	  				
	  				 for (int r = 0; r < number_of_scenarios; r += 1)
	  				 {
	  					ArrayList<ArrayList<Double>> data_set_f_r = Instance.GenerateNewDataSet(I, data_set_f, I.number_of_items, k_lf, max_noise, shifts, missed);
	  					ArrayList<ArrayList<Double>> delta_f_r = Instance.ComputeDelta(I, data_set_f_r, support_constraints);  
	  					
	  					ArrayList<Double> average_f_r = new ArrayList<Double>();
	  					
	  					for (int i = 0; i < I.number_of_items; i++)
	  				    {
	  				    	    average_f_r.add(Instance.AverageOfArray(data_set_f_r.get(i)));	
	  				    }
	  					
	  					averages_f.put(r, average_f_r);
	  					data_sets_f.put(r, data_set_f_r);
	  					deltas_f.put(r, delta_f_r);
	  				 }
		  				  
		  				
		  		     start = System.currentTimeMillis();
									        
				  	 ArrayList<Double> sol_sp_5 =  Instance.SolveLeadersProblemScenarios(cplex, I, number_of_scenarios/2, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_sets_f,
			            			support_constraints, "Risk-averse", "Risk-averse", delta_l,  deltas_f, average_l, averages_f, M1);
					        
				  	 finish = System.currentTimeMillis();
				  	
				  	 if (sol_sp_5.size() == 0)
		  			 {
		  				 max_index = l;
		  				 break;
		  			 }
				  	 
				  	 times_sp_5.put(instance, (finish-start)*Math.pow(10, -3));
				  	 
				  	 start = System.currentTimeMillis();
			        
				  	 ArrayList<Double> sol_sp_10 =  Instance.SolveLeadersProblemScenarios(cplex, I, number_of_scenarios, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_sets_f,
			            			support_constraints, "Risk-averse", "Risk-averse", delta_l,  deltas_f, average_l, averages_f, M1);
					        
				  	 finish = System.currentTimeMillis();
                     
				  	 if (sol_sp_10.size() == 0)
		  			 {
		  				 max_index = l;
		  				 break;
		  			 }
				  	 
				  	 times_sp_10.put(instance, (finish-start)*Math.pow(10, -3));    
		        }
		     }
		              
	      
	        //Compute mean and variance
	        ArrayList<Double> mean_time_b = new ArrayList<Double>();
	        ArrayList<Double> var_time_b = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_time_sp_5 = new ArrayList<Double>();
	        ArrayList<Double> var_time_sp_5 = new ArrayList<Double>();
	        
	        ArrayList<Double> mean_time_sp_10 = new ArrayList<Double>();
	        ArrayList<Double> var_time_sp_10 = new ArrayList<Double>();
	        
	        
	        int num_param = 1;
	        for (int l = 0; l < max_index; l++)
	        {
	        	
	        	ArrayList<Double> times_param_b = new ArrayList<Double>();
	        	ArrayList<Double> times_param_sp_5 = new ArrayList<Double>();
	        	ArrayList<Double> times_param_sp_10 = new ArrayList<Double>();
	        	
	        	for (int t = 1; t <= number_of_instances_1; t++)
	        	{
	        		ArrayList<Integer> instance = new ArrayList<Integer>();
		            instance.add(t);
		            instance.add(num_param);
		            
		            times_param_b.add(times_b.get(instance));  
		            times_param_sp_5.add(times_sp_5.get(instance));  
		            times_param_sp_10.add(times_sp_10.get(instance));  
		             
	        	}
	        	
	        	
	        	double mean_b = Instance.AverageOfArray(times_param_b);
	        	double var_b = Instance.VarianceOfArray(times_param_b);
	        	
	        	double mean_sp_5 = Instance.AverageOfArray(times_param_sp_5);
	        	double var_sp_5 = Instance.VarianceOfArray(times_param_sp_5);
	        	
	        	double mean_sp_10 = Instance.AverageOfArray(times_param_sp_10);
	        	double var_sp_10 = Instance.VarianceOfArray(times_param_sp_10);
	        
	        	mean_time_b.add(mean_b);
	        	var_time_b.add(var_b);
	        	
	        	mean_time_sp_5.add(mean_sp_5);
	        	var_time_sp_5.add(var_sp_5);
	        	
	        	mean_time_sp_10.add(mean_sp_10);
	        	var_time_sp_10.add(var_sp_10);
	        	
	        	num_param ++;
	        }
	        
	        FileWriter fileWriter = new FileWriter("Experiment_9.txt");
		       
	        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
	        {
	        	printWriter.printf("%d distributions \n", number_of_instances_1); 
	        	printWriter.println();
	        	
	        	printWriter.printf("Parameters: n = %d, max_noise = %.2f, const_l = %.2f, const_f = %.2f \n", n, max_noise, const_l, const_f); 
	        	
	        	printWriter.println();
	        	
	        	num_param = 1;
				
	        	for (int l = 0; l < max_index; l++)
			     {
			        k_l = sample_sizes[l];
			      
			        printWriter.printf("k_l = k_f = %d : %.2f (%.2f) | %.2f (%.2f) | %.2f (%.2f) \n", k_l, mean_time_b.get(num_param - 1), var_time_b.get(num_param - 1), mean_time_sp_5.get(num_param - 1), var_time_sp_5.get(num_param - 1), mean_time_sp_10.get(num_param - 1), var_time_sp_10.get(num_param - 1));
			        num_param++;
				 }
				  
				printWriter.flush();
				printWriter.close();
			}
        }
	  
//***************************************************************************************	        
	    	//Experiment 10 (solution times n)        
		    //num_of_experiment = 10;
	        if (num_of_experiment == 10)
	        { 
	        	System.out.println("Experiment 10: \n");
	        	
	        	int[] sizes = {5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 33, 39, 46, 50};
	        	//int[] sizes = {5};
	        	int max_index = sizes.length;
	        	
	        	Map<Integer, Integer> ind1 = new HashMap<Integer, Integer>();
	        	Map<Integer, Integer> ind2 = new HashMap<Integer, Integer>();
	        	Map<Integer, Integer> ind3 = new HashMap<Integer, Integer>();
	        	
	        	for (int l = 0; l < max_index; l++)
	        	{
	        		ind1.put(l, 1);
	        		ind2.put(l, 1);
	        		ind3.put(l, 1);
	        	}
	        	
		        k_l = 30;
	        	k_f = 30;
		        
		        int number_of_scenarios = 10;
	        	double max_noise = 0.2;
	        	
	        	const_l = 0.1;
	        	const_f = 0.1; 
	         
	        	epsilon_l = const_l/Math.sqrt(k_l);
	        	epsilon_f = const_l/Math.sqrt(k_l);
	        	
		        Map<ArrayList<Integer>, Double> times_b = new HashMap<ArrayList<Integer>, Double>();
		        Map<ArrayList<Integer>, Double> times_sp_5 = new HashMap<ArrayList<Integer>, Double>();
		        Map<ArrayList<Integer>, Double> times_sp_10 = new HashMap<ArrayList<Integer>, Double>();
		        
		        for (int t1 = 0; t1 < number_of_instances_1; t1++)
		        {
		        	int number_of_instance = t1 + 1;
		        	System.out.println(number_of_instance);	
		        
			        int num_param = 1;
		        	for (int l = 0; l < max_index; l++)
		        	{
		        		long instanceSeed = globalRandom.nextLong();
			        	Random random = new Random(instanceSeed);
		        		Instance I = new Instance(sizes[l], p, d_f, w_l, w_f, 0, random); 
			        
			            ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
			            int number_of_sc = support_constraints.size() - 1; 
			            
			            ArrayList<ArrayList<Double>> shifts = new ArrayList<ArrayList<Double>>();
		  		        ArrayList<ArrayList<Double>> missed = new ArrayList<ArrayList<Double>>();
		  		        
		  		        for (int i = 0; i < I.number_of_items; i++)
		  		        {
		  		        	shifts.add(new ArrayList<Double>());
		  		        	missed.add(new ArrayList<Double>());
		  		        	
		  		            for (int j = 0; j < k_f; j++)
		  		            {
		  		            	shifts.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
		  		            	missed.get(i).add(ThreadLocalRandom.current().nextDouble(0, 1));
		  		            }
		  		        }
			            
			            
			        	ArrayList<Integer> instance = new ArrayList<Integer>();
			            instance.add(number_of_instance);
			            instance.add(num_param);
			             
			            System.out.println(num_param);
			            num_param ++; 
			            
			            int k_lf = (int)(2./3.*k_f);
    		        	 
    		        	ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
 				        ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I, data_set_l, support_constraints);
   				        
			            ArrayList<ArrayList<Double>> data_set_f = Instance.TakeSubsetColumns(data_set_l, I.number_of_items, k_f); 
					    ArrayList<ArrayList<Double>> delta_f = Instance.TakeSubsetRows(delta_l, k_f, number_of_sc);  
				              
					    ArrayList<Double> average_l = new ArrayList<Double>();    
					    ArrayList<Double> average_f = new ArrayList<Double>();
					     
					    for (int i = 0; i < I.number_of_items; i++)
					    	{
					    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
					    	    average_f.add(Instance.AverageOfArray(data_set_f.get(i)));	
					    	} 
			             
			             M1 = Math.max(I.number_of_items, 1 + epsilon_f/(1 - alpha_f));
			             
			             double start = 0;
			             double finish = -1;
			             
			             if (ind1.get(l) == 1)
			             {
				             start = System.currentTimeMillis();
								   
				             ArrayList<ArrayList<Double>> sol1 = Instance.SolveLeadersProblem(cplex, I, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_set_f,
						    			support_constraints, "Risk-averse", "Risk-averse", delta_l, delta_f, average_l, average_f, M1);
						     
				             finish = System.currentTimeMillis();
				             
				             if (sol1.size() == 0)
				  			 {
				            	 for (int l1 = l; l1 < max_index; l1++)
				            		 ind1.put(l1, 0);
				            	 
				            	 start = 0;
				            	 finish = -1;
				  			 }
				             
			             }
						 
                         times_b.put(instance, (finish-start)*Math.pow(10, -3));
			             
                         //Solve semi-pessimistic formulation
		  			     Map<Integer, ArrayList<ArrayList<Double>>> data_sets_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
		  				 Map<Integer, ArrayList<ArrayList<Double>>> deltas_f = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
		  				 Map<Integer, ArrayList<Double>> averages_f = new HashMap<Integer, ArrayList<Double>>();
		  				
		  				 for (int r = 0; r < number_of_scenarios; r += 1)
		  				 {
		  					ArrayList<ArrayList<Double>> data_set_f_r = Instance.GenerateNewDataSet(I, data_set_f, I.number_of_items, k_lf, max_noise, shifts, missed);
		  					ArrayList<ArrayList<Double>> delta_f_r = Instance.ComputeDelta(I, data_set_f_r, support_constraints);  
		  					
		  					ArrayList<Double> average_f_r = new ArrayList<Double>();
		  					
		  					for (int i = 0; i < I.number_of_items; i++)
		  				    {
		  				    	    average_f_r.add(Instance.AverageOfArray(data_set_f_r.get(i)));	
		  				    }
		  					
		  					averages_f.put(r, average_f_r);
		  					data_sets_f.put(r, data_set_f_r);
		  					deltas_f.put(r, delta_f_r);
		  				 }
			  				  
			  				
		  				 start = 0;
			             finish = -1;
			             
			             if (ind2.get(l) == 1)
			             {
				             start = System.currentTimeMillis();
								   			        
				             ArrayList<Double> sol_sp_5 =  Instance.SolveLeadersProblemScenarios(cplex, I, number_of_scenarios/2, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_sets_f,
				            			support_constraints, "Risk-averse", "Risk-averse", delta_l,  deltas_f, average_l, averages_f, M1);
						        
					  	     finish = System.currentTimeMillis();
					  	
						  	 if (sol_sp_5.size() == 0)
						  	 {
						  		 for (int l1 = l; l1 < max_index; l1++)
				            		 ind2.put(l1, 0);
						  		 
						  		 start = 0;
				            	 finish = -1;
					  		 }
			             }
			             
					  	 times_sp_5.put(instance, (finish-start)*Math.pow(10, -3));
					  	 
					  	 start = 0;
			             finish = -1;
					  	 
					  	 if (ind3.get(l) == 1)
			             {
				             start = System.currentTimeMillis();
								   			        
				             ArrayList<Double> sol_sp_10 =  Instance.SolveLeadersProblemScenarios(cplex, I, number_of_scenarios, epsilon_l, epsilon_f, alpha_l, alpha_f, k_l, k_f, data_set_l, data_sets_f,
				            			support_constraints, "Risk-averse", "Risk-averse", delta_l,  deltas_f, average_l, averages_f, M1);
						        
					  	     finish = System.currentTimeMillis();
					  	
						  	 if (sol_sp_10.size() == 0)
						  	 {
						  		 for (int l1 = l; l1 < max_index; l1++)
				            		 ind3.put(l1, 0);
						  		 
						  		 start = 0;
				            	 finish = -1;
					  		 }
			             }
			             
					  	 times_sp_10.put(instance, (finish-start)*Math.pow(10, -3));
				        }
		        	
		        	    /*System.out.println(ind1);
		        	    System.out.println(ind2);
		        	    System.out.println(ind3);*/
			        }
			           
		      
		        //Compute mean and variance
		        ArrayList<Double> mean_time_b = new ArrayList<Double>();
		        ArrayList<Double> var_time_b = new ArrayList<Double>();
		        
		        ArrayList<Double> mean_time_sp_5 = new ArrayList<Double>();
		        ArrayList<Double> var_time_sp_5 = new ArrayList<Double>();
		        
		        ArrayList<Double> mean_time_sp_10 = new ArrayList<Double>();
		        ArrayList<Double> var_time_sp_10 = new ArrayList<Double>();
		        
		        
		        int num_param = 1;
		        for (int l = 0; l < max_index; l++)
		        {
		        	
		        	ArrayList<Double> times_param_b = new ArrayList<Double>();
		        	ArrayList<Double> times_param_sp_5 = new ArrayList<Double>();
		        	ArrayList<Double> times_param_sp_10 = new ArrayList<Double>();
		        	
		        	for (int t = 1; t <= number_of_instances_1; t++)
		        	{
		        		ArrayList<Integer> instance = new ArrayList<Integer>();
			            instance.add(t);
			            instance.add(num_param);
			            
			            times_param_b.add(times_b.get(instance));  
			            times_param_sp_5.add(times_sp_5.get(instance));  
			            times_param_sp_10.add(times_sp_10.get(instance));  
			             
		        	}
		        	
		        	
		        	double mean_b = Instance.AverageOfArray(times_param_b);
		        	double var_b = Instance.VarianceOfArray(times_param_b);
		        	
		        	double mean_sp_5 = Instance.AverageOfArray(times_param_sp_5);
		        	double var_sp_5 = Instance.VarianceOfArray(times_param_sp_5);
		        	
		        	double mean_sp_10 = Instance.AverageOfArray(times_param_sp_10);
		        	double var_sp_10 = Instance.VarianceOfArray(times_param_sp_10);
		        
		        	mean_time_b.add(mean_b);
		        	var_time_b.add(var_b);
		        	
		        	mean_time_sp_5.add(mean_sp_5);
		        	var_time_sp_5.add(var_sp_5);
		        	
		        	mean_time_sp_10.add(mean_sp_10);
		        	var_time_sp_10.add(var_sp_10);
		        	
		        	num_param ++;
		        }
		        
		        FileWriter fileWriter = new FileWriter("Experiment_10.txt");
			       
		        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
		        {
		        	printWriter.printf("%d distributions \n", number_of_instances_1); 
		        	printWriter.println();
		        	
		        	printWriter.printf("Parameters: n = %d, max_noise = %.2f, const_l = %.2f, const_f = %.2f \n", n, max_noise, const_l, const_f); 
		        	
		        	printWriter.println();
		        	
		        	num_param = 1;
					
		        	for (int l = 0; l < max_index; l++)
				     {
				        int size = sizes[l];
				      
				        printWriter.printf("n = %d : %.2f (%.2f) | %.2f (%.2f) | %.2f (%.2f) \n", size, mean_time_b.get(num_param - 1)*ind1.get(l), var_time_b.get(num_param - 1)*ind1.get(l), mean_time_sp_5.get(num_param - 1)*ind2.get(l), var_time_sp_5.get(num_param - 1)*ind2.get(l), mean_time_sp_10.get(num_param - 1)*ind3.get(l), var_time_sp_10.get(num_param - 1)*ind3.get(l));
				        num_param++;
					 }
					  
					printWriter.flush();
					printWriter.close();
				}
	        }
	     
//***************************************************************************************
			  	 //Experiment 11 (solution times pessimistic k_l, k_f)        
			     //num_of_experiment = 11;     
			       if (num_of_experiment == 11)
			        { 
			        	System.out.println("Experiment 11: \n");
			         
				        int[] sample_sizes = {10, 15, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200, 250, 300};
				        int max_index = sample_sizes.length;
			        	
			        	Map<Integer, Integer> ind1 = new HashMap<Integer, Integer>();
			        	Map<Integer, Integer> ind2 = new HashMap<Integer, Integer>();
			        	Map<Integer, Integer> ind3 = new HashMap<Integer, Integer>();
			        	
			        	for (int l = 0; l < max_index; l++)
			        	{
			        		ind1.put(l, 1);
			        		ind2.put(l, 1);
			        		ind3.put(l, 1);
			        	}
			        	
				        double alpha_l_mod = 0.9;
		  		        
				        Map<ArrayList<Integer>, Double> times_paf = new HashMap<ArrayList<Integer>, Double>();
				        Map<ArrayList<Integer>, Double> times_prn = new HashMap<ArrayList<Integer>, Double>();
				       
				        for (int t1 = 0; t1 < number_of_instances_1; t1++)
				        {  
				        	int number_of_instance = t1 + 1;
					        System.out.println(number_of_instance);	
					        
					        Instance I = random_instances.get(t1);
					        Instance I_u = random_instances_u.get(t1);
					        
					        ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
					        int number_of_sc = support_constraints.size() - 1; 
					        
					    
					        //Main parameter loop
		                    int num_param = 1;
		                    for (int l = 0; l < max_index; l++)
		    		        {
		    		        	
					        	 ArrayList<Integer> instance = new ArrayList<Integer>();
					             instance.add(number_of_instance);
					             instance.add(num_param);
					             
					             System.out.println(num_param);
					             num_param ++; 
					            
					             k_l = sample_sizes[l];
					           
					             const_l = 0.1;
							     epsilon_l = const_l/Math.sqrt(k_l);
					  			 
		    		        	 ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
		 				         
		 				         ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I_u, data_set_l, support_constraints);
		   				        
							     ArrayList<Double> average_l = new ArrayList<Double>();    
							     
							       for (int i = 0; i < I.number_of_items; i++)
							    	{
							    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
							    	} 
					             
							     
							     //Solve pessimistic formulation
					  		     double lower_bound = 0.;
					  			 double upper_bound = I.number_of_items; 
					  		     
					  			 M1 = Math.max(I.number_of_items, 1.);
					  			 
					  			 double start = 0;
					             double finish = -1;
							  	 
							  	 if (ind1.get(l) == 1)
					             {
						             start = System.currentTimeMillis();
										   			        
						             ArrayList<Double> sol_paf = Instance.SolveMasterProblemAmbiguityFree(cplex, cplex2, I, lower_bound, upper_bound, alpha_l_mod, k_l, data_set_l, M1, start);
							  			   
							  	     finish = System.currentTimeMillis();
							  	
								  	 if (sol_paf.size() == 0)
								  	 {
								  		 for (int l1 = l; l1 < max_index; l1++)
						            		 ind1.put(l1, 0);
								  		 
								  		 start = 0;
						            	 finish = -1;
							  		 }
					             }
					  			 	 
					  			 times_paf.put(instance, (finish-start)*Math.pow(10, -3));
					             
							     //Solve pessimistic formulation
					  		     double M_objective = I.number_of_items;
					  			  
					  		     start = 0;
					             finish = -1;
							  	 
							  	 if (ind2.get(l) == 1)
					             {
						             start = System.currentTimeMillis();
										   			        
						             ArrayList<Double> sol_prn = Instance.SolveMasterProblemRiskNeutral(cplex, cplex2, I_u, lower_bound, upper_bound, support_constraints, epsilon_l, k_l, data_set_l, delta_l, average_l, M_objective, binary_indicator, start);
							  				   
							  	     finish = System.currentTimeMillis();
							  	
								  	 if (sol_prn.size() == 0)
								  	 {
								  		 for (int l1 = l; l1 < max_index; l1++)
						            		 ind2.put(l1, 0);
								  		 
								  		 start = 0;
						            	 finish = -1;
							  		 }
					             }
					  		     
		                         times_prn.put(instance, (finish-start)*Math.pow(10, -3));
					             
					        }
					     }
					              
				      
				        //Compute mean and variance
				        ArrayList<Double> mean_time_paf = new ArrayList<Double>();
				        ArrayList<Double> var_time_paf = new ArrayList<Double>();
				        
				        ArrayList<Double> mean_time_prn = new ArrayList<Double>();
				        ArrayList<Double> var_time_prn = new ArrayList<Double>();
				         
				        int num_param = 1;
				        for (int l = 0; l < max_index; l++)
				        {
				        	
				        	ArrayList<Double> times_param_paf = new ArrayList<Double>();
				        	ArrayList<Double> times_param_prn = new ArrayList<Double>();
				        	
				        	for (int t = 1; t <= number_of_instances_1; t++)
				        	{
				        		ArrayList<Integer> instance = new ArrayList<Integer>();
					            instance.add(t);
					            instance.add(num_param);
					            
					            times_param_paf.add(times_paf.get(instance));  
					            times_param_prn.add(times_prn.get(instance));  
					             
				        	}
				        	
				        	
				        	double mean_paf = Instance.AverageOfArray(times_param_paf);
				        	double var_paf = Instance.VarianceOfArray(times_param_paf);
				        	
				        	double mean_prn = Instance.AverageOfArray(times_param_prn);
				        	double var_prn = Instance.VarianceOfArray(times_param_prn);
				        	
				        	mean_time_paf.add(mean_paf);
				        	var_time_paf.add(var_paf);
				        	
				        	mean_time_prn.add(mean_prn);
				        	var_time_prn.add(var_prn);
				   
				        	num_param ++;
				        }
				        
				        FileWriter fileWriter = new FileWriter("Experiment_11.txt");
					       
				        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
				        {
				        	printWriter.printf("%d distributions \n", number_of_instances_1); 
				        	printWriter.println();
				        	
				        	printWriter.printf("Parameters: n = %d, k_f = %d, const_l = %.2f \n", n, k_f, const_l); 
				        	
				        	printWriter.println();
				        	
				        	num_param = 1;
							
				        	for (int l = 0; l < max_index; l++)
						     {
						        k_l = sample_sizes[l];
						      
						        printWriter.printf("k_l = %d : %.2f (%.2f) | %.2f (%.2f) \n", k_l, mean_time_paf.get(num_param - 1)*ind1.get(l), var_time_paf.get(num_param - 1)*ind1.get(l), mean_time_prn.get(num_param - 1)*ind2.get(l), var_time_prn.get(num_param - 1)*ind2.get(l));
						        num_param++;
							 }
							  
							printWriter.flush();
							printWriter.close();
						}
			        }
			     
//***************************************************************************************
			  	 //Experiment 12 (solution times pessimistic k_l, k_f)        
			     //num_of_experiment = 12;     
			        if (num_of_experiment == 12)
			        { 
			        	System.out.println("Experiment 12: \n");
			         
			        	int[] sizes = {5, 6, 7, 8, 10, 12, 14, 17, 20, 24, 28, 33, 39, 46, 50};
			        	int max_index = sizes.length;
				        
			        	Map<Integer, Integer> ind1 = new HashMap<Integer, Integer>();
			        	Map<Integer, Integer> ind2 = new HashMap<Integer, Integer>();
			        	Map<Integer, Integer> ind3 = new HashMap<Integer, Integer>();
			        	
			        	for (int l = 0; l < max_index; l++)
			        	{
			        		ind1.put(l, 1);
			        		ind2.put(l, 1);
			        		ind3.put(l, 1);
			        	}
			        	
			        	k_l = 30;
			        	k_f = 30;
			        	
			        	const_l = 0.1;
					    epsilon_l = const_l/Math.sqrt(k_l);
				        
			        	double alpha_l_mod = 0.9;
		  		  
				        Map<ArrayList<Integer>, Double> times_paf = new HashMap<ArrayList<Integer>, Double>();
				        Map<ArrayList<Integer>, Double> times_prn = new HashMap<ArrayList<Integer>, Double>();
				       
				        for (int t1 = 0; t1 < number_of_instances_1; t1++)
				        {
				        	int number_of_instance = t1 + 1;
				        	System.out.println(number_of_instance);	
				        
					        int num_param = 1;
				        	for (int l = 0; l < max_index; l++)
				        	{
				        		long instanceSeed = globalRandom.nextLong();
					        	Random random = new Random(instanceSeed);
				        		Instance I = new Instance(sizes[l], p, d_f, w_l, w_f, 0, random); 
				        		Instance I_u = new Instance(sizes[l], p, w_l, w_f, random);
					           
				        		ArrayList<ArrayList<Double>> support_constraints = Instance.ConstructSupportConstraints(I);
					            int number_of_sc = support_constraints.size() - 1; 
				        
					            ArrayList<Integer> instance = new ArrayList<Integer>();
					            instance.add(number_of_instance);
					            instance.add(num_param);
					             
					            System.out.println(num_param);
					            num_param ++; 
					        
					            ArrayList<ArrayList<Double>> data_set_l = Instance.GenerateData(I, k_l); 
		 				        ArrayList<ArrayList<Double>> delta_l = Instance.ComputeDelta(I_u, data_set_l, support_constraints);
		   				        
						        ArrayList<Double> average_l = new ArrayList<Double>();    
						     
						        for (int i = 0; i < I.number_of_items; i++)
						    	{
						    	    average_l.add(Instance.AverageOfArray(data_set_l.get(i)));	
						    	} 
					             
							     
							    //Solve pessimistic formulation
					  		    double lower_bound = 0.;
					  			double upper_bound = I.number_of_items; 
					  		     
					  			M1 = Math.max(I.number_of_items, 1.);
					  			 
					  			double start = 0;
					            double finish = -1;
							  	 
							  	if (ind1.get(l) == 1)
					            {
						             start = System.currentTimeMillis();
										   			        
						             ArrayList<Double> sol_paf = Instance.SolveMasterProblemAmbiguityFree(cplex, cplex2, I, lower_bound, upper_bound, alpha_l_mod, k_l, data_set_l, M1, start);
							  				   
							  	     finish = System.currentTimeMillis();
							  	
								  	 if (sol_paf.size() == 0)
								  	 {
								  		 for (int l1 = l; l1 < max_index; l1++)
						            		 ind1.put(l1, 0);
								  		 
								  		 start = 0;
						            	 finish = -1;
							  		 }
					            }
					  			
					  			times_paf.put(instance, (finish-start)*Math.pow(10, -3));
					             
							    //Solve pessimistic formulation
					  		    double M_objective = I.number_of_items;
					  			
					  		    start = 0;
					            finish = -1;
							  	 
							  	if (ind2.get(l) == 1)
					            {
						             start = System.currentTimeMillis();
										   			        
						             ArrayList<Double> sol_prn = Instance.SolveMasterProblemRiskNeutral(cplex, cplex2, I_u, lower_bound, upper_bound, support_constraints, epsilon_l, k_l, data_set_l, delta_l, average_l, M_objective, binary_indicator, start);
							  					   
							  	     finish = System.currentTimeMillis();
							  	
								  	 if (sol_prn.size() == 0)
								  	 {
								  		 for (int l1 = l; l1 < max_index; l1++)
						            		 ind2.put(l1, 0);
								  		 
								  		 start = 0;
						            	 finish = -1;
							  		 }
					            }
					  		     
		                        times_prn.put(instance, (finish-start)*Math.pow(10, -3));    
					        }
					     }
					              
				      
				        //Compute mean and variance
				        ArrayList<Double> mean_time_paf = new ArrayList<Double>();
				        ArrayList<Double> var_time_paf = new ArrayList<Double>();
				        
				        ArrayList<Double> mean_time_prn = new ArrayList<Double>();
				        ArrayList<Double> var_time_prn = new ArrayList<Double>();
				         
				        int num_param = 1;
				        for (int l = 0; l < max_index; l++)
				        {
				        	
				        	ArrayList<Double> times_param_paf = new ArrayList<Double>();
				        	ArrayList<Double> times_param_prn = new ArrayList<Double>();
				        	
				        	for (int t = 1; t <= number_of_instances_1; t++)
				        	{
				        		ArrayList<Integer> instance = new ArrayList<Integer>();
					            instance.add(t);
					            instance.add(num_param);
					            
					            times_param_paf.add(times_paf.get(instance));  
					            times_param_prn.add(times_prn.get(instance));  
					             
				        	}
				        	
				        	
				        	double mean_paf = Instance.AverageOfArray(times_param_paf);
				        	double var_paf = Instance.VarianceOfArray(times_param_paf);
				        	
				        	double mean_prn = Instance.AverageOfArray(times_param_prn);
				        	double var_prn = Instance.VarianceOfArray(times_param_prn);
				        	
				        	mean_time_paf.add(mean_paf);
				        	var_time_paf.add(var_paf);
				        	
				        	mean_time_prn.add(mean_prn);
				        	var_time_prn.add(var_prn);
				   
				        	num_param ++;
				        }
				        
				        FileWriter fileWriter = new FileWriter("Experiment_12.txt");
					       
				        try (PrintWriter printWriter = new PrintWriter(fileWriter)) 
				        {
				        	printWriter.printf("%d distributions \n", number_of_instances_1); 
				        	printWriter.println();
				        	
				        	printWriter.printf("Parameters: k_l = %d, k_f = %d, const_l = %.2f \n", k_l, k_f, const_l); 
				        	
				        	printWriter.println();
				        	
				        	num_param = 1;
							
				        	for (int l = 0; l < max_index; l++)
						     {
				        		int size = sizes[l];
						      
						        printWriter.printf("n = %d : %.2f (%.2f) | %.2f (%.2f) \n", size, mean_time_paf.get(num_param - 1)*ind1.get(l), var_time_paf.get(num_param - 1)*ind1.get(l), mean_time_prn.get(num_param - 1)*ind2.get(l), var_time_prn.get(num_param - 1)*ind2.get(l));
						        num_param++;
							 }
							  
							printWriter.flush();
							printWriter.close();
						}
			        }			  
	        
	  }
	  }
