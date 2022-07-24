using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;

namespace multiObjectiveSearch
{
	public class NSGA2
	{
		public NSGA2(string rs)
		{
			settings = new RunSettings(rs);
			
			RandomGenerator = new Random();
			pointRandomGenerator = new Random();
			crossoverRandomGenerator = new Random();
			mutationRandomGenerator = new Random();
			mutation_2_RandomGenerator = new Random();
			mutation_3_RandomGenerator = new Random();
			crossover_2_RandomGenerator = new Random();
		}
		private Random RandomGenerator = null, pointRandomGenerator = null,
				mutationRandomGenerator = null ,mutation_2_RandomGenerator = null, mutation_3_RandomGenerator = null, 
				crossoverRandomGenerator = null, crossover_2_RandomGenerator = null;
		private RunSettings settings = null;
		#region NSGA2 functions
		private List<chromosome> FastNondominatedSort(List<chromosome> population)
		{
			List<chromosome> z = new List<chromosome>();
			
			validatePopulation(ref population);
			//EvaluateAgainstObjectiveFunctions(population);
			
			int size = population.Count, M = settings.dimentions;
			int[] n_ = new int[size];
			List<int> [] p_ = new List<int> [size];
			List<List<int>> fronts = new List<List<int>>(0);
			fronts.Add(new List<int>(0));
			int front = 0;
			
			#region first scan
			for(int i = 0; i < size; i++)
			{
				n_[i] = 0;
				p_[i] = new List<int> (0);
				for(int j = 0; j < size; j++)
				{
					if(i == j) continue;
					
					int less = 0, more = 0, equal = 0;
					for(int k = 0; k < M; k++)
					{
						if(population[i].rank[k] > population[j].rank[k])
							less++;
						else if(population[i].rank[k] < population[j].rank[k])
							more++;
						else equal++;
					}
					if(less == 0 && equal != M)
						n_[i] ++;
					else if(more == 0 && equal != M)
						p_[i].Add(j);
				}
				if(n_[i] == 0)
				{
					population[i].totalRank = 1;
					fronts[front].Add(i);
				}
			}
			#endregion
			
			#region find subsequent fronts
			while(fronts[front].Count != 0)
			{
				List<int> Q = new List<int>(0);
				for(int i = 0; i < fronts[front].Count; i++)
				{
					if(p_[fronts[front][i]].Count > 0)
					{
						for(int j = 0; j < p_[fronts[front][i]].Count; j++)
						{
							n_[p_[fronts[front][i]][j]] = n_[p_[fronts[front][i]][j]] - 1;
							if(n_[p_[fronts[front][i]][j]] == 0)
							{
								population[p_[fronts[front][i]][j]].totalRank = front + 1;
								Q.Add(p_[fronts[front][i]][j]);
							}
						}
					}
				}
				front ++;
				if(Q.Count == 0)
					break;
				fronts.Add(new List<int>(0));
				fronts[front] = Q;
			}
			
			#endregion
			
			#region sorting population based on fronts on Ascented format
			int[] fronts_sorted_index = new int[size];
			for(int i = 0; i < size; i++)
				fronts_sorted_index[i] = i;
			for(int i = 0; i < size; i++)
			{
				for(int j = 1; j < size; j++)
				{
					if(population[j].totalRank > population[j - 1].totalRank)
					{
						int temp = fronts_sorted_index[j];
						fronts_sorted_index[j] = fronts_sorted_index[j - 1];
						fronts_sorted_index[j - 1] = temp;
					}
				}
			}
			#endregion
			
			#region Crowding Distance
			int current_index = 0;
			for(front = 0; front < fronts.Count; front++)
			{
				List<chromosome> y = new List<chromosome>(0);
				for(int i = 0; i < fronts[front].Count; i++)
					y.Add(population[fronts_sorted_index[i + current_index]]);
				current_index += fronts[front].Count;
				for(int i = 0; i < settings.dimentions; i++)
				{
					#region sort y based on objective[i]
					int[] obj_sorted_index = new int[y.Count];
					for(int a = 0; a < y.Count; a++)
						obj_sorted_index[a] = a;
					for(int a = 0; a < y.Count; a++)
					{
						for(int b = 1; b < y.Count; b++)
						{
							if(y[obj_sorted_index[b]].rank[i] < y[obj_sorted_index[b - 1]].rank[i])
							{
								int temp = obj_sorted_index[b];
								obj_sorted_index[b] = obj_sorted_index[b - 1];
								obj_sorted_index[b - 1] = temp;
							}
						}
					}
					List<chromosome> y2 = new List<chromosome>(0);
					for(int a = 0; a < fronts[front].Count; a ++)
						y2.Add(y[obj_sorted_index[a]]);
					#endregion
					
					double fmax = y2[y2.Count - 1].rank[i];
					double fmin = y2[0].rank[i];
					y[obj_sorted_index[y.Count -1]].CrowdingDistance[i] = int.MaxValue;
					y[obj_sorted_index[0]].CrowdingDistance[i] = int.MaxValue;
					
					for(int j = 1; j < y.Count - 1; j ++)
					{
						double next_obj = y2[j + 1].rank[i];
						double prev_obj = y2[j - 1].rank[i];
						if(fmax - fmin == 0)
							y[obj_sorted_index[j]].CrowdingDistance[i] = int.MaxValue;
						else
							y[obj_sorted_index[j]].CrowdingDistance[i] = (next_obj - prev_obj) / (fmax - fmin);
					}
				}
				
				for(int i = 0; i < y.Count; i++)
				{
					//double tmp = 0;
					double maxd = 0;
					for(int j = 0; j < settings.dimentions; j++)
					{
						//tmp += Math.Pow (y[i].CrowdingDistance[j], 2);
						y[i].Distance += y[i].CrowdingDistance[j];
						if(y[i].CrowdingDistance[j] > maxd) maxd = y[i].CrowdingDistance[j];
					}
					//y[i].Distance = Math.Pow (tmp, 0.5);
					//y[i].Distance = Median(ref y[i].CrowdingDistance);
					//y[i].Distance /= settings.dimentions;
					y[i].Distance = maxd;
				}
				for(int i = 0; i < y.Count; i++)
					z.Add(y[i]);
			}
			#endregion
			
			return z;
		}
		private List<chromosome> TournomentSelection(List<chromosome> population)
		{				
			int population_size = population.Count;
			int poolsize = population_size / 2;
			int toursize = 2;
			List<chromosome> pool = new List<chromosome>();
			int[] selectedmask = new int[population_size];
			Random r = new Random();
			
			for(int i = 0; i < poolsize; i++)
			{
				int[] candidate= new int[toursize];
				double[] objrank = new double[toursize];
				double[] objdist = new double[toursize];
				
				//double minrank = double.MaxValue;
				double minrank = 0;
				int minidx = -1;
				for(int j = 0; j < toursize; j++)
				{
					candidate[j] = Convert.ToInt32 (r.NextDouble() * population_size - 1);
					if(candidate[j] == -1) candidate[j] = 0;
					while(selectedmask[candidate[j]] == 1)
						candidate[j] = (candidate[j] + 1) % population_size;
					selectedmask[candidate[j]] = 1;
					
					objrank[j] = population[candidate[j]].totalRank;
					objdist[j] = population[candidate[j]].Distance;
					
					if(objrank[j] > minrank)
					{
						minrank = objrank[j];
						minidx = j;
					}
				}
				
				//if there is more than 1 obj with min rank, get the 1st max dist
				int mincnt = 0;
				for(int j = 0; j < toursize; j++)
					if(objrank[j] == minrank)
						mincnt++;
				if(mincnt > 1)
				{
					int maxidx = -1;
					double maxdist = double.MinValue;
					for(int j = 0; j < toursize; j++)
					{
						if(objrank[j] == minrank && objdist[j] > maxdist)
						{
							maxdist = objdist[j];
							maxidx = j;
						}
					}
					if(maxdist == -1)
						throw new Exception("این چه وضعشه آخه!!؟؟؟");
					pool.Add(population[candidate[maxidx]]);
				}
				else
					pool.Add(population[candidate[minidx]]);
			}
			
			return pool;
		}
		private List<chromosome> CrossowerAndMutation(List<chromosome> population, int population_size,
		                                           double p_crossower, double p_mutation, int Mutation_MAXDIST)
		{
			//int needed_pop = population_size - population.Count;
			
			//crate 'population_size' childrens
			List <chromosome> children = new List<chromosome>();
			if(population.Count > 2)
			{
				for(int i = 0; i < population.Count; i++)
				{
					double r2 = crossoverRandomGenerator.NextDouble ();
					if(r2 > p_crossower)
						continue;
					int nextp = 0;
					do { nextp = crossover_2_RandomGenerator.Next (0, population.Count); }
					while(nextp == i );
					chromosome c1 = new chromosome(population[i], ref settings);
					chromosome c2 = new chromosome(population[nextp], ref settings);
					
					c1.y = population[nextp].y;
					c2.y = population[i].y;
					
					children.Add (c1);
					if(children.Count >= population_size)
						break;
					children.Add (c2);
					if(children.Count >= population_size)
						break;
				}
			}
			//mutation
			int size2 = children.Count;
			for(int i = 0; i < size2; i++)
			{
				double R = mutationRandomGenerator.NextDouble();
				if(R < p_mutation)
				{
					chromosome tmp = new chromosome(children[i], ref settings);
					int difdist = mutation_2_RandomGenerator.Next (0, Mutation_MAXDIST);
					double d2 = mutation_3_RandomGenerator.NextDouble();
					if(d2 < 0.25)
						tmp.x += difdist;
					else if(d2 < 0.5)
						tmp.x -= difdist;
					if(d2 < 0.75)
						tmp.y += difdist;
					else
						tmp.y -= difdist;
					if(tmp != null)
						children.Add(tmp);
				}
			}
			
			return removerepetitiveanswers(children);
			//return removeWeeklyDivergAnswers(ref childs, benchmarkindex);
		}
		private List<chromosome> Merge(List<chromosome> children, List<chromosome> population)
		{
			int childsize = children.Count;
			int newsize = childsize + population.Count;
			List<chromosome> union = new List<chromosome>();
			for(int i = 0; i < childsize; i++)
				union.Add(children[i]);
			for(int i = childsize; i < newsize; i++)
				union.Add(population[i - childsize]);
			return union;
		}
		private List<chromosome> removerepetitiveanswers(List<chromosome> intermediate_chromosome)
		{
			if(intermediate_chromosome == null || intermediate_chromosome.Count == 0)
				return null;
			
			int size = intermediate_chromosome.Count;
			List<string> chrstr = new List<string>(0);
			chrstr.Add(intermediate_chromosome[0].PrintRawString(","));
			int reps = 0;
			int []repsa = new int[size];
			for(int i = 1; i < size; i++)
			{
				string str = intermediate_chromosome[i].PrintRawString(",");
				for(int j = 0; j < chrstr.Count; j++)
				{
					if(chrstr[j] == str)
					{
						reps++;
						repsa[i] = 1;
						break;
					}
				}
				if(repsa[i] != 1)
					chrstr.Add(str);
			}
			//Console.Write(" {0} repetetive answers found.", reps);
			List<chromosome> tmppopnr = new List<chromosome>(0);
			for(int i = 0; i < size; i++) 
				if(repsa[i] != 1) 
					tmppopnr.Add (intermediate_chromosome[i]);
			size -= reps;
			List<chromosome> intermediate_chromosome2 = new List<chromosome>();
			for(int i = 0; i < size; i++) 
				intermediate_chromosome2.Add(tmppopnr[i]);
			return intermediate_chromosome2;
		}
		
		private List<chromosome> GenerateExtraPopulation(int PopulationSize)
		{
			List<chromosome> population = new List<chromosome>();
			
			for(int i = 0; i < PopulationSize; i++)
			{
				double x = pointRandomGenerator.Next(Convert.ToInt32(Math.Round(settings.region.Left)), 
				                                     Convert.ToInt32(Math.Round(settings.region.Right)));
				double y = pointRandomGenerator.Next(Convert.ToInt32(Math.Round(settings.region.Top)),
				                                     Convert.ToInt32(Math.Round(settings.region.Bottom)));
				population.Add(new chromosome(x, y, ref settings));
			}
			
			return population;
		}
		private void validatePopulation(ref List<chromosome> pop)
		{
			for(int i = 0; i < pop.Count; i++)
			{
				if(i > pop.Count - 1)
					return;
				int d = pop[i].ValidatePoint();
				if(d < 0)
				{
					pop.RemoveAt(i);
					i--;
					continue;
				}
				pop[i].EvaluatePoint();
			}
		}
		#endregion
		
		public List<chromosome> SearchDesignSpace()
		{
			#region reading settings
			int t = 0;
			int Maxt = settings.ReadValueOfInteger(settings.RootAddress, "NSGA2_InnerLoopMaxRuns");
			int PopulationSize = settings.ReadValueOfInteger(settings.RootAddress, "NSGA2_population_size");
			double P_Crossover = settings.ReadValueOfDouble(settings.RootAddress, "NSGA2_P_Crossover");
			double P_Mutation = settings.ReadValueOfDouble(settings.RootAddress, "NSGA2_P_Mutation");
			int MAX_NSGA2_Runs = settings.ReadValueOfInteger (settings.RootAddress, "NSGA2_OuterLoopMaxRuns") - 1;
			bool initial_import_type = settings.ReadValueOfString(settings.RootAddress, "initial_import_type")
				== "add";
			#endregion
			
			List<chromosome> pop_collection = new List<chromosome>();
			for(int NSGA2_Counter = 0; NSGA2_Counter  < MAX_NSGA2_Runs; NSGA2_Counter ++)
			{
				Console.WriteLine("\n\nrunning NSGA2 for " + (NSGA2_Counter  + 1).ToString() + "th times :");
				
				//initialize population
				#region main loop
				//List<chromosome> parents = GenerateExtraPopulation(PopulationSize);
				List<chromosome> parents = settings.ReadInitialPoints(settings.RootAddress, "initial_points", ref settings);
				if(parents == null)
					parents.AddRange(GenerateExtraPopulation(PopulationSize));
				else if(initial_import_type)
					parents.AddRange(GenerateExtraPopulation(PopulationSize - parents.Count));
				validatePopulation(ref parents);
				// running for a fixed iterations
				t = 0;
				while(t < Maxt)
				{
					if(parents == null || parents.Count == 0)
						break;
					List<chromosome>selected = FastNondominatedSort(parents);
					if(selected == null || selected.Count == 0)
						break;
					selected = TournomentSelection(selected);
					if(selected == null || selected.Count == 0)
						break;
                    t++;
					List<chromosome> ltmp = FastNondominatedSort(selected);
					Console.WriteLine(
						  "iteration {0}: population-size = {1}, best answer : [{2},{3},{4},{5},{6}]",
                          t, PopulationSize,
		                  Convert.ToInt32(Math.Round(ltmp[0].rank[0])).ToString(),
		                  Convert.ToInt32(Math.Round(ltmp[0].rank[1])).ToString(),
		                  Convert.ToInt32(Math.Round(ltmp[0].rank[2])).ToString(),
		                  Convert.ToInt32(Math.Round(ltmp[0].rank[3])).ToString(),
		                  Convert.ToInt32(Math.Round(ltmp[0].rank[4])).ToString());
					
					List<chromosome> children = CrossowerAndMutation(selected, PopulationSize / 2, P_Crossover, P_Mutation, 3000000);
					List<chromosome> randomPop = new List<chromosome>();
					if(children == null || children.Count == 0)
						randomPop = GenerateExtraPopulation(PopulationSize - selected.Count);
					else
						randomPop = GenerateExtraPopulation(PopulationSize - selected.Count - children.Count);
					//List<chromosome> randomPop = GenerateExtraPopulation(children.Count / 2);
					List<chromosome> union = new List<chromosome>();
					if(randomPop == null || randomPop.Count == 0)
						continue;
					else
						union = Merge(Merge(children, selected), randomPop);
					parents = removerepetitiveanswers(union);
					if(union == null || union.Count == 0)
						break;
					validatePopulation(ref parents);
					     
					
				}
				#endregion
				pop_collection.AddRange(parents);
				pop_collection = removerepetitiveanswers(pop_collection);
			}
			
			//select best answer from all iteration
			return TournomentSelection(FastNondominatedSort(pop_collection));
			
			//return pop_collection;
		}
	}
	
}

