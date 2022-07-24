using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using Catfood.Shapefile;

namespace multiObjectiveSearch
{
	public class HypE
	{
		#region basic area
		public HypE(string rs)
		{
			HyPE_drand_RandomGenerator = new Random();
			pointRandomGenerator = new Random();
			mutationRandomGenerator = new Random();
			mutation_2_RandomGenerator = new Random();
			mutation_3_RandomGenerator = new Random();
			crossoverRandomGenerator = new Random();
			crossover_2_RandomGenerator = new Random();
			
			settings = new RunSettings(rs);
		}
		private Random HyPE_drand_RandomGenerator = null, pointRandomGenerator = null,
				mutationRandomGenerator = null, mutation_2_RandomGenerator = null, mutation_3_RandomGenerator = null,
				crossoverRandomGenerator = null, crossover_2_RandomGenerator = null;
		private RunSettings settings = null;
		#endregion 
		
		private enum comp {a_better_b, b_better_a, incomparable, indifferent};
			
		#region HypE functions
		List<List<chromosome>> generateFrontPartition(List<chromosome> pop_all)
		{
			List<List<chromosome>> front_part = new List<List<chromosome>>();
			
			int actFront = 0;
			int notDominated;
			int[] checkeda = new int[pop_all.Count];
			int[] actcheckeda = new int[pop_all.Count];
			int added;
		
			for(int i = 0; i < pop_all.Count; i++)
				checkeda[i] = 0;
		
			int j = 0;
			added = 0;
			actFront = 0;
			while( added < pop_all.Count )
			{
				for(int i = 0; i < pop_all.Count; i++ )
				{
					if( checkeda[i] == 1 )
						continue;
					actcheckeda[i] = 0;
					notDominated = 1;
					j = 0;
					while( notDominated == 1 && j < pop_all.Count )
					{
						if( i != j && checkeda[j] == 0 && dominates(j, i, pop_all) )
							notDominated = 0;
						j++;
					}
					if( notDominated == 1)
					{
						actcheckeda[i] = 1;
						if(front_part.Count < actFront + 1) front_part.Add (new List<chromosome>());
						front_part[actFront].Add(pop_all[i]);
						added++;
					}
				}
				for( j = 0; j < pop_all.Count; j++ )
					if( actcheckeda[j] == 1 )
						checkeda[j] = 1;
				actFront++;
			}
			
			return front_part;
		}
		bool dominates(int a, int b, List<chromosome>pop_all)
		{
			int i;
			bool a_is_worse = false;
			bool equal = true;
			int dim = settings.dimentions;
			
			for (i = 0; i < dim && !a_is_worse; i++)
			{
				a_is_worse = pop_all[a].rank[i] > pop_all[b].rank[i];
				equal = (pop_all[a].rank[i] == pop_all[b].rank[i]) && equal;
			}
		
			return (!equal && !a_is_worse);
		}
		List<List<chromosome>> environmentalSelection(List<List<chromosome>> front_part, int populationsize, double bound, int nrOfSamples)
		{
			/** Start with front wise reduction */
			for(int i = front_part.Count - 1; i >= 0; i--)
			{
				if(front_part.Count - front_part[i].Count >= populationsize)
					front_part.RemoveAt(i);
				else
					break;
			}
			/** Then remove from worst front */
			if(front_part.Count > populationsize)
				front_part[front_part.Count - 1] = hypeReduction(ref front_part, front_part[front_part.Count - 1], populationsize, bound, nrOfSamples);
			
			int size = 0;
			for(int i = 0; i < front_part.Count; i++) size += front_part[i].Count;
			//if(size != populationsize)
			//	throw new Exception("err in function 'environmentalSelection(List<List<chromosome>> front_part, int populationsize)' about size : \n\tfront_part.Count != populationsize");
			
			return front_part;
		}
		//Iteratively remove individuals from front based on sampled hypeIndicator value
		List<chromosome> hypeReduction(ref List<List<chromosome>> part_p, List<chromosome> fp, int alpha, double bound, int nrOfSamples)
		{
			double[] val = new double[fp.Count];
			int dim = settings.dimentions;
			double[] points = new double[fp.Count * dim];
			double checkRate = 0;
			double minRate = -1;
			int sel = -1;
			int i;
		
			if(fp.Count == 0) throw new Exception("err in function 'hypeReduction' about 'fp.Count == 0'");
			
			getObjectiveArray(fp, ref points);
		
			int size = 0;
			for(int a = 0; a < part_p.Count; a++) size += part_p[a].Count;
			while(size > alpha)
			{
				hypeIndicator(ref val, fp.Count, 0.0, bound, nrOfSamples, size - alpha, ref points, bound);
		
				sel = -1;
				for( i = 0; i < fp.Count; i++)
				{
					checkRate = val[i];
					if(sel == -1  || checkRate < minRate) 
					{
						minRate = checkRate;
						sel = i;
					}
				}
				
				assert(sel >= 0, "hypeReduction", "sel >= 0");
				assert(sel < fp.Count, "hypeReduction", "sel < fp.Count");
				
				//memcpy( points + dim*sel, points + dim*(fp->size - 1), sizeof(double)*dim );
				int stratindex = dim * (fp.Count - 1);
				for(int a = stratindex; a < points.Length; a++)
					points[dim * sel + a - stratindex] = points[a];
				double[] points_ = new double[points.Length - stratindex + 1 + dim * sel];
				for(int a = 0; a < points_.Length; a++) points_[a] = points[a];
				points = points_;
				
				fp.RemoveAt (sel);
			}
			return fp;
		}
		void getObjectiveArray(List<chromosome> A, ref double[] pointArray)
		{
			int dim = settings.dimentions;
			for(int i = 0; i < A.Count; i++) for(int k = 0; k < dim; k++) pointArray[i * dim + k] = A[i].rank[k];
		}
		void getObjectiveArray(ref int[] indices, List<chromosome> A, ref double[] pointArray)
		{
			int dim = settings.dimentions;
			for(int i = 0; i < A.Count; i++) for(int k = 0; k < dim; k++) pointArray[i * dim + k] = A[indices[i]].rank[k];
		}
		void hypeIndicator(ref double[] val, int popsize, double lowerbound, double upperbound, int nrOfSamples, int param_k, ref double[] points, double bound)
		{
			int i,j;
			double[] rho = new double[param_k + 1];
			
			rho[0] = 0;
			for(i = 1; i <= param_k; i++)
			{
				rho[i] = 1.0 / (double)i;
				for(j = 1; j <= i - 1; j++) rho[i] *= (double)(param_k - j ) / (double)( popsize - j );
			}
			for(i = 0; i < popsize; i++) val[i] = 0.0;
		
			if(nrOfSamples < 0)
				hypeExact(ref val, popsize, lowerbound, upperbound, param_k, ref points, ref rho, bound);
			else
				hypeSampling(ref val, popsize, lowerbound, upperbound, nrOfSamples, param_k, ref points, ref rho);
		}
		void hypeExact(ref double[] val, int popsize, double lowerbound, double upperbound, int param_k, ref double[] points, ref double[] rho, double bound)
		{
			int i;
			int dim = settings.dimentions;
			double[] boundsVec = new double[dim];
			int[] indices = new int[popsize];
			for( i = 0; i < dim; i++ )
				boundsVec[i] = bound;
			for( i = 0; i < popsize; i++  )
				indices[i] = i;
		
			hypeExactRecursive(ref points, popsize, dim, popsize, dim - 1, ref boundsVec, ref indices, ref val, ref rho, param_k);
		}
		void hypeExactRecursive(ref double[] input_p, int pnts, int dim, int nrOfPnts, int actDim, ref double[] bounds, ref int[] input_pvec,
		                        ref double[] fitness, ref double[] rho, int param_k)
		{
			int i, j;
			double extrusion;
			int[] pvec = new int[pnts];
			double[] p = new double[pnts * dim];
			for(i = 0; i < pnts; i++)
			{
				fitness[i] = 0;
				pvec[i] = input_pvec[i];
			}
			for(i = 0; i < pnts * dim; i++)
				p[i] = input_p[i];
		
			rearrangeIndicesByColumn(ref p, nrOfPnts, dim, actDim, ref pvec);
		
			for(i = 0; i < nrOfPnts; i++)
			{
				if(i < nrOfPnts - 1)
					extrusion = p[ (pvec[i + 1]) * dim + actDim ] - p[ pvec[i] * dim + actDim ];
				else
					extrusion = bounds[actDim] - p[ pvec[i] * dim + actDim ];
		
				if(actDim == 0) 
				{
					if(i + 1 <= param_k)
						for(j = 0; j <= i; j++)
							fitness[pvec[j]] = fitness[pvec[j]] + extrusion * rho[i + 1];
				}
				else if(extrusion > 0)
				{
					double[] tmpfit = new double[pnts];
					hypeExactRecursive(ref p, pnts, dim, i+1, actDim-1, ref bounds, ref pvec, ref tmpfit, ref rho, param_k);
					for(j = 0; j < pnts; j++)
						fitness[j] += extrusion * tmpfit[j];
				}
			}
		}
		void rearrangeIndicesByColumn(ref double[] mat, int rows, int columns, int col, ref int[] ind)
		{
			const int  MAX_LEVELS  = 300;
			int[] beg = new int[MAX_LEVELS];
			int[] end = new int[MAX_LEVELS];
			int i = 0, L = 0, R = 0, swap = 0;
			double pref, pind;
			double[] ref_ = new double[rows];
				ref_[i] = mat[col + ind[i] * columns];
			i = 0;
			
			beg[0] = 0; end[0] = rows;
			while (i >= 0) 
			{
				L = beg[i];
				R = end[i]-1;
				if(L < R)
				{
					pref = ref_[L];
					pind = ind[L];
					while(L < R) 
					{
						while(ref_[ R ] >= pref && L < R)
							R--;
						if( L < R ) 
						{
							ref_[ L ] = ref_[ R ];
							ind[ L++] = ind[R];
						}
						while( ref_[L] <= pref && L < R )
							L++;
						if( L < R) 
						{
							ref_[ R ] = ref_[ L ];
							ind[ R--] = ind[L];
						}
					}
					ref_[ L ] = pref;
					ind[L] = (int)pind;
					beg[i + 1] = L + 1;
					end[i + 1] = end[i];
					end[i++] = L;
					if(end[i] - beg[i] > end[i-1] - beg[i-1])
					{
						swap = beg[i]; beg[i] = beg[i-1]; beg[i-1] = swap;
						swap = end[i]; end[i] = end[i-1]; end[i-1] = swap;
					}
				}
				else 
					i--;
			}
		}
		void hypeSampling(ref double[] val, int popsize, double lowerbound, double upperbound, int nrOfSamples, int param_k, ref double[] points, ref double[] rho )
		{
			assert(popsize >= 0, "hypeSampling", "popsize >= 0");
			assert(lowerbound <= upperbound, "hypeSampling", "lowerbound <= upperbound");
			assert(param_k >= 1, "hypeSampling", "param_k >= 1");
			assert(param_k <= popsize, "hypeSampling", "param_k <= popsize");
		
			int i, s, k;
			int[] hitstat = new int[popsize];
			int domCount;
			
			int dim = settings.dimentions;
			double[] sample = new double[dim];
			for(s = 0; s < nrOfSamples; s++)
			{
				for(k = 0; k < dim; k++) sample[ k ] = drand(lowerbound, upperbound);
		
				domCount = 0;
				for(i = 0; i < popsize; i++)
				{
					int startindex = i * dim;
					double[] p_ = new double[points.Length - startindex + 1];
					for(int q = startindex; q < points.Length; q++) p_[q - startindex] = points[q];
					if(weaklyDominates(ref p_, ref sample, dim) == 1)
					{
						domCount++;
						if( domCount > param_k ) break;
						hitstat[i] = 1;
					}
					else
						hitstat[i] = 0;
				}
				if(domCount > 0 && domCount <= param_k)
				{
					for(i = 0; i < popsize; i++)
						if(hitstat[i] == 1)
							val[i] += rho[domCount];
				}
			}
			for(i = 0; i < popsize; i++)
			{
				val[i] = val[i] * Math.Pow( (upperbound-lowerbound), dim ) / (double)nrOfSamples;
			}
		}
		double RAND_MAX = 0.9;
		double drand(double from_, double to)
		{
			double j;
			j = from_ + (double)( (to - from_) * HyPE_drand_RandomGenerator.NextDouble() / (RAND_MAX + 1.0));
			return (j);
		}
		int weaklyDominates(ref double[] point1, ref double[] point2, int no_objectives)
		{
			bool better = true;
			int i = 0;
		
			while(i < no_objectives && better)
			{
				better = point1[i] <= point2[i];
				i++;
			}
			return better == true ? 1 : 0;
		}		
		List<chromosome> cleanUpArchive(List<List<chromosome>> partp, List<chromosome> pop_all)
		{
			List<chromosome> new_pop = new List<chromosome>();
			int i, j;
		
			//List<chromosome> new_pop = create_pop(alpha + lambda, dim);
			for (i = partp.Count - 1; i >= 0; i--)
			{
				for (j = partp[i].Count - 1; j >= 0; j--)
				{
					chromosome c = partp[i][j];
					int index = -1;
					for(int q = 0; q < pop_all.Count; q++)
					{ 
						bool iss = true;
						if(c.x != pop_all[q].x || c.y != pop_all[q].y)
							iss = false;
						if(iss == true)
						{
							index = q;
							break;
						}
					}		
					new_pop.Add(pop_all[index]);
					pop_all.RemoveAt(index);
				}
			}
			pop_all.Clear ();
			pop_all = new_pop;
			
			return pop_all;
		}
		void matingSelection(List<List<chromosome>>front_part, ref List<chromosome> pp_all, ref List<chromosome> pp_sel, double bound, int mu,
		                     int tournament, int mating, int nrOfSamples)
		{
			int winner = 0, winnerFront = 0;
			int opponent = 0, opponentFront = 0;
			int i,j;
		
			if(mating == 1)
				hypeFitnessMating(bound, nrOfSamples, pp_all);
			else
				for(i = 0; i < pp_all.Count; i++)
					pp_all[i].totalRank = 0.0;
		
			pp_sel = new List<chromosome>();
			for (i = 0; i < mu; i++)
			{
				int sort = mating == 2 ? 0 : 1;
				determineIndexAndFront(front_part, pp_all, irand(pp_all.Count), ref winner, ref winnerFront, sort);
				assert(winner < pp_all.Count , "matingselection", "winner < pp_all->size");
				assert(winner >= 0, "matingselection", "winner >= 0");
				for (j = 1; j < tournament; j++)
				{
					sort = mating == 2 ? 0 : 1;
					determineIndexAndFront(front_part, pp_all, irand(front_part.Count), ref opponent, ref opponentFront, sort);
					assert(opponent < pp_all.Count, "matingSelection", "opponent < pp_all->size");
					assert(opponent >= 0, "matingSelection", "opponent < pp_all->size");
					if (opponentFront < winnerFront || (opponentFront == winnerFront && pp_all[opponent].totalRank > pp_all[winner].totalRank ))
					{
						winner = opponent;
						winnerFront = opponentFront;
					}
				}
				addToSelection(i, winner, ref pp_all, ref pp_sel);
			}
		}
		//Calculates the fitness of all individuals in pp_all based on the hypE indicator
		void hypeFitnessMating(double bound, int nrOfSamples, List<chromosome> pp_all)
		{
			int i;
			double[] val = new double[pp_all.Count];
			int dim = settings.dimentions;
			double[] points = new double[pp_all.Count * dim];
			int[] indices = new int[pp_all.Count];
			for(i = 0; i < pp_all.Count; i++) indices[i] = i;
			getObjectiveArray(ref indices, pp_all, ref points);
			hypeIndicator(ref val, pp_all.Count, 0, bound, nrOfSamples, pp_all.Count, ref points, bound);
			for(i = 0; i < pp_all.Count; i++)
				pp_all[i].totalRank = val[i];
		}
		// Determine the front and the referencing index of the nth individual
		void determineIndexAndFront(List<List<chromosome>> partp, List<chromosome> pp_all, int n, ref int index, ref int front, int sort)
		{
			assert(n >= 0, "determineIndexAndFront", "n >= 0");
			int size = 0;
			for(int q = 0; q < partp.Count; q++) size += partp[q].Count;
			assert(n <  size, "determineIndexAndFront", "n <  partp->size");
			if(sort == 0)
				index = n;
		
			int i = partp.Count - 1;
			while (i >= 0 && (n - partp[i].Count) >= 0)
			{
				n -= partp[i].Count;
				i--;
			}
			assert(i >= 0 ,"determineIndexAndFront", "i >= 0");
			assert(n >= 0 ,"determineIndexAndFront", "n >= 0");
			assert(n < partp[i].Count, "determineIndexAndFront", "n < partp->front_array[i].size");
			
			index = -1;
			for(int q = 0; q < pp_all.Count; q++)
			{
				bool iss = true;
				if(pp_all[q].x != partp[i][n].x || pp_all[q].y != partp[i][n].y)
					iss = false;
				if(iss == true)
				{
					index = q;
					break;
				}
			}
			//index = partp[i][n];
			front = i;
		}
		int irand(int range)
		{
			int j;
			j=(int) ((double)range * (double) HyPE_drand_RandomGenerator.NextDouble () / (RAND_MAX + 1.0));
			return (j);
		}
		void addToSelection(int i, int c, ref List<chromosome> pp_all, ref List<chromosome> pop_sel)
		{
			assert(0 <= i, "addToSelection", "0 <= i");
			assert(i <= pop_sel.Count, "addToSelection", "i < pp_sel.Count");
			assert(0 <= c, "addToSelection", "0 <= c");
			assert(c < pp_all.Count, "addToSelection", "c < pp_all.Count");
			if(pop_sel.Count <= i) 
			{
				pop_sel.Add(pp_all[c]);
				pop_sel[pop_sel.Count - 1].totalRank = pp_all[c].totalRank;
			}
			else
			{
				pop_sel[i] = pp_all[c];
				pop_sel[i].totalRank = pp_all[c].totalRank;
			}
		}
		#endregion
		
		#region genetic functions
		private List<chromosome> CrossowerAndMutation(ref List<chromosome> population,
		                                              double p_crossower, double p_mutation, double MUTATION_CHANGE)
		{
			//crate 'population_size' childrens
			List <chromosome> children = new List<chromosome>();
			int[] selectedmask = new int[population.Count];
			int selecteds = 0;
			for(int i = 0; i < population.Count; i++)
			{
				if(crossoverRandomGenerator.NextDouble() > p_crossower)
					continue;
				
				if(population[i] == null) continue;
				if(selectedmask[i] == 1) continue;
				if(selecteds >= population.Count - 1) continue;
				
				//crossower
				int nextp = 0;
				do { nextp = crossover_2_RandomGenerator.Next (0, population.Count); }
				while(nextp == i || selectedmask[nextp] == 1);
				
				selectedmask[i] = 1;
				selectedmask[nextp] = 1;
				selecteds += 2;
				
				chromosome c1 = population[i];
				chromosome c2 = population[nextp];
				chromosome c3 = new chromosome(c1.x, c2.y, ref settings);
				chromosome c4 = new chromosome(c2.x, c1.y, ref settings);
				children.Add(c3);
				children.Add(c4);
			}
			
			//mutation
			for(int i = 0; i < children.Count; i++)
			{
				if(mutationRandomGenerator.NextDouble() > p_mutation)
					continue;
				
				double r = mutation_2_RandomGenerator.NextDouble();
				if(r < 0.5)
				{
					if(r < 0.25)
						children[i].x += pointRandomGenerator.NextDouble() * MUTATION_CHANGE;
					else
						children[i].x -= pointRandomGenerator.NextDouble() * MUTATION_CHANGE;
				}
				else
				{
					if(r < 0.75)
						children[i].y += pointRandomGenerator.NextDouble() * MUTATION_CHANGE;
					else
						children[i].y -= pointRandomGenerator.NextDouble() * MUTATION_CHANGE;
				}
			}
			
			
			if(children.Count != 0) children = removerepetitiveanswers(ref children);
			return children;
		}
		#endregion
			
		#region tools
		
		private List<chromosome> GenerateNewPopulation(int PopulationSize)
		{
			List<chromosome> population = new List<chromosome>();
			double scalex = settings.region.Right - settings.region.Left;
			double scaley = settings.region.Bottom - settings.region.Top;
			for(int i = 0; i < PopulationSize; i++)
			{
				population.Add(new chromosome(
					pointRandomGenerator.NextDouble() * scalex + settings.region.Left,
					pointRandomGenerator.NextDouble() * scaley + settings.region.Top,
					ref settings));
			}
			return population;
		}
		private List<chromosome> Merge(List<chromosome> l1, List<chromosome> l2)
		{
			List<chromosome> l3 = new List<chromosome>();
			for(int i = 0; i < l1.Count; i++) l3.Add (l1[i]);
			for(int i = 0; i < l2.Count; i++) l3.Add (l2[i]);
			return l3;
		}
		private void assert(bool ap, string functionname, string app)
		{
			if(! ap) throw new Exception("Err on function'"+ functionname + "' about '" + app + "'");
		}
		private List<chromosome> removerepetitiveanswers(ref List<chromosome> intermediate_point)
		{
			int size = intermediate_point.Count;
			List<string> chrstr = new List<string>(0);
			chrstr.Add(intermediate_point[0].PrintCMDString());
			int reps = 0;
			int []repsa = new int[size];
			for(int i = 1; i < size; i++)
			{
				string str = intermediate_point[i].PrintCMDString();
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
			for(int i = 0; i < size; i++) if(repsa[i] != 1) tmppopnr.Add (intermediate_point[i]);
			return tmppopnr;
		}
		private List<chromosome> MortezaSelection(ref List<chromosome> pop_sel, ref double[] baserank, int minselesctsize, int maxselectsize)
		{
			double[,] v = new double[pop_sel.Count, settings.dimentions];
			for(int i = 0; i < pop_sel.Count; i++) for(int j = 0; j < settings.dimentions; j++) v[i, j] = pop_sel[i].rank[j] / baserank[j];
			
			List<chromosome> morethan1 = new List<chromosome>();
			int[] selectedmask = new int[pop_sel.Count];
			int seletedcounts = 0;
			double[] maxv = new double[pop_sel.Count];
			double[] minv = new double[pop_sel.Count];
			for(int i = 0; i < pop_sel.Count; i++)
			{
				double min = 100000, max = 0;
				for(int j = 0; j < settings.dimentions; j++)
				{
					if(v[i, j] < min) min = v[i, j];
					if(v[i, j] > max) max = v[i, j];
				}
				maxv[i] = max;
				minv[i] = min;
				if(min >= 1)
				{
					morethan1.Add(pop_sel[i]);
					selectedmask[i] = 1;
					seletedcounts ++;
					if(seletedcounts == maxselectsize)
						return morethan1;
				}
			}
			
			if(morethan1.Count < minselesctsize)
			{
				//bubble-sort based on max speedup
				int[] selectedmaxvindex = new int[pop_sel.Count];
				for(int i = 0; i < pop_sel.Count; i++) selectedmaxvindex[i] = i;
				for(int i = 0; i < pop_sel.Count; i++)
				{
					for(int j = 1; j < pop_sel.Count; j++)
					{
						if(i != j && maxv[selectedmaxvindex[j]] > maxv[selectedmaxvindex[i]])
						{
										 int tmp = selectedmaxvindex[j];
							selectedmaxvindex[j] = selectedmaxvindex[i];
							selectedmaxvindex[i] = tmp;
						}
					}
				}
				
				for(int i = 0; i < pop_sel.Count; i++)
				{
					if(selectedmaxvindex[i] == 0)
						morethan1.Add(pop_sel[selectedmaxvindex[i]]);
					if(morethan1.Count > minselesctsize)
						break;
				}
			}
			
			return morethan1;
		}
		public void validatePopulation(ref List<chromosome> pop)
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
		//public List<List<chromosome>> SearchDesignSpace3()
		{
			#region reading settings
			int maxt = settings.ReadValueOfInteger(settings.RootAddress, "HyPE_InnerLoopMaxRuns");
			int PopulationSize = settings.ReadValueOfInteger(settings.RootAddress, "HypE_populationsize");
			int toursize = settings.ReadValueOfInteger (settings.RootAddress, "HypE_tournament");
			int mating = settings.ReadValueOfInteger(settings.RootAddress, "HypE_mating");
			int SelectionSize = settings.ReadValueOfInteger(settings.RootAddress, "HypE_SelectionSize");
			double bound = settings.ReadValueOfDouble(settings.RootAddress, "HypE_bound");
			int numberofSamples = settings.ReadValueOfInteger (settings.RootAddress, "HypE_nrOfSamples");
			int MAX_HYpE_Runs = settings.ReadValueOfInteger (settings.RootAddress, "HypE_OuterLoopMaxRuns") - 1;
			double P_Crossover = settings.ReadValueOfDouble(settings.RootAddress, "NSGA2_P_Crossover");
			double P_Mutation = settings.ReadValueOfDouble(settings.RootAddress, "NSGA2_P_Mutation");
			//double P_Mutation_maxDist = settings.ReadValueOfDouble(settings.RootAddress, "NSGA2_Mutation_maxDist");
			bool initial_import_type = settings.ReadValueOfString(settings.RootAddress, "initial_import_type")
				== "add";
			#endregion
			
			#region preparing
			//initialize population
			List<chromosome> pop_collection = new List<chromosome>();
			List<chromosome> pop_all = new List<chromosome>();
			List<chromosome> pop_sel = new List<chromosome>();
			//fronts
			List<List<chromosome>> front_part = new List<List<chromosome>>();
			//collections
			#endregion
			
			for(int HYpE_Counter = 0; HYpE_Counter < MAX_HYpE_Runs; HYpE_Counter++)
			{
				Console.WriteLine("\n\nrunning HYpE for " + (HYpE_Counter + 1).ToString() + "th times :");
				
				#region main loop
				//initializing the first generation is diiferent from others, because of 
				//  improving efficiancy of answers
				//List<chromosome> pop_new = GenerateInitialPopulation (PopulationSize);
				List<chromosome> pop_new = settings.ReadInitialPoints(settings.RootAddress, "initial_points", ref settings);
				if(pop_new == null)
					pop_new = GenerateNewPopulation(PopulationSize);
				if(initial_import_type)
					pop_new.AddRange(GenerateNewPopulation(PopulationSize - pop_new.Count));
				validatePopulation(ref pop_new);
				int t = 0;
				//double bestanswerminipc_prev = double.MaxValue;
				double bestanswerminipc = double.MaxValue;
				while(t < maxt)
				{
					double P_Mutation_maxDist = mutation_3_RandomGenerator.Next(0, Convert.ToInt32(Math.Round(Math.Min(
						(settings.region.Bottom - settings.region.Top) / 2,
						(settings.region.Right - settings.region.Left)/ 2))));
					pop_all = CrossowerAndMutation(ref pop_new, P_Crossover, P_Mutation, P_Mutation_maxDist);
					List<chromosome> tmp = GenerateNewPopulation(PopulationSize - pop_all.Count - pop_new.Count + 1);
					pop_all.AddRange(tmp);
					pop_all.AddRange(pop_new);
					pop_all = removerepetitiveanswers(ref pop_all);
					validatePopulation(ref pop_all);
					
					//assert
					front_part = generateFrontPartition(pop_all);
					front_part = environmentalSelection(front_part, PopulationSize, bound, numberofSamples);
					pop_all = cleanUpArchive(front_part, pop_all);
					
					pop_sel.Clear();
					matingSelection(front_part, ref pop_all, ref pop_sel, bound, SelectionSize, toursize, mating, pop_all.Count);
					pop_new = pop_sel;
					
					#region reporting
					double min = double.MaxValue;
					int index = -1;
					for(int i = 0; i < pop_new.Count; i++)
					{
						for(int j = 0; j < settings.dimentions; j++)
						{
							double v = pop_new[i].rank[j];
							if(v < min) min = v;
							index = i;
						}
						if(min < bestanswerminipc)
						{
							//index = i;
							bestanswerminipc = min;
						}
					}
                    t++;
					Console.WriteLine("iteration{0}: population={1}, minr={2},\tranks=[{3},{4},{5},{6},{7}]",
                                      t, PopulationSize, Math.Round(bestanswerminipc),
				                  Math.Round(pop_new[index].rank[0]),
				                  Math.Round(pop_new[index].rank[1]),
				                  Math.Round(pop_new[index].rank[2]),
				                  Math.Round(pop_new[index].rank[3]),
				                  Math.Round(pop_new[index].rank[4]));
					
					
					#endregion
				}
				#endregion
				
				pop_collection.AddRange(pop_new);
				pop_collection = removerepetitiveanswers(ref pop_collection);
			}
			
			//select best answers from all iterations
			front_part = generateFrontPartition(pop_collection);
			front_part = environmentalSelection(front_part, PopulationSize, bound, numberofSamples);
			pop_collection = cleanUpArchive(front_part, pop_collection);
			pop_sel.Clear();
			matingSelection(front_part, ref pop_collection, ref pop_sel, bound, SelectionSize, toursize, mating, pop_collection.Count);
			return pop_sel;
					
			//return pop_collection;
			//return generateFrontPartition(pop_collection);
		}
		
	}
}

