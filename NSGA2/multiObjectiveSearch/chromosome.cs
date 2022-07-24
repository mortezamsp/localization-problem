/*
 * Created by SharpDevelop.
 * User: pc
 * Date: 11/30/2017
 * Time: 11:44 PM
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.Collections.Generic;
using Catfood.Shapefile;

namespace multiObjectiveSearch
{
	/// <summary>
	/// Description of point.
	/// </summary>
	public class chromosome
	{
		/// <summary>
		/// x and y showing a point, that is the goal we want to optimize
		/// </summary>
		public double x = 0, y = 0;
		/// <summary>
		/// dimentions is number of objectives, number of fittness functions we should satisfy
		/// </summary>
		private int dimentions = 0;
		/// <summary>
		/// rank[] is an array, indicating output of each fitness function.
		/// 	length of rank[] array is equal to dimentions.
		/// 	ranks[i] is value of this answer for fitness functions i.
		/// </summary>
		public double[] rank = null;
		/// <summary>
		/// to integrate all values of rank[] into one value, we use totalRank.
		/// 	it shows total value of an answer around all dimentions.
		/// 	it is used when i want to compare two answers, without using multi-objective comparing functions.
		/// 	it could be average of rank[] values.
		/// </summary>
		public double totalRank = 0;
		/// <summary>
		/// settings contaning objective dimentions, region, ...
		/// </summary>
		private RunSettings settings;
		/// <summary>
		/// random generator used for creating new points
		/// </summary>
		private Random pointRandomGenerator = null;
		/// <summary>
		/// indicating which river and which sub-river is containg the point
		/// </summary>
		public double lineIndex = 0, riverIndex = 0;
		/// <summary>
		/// used for NSGA2 algorithm
		/// </summary>
		public double[] CrowdingDistance = null;
		public double Distance = 0;
		
		public chromosome(double x_, double y_, ref RunSettings settings_)
		{
			x = x_;
			y = y_;
			dimentions = settings_.dimentions;
			settings = settings_;
			rank = new double[dimentions];
			totalRank = 0;
			pointRandomGenerator = new Random();
			CrowdingDistance = new double[dimentions];
		}
		public chromosome(chromosome c, ref RunSettings settings_)
		{
			x = c.x;
			y = c.y;
			dimentions = settings_.dimentions;
			settings = settings_;
			rank = new double[dimentions];
			totalRank = 0;
			pointRandomGenerator = new Random();
			CrowdingDistance = new double[dimentions];
		}
		
		public string PrintCMDString()
		{
			string s = "";
			s = s + "x="+ x.ToString() + ",y=" + y.ToString() + ",ranks={";
			for(int i=0; i<dimentions - 1; i++)
				s = s + rank[i].ToString() + ",";
			s = s + rank[dimentions - 1].ToString() + "}";
			return s;
		}
		
		public string PrintRawStringComplete(string delim)
		{
			string s = "";
			s = s + x.ToString() + delim + y.ToString() + delim;
			for(int i=0; i<dimentions - 1; i++)
				s = s + rank[i].ToString() + delim;
			s = s + rank[dimentions - 1].ToString();
			return s;
		}
		
		public string PrintRawString(string delim)
		{
            //string s = "";
            //s = s + x.ToString() + delim + y.ToString() + delim;
            //return s;
            string s = "";
            s = s + x.ToString() + delim + y.ToString() + delim;
            for (int i = 0; i < dimentions - 1; i++)
                s = s + rank[i].ToString() + delim;
            s = s + rank[dimentions - 1].ToString();
            return s;
		}
		
		
		/// <summary>
		/// determining if a point is placed on a line, and if not, put it on the nearest line.
		/// </summary>
		/// <param name="pop"></param>
		/// <returns>-1 = point is not valid, 0 = point is valid</returns>
		public int ValidatePoint()
		{
			double mind = double.MaxValue;
			int mindline = 0, mindshape = 0;
			//int shapec = 0;
			
			//Shapefile sh = settings.GetShape(0);
			//foreach(Catfood.Shapefile.ShapePolyLine s in sh)
			for(int i = 0; i < settings.rivers.Count; i++)
			{
				//points are on rivers
				for(int j = 1; j <settings.rivers[i].Count; j++)
				{
					double d = Math.Abs(Geometry.DistanceOfPointToLine(settings.rivers[i][j],
					                                                   settings.rivers[i][j - 1],
					                                                   new PointD(x, y)));
					if(d < mind)
					{
						mind = d;
						mindline = j;
						mindshape = i;
						if(mind == 0)
							break;
					}
				}
				//shapec++;
			}
			//sh.Dispose();
			
			//if mind == 0, it means that this point is placed on one of lines
			if(mind == 0)
				return 0;
			
			//fixing point
			//change position of points, lie it on lines
			//shapec = 0;
			//sh = settings.GetShape(0);
			//foreach(Catfood.Shapefile.ShapePolyLine s in sh)
			//for(int i = 0; i < settings.rivers.Count; i++)
			//{
				//if(shapec < mindshape)
				//if(i < mindshape)
				//{
					//shapec++;
				//	continue;
				//}
				
				double dx = settings.rivers[mindshape][mindline].X - settings.rivers[mindshape][mindline - 1].X;
				double dy = settings.rivers[mindshape][mindline].Y - settings.rivers[mindshape][mindline - 1].Y;
				double b = settings.rivers[mindshape][mindline].Y - (dy / dx)*settings.rivers[mindshape][mindline].X;
				
				if(dx == 0)
				{
					x = settings.rivers[mindshape][mindline].X;
					if(settings.rivers[mindshape][mindline - 1].Y < settings.rivers[mindshape][mindline].Y)
						y = pointRandomGenerator.NextDouble()*Math.Abs(dy) + settings.rivers[mindshape][mindline - 1].Y;
					else
						y = pointRandomGenerator.NextDouble()*Math.Abs(dy) + settings.rivers[mindshape][mindline].Y;
				}
				else if(dy == 0)
				{
					if(settings.rivers[mindshape][mindline - 1].X < settings.rivers[mindshape][mindline].X)
						x = pointRandomGenerator.NextDouble()*Math.Abs(dx) + settings.rivers[mindshape][mindline - 1].X;
					else
						x = pointRandomGenerator.NextDouble()*Math.Abs(dx) + settings.rivers[mindshape][mindline].X;
					y = settings.rivers[mindshape][mindline].Y;
				}
				else
				{
					if(settings.rivers[mindshape][mindline - 1].X < settings.rivers[mindshape][mindline].X)
						x = pointRandomGenerator.NextDouble()*Math.Abs(dx) + settings.rivers[mindshape][mindline - 1].X;
					else
						x = pointRandomGenerator.NextDouble()*Math.Abs(dx) + settings.rivers[mindshape][mindline].X;
					y = (dy / dx)*x + b;
				}
				
				double d2 = Geometry.DistanceOfPointToLine(settings.rivers[mindshape][mindline],
				                                           settings.rivers[mindshape][mindline - 1],
				                                           new PointD(x, y));
				if(Math.Abs(d2) > 1E-5)
				{
					Console.WriteLine("error on fixing points");
					return -1;
				}
				//break;
				
				lineIndex = mindline;
				riverIndex = mindshape;
			//}
			//sh.Dispose();
			return 0;
		}
		
		/// <summary>
		/// calculates min distance from current point with important stations or points of map.
		/// </summary>
		/// <param name="answer"></param>
		/// <returns>min distances from important points or stations.</returns>
		public void EvaluatePoint()
		{
			
			//objective 1 : distance from water supply stations
			double mindists = double.MaxValue;
			Shapefile sh = settings.GetShape(1);
			foreach(Catfood.Shapefile.ShapePolygon s in sh)
			{
				double mind = double.MaxValue;
				for(int i = 1; i <s.Parts[0].Length; i++)
				{
					double d = Math.Abs(Geometry.DistanceOfPointToLine(
						new PointD(s.Parts[0][i].X, s.Parts[0][i].Y),
						new PointD(s.Parts[0][i - 1].X, s.Parts[0][i - 1].Y),
						new PointD(x, y)));
					if(d < mind)
						mind = d;
				}
				double d2 = Math.Abs(Geometry.DistanceOfPointToLine(
						new PointD(s.Parts[0][0].X, s.Parts[0][0].Y),
						new PointD(s.Parts[0][s.Parts[0].Length - 1].X, s.Parts[0][s.Parts[0].Length - 1].Y),
						new PointD(x, y)));
				if(d2 < mind)
					mind = d2;
				
				if(mind < mindists)
					mindists = mind;
			}
			sh.Dispose();
			
			//objective 2 : distance from water polution stations
			/*
			//double sumdistu = 0;
			double mindistu = double.MaxValue;
			sh = settings.GetShape(2);
			foreach(Catfood.Shapefile.ShapePolygon s in sh)
			{
				//double mind = double.MaxValue;
				for(int i = 1; i <s.Parts[0].Length; i++)
				{
					//ignore polution center, if it is after sencor
					//if(s.Parts[0][i].Y > y)
					//	continue;
					
					double d = Math.Abs(Geometry.DistanceOfPointToLine(
						new PointD(s.Parts[0][i].X, s.Parts[0][i].Y),
						new PointD(s.Parts[0][i - 1].X, s.Parts[0][i - 1].Y),
						new PointD(x, y)));
					//sumdistu += d;
					if(d < mindistu)
						mindistu = d;
				}
				double d2 = Math.Abs(Geometry.DistanceOfPointToLine(
						new PointD(s.Parts[0][0].X, s.Parts[0][0].Y),
						new PointD(s.Parts[0][s.Parts[0].Length - 1].X, s.Parts[0][s.Parts[0].Length - 1].Y),
						new PointD(x, y)));
				//sumdistu += d2;
				if(d2 < mindistu)
					mindistu = d2;
			}
			sh.Dispose();
			*/
			double mindistu = double.MaxValue;
			double sumdistu = 0;
			sh = settings.GetShape(2);
			foreach(Catfood.Shapefile.ShapePolygon s in sh)
			{
				//double mind = double.MaxValue;
				for(int i = 1; i <s.Parts[0].Length; i++)
				{
					//ignore polution center, if it is after sencor
					if(settings.flow1 == "top" && s.Parts[0][i].Y > y)
						continue;
					if(settings.flow1 == "buttom" && s.Parts[0][i].Y < y)
						continue;
					if(settings.flow2 == "left" && s.Parts[0][i].X > x)
						continue;
					if(settings.flow2 == "right" && s.Parts[0][i].X < x)
						continue;
					
					double d = Math.Abs(Geometry.DistanceOfPointToLine(
						new PointD(s.Parts[0][i].X, s.Parts[0][i].Y),
						new PointD(s.Parts[0][i - 1].X, s.Parts[0][i - 1].Y),
						new PointD(x, y)));
					sumdistu += 1 / Math.Pow(Math.E, d);
					//sumdistu += d;
					if(d < mindistu)
						mindistu = d;
				}
				double d2 = Math.Abs(Geometry.DistanceOfPointToLine(
						new PointD(s.Parts[0][0].X, s.Parts[0][0].Y),
						new PointD(s.Parts[0][s.Parts[0].Length - 1].X, s.Parts[0][s.Parts[0].Length - 1].Y),
						new PointD(x, y)));
				//sumdistu += d2;
				if(d2 < mindistu)
					mindistu = d2;
				sumdistu += 1 / Math.Pow(Math.E, d2);
			}
			sh.Dispose();
			
			//objective 3 : being placed after river joints
			/*
			// find all river joints before point
			//sh = settings.GetShape(3);
			//double mindj = double.MaxValue;
			//foreach(Catfood.Shapefile.ShapePoint s in sh)
			//{
				int numjoint = 0;
				//int[] pos = settings.getLineOfPointIndex(s.Point);
				int[] pos = settings.getLineOfPointIndex(new PointD(x, y));
				if(pos != null)
				{
					Queue<int>qlines = new Queue<int>();
					qlines.Enqueue(pos[0]);
					while(qlines.Count != 0)
					{
						int l = qlines.Dequeue();
						for(int i = 0; i < settings.rivers.Count; i++)
						{
							if(settings.riverJoints[i, l] == 1)
							{
								//if(settings.
								qlines.Enqueue(i);
								numjoint++;
							}
						}
					}
				}
				//double d = Math.Abs(Geometry.DistanceOfPointToPoint(new PointD(x, y), s.Point));
				//if(d < mindj)
				//	mindj = d;
			//}
			//sh.Dispose();
			*/
			//finding nearest joint point to river containing this point
			int[] pos = settings.getLineOfPointIndex(new PointD(x, y));
			double mindistj = double.MaxValue;
			for(int i = 0; i < settings.numJoints; i++)
			{
				double d = Math.Abs(Geometry.DistanceOfPointToLine(settings.rivers[pos[0]][pos[1]],
				                                                   settings.rivers[pos[0]][pos[1]-1],
				                                                   new PointD(x, y)));
				if(d < mindistj)
					mindistj = d;
			}
			
			//objective 4 : distance to power lines
			double mindistp = double.MaxValue;
			sh = settings.GetShape(4);
			foreach(Catfood.Shapefile.ShapePolyLine s in sh)
			{
				double mind = double.MaxValue;
				for(int i = 1; i <s.Parts[0].Length; i++)
				{
					//distance from line
					double d = Math.Abs(Geometry.DistanceOfPointToLine(
						new PointD(s.Parts[0][i].X, s.Parts[0][i].Y),
						new PointD(s.Parts[0][i - 1].X, s.Parts[0][i - 1].Y),
						new PointD(x, y)));
					if(d < mind)
						mind = d;
				}
				
				if(mind < mindistp)
					mindistp = mind;
			}
			sh.Dispose();
			
			//objective 5 : distance to roads
			double mindistr = double.MaxValue;
			sh = settings.GetShape(5);
			foreach(Catfood.Shapefile.ShapePolyLine s in sh)
			{
				double mind = double.MaxValue;
				for(int i = 1; i <s.Parts[0].Length; i++)
				{
					double d = Math.Abs(Geometry.DistanceOfPointToLine(
						new PointD(s.Parts[0][i].X, s.Parts[0][i].Y),
						new PointD(s.Parts[0][i - 1].X, s.Parts[0][i - 1].Y),
						new PointD(x, y)));
					if(d < mind)
						mind = d;
				}
				
				if(mind < mindistr)
					mindistr = mind;
			}
			sh.Dispose();
			
			//return Math.Sqrt(Math.Pow(x-points[objectiveFunctionIndex].x,2)+
			//                 Math.Pow(y-points[objectiveFunctionIndex].y,2));
			rank = new double[]{mindists, sumdistu, mindistj, mindistp, mindistr};
			//return new double[]{mindists+ sumdistu+ mindistj+ mindistp+ mindistr};
		}
	}
	
	
	
}
