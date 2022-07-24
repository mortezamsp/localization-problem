using System;
using System.IO;
using System.Collections.Generic;
using Catfood.Shapefile;

namespace multiObjectiveSearch
{
	/// <summary>
	/// Description of runsettings.
	/// </summary>
	public class RunSettings
	{
		public RunSettings(string root_address_)
		{
			RootAddress = root_address_;
			dimentions = ReadValueOfInteger(RootAddress, "HypE_dimentions");
			ReadPoints();
			ReadRegion();
		}
		
		public int dimentions;
		public Shapefile[] shapes = null;
		public string RootAddress = null;
		public RectangleD region;
		public double[,] riverJoints = null;
		public List<List<PointD>> rivers = null;
		public int numJoints = 0;
		public List<PointD> jointPoints =  null;
		public string flow1 = null, flow2 = null;
		
		/// <summary>
		/// reads 'key's value from text file 'rootaddress'
		/// </summary>
		/// <param name="rootaddress"></param>
		/// <param name="key"></param>
		/// <returns>nit.MinValue : not found</returns>
		public int ReadValueOfInteger(string rootaddress, string key)
		{
			StreamReader sr = null;
			string s = null;
			try
			{
				sr = new StreamReader(rootaddress);
				s = sr.ReadToEnd();
			}
			catch
			{
				if(sr != null)
					sr.Close();
				//throw new Exception("FILEERR");
				return int.MinValue;
			}
			finally
			{
				if(sr!= null)
					sr.Close();
			}
			string[] lines = s.Split(new char[]{'\n','\r'}, StringSplitOptions.RemoveEmptyEntries);
			bool found = false;
			for(int i = 0; i < lines.Length; i++)
			{
				if(lines[i].StartsWith("//"))
					continue;
				string[] parms = lines[i].Split(new char[]{' ','\t',',',';',':'}, StringSplitOptions.RemoveEmptyEntries);
				if(parms[0] == key)
				{
					int ret = int.MinValue;
					try
					{
						ret = int.Parse(parms[1]);
						found = true;
					}
					catch
					{
						//
					}
					if(found)
						return ret;
				}
			}
			return int.MinValue;
		}
		
		public double ReadValueOfDouble(string rootaddress, string key)
		{
			StreamReader sr = null;
			string s = null;
			try
			{
				sr = new StreamReader(rootaddress);
				s = sr.ReadToEnd();
			}
			catch
			{
				if(sr != null)
					sr.Close();
				//throw new Exception("FILEERR");
				return double.MinValue;
			}
			finally
			{
				if(sr!= null)
					sr.Close();
			}
			string[] lines = s.Split(new char[]{'\n','\r'}, StringSplitOptions.RemoveEmptyEntries);
			bool found = false;
			for(int i = 0; i < lines.Length; i++)
			{
				if(lines[i].StartsWith("//"))
					continue;
				string[] parms = lines[i].Split(new char[]{' ','\t',',',';',':'}, StringSplitOptions.RemoveEmptyEntries);
				if(parms[0] == key)
				{
					double ret = double.MinValue;
					try
					{
						ret = double.Parse(parms[1]);
						found = true;
					}
					catch
					{
						//
					}
					if(found)
						return ret;
				}
			}
			return double.MinValue;
		}
		
		public List<chromosome> ReadInitialPoints(string rootaddress, string key, ref RunSettings settings)
		{
			string address = ReadValueOfString(rootaddress, key);
			if(address == null)
				return null;
			
			List<chromosome> population = new List<chromosome>();
			StreamReader sr = null;
			string s = null;
			try
			{
				sr = new StreamReader(address);
				s = sr.ReadToEnd();
			}
			catch
			{
				if(sr != null)
					sr.Close();
				//throw new Exception("FILEERR");
				return null;
			}
			finally
			{
				if(sr!= null)
					sr.Close();
			}
			string[] lines = s.Split(new char[]{'\n','\r'}, StringSplitOptions.RemoveEmptyEntries);
			bool found = false;
			for(int i = 0; i < lines.Length; i++)
			{
				if(lines[i].StartsWith("//"))
						continue;
				string[] parms = lines[i].Split(new char[]{' ','\t',',',';',':'}, StringSplitOptions.RemoveEmptyEntries);
				double p1=0, p2=0;
				try
				{
					p1 = double.Parse(parms[0]);
					p2 = double.Parse(parms[1]);
					
					found = true;
				}
				catch
				{
					//
				}
				if(found)
					//break;
					population.Add(new chromosome(p1,p2,ref settings));
			}
			return population;
		}
		
		public string ReadValueOfString(string rootaddress, string key)
		{
			StreamReader sr = null;
			string s = null;
			try
			{
				sr = new StreamReader(rootaddress);
				s = sr.ReadToEnd();
			}
			catch
			{
				if(sr != null)
					sr.Close();
				//throw new Exception("FILEERR");
				return null;
			}
			finally
			{
				if(sr!= null)
					sr.Close();
			}
			string[] lines = s.Split(new char[]{'\n','\r'}, StringSplitOptions.RemoveEmptyEntries);
			bool found = false;
			for(int i = 0; i < lines.Length; i++)
			{
				string[] parms = lines[i].Split(new char[]{' ','\t',',',';',':'}, StringSplitOptions.RemoveEmptyEntries);
				if(parms[0] == key)
				{
					if(lines[i].StartsWith("//"))
						continue;
					string ret = null;
					try
					{
						ret = parms[1];
						found = true;
					}
					catch
					{
						//
					}
					if(found)
						return ret;
				}
			}
			return null;
		}
		
		public void ReadPoints()
		{
			shapes = new Shapefile[6];
			
			//save rivers onto a list
			shapes[0] = new Shapefile("points\\river.sh");
			int shapec = 0;
			rivers = new List<List<PointD>>();
			foreach(ShapePolyLine l in shapes[0])
			{
				rivers.Add(new List<PointD>());
				for(int i=0; i<l.Parts[0].Length; i++)
					rivers[shapec].Add(l.Parts[0][i]);
				shapec++;
			}
			shapes[0].Dispose();
			
			//finding river joints count
			shapes[3] = new Shapefile("points\\f3.sh");
			jointPoints = new List<PointD>();
			foreach(Catfood.Shapefile.ShapePoint s in shapes[3])
			{
				jointPoints.Add(s.Point);
				numJoints ++;
			}
			/*
			//save riverjoints in a graph. graph[i,j] shows that river i  is connected to river j
			riverJoints = new double[shapec, shapec];
			for(int A = 0; A < numJoints; A++)
			{
				double mind = double.MaxValue;
				int[] pos = getLineOfPointIndex(jointPoints[A]);
				
				//for(int i = 0; i < rivers.Count; i++)
				//{
					for(int j = 0; j < rivers.Count; j++)
					{
						//if(i == j) continue;
						if(j == pos[0]) continue;
						for(int k = 1; k < rivers[j].Count; k++)
						{
							double d = Geometry.DistanceOfPointToLine(rivers[j][k], rivers[j][k-1], jointPoints[A]);//rivers[pos[0]][rivers[pos[0]].Count - 1]);
							if(Math.Abs(d) < 1E-5)
							{
								riverJoints[pos[0], j] = 1;
								break;
							}
							if(d < mind)mind=d;
						}
					}
					mind=mind;
				//}
			}
			*/
			shapes[3].Dispose();
			
			//rivers
			shapes[0] = new Shapefile("points\\river.sh");
			//water supply
			shapes[1] = new Shapefile("points\\f1.sh");
			//water polution
			shapes[2] = new Shapefile("points\\f2.sh");
			//river joints
			shapes[3] = new Shapefile("points\\f3.sh");
			//power lines
			shapes[4] = new Shapefile("points\\f4.sh");
			//roads
			shapes[5] = new Shapefile("points\\f5.sh");
			
			/*
			Shapefile s1 = new Shapefile("points\\river.sh");
			foreach (ShapePolyLine shape in s1)
		    {
				PointD[] p = shape.Parts[0];
				for(int i = 0; i < p.Length; i++)
					
		    }
			s1.Close();
			*/
			
			flow1 = ReadValueOfString(RootAddress, "flow1");
			flow2 = ReadValueOfString(RootAddress, "flow2");
			
			//return null;
		}
		
		private void ReadRegion()
		{
			double top = Math.Min(shapes[0].BoundingBox.Top,
			               Math.Min(shapes[1].BoundingBox.Top,
			               Math.Min(shapes[2].BoundingBox.Top,
			               Math.Min(shapes[3].BoundingBox.Top,
			               Math.Min(shapes[4].BoundingBox.Top,
			                                                   shapes[5].BoundingBox.Top)))));
			double buttom = Math.Max(shapes[0].BoundingBox.Bottom,
			               Math.Max(shapes[1].BoundingBox.Bottom,
			               Math.Max(shapes[2].BoundingBox.Bottom,
			               Math.Max(shapes[3].BoundingBox.Bottom,
			               Math.Max(shapes[4].BoundingBox.Bottom,
			                                                   shapes[5].BoundingBox.Bottom)))));
			double left = Math.Min(shapes[0].BoundingBox.Left,
			               Math.Min(shapes[1].BoundingBox.Left,
			               Math.Min(shapes[2].BoundingBox.Left,
			               Math.Min(shapes[3].BoundingBox.Left,
			               Math.Min(shapes[4].BoundingBox.Left,
			                                                   shapes[5].BoundingBox.Bottom)))));
			double right = Math.Max(shapes[0].BoundingBox.Right,
			               Math.Max(shapes[1].BoundingBox.Right,
			               Math.Max(shapes[2].BoundingBox.Right,
			               Math.Max(shapes[3].BoundingBox.Right,
			               Math.Max(shapes[4].BoundingBox.Right,
			                                                   shapes[5].BoundingBox.Right)))));
			region = new RectangleD(left,top,right,buttom);
		}
		
		public Shapefile GetShape(int shapeIndex)
		{
			switch(shapeIndex)
			{
				case 0: return new Shapefile("points\\river.sh"); 
				case 1: return new Shapefile("points\\f1.sh"); 
				case 2: return new Shapefile("points\\f2.sh"); 
				case 3: return new Shapefile("points\\f3.sh"); 
				case 4: return new Shapefile("points\\f4.sh"); 
				case 5: return new Shapefile("points\\f5.sh"); 
				default : return null;
			}
		}
		
		public int[] getLineOfPointIndex(PointD p)
		{
			for(int i = 0; i < rivers.Count; i++)
			{
				for(int j = 1; j < rivers[i].Count; j++)
				{
					double d = Geometry.DistanceOfPointToLine(rivers[i][j], rivers[i][j-1], p);
					if(Math.Abs(d) <= 1E-5)
						return new int[]{i, j};
				}
			}
			return null;
		}
	}
	
	public class point
	{
		public int x, y;
		point(int x_, int y_)
		{
			x = x_;
			y = y_;
		}
	}

}
